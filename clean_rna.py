#!/usr/bin/env python3

import sys
import os
import tempfile
import subprocess
import shutil
import copy
from io import StringIO
import argparse
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.motifs import Motif
from sklearn.mixture import GaussianMixture
import numpy as np
import shlex


def parse_args():

    # Required arguments

    parser = argparse.ArgumentParser(description="Clean and align fasta sequences")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--check", type=str, help="Output fasta for rejected sequences")
    parser.add_argument("-g", "--good", type=str, help="Output fasta for good sequences")

    # Optional arguments

    # Additional output
    parser.add_argument("-s", "--save", type=str, help="File to save initial alignment with all sequencs, before filtering")

    # Use profile for improved alignment. Sequences should have PROFILE in fasta ID
    parser.add_argument("-p", "--profile", type=str, help="Profile fasta")

    # Specify thresholds
    parser.add_argument("-ct", "--consensus_threshold", type=float, help="Minimum proportion of characters to accept as consensus character. Default 0.7")
    parser.add_argument("-gt", "--gap_threshold", type=float, help="Maximum proportion gaps to accept. Default 0.95")
    parser.add_argument("-ml", "--min_length", type=float, help="Minimum sequence length. Default 100 bps")

    # Ignore warning if no sequence passes match threshold
    parser.add_argument("-w", "--warning", action='store_true', help="Ignore warnings, try to proceed with script")

    # MAFFT options for final alignment, instead of default fast alignment
    parser.add_argument('-m', "--mafft_command", type=str, help="Custom MAFFT command")
    # Default with profile: 'mafft --addfragments input --retree 1 --maxiterate 0 --nofft --parttree --adjustdirection --anysymbol --thread autodetect profile'
    # Default without profile: 'mafft --retree 1 --maxiterate 0 --nofft --parttree --adjustdirection --anysymbol --thread autodetect input'
    # Example custom input: 'mafft --addfragments input --globalpair --maxiterate 1000 --adjustdirection --anysymbol --thread autodetect profile'

    return parser.parse_args()


def remove_empty_sequences(records):
    output = [rec for rec in records if rec.seq.count('-') < len(rec.seq)]
    if len(output) < len(records):
        print(f'Removed {len(records) - len(output)} empty sequences')
    return output


def align_sequences(records, profile=False, command=False):
    print('Aligning sequences with MAFFT')
    try:
        # Set temp directory for mafft
        temp_scratch_dir = os.path.expanduser("~/scratch/tmp")
        os.makedirs(temp_scratch_dir, exist_ok=True)
        temp_dir = tempfile.mkdtemp(prefix="mafft_temp_", dir=temp_scratch_dir)
        input_file_path = os.path.join(temp_dir, "input.fasta")
        with open(input_file_path, "w") as temp_input:
            SeqIO.write(records, temp_input, "fasta")
        # Get thread count
        threads = str(os.cpu_count())
        threads = 1 if threads is None else threads

        # Call MAFFT
        if command:
            command = command.replace('autodetect', threads).replace('input', input_file_path).replace('profile', profile)#.split(' ')
            command = shlex.split(command)
        else:
            if profile:
                command = ['mafft', '--addfragments', input_file_path, '--retree', '1', '--maxiterate', '0', '--nofft', '--parttree', '--adjustdirection', '--anysymbol', '--thread', str(threads), profile]
            else:
                command = ['mafft', '--retree', '1', '--maxiterate', '0', '--nofft', '--parttree', '--adjustdirection', '--anysymbol', '--thread', str(threads), input_file_path]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout.strip():
            print("MAFFT failed:")
            print("stderr:", result.stderr)
            return []
        aligned_io = StringIO(result.stdout)
        results = list(SeqIO.parse(aligned_io, "fasta"))
        # for rec in results:
        #     rec.id = rec.id.replace('_R_', '')

    except Exception as e:
        print(f"Error: {e}")
        return []

    finally:
        # Remove tmp files
        if temp_dir and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    return results


def find_outliers(records, consensus_threshold):
    seq_list = list(records)
    for rec in seq_list:
        rec.seq = rec.seq.upper()
    alignment = MultipleSeqAlignment(seq_list)
    alignment_length = alignment.get_alignment_length()
    sequences = [record.seq for record in alignment]
    num_seqs = len(sequences)
    uncertain = 'N'
    profile = False
    # Check for profile
    if any('PROFILE' in rec.id for rec in alignment):
        profile_seqs = [rec for rec in alignment if 'PROFILE' in rec.id]
        if len(profile_seqs) > 50:
            profile = True
            print('Using profile sequences to calculate consensus')
    else:
        profile_seqs = [rec for rec in alignment]

    # Calculate consensus sequence
    consensus = []
    seq_len = len(profile_seqs[0])
    for i in range(seq_len):
        column = [seq[i] for seq in profile_seqs]
        counts = Counter(column)
        base, count = counts.most_common(1)[0]
        identity = count / len(column)
        consensus.append(base if identity >= consensus_threshold else uncertain)
    consensus = ''.join(consensus)

    good = []
    check = []
    for rec in records:
        rec.seq = rec.seq.replace(uncertain, '-')

    # Calculate sequence variation within profile
    profile_consensus_range = []
    sequence_consensus_range = []
    if profile == True:
        for seq in profile_seqs:
            # Count bases matching consensus
            pairs = [(s, c) for s, c in zip(seq, consensus)]
            match_count = sum(1 for base, cons in pairs if base == cons)
            profile_consensus_range.append(match_count / len(pairs))
        similarity_threshold = min(profile_consensus_range)# - 0.1 if min(profile_consensus_range) >= 0.6 else 0.5
                
    for rec in records:
        outlier = False
        seq = str(rec.seq)
        # Count bases matching consensus
        start = next((i for i, c in enumerate(seq) if c != '-'), None)
        if start == None:
            continue
        end = len(seq) - next(i for i, c in enumerate(reversed(seq)) if c != '-')
        pairs = [(s, c) for s, c in zip(seq[start:end], consensus[start:end])]
        match_count = sum(1 for base, cons in pairs if base == cons)
        match = match_count / len(pairs)
        sequence_consensus_range.append(match)
        if profile:
            if 'PROFILE' not in rec.id:
                good.append(rec) if match >= similarity_threshold else check.append(rec)
    if not profile:
        # Clusering approach if no profile
        print('No profile provided; using clustering to determine acceptance threshold')
        match_scores = np.array(sequence_consensus_range).reshape(-1, 1)
        gmm = GaussianMixture(n_components=2).fit(match_scores)
        labels = gmm.predict(match_scores)
        means = gmm.means_.flatten()
        high_match_label = np.argmax(means)
        high_match_indices = np.where(labels == high_match_label)[0]
        good = [records[i] for i in high_match_indices]
        check = [records[i] for i in range(len(records)) if i not in high_match_indices]

    print(f'Sequence similarity to consensus ranges from {min(sequence_consensus_range):.3f} to {max(sequence_consensus_range):.3f}')
    if profile:
        print(f'Profile similarity to consensus ranges from {min(profile_consensus_range):.3f} to {max(profile_consensus_range):.3f}')
        print(f'{len(good)} sequences passed match score threshold of {(similarity_threshold):.3f}')
    else:
        print(f'{len(good)} sequences in high match cluster')

    return check, good


def delete_gappy_columns(records, gap_threshold):
    # Get gap percentage for each position in alignment
    gap_threshold = gap_threshold if gap_threshold is not None else 1
    print(f'Deletiing columns >= {gap_threshold * 100:.0f}% gaps')
    gap_counts = [sum(1 for rec in records if rec.seq[i] == '-') for i in range(len(records[0]))]
    gap_freq = [count / len(records) for count in gap_counts]
    keep = [i for i in range(len(gap_freq)) if gap_freq[i] < gap_threshold]

    good = []
    # Find outlier sequences
    x = 0
    for rec in records:
        # Remove gappy columns
        rec.seq = Seq(''.join(rec.seq[i] for i in keep))
        good.append(rec)
    print(f'{len(records[0].seq) - len(keep)} columns removed')
    return good


def main():
    print('Running clean_and_align.py')
    args = parse_args()

    # Check sequences
    records = list(SeqIO.parse(args.input, 'fasta'))
    all_nt_records = [copy.deepcopy(rec) for rec in records]
    print(f'Found {len(records)} sequences in {args.input}')
    ids = set(rec.id for rec in records)
    if len(ids) < len(records):
        print('WARNING: Sequences IDs not all unique. This may cause errors in filtering.\n')
        if not args.warning:
            sys.exit()

    records = remove_empty_sequences(records)
    shortest = min(len(rec.seq.replace('-', '')) for rec in records)
    longest = max(len(rec.seq.replace('-', '')) for rec in records)
    print(f'Sequence length ranges from {shortest} to {longest} characters')
    min_length = int(args.min_length) if args.min_length else 100
    records = [rec for rec in records if len(rec.seq.replace('-', '')) >= min_length]
    print(f'Removed {len(all_nt_records) - len(records)} sequences shorter than {min_length} bps: {len(records)} sequences remaining')

    print('\nFiltering seqeunces')
    
    # Align
    aligned = align_sequences(records, args.profile)
    if args.save:
        SeqIO.write(aligned, args.save, 'fasta')
    for rec in aligned:
        rec.id = rec.id.replace('_R_', '')
        rec.description = ''

    # Find outliers
    consensus_threshold = args.consensus_threshold if args.consensus_threshold else 0.7
    gap_threshold = args.gap_threshold if args.gap_threshold else 0.95
    check, good = find_outliers(aligned, consensus_threshold)
    if good == []:
        print('WARNING: No sequences passed filtering criteria\n'
                'Consider checking profile and alignment and adjusting thresholds')
        if args.warning:
            print('Will attempt to clean all sequneces\n')
            good = [rec for rec in aligned]
        else:
            sys.exit()


    good_ids = [rec.id for rec in good]
    check_ids = [rec.id for rec in check]
    good_nt = [rec for rec in all_nt_records if rec.id in good_ids]


    # Align
    if good_nt != []:
        print('\nFinal nucleotide alignment:')
        if args.mafft_command:
            good_nt = align_sequences(good_nt, args.profile, command=args.mafft_command)
        else:
            good_nt = align_sequences(good_nt, args.profile)

        for rec in good_nt:
            rec.id = rec.id.replace('_R_', '')
            rec.description = ''
        good_nt = delete_gappy_columns(good_nt, gap_threshold=gap_threshold)
        good_nt = remove_empty_sequences(good_nt)
        good_nt  = [rec for rec in good_nt if 'PROFILE' not in rec.id]

        with open(args.good, 'w') as file:
            SeqIO.write(good_nt, file, 'fasta')
            print(f'\n{len(good_nt)} sequences with {len(good_nt[0].seq)} columns written to {args.good}')
    
    # Save sequences that did not pass filters
    if check != []:
        print('\nSaving seqeunces that did not pass filters')
        good_ids = [rec.id for rec in good_nt]
        check = [rec for rec in all_nt_records if rec.id not in good_ids]
        check = align_sequences(check, args.profile)
        check  = [rec for rec in check if 'PROFILE' not in rec.id]
        for rec in check:
            rec.id = rec.id.replace('_R_', '')
            rec.description = ''
        with open(args.check, 'w') as file:
            SeqIO.write(check, file, 'fasta')
            print(f'{len(check)} sequences written to {args.check}')
    else:
        print('All sequences passed cleaning filters')


if __name__ == "__main__":
    main()
