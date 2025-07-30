import sys
import os
import tempfile
import subprocess
from io import StringIO
import argparse, argcomplete
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


def parse_args():

    # Required arguments

    parser = argparse.ArgumentParser(description="Clean and align fasta sequences")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--check", type=str, help="Output fasta for rejected sequences")
    parser.add_argument("-g", "--good", type=str, help="Output fasta for good sequences as nucleotide ")
    parser.add_argument("-l", "--locus", choices=['mito', 'nuc', 'rna'], default='mito', help="Specify if sequences are invertebrate mitochondrial coding or nuclear coding, or RNA")

    # Optional arguments

    # Additional output
    parser.add_argument("-s", "--save", type=str, help="File to save initial alignment with all sequencs, before filtering")
    parser.add_argument("-a", "--aa_output", type=str, help="Output fasta for good sequences as protein")

    # Use profile for improved alignment. Sequences should have PROFILE in fasta ID
    parser.add_argument("-np", "--nt_profile", type=str, help="Nucleotide profile fasta")
    parser.add_argument("-ap", "--aa_profile", type=str, help="Amino acid profile fasta, if using translation option")

    # Specify threshold for consensus seqeunce
    parser.add_argument("-ct", "--consensus_threshold", type=float, help="Minimum proportion of characters to accept as consensus character. Default 0.7")
    # parser.add_argument("-mt", "--match_threshold", type=float, help="Mimumum similarity to consensus"
    #                                                                  "For coding sequences: acts on chunks of 20bp, default 0.5"
    #                                                                  "For RNA, works on whole sequence, default 0.7")
    parser.add_argument("-gt", "--gap_threshold", type=float, help="Maximum proportion gaps to accept. Default 0.95")

    # Ignore warning if no sequence passes match threshold
    parser.add_argument("-w", "--warning", action='store_true', help="If no sequences pass initial threshold checks, proceed with cleaning anyway. Default exit script.")

    # MAFFT options for final alignment, instead of default fast alignment
    parser.add_argument('-m', "--mafft_command", type=str, help="Custom MAFFT command")
    # Default with profile: 'mafft --addfragments input --retree 1 --maxiterate 0 --adjustdirection --anysymbol --thread autodetect profile
    # Default without profile: 'mafft --retree 1 --maxiterate 0 --adjustdirection --anysymbol --thread autodetect input'
    # Example custom input: 'mafft --addfragments input --globalpair --maxiterate 1000 --adjustdirection --anysymbol --thread autodetect profile'

    argcomplete.autocomplete(parser)
    return parser.parse_args()

#args = argparse.Namespace(input='test.fasta', check='check.fasta', good='good.fasta', translate=True, locus='mito', nt_profile='profiles/0_NT_profiles/COX1.fasta', aa_profile='profiles/0_AA_profiles/COX1.fasta')
def remove_empty_sequences(records):
    output = [rec for rec in records if rec.seq.count('-') < len(rec.seq)]
    if len(output) < len(records):
        print(f'Removed {len(records) - len(output)} empty sequences')
    return output


def translate(records, trans_table=5):
    print('Translating sequences to protein')
    aa_records = []
    stop_codons = ['TAA', 'TAG'] if trans_table == 5 else ['TAA', 'TAG', 'TGA']    # Add extra stop codon for nuclear DNA
    reading_frames = {}
    for rec in records:
        # Remove gaps
        seq = rec.seq.upper().replace('-', '')
        frame_totals = [0, 0, 0]
        # Find reading frame
        for f in range(3):
            x = 0
            codons = [seq[i:i+3] for i in range(f, len(seq), 3)]
            for codon in codons:
                if codon in stop_codons:
                    x += 1
            frame_totals[f] += x
        # Choose reading frame with fewest stop codons
        frame = frame_totals.index(min(frame_totals))
        reading_frames[rec.id] = frame
        # Add Ns to complete last codon
        leftover = len(seq[frame:]) % 3
        seq = seq[frame:] + ('N' * (3 - leftover))
        # Translate
        aa_rec = rec
        aa_rec.seq = seq.translate(table=trans_table)
        aa_records.append(aa_rec)
    return aa_records, reading_frames


def align_sequences(records, profile=False, command=False):
    # Write temporary input file for mafft
    print('Aligning sequences with MAFFT')
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        SeqIO.write(records, temp_input, "fasta")
        temp_input.flush()
        input = temp_input.name
    # Get thread count
    threads = str(os.cpu_count())
    threads = 1 if threads is None else threads
    # Call MAFFT
    if command:
        command = command.replace('autodetect', threads).replace('input', input).replace('profile', profile).split(' ')
    else:
        if profile:
            command = ['mafft', '--addfragments', input, '--retree', '1', '--maxiterate', '0', '--adjustdirection', '--anysymbol', '--thread', str(threads), profile]
        else:
            command = ['mafft', '--retree', '1', '--maxiterate', '0', '--adjustdirection', '--anysymbol', '--thread', str(threads), input]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        print("MAFFT failed:")
        print("stderr:", result.stderr)
        return []
    aligned_io = StringIO(result.stdout)
    results = list(SeqIO.parse(aligned_io, "fasta"))
    for rec in results:
        rec.id = rec.id.replace('_R_', '')
    # # Remove profile seqeunces
    # if profile:
    #     rec_ids = set(rec.id for rec in records)
    #     results = [rec for rec in results if rec.id in rec_ids]
    return results


def find_outliers(records, consensus_threshold, data, locus, chunk=False):
    print('Finding outliers')
    seq_list = list(records)
    for rec in seq_list:
        rec.seq = rec.seq.upper()
    alignment = MultipleSeqAlignment(seq_list)
    alignment_length = alignment.get_alignment_length()
    sequences = [record.seq for record in alignment]
    num_seqs = len(sequences)
    uncertain = 'N' if data == 'nt' else 'X'
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

    # Break in to chunks to find local errors
    if chunk:
        for rec in records:
            if 'PROFILE' not in rec.id:
                outlier = False
                seq = str(rec.seq)
                chunk_size = 10
                for i in range(0, len(pairs), chunk_size):
                    chunk = pairs[i:i+chunk_size]
                    # Get last bit
                    if len(chunk) < chunk_size:
                        chunk = pairs[-chunk_size:]
                    match_count = sum(1 for base, cons in chunk if base == cons)
                    match = match_count / 10
                    if match < 0.5:
                        check.append(rec)
                        outlier = True
                        break
            if outlier == False:
                good.append(rec)

    # Consider whole sequence at once
    else:
        # Calculate sequence variation within profile
        profile_similarity_consensus = []
        sequence_similarity_consensus = []
        if profile == True:
            for seq in profile_seqs:
                # Count bases matching consensus
                pairs = [(s, c) for s, c in zip(seq, consensus)]
                match_count = sum(1 for base, cons in pairs if base == cons)
                profile_similarity_consensus.append(match_count / len(pairs))
                 
        for rec in records:
            outlier = False
            seq = str(rec.seq)
            # Count bases matching consensus
            start = next(i for i, c in enumerate(seq) if c != '-')
            end = len(seq) - next(i for i, c in enumerate(reversed(seq)) if c != '-')
            pairs = [(s, c) for s, c in zip(seq[start:end], consensus[start:end])]
            match_count = sum(1 for base, cons in pairs if base == cons)
            match = match_count / len(pairs)
            sequence_similarity_consensus.append(match_count / len(pairs))
            if profile:
                if 'PROFILE' not in rec.id:
                    good.append(rec) if match_count/len(pairs) >= min(profile_similarity_consensus) else check.append(rec)
        if not profile:
            # Clusering approach if no profile
            print('No profile provided; using clustering to determine acceptance threshold')
            match_scores = np.array(sequence_similarity_consensus).reshape(-1, 1)
            gmm = GaussianMixture(n_components=2).fit(match_scores)
            labels = gmm.predict(match_scores)
            means = gmm.means_.flatten()
            high_match_label = np.argmax(means)
            high_match_indices = np.where(labels == high_match_label)[0]
            good = [records[i] for i in high_match_indices]
            check = [records[i] for i in range(len(records)) if i not in high_match_indices]

    print(f'Sequence similarity to consensus ranges from {min(sequence_similarity_consensus):.3f} to {max(sequence_similarity_consensus):.3f}')
    if profile:
        print(f'Profile similarity to consensus ranges from {min(profile_similarity_consensus):.3f} to {max(profile_similarity_consensus):.3f}')
        print(f'{len(good)} sequences passed match score threshold of {min(profile_similarity_consensus):.3f}')
    else:
        print(f'{len(good)} sequences in high match cluster')

    return check, good

# def gap_ratio_filter(records, consensus_threshold, match_thresold, data):
#     # Find sequneces with atypical gap ratio
#     print('Finding atypical gaps')
#     # Check for profile
#     if any('PROFILE' in rec.id for rec in records):
#         profile_seqs = [rec for rec in records if 'PROFILE' in rec.id]
#     else:
#         profile_seqs = [rec for rec in records]
#     # Calculate gap consensus
#     gap_counts = [sum(1 for seq in profile_seqs if seq[i] == '-') for i in range(len(records[0]))]
#     gap_freq = [count / len(records) for count in gap_counts]
#     gap_consensus = ''.join(['-' if i >= consensus_threshold else 'N' for i in gap_freq ])
#     gap_count_consensus = sum(1 for i in gap_freq if i >= consensus_threshold)
#     gap_ratio = gap_count_consensus / len(gap_freq)
    
#     # Find outlier sequences
#     good = []
#     check = []
#     uncertain = 'N' if data == 'nt' else 'X'
#     for rec in records:
#         if 'PROFILE' in rec.id:
#             good.append(rec)
#             continue
#         # Find sequence start and end
#         seq = rec.seq.replace(uncertain, '-')
#         start = next(i for i, c in enumerate(seq) if c != '-')
#         end = len(seq) - next(i for i, c in enumerate(reversed(seq)) if c != '-')
#         seq = seq[start:end]
#         con = gap_consensus[start:end]
#         # Match gaps to consensus
#         gap_matches = sum(1 for a, b in zip(seq, con) if a == '-' and b == '-')
#         char_matches = sum(1 for a, b in zip(seq, con) if a != '-' and b != '-')
#         match_ratio = (gap_matches + char_matches )/ (end - start)

#         if match_ratio >= match_thresold:
#             good.append(rec)
#         else:
#             check.append(rec)

#     print(f'{len(good)} sequences passed gap ratio filter')
#     return check, good


def find_internal_stop_codons(records, data, locus, reading_frames):
    print('Finding internal stop codons')
    good = []
    check = []
    for rec in records:
        end = len(rec.seq) - next(i for i, c in enumerate(reversed(rec.seq)) if c != '-')
        seq = rec.seq[:end]
        if data == 'nt':
            stop_codons = ['TAA', 'TAG'] if locus == 'mito' else ['TAA', 'TAG', 'TGA']
            frame = reading_frames[rec.id] - 1
            codons = [seq[i:i+3] for i in range(frame, len(seq), 3)]
            check.append(rec) if any(stop for stop in stop_codons in codons[:-1]) else good.append(rec)
        if data == 'aa':
            for rec in records:
                check.append(rec) if '*' in seq[:-1] else good.append(rec)
    print(f'{len(good)} sequences without internal stop codons')
    return check, good


def replace_partial_codons(records, reading_frames, trans_table=5):
    print('Searching for partial codons')
    check = []
    good = []
    x = 0
    y = 0
    for rec in records:
        if 'PROFILE' not in rec.id:
            frame = reading_frames[rec.id] - 1
            seq = ''
            has_stop_codon = False
            sequence = str(rec.seq).upper()
            # Split into codons according to reading frame
            codons = [sequence[i:i+3] for i in range(frame, len(sequence), 3)]
            # Check for partial codons
            for codon in codons:
                if '-' in codon and any(c.isalpha() for c in codon):
                    codon = codon.replace('-', 'N')
                    y += 1
                seq = seq + codon
            rec.seq = Seq(seq)
        check.append(rec)
    print(f'{y} partial codons corrected')
    return check


def delete_gappy_columns(records, gap_threshold):
    # Get gap percentage for each position in alignment
    print('Finding gappy columns')
    gap_threshold = gap_threshold if gap_threshold is not None else 1
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
    print(f'{len(records[0].seq) - len(keep)} gappy alignment positions removed')
    return good


def trim_to_profile(records):
    # Trim start to match profile and trim end to longest sequence
    print('Trimming alignment to match profile')
    length = len(records[0].seq)
    start = next((i for i, c in enumerate(records[0].seq) if c != '-'), None)
    if start != 0:
        for rec in records:
            rec.seq = rec.seq[start:]
    records = [rec for rec in records if 'PROFILE' not in rec.id]
    gap_counts = [sum(1 for rec in records if rec.seq[i] == '-') for i in range(len(records[0]))]
    # Find last column with characters
    for i in reversed(range(len(records[0].seq))):
        if gap_counts[i] < len(records):
            end = i + 1
            break
    # Remove gap only columns from end
    if end < len(records[0].seq):
        for rec in records:
            rec.seq = rec.seq[:end]
    if len(records[0].seq) < length:
        print(f'{length - len(records[0].seq)} alignment positions removed')
    return records


def main():
    print('Running clean_and_align.py')
    args = parse_args()

    # Check sequences
    nt_records= list(SeqIO.parse(args.input, 'fasta'))
    records = list(SeqIO.parse(args.input, 'fasta'))
    print(f'Found {len(records)} sequences in {args.input}')
    ids = set(rec.id for rec in records)
    if len(ids) < len(records):
        print('WARNING: Sequences IDs not all unique. This may cause errors in filtering.\n')

    records = remove_empty_sequences(records)
    shortest = min(len(rec.seq.replace('-', '')) for rec in records)
    longest = max(len(rec.seq.replace('-', '')) for rec in records)
    print(f'Sequence length ranges from {shortest} to {longest} characters')

    print('\nFiltering seqeunces')
    # Translate and align
    if args.locus == 'rna':
        aligned = align_sequences(records, args.nt_profile)
        data = 'nt'
        translation=False
    else:
        trans_table = 5 if args.locus == 'mito' else 1
        aa_recs, reading_frames = translate(records, trans_table)
        aligned = align_sequences(aa_recs, args.aa_profile)
        data = 'aa'
        translation=True

    if args.save:
        SeqIO.write(aligned, args.save, 'fasta')

    # Find outliers
    consensus_threshold = args.consensus_threshold if args.consensus_threshold else 0.7
    gap_threshold = args.gap_threshold if args.gap_threshold else 0.95
    # match_threshold = args.match_threshold if args.match_threshold else 0.7 if args.locus == 'rna' else 0.5

    check, good = find_outliers(aligned, consensus_threshold, data=data, locus=args.locus)
    if good == []:
        print('WARNING: No sequences passed filtering criteria\n'
                'Consider checking profile and alignment and adjusting thresholds')
        if args.warning:
            print('Will attempt to clean all sequneces\n')
            good = [rec for rec in aligned]
        else:
            sys.exit()

    # Skip stop codon filter for RNA
    if args.locus != 'rna':
        # Stop codon filter

        check_stop, good = find_internal_stop_codons(good, data=data, locus=args.locus, reading_frames=reading_frames)
        check.extend(check_stop)

    good_ids = [rec.id for rec in good]
    good_nt = [rec for rec in nt_records if rec.id in good_ids]

    if translation:
        if (len(check) > 0):
            print(f'{len(check)} sequences failed filtering criteria')


            # Realign 'check' sequences as nucleotide
            print('\nRealigning sequences')
            input_ids = [rec.id for rec in nt_records]
            check_ids = [rec.id for rec in check if rec.id in input_ids]
            check = [rec for rec in nt_records if rec.id in check_ids]
            check = align_sequences(check, args.nt_profile)

            # Remove gappy columns
            check = delete_gappy_columns(check, gap_threshold=args.gap_threshold)

            print('\nCleaning sequences')
            # Fill partial codons
            check_nt = replace_partial_codons(check, reading_frames, trans_table)

            # Align as AA again
            aa_recs, reading_frames = translate(check_nt, trans_table)
            if good != []:
                for rec in good:
                    rec.id = 'PROFILE::' + rec.id
                    aa_recs.append(rec)
            check = align_sequences(aa_recs, args.aa_profile)

            # Filter outliers again
            check, checked_good = find_outliers(check, consensus_threshold, data='aa', locus=args.locus)
            good_add_ids = [rec.id for rec in checked_good]
            good.extend(checked_good)
            good_ids = good_ids + [rec for rec in check_nt if rec.id in good_add_ids]

            # Align again, if more sequence were added to 'good' sequences
            if checked_good != []:
                print('\nFinal protein alignment')
                # AA output
                if args.mafft_command:
                    good_aa = align_sequences(good, args.aa_profile, command=args.mafft_command)
                else:
                    good_aa = align_sequences(good, args.aa_profile)
            else:
                good_aa = [rec for rec in good]
            good_aa = delete_gappy_columns(good_aa, gap_threshold)
            if args.aa_profile:
                good_aa = trim_to_profile(good_aa)

        else:
            good_aa = [rec for rec in aligned]

        if args.aa_output:
            with open(args.aa_output, 'w') as file:
                SeqIO.write(good_aa, file, 'fasta')
                print(f'{len(good_aa)} protein sequences with {len(good_aa[0].seq)} columns written to {args.aa_output}')
        # NT output
        good_ids = [rec.id for rec in good_aa if 'PROFILE' not in rec.id]
        good_nt = [rec for rec in nt_records if rec.id in good_ids]
   

    # Align as nucleotide

    if good_nt != []:
        print('\nFinal nucleotide alignment:')
        if args.mafft_command:
            good_nt = align_sequences(good_nt, args.nt_profile, command=args.mafft_command)
        else:
            good_nt = align_sequences(good_nt, args.nt_profile)

        good_nt = delete_gappy_columns(good_nt, gap_threshold=args.gap_threshold)
        if args.nt_profile:
            good_nt = trim_to_profile(good_nt)
        with open(args.good, 'w') as file:
            SeqIO.write(good_nt, file, 'fasta')
            print(f'\n{len(good_nt)} sequences with {len(good_nt[0].seq)} columns written to {args.good}')
    
    # Save sequences that did not pass filters
    if check != []:
        print('\nSaving seqeunces that did not pass filters')
        check_ids = [rec.id for rec in check]
        check = [rec for rec in nt_records if rec.id in check_ids]
        check = align_sequences(check, args.nt_profile)

        with open(args.check, 'w') as file:
            SeqIO.write(check, file, 'fasta')
            print(f'{len(check)} sequences written to {args.check}')
    else:
        print('All sequences passed cleaning filters')


if __name__ == "__main__":
    main()


# Add min sequence length?