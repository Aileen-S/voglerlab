from Bio import SeqIO
from Bio import AlignIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.motifs import Motif
import argparse, argcomplete


def parse_args():
    parser = argparse.ArgumentParser(description="Filter out outlier sequences in alignmnet")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--check", type=str, help="Output fasta for sequences to be checked")
    parser.add_argument("-g", "--good", type=str, help="Output fasta for good sequences")

    # Find outlier sequences
    parser.add_argument("-ct", "--consensus_threshold", type=float, help="Minimum proportion of characters to accept as consensus character. Default 0.7")
    parser.add_argument("-mt", "--match_threshold", type=float, help="Mimumum similarity to consensus to accept. Default 0.7")

    parser.add_argument("-gt", "--gap_threshold", type=float, help="Maximum proportion gaps to accept" 
                        "(eg 0.95 would filter out sequences that have characters in a position where >95% of seqeunces have a gap)")


    argcomplete.autocomplete(parser)
    return parser.parse_args()


def find_outliers(fasta, consensus_threshold=0.7, match_threshold=0.7):
    alignment = AlignIO.read(fasta, "fasta")
    alignment_length = alignment.get_alignment_length()
    print(f'Alignment length: {alignment_length}')
    sequences = [str(record.seq) for record in alignment]
    num_sequences = len(sequences)

    columns_to_keep = []

    # Iterate through each column
    for col_index in range(alignment_length):
        column = [seq[col_index] for seq in sequences]
        # Check if the column contains only gap characters
        gap_only = all(char == '-' or char == 'N' for char in column)
        if not gap_only:
            columns_to_keep.append(col_index)

    # Reconstruct the alignment with only the columns to keep
    trimmed_seqs = []
    for record in alignment:
        trimmed_seq = ''.join(record.seq[i] for i in columns_to_keep)
        trimmed_seqs.append(SeqRecord(Seq(trimmed_seq), id=record.id, description=record.description))
    print(f'Trimmed {alignment_length - len(columns_to_keep)} gap-only columns from alignment')

    # msa = MultipleSeqAlignment(trimmed_seqs)
    # motif = Motif(alignment=msa)
    motif = motifs.create([record.seq for record in trimmed_seqs])
    counts = motif.counts
    consensus = counts.calculate_consensus(identity=consensus_threshold)
    consensus_seq = str(consensus)

    # Store distances from consensus
    good = []
    check = []
    for record in trimmed_seqs:
        seq = str(record.seq)
        
        # Count bases matching consensus
        match_count = 0
        positions = 0
        for i, base in enumerate(seq):
            if base in ('-', 'N') or consensus[i] in ('-', 'N'):
                continue
            positions += 1
            if base == consensus[i]:
                match_count += 1
            # Calculate percentage of non-gap, non-ambiguous characters matching consensus
            match = match_count / positions
        if match >= match_threshold:
            good.append(record)
        else:
            check.append(record)
    print(f'{len(good)} sequences passed {match_threshold * 100}% consensus threshold')

    return check, good



def find_atypical_gaps(records, gap_threshold):
    # Get gap percentage for each position in alignmentz
    gap_counts = [sum(1 for rec in records if rec.seq[i] == '-') for i in range(len(records[0]))]
    gap_freq = [count / len(records) for count in gap_counts]
    good = []
    check = []
    # Find outlier sequences
    x = 0
    for rec in records:
        # Replace N with gap in positions where most have gaps
        for i in range(len(rec.seq)):
            if rec.seq[i] == 'N' and gap_freq[i] > gap_threshold:
                rec.seq = rec.seq[:i] + '-' + rec.seq[i+1:]

        if any(rec.seq[i] != '-' and gap_freq[i] > gap_threshold for i in range(len(rec))):
            check.append(rec)
        else:
            good.append(rec)
        x += 1
    print(f'{len(good)} sequences passed gap filter')
    return check, good


def main():
    args = parse_args()

    # Check sequences
    records = list(SeqIO.parse(args.input, 'fasta'))
    print(f'Found {len(records)} sequences in {args.input}')
    lengths = set(len(rec) for rec in records)
    if len(lengths) > 1:
        print('Warning: Sequences are not all the same length. Please make sure sequences are aligned')

    # Find outliers
    consensus_threshold = args.consensus_threshold if args.consensus_threshold else 0.7
    match_threshold = args.match_threshold if args.match_threshold else 0.7
    check, good = find_outliers(args.input, consensus_threshold, match_threshold)

    # Gap threshold filter
    if args.gap_threshold:
        check1, good = find_atypical_gaps(good, gap_threshold=args.gap_threshold)
        check.extend(check1)

    # Save sequences without internal stop codons
    with open(args.good, 'w') as file:
        SeqIO.write(good, file, 'fasta')
        print(f'{len(good)} sequences saved to {args.good}')
        
    # Save sequences containing internal stop codons
    with open(args.check, 'w') as file:
        SeqIO.write(check, file, 'fasta')
        print(f'{len(check)} sequences written to {args.check}')


if __name__ == "__main__":
    main()
