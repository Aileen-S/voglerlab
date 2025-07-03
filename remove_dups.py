
import argparse
from Bio import SeqIO


# Argument parserremoved and 
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Remove duplicate sequences from fasta.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument('-o', '--output', type=str, help="Output fasta with duplicates removed")
parser.add_argument('-d', '--dups', type=str, help="Optional output file for removed duplicate sequences")
parser.add_argument('-m', '--min_overlap', type=int, help="Minumum overlap length for sequences to be compared and removed if identical")


args = parser.parse_args()         

sequences = []
check = set()
seq_ids = []
dup_ids = []

initial_count = 0
# Remove identical sequences
with open(args.input) as file:
    record = SeqIO.parse(file, "fasta")
    for rec in record:
        initial_count += 1
        if rec.seq not in sequences:
            sequences.append(rec.seq)
            check.add(rec.seq)
            seq_ids.append(rec.id)
        else:
            dup_ids.append(rec.id)
records = (r for r in SeqIO.parse(args.input, "fasta") if r.id in seq_ids)

# Check for sequences with identical overlapping potrion
# New dict for sequences
final_seqs = {}
x = 0

# Threshold for minimum overlap length 
if args.min_overlap:
    min_overlap_length = args.min_overlap
else:
    min_overlap_length = 200

# Interrate throguh records
for rec in records:
    if x == 0:
        # Add first sequence to final_seqs
        final_seqs[rec.id] = rec.seq
        x += 1
        continue
    else:
        x += 1
        for seq_id in list(final_seqs.keys()):  # Use list to avoid modifying dict while iterating
            sequence_match = True
            overlap_count = 0
            for i, char in enumerate(rec.seq):
                if char == '-' or final_seqs[seq_id][i] == '-':
                    continue  # Skip gaps
                elif char != final_seqs[seq_id][i]:
                    sequence_match = False
                    break  # Exit the loop on mismatch
                else: 
                    overlap_count += 1 

            # Check if sequences have significant overlap
            if overlap_count < min_overlap_length: 
                sequence_match = False


            if sequence_match:
                # Count non-gap characters
                print(f'{seq_id} identical to {rec.id}')
                rec_non_gap_count = sum(c != '-' for c in rec.seq)
                final_non_gap_count = sum(c != '-' for c in final_seqs[seq_id])
                
                if rec_non_gap_count > final_non_gap_count:
                    # New sequence has fewer gaps, keep it
                    dup_ids.append(seq_id)
                    final_seqs[rec.id] = rec.seq
                    del final_seqs[seq_id]
                else:
                    # Existing sequence has fewer gaps, add the new sequence ID to duplicates
                    dup_ids.append(rec.id)
                break  # Exit loop after handling match
            else:
                # No match found, add the sequence to final_seqs
                final_seqs[rec.id] = rec.seq


count = len(final_seqs)
with open(args.output, 'w') as output:
    for rec_id, seq in final_seqs.items():
        output.write(f'>{rec_id}\n{seq}\n')
print(f"Saved {count} unique records of {initial_count} from {args.input} to {args.output}")

if args.dups:
    dups = (r for r in SeqIO.parse(args.input, "fasta") if r.id in dup_ids)
    count = SeqIO.write(dups, args.dups, "fasta")
    print(f"{count} duplicate records saved to {args.dups}")
else:
    print(f"{len(dup_ids)} duplicate records removed")
