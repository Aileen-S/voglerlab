#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re


def is_nucleotide(records):
    acgt  = ['A', 'C', 'G', 'T']
    nuc = 0
    other = 0
    for rec in records:
        seq = rec.seq.upper().replace('-', '')
        for base in seq:
            if base in acgt:
                nuc += 1
            else:
                other += 1
    return True if nuc > other else False


def replace_ambiguous(records, nuc):
    ambiguous = '[RYSWKMBDHVN]' if nuc else '[BZJX]'
    for rec in records:
        seq = str(rec.seq.upper())
        seq = re.sub(ambiguous, '-', seq)
        rec.seq = Seq(seq)
    return records


def remove_identical(records):
    sequences = []
    seq_ids = []

    for rec in records:
        if rec.seq not in sequences:
            sequences.append(rec.seq)
            seq_ids.append(rec.id)

    return seq_ids


def is_identical_overlap(seq_a, seq_b, min):
    overlap = [
        (a, b) for a, b in zip(seq_a, seq_b)
        if a != '-' and b != '-'
    ]
    if len(overlap) > min:
        return all(a == b for a, b in overlap)


# Argument parserremoved and 
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Remove duplicate sequences from fasta.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument('-o', '--output', type=str, help="Output fasta with duplicates removed")
parser.add_argument('-d', '--dups', type=str, help="Optional output file for removed duplicate sequences")
parser.add_argument('-m', '--min_overlap', type=int, help="Minumum overlap length for sequences to be compared and removed if identical (default 200 characters)")

args = parser.parse_args()         


def main():
    # Remove identical sequences
    with open(args.input) as file:
        records = list(SeqIO.parse(file, "fasta"))

    nuc = is_nucleotide(records)
    clean = replace_ambiguous(records, nuc)
    dedup = [rec for rec in records if rec.id in remove_identical(clean)]

    final_ids = []
    final_seqs = []

    min_overlap = args.min_overlap if args.min_overlap else 200
    for rec in dedup:
        for s in final_seqs:
            if is_identical_overlap(rec.seq, s, min_overlap):
                break
        else:
            final_ids.append(rec.id)
            final_seqs.append(rec.seq)

    output_records = [rec for rec in records if rec.id in final_ids]

    with open(args.output, 'w') as output:
        SeqIO.write(output_records, output, "fasta")
    print(f"Saved {len(output_records)} unique records of {len(records)} from {args.input} to {args.output}")

    if args.dups:
        dups = [rec for rec in records if rec.id not in final_ids]
        with open(args.dups, 'w') as output:
            SeqIO.write(dups, output, "fasta")
        print(f"{len(dups)} duplicate records saved to {args.dups}")
    else:
        print(f"{len(records) - len(output_records)} duplicate records removed")

if __name__ == "__main__":
    main()
