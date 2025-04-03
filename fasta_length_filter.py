#!/usr/bin/env python3

from Bio import SeqIO
import argparse, argcomplete
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output fasta file")
parser.add_argument("-l", "--length", type=str, help="Minimum sequence length")

argcomplete.autocomplete(parser)
args = parser.parse_args()

x = 0
y = 0
recs = {}

records = SeqIO.parse(args.input, "fasta")
with open(args.output, 'w') as output:
    for rec in records:
            rec_seq = rec.seq.upper()
            gaps = rec_seq.count('-') + rec_seq.count('N') + rec_seq.count('X') + rec_seq.count('*')
            length = len(rec_seq) - gaps
            if length >= int(args.length):
                x += 1
                recs[rec.id] = rec_seq
                output.write(f'>{rec.id}\n{rec_seq}\n')
            else:
                y += 1

print(f'{x} sequence(s) longer than {args.length} characters: written to {args.output}')
print(f'{y} sequence(s) shorter than {args.length} characters: removed')
