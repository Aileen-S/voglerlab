#!/usr/bin/env python3
import csv
import argparse
from Bio import SeqIO

# Argument parser
parser = argparse.ArgumentParser(description="Select longest sequences for PTP delimited species")
parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
parser.add_argument("-s", "--supermatrix", type=str, help="fasta with the tree's supermatrix")
parser.add_argument("-k", "--keep", type=str, help="file with list of taxa to keep (eg constraint, outgroup)")
parser.add_argument("-n", "--named", action='store_true', help="Keep all taxa with binomial (string after last underscore all lower case))")
parser.add_argument("-o", "--output", type=str, help="output file with filtered fasta IDs")

args = parser.parse_args()

# Get nucleotide ocunt for each taxon
count = {}
records = SeqIO.parse(args.supermatrix, "fasta")
x = 0
for rec in records:
    x += 1
    # Count gaps
    total = len(rec.seq)
    rec.seq = rec.seq.upper()
    gaps = rec.seq.count('-') + rec.seq.count('N') + rec.seq.count('X') + rec.seq.count('*')
    length = total - gaps
    count[rec.id] = length
print(f'{x} taxa in supermatrix')

keep = []
if args.keep:
    file = open(args.keep)
    lines = file.readlines()
    for line in lines:
        keep.append(line.strip)

file = open(args.input)

# Save mPTP delimited species lists
species_lists = []
temp = []
x = 0
lines = file.readlines()
for line in lines:
    if 'Number of delimited species' in line:
        print(line)
    x += 1
    if x < 10:
        continue
    if 'Species ' in line:
        continue
    line = line.strip()
    if line != '':
        temp.append(line)

        # Keep taxa with species/subspecies at end of id string.
        if args.named:
            parts = line.rsplit('_', 1)
            try:
                if parts[1].islower():
                    keep.append(line)
            except IndexError:
                pass
    else:
        species_lists.append(temp)
        temp = []

# Get taxon with most genes in each species list
chosen = []
for species in species_lists:
    nt = 0   # nucleotide count
    ch = ''  # chosen taxon
    lab = 0
    for s in species:
        if s in keep:
            chosen.append(s)
        else:
            try:
                if count[s] > nt:
                    nt = count[s]
                    ch = s
            except KeyError:
                print(f'{s} is not in metadata')

    if lab == 0:
        chosen.append(ch)

output = open(args.output, 'w')
for m in chosen:
    if m != '':
        output.write(f'{m}\n')

chosen = list(set(chosen))

print(f'{len(chosen)} taxa saved to {args.output}')
