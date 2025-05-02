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

keep = set()
if args.keep:
    file = open(args.keep)
    lines = file.readlines()
    for line in lines:
        keep.add(line.strip)

file = open(args.input)

# Save mPTP delimited species lists
species_lists = {}
temp = []
start_processing = False

lines = file.readlines()
for line in lines:
    if 'Number of delimited species' in line:
        print(line)

    # skip initial lines
    if 'Species ' in line:
        if temp:
            species_lists[species_no] = temp
        species_no = line.strip().replace('Species ', '').replace(':', '')
        temp = []
        start_processing = True
    else:
        # process non-empty lines
        if start_processing and line:

            line = line.strip()
            if line != '':
                temp.append(line)

                # Keep taxa with species/subspecies at end of id string.
                if args.named:
                    parts = line.rsplit('_', 1)
                    try:
                        if parts[1].islower():
                            keep.add(line)
                    except IndexError:
                        #print(f'IndexError: {line}')
                        pass
if temp:
    species_lists[species_no] = temp

# Get taxon with most genes in each species list
for species_no, species_list in species_lists.items():
    nt = 0   # nucleotide count
    ch = ''  # chosen taxon
    keepers = []
    others = []
    for spec in species_list:
        if spec in keep:
            keepers.append(spec)
        else:
            others.append(spec)
    # Find longest sequence in keep list (if used)
    for spec in keepers:
        try:
            if count[spec] > nt:
                nt = count[spec]
                longest = spec
        except KeyError:
            print(f'{s} is not in metadata')
    # Find longest sequence in others list
    for spec in others:
        try:
            if count[spec] > nt:
                nt = count[spec]
                longest = spec
        except KeyError:
            print(f'{s} is not in metadata')
    # Keep other if longer than longest in keep
    keep.add(longest)

output = open(args.output, 'w')
for m in keep:
    if m != '':
        output.write(f'{m}\n')

#chosen = list(set(chosen))
print(f'{len(keep)} taxa saved to {args.output}')
