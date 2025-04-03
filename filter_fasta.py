#!/usr/bin/env python3

from Bio import SeqIO, Phylo
import argparse, argcomplete
import re

parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
parser.add_argument("-f", "--found", type=str, help="Output file for sequences found in search filter")
parser.add_argument("-n", "--notfound", type=str, help="Output file for input sequences not present in search filter")
parser.add_argument("-s", "--search", type=str, help="File with IDs to search for. Specify file type with -t.")

parser.add_argument('-t', '--type', type=str, choices=['list', 'fasta', 'hmmer', 'tree'],
                    help="Search file type:\n"
                         "list: file with list of complete IDs\n" #or partial fasta 
                         "fasta: fasta file\n"
                         "hmmer: file with HMMER output from --tblout option\n")

argcomplete.autocomplete(parser)
args = parser.parse_args()         # Process input args from command line

filter = []
with open(args.search) as infile:
    if args.type == 'list':
        lines = infile.readlines()
        for line in lines:
            if line != '\n':
                filter.append(line.strip())

    elif args.type == 'fasta':
        recs = SeqIO.parse(args.search, "fasta")
        for rec in recs:
            filter.append(rec.id)

    elif args.type == 'hmmer':
        lines = infile.readlines()
        for line in lines:
            if '#' in line:
                continue
            line = line.split(' ')
            filter.append(line[0])

    elif args.type == 'tree':
        tree = Phylo.read(infile, "newick")
        # Extract taxa (terminal nodes)
        taxa = tree.get_terminals()
        filter = [taxon.name for taxon in taxa]

# Remove frame tags
filter = [re.sub(r';frame=\d', '', f) for f in filter]
print(f'{len(filter)} IDs in {args.search}')

records = SeqIO.parse(args.input, "fasta")
print(f'{len(list(records))} records in {args.input}')

# Define frame tag
frametag = re.compile(r';frame=\d')

# For each record, check if it is present in filter file (without frame tag)
# Write to output
if args.found:
    records = (r for r in SeqIO.parse(args.input, "fasta") if frametag.sub('', r.id) in filter)
    count = SeqIO.write(records, args.found, "fasta")
    print(f"{count} of {len(filter)} records from {args.search} found in {args.input}, saved to {args.found}")

if args.notfound:
    records = (r for r in SeqIO.parse(args.input, "fasta") if frametag.sub('', r.id) not in filter)
    count = SeqIO.write(records, args.notfound, "fasta")
    print(f"{count} records from {args.input} not present in {args.search}, saved to {args.notfound}")

print('\n')
