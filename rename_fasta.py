#!/usr/bin/env python3

from Bio import SeqIO
import argparse, argcomplete
import csv

parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-c", "--csv", type=str, help="CSV file, new names in first column, old names in second")
parser.add_argument("-l", "--list", action='store_true', help='Specify CSV second column is a comma delimited list of IDs')
parser.add_argument("-d", "--dups", action='store_true', help='Keep duplicates: ie, if multiple sequences have the same new ID, only change the name of the longest sequence.')
parser.add_argument("-o", "--output", type=str, help="Output fasta file")
parser.add_argument("-r", "--renamed", type=str, help="Optional output csv with old and new names")

argcomplete.autocomplete(parser)
args = parser.parse_args()


recs = {}

if args.csv:
    # new names in first column, old names in second
    meta = {}
    with open(args.csv) as file:
        metadata = csv.reader(file)
        for row in metadata:
            if row [0] != '' and row[0] != 'NA':
                if args.list:
                    meta[row[0]] = row[1].split(',')
                else:
                    meta[row[1]] = row[0]

    records = SeqIO.parse(args.input, "fasta")
    for rec in records:
        new_id = rec.id
        if ';frame=' in rec.id:
            r_id, frame = rec.id.split(';')
            if args.list:
                # For list mode, check if r_id is in meta keys
                for k, v in meta.items():
                    if r_id in v:
                        new_id = f'{k};{frame}'
                        break  # Break once found
            else:
                # For non-list mode, check exact match
                if r_id in meta:
                    new_id = f'{meta[r_id]};{frame}'
        else:
            if args.list:
                for k, v in meta.items():
                    if rec.id in v:
                        new_id = k
            else:
                if rec.id in meta:
                    new_id = meta[rec.id]
        if new_id in recs:
            recs[new_id][rec.id] = rec.seq
        else:
            recs[new_id] = {rec.id: rec.seq}


# Check for duplicate names
# Save longest sequence for duplicates
selected = {}
removed = {}

if args.dups:
    print("Renaming longest sequence for each duplicate new name; keeping shorter sequences under old names")
else:
    print("Renaming longest sequence for each duplicate new name; removing shorter sequences")

for new, old in recs.items():
    x = 0
    for rec_id, rec_seq in old.items():
        # Get sequence length
        rec_seq = rec_seq.upper()
        gaps = rec_seq.count('-') + rec.seq.count('N') + rec.seq.count('XX')
        length = len(rec_seq) - gaps

        if x == 0:
            max_rec = {rec_id: rec_seq}
            max_len = length
            x += 1
        else:
            # Find longest sequence
            if length > max_len:
                # Keep short records under old ID if dups option is used
                if args.dups:
                    selected[list(max_rec.keys())[0]] = max_rec
                else:
                    removed[new] = max_rec
                max_rec = {rec_id: rec_seq}
                max_len = length
            else:
                if args.dups:
                    selected[rec_id] = {rec_id: rec_seq}
                else:
                    removed[new] = {rec_id: rec_seq}
    selected[new] = max_rec


if args.renamed:
    with open(args.renamed, 'w') as id_list:
        writer = csv.writer(id_list)
        writer.writerow(['Old ID', 'New ID', 'Removed'])
        for new, rec in selected.items():
            #print(f'{new}:     {rec}')
            for old, seq in rec.items():
                writer.writerow([old, new])
        for new, old in removed.items():
            writer.writerow([old, new, 'yes'])
    print(f'Saved original IDs to {args.renamed}')


with open(args.output, 'w') as output:
    for new, rec in selected.items():
        for old, seq in rec.items():
            output.write(f'>{new}\n{seq}\n')

    print(f'Saved renamed fasta to {args.output}')
