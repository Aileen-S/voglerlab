#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument("-o", "--output", type=str, help="Output fasta")

parser.add_argument("-c", "--csv", type=str, help="CSV file: default new names in first column, old names in second, else use -t and -l flags")
parser.add_argument("-t", "--tips", type=str, help="Column with fasta IDs, if not second column")
parser.add_argument("-l", "--label", type=str, help="Column names required in labels, if not first column. Comma delimited list.")

#parser.add_argument("-k", "--list", action='store_true', help='Specify CSV second column is a comma delimited list of IDs')
parser.add_argument("-d", "--dups", action='store_true', help='Keep duplicate new IDs: ie, if multiple sequences have the same new ID, only change the name of the longest sequence.')
parser.add_argument("-r", "--renamed", type=str, help="Optional output csv with old and new names")

args = parser.parse_args()


#args = argparse.Namespace(input='test.fasta',csv='test.csv', tips='a', label='b,c,d') # This is how I step through the script interactively

with open(args.csv) as file:
    # Get column with old IDs
    tips = args.tips if args.tips else 1
    # Get column(s) with new IDs
    cols = args.label.split(',') if args.label else 0
    if cols:
        cols = [tips] + [col for col in cols if col != tips]
    meta = pd.read_csv(file, index_col=tips, usecols=cols)
    meta = meta.drop_duplicates()

    # Check for duplicate values in the specified tips column
    if meta.index.duplicated().any():
        duplicates = meta.index[meta.index.duplicated()]
        # Keep unique values
        duplicates = set(duplicates)
        print(f"Duplicate IDs in {tips} column: {list(duplicates)}")
print(meta)

if tips in cols:
    cols.remove(tips)

recs = {}
records = SeqIO.parse(args.input, "fasta")
for rec in records:
    print(rec.id)
    new = []
    new_id = rec.id # Default keep old ID
    # Check if ID is in metadata
    # Join values in chosen column(s) to get new ID
    r_id = rec.id.split(';')[0]
    if r_id in meta.index:
        for col in cols:
            new.append(str(meta.loc[r_id][col]))
        new_id = str('_'.join(new))
        # Append new ID to old ID, preserving frame tags
        if ';frame=' in rec.id:
            r_id, frame = rec.id.split(';')
            new_id = f'{r_id}_{new_id};{frame}'
        else:
            new_id = f'{r_id}_{new_id}'


    if new_id in recs:
        recs[new_id][rec.id] = rec.seq
    else:
        recs[new_id] = {rec.id: rec.seq}

# for k, v in recs.items():
#     print(k)
#     print(v)

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
        gaps = rec_seq.count('-') + rec_seq.count('N') + rec_seq.count('X') + rec_seq.count('*')
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
