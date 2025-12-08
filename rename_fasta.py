#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Rename sequences in fasta file from CSV")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--csv", type=str, help="CSV file, new names in first column, old names in second")
    parser.add_argument("-l", "--list", action='store_true', help='Specify CSV second column is a comma delimited list of IDs')
    parser.add_argument("-d", "--dups", action='store_true', help='Keep duplicate new IDs: '
                        'ie, if multiple sequences have the same new ID, only change the name of the longest sequence. '
                        'Default is to remove duplicates from fasta (keep longest sequence for each new ID)')
    parser.add_argument("-o", "--output", type=str, help="Output fasta file")
    parser.add_argument("-r", "--renamed", type=str, help="Optional output csv with old and new names")
    parser.add_argument("-a", "--amino_acid", type=str, help="Sequence is amino acid (default nucleotide). Needed to count gap characters in alignments.")

    return parser.parse_args()


recs = {}

def process_metadata(input_csv, list=False):
    # new names in first column, old names in second
    meta = {}
    with open(input_csv) as file:
        metadata = csv.reader(file)
        for row in metadata:
            if len(row) >= 2:
                if row [0] != '' and row[0] != 'NA':
                    if list:
                        meta[row[0]] = row[1].split(',')
                    else:
                        meta[row[1]] = row[0]

    return meta

def process_fasta_ids(input_fasta, meta, list=False):
    records = SeqIO.parse(input_fasta, "fasta")
    for rec in records:
        new_id = rec.id
        if ';frame=' in rec.id:
            r_id, frame = rec.id.split(';')
            if list:
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
            if list:
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
    return recs

def process_duplicate_ids(recs, dups=False, aa=False):
    # Check for duplicate new IDs
    selected = {}
    removed = {}
    if dups:
        print("Renaming longest sequence for each duplicate new name; keeping shorter sequences under old names")
    else:
        print("Renaming longest sequence for each duplicate new name; removing shorter sequences")
    for new, old in recs.items():
        max_len = -1
        max_rec = {}
        for rec_id, rec_seq in old.items():
            # Get sequence length
            rec_seq = rec_seq.upper()
            if aa:
                gaps = rec_seq.count('-') + rec_seq.count('X') + rec_seq.count('*')
            else:
                gaps = rec_seq.count('-') + rec_seq.count('N')
            length = len(rec_seq) - gaps
            # Find longest sequence
            if length > max_len:
                # Process previous max_rec if a longer match has been found
                if max_rec:
                    if dups:
                        selected[list(max_rec.keys())[0]] = max_rec
                    removed[rec_id] = new
                max_rec = {rec_id: rec_seq}
                max_len = length
            else:
                if dups:
                    selected[rec_id] = {rec_id: rec_seq}
                else:
                    removed[rec_id] = new

        selected[new] = max_rec
    return selected, removed

def write_fasta(selected, output_fasta):
    with open(output_fasta, 'w') as output:
        for new, rec in selected.items():
            for old, seq in rec.items():
                output.write(f'>{new}\n{seq}\n')

        print(f'Saved renamed fasta to {output_fasta}')

def write_csv(selected, removed, output_csv):
    with open(output_csv, 'w') as id_list:
        writer = csv.writer(id_list)
        writer.writerow(['old_id', 'new_id', 'removed'])
        for new, rec in selected.items():
            #print(f'{new}:     {rec}')
            for old, seq in rec.items():
                writer.writerow([old, new, ''])

        for old, new in removed.items():
            writer.writerow([old, new, 'yes'])
    print(f'Saved old and new fasta IDs to {output_csv}')


def main():
    args = parse_args()
    meta = process_metadata(args.csv, args.list)
    recs = process_fasta_ids(args.input, meta, args.list)
    selected, removed = process_duplicate_ids(recs, args.dups, args.amino_acid)
    write_fasta(selected, args.output)
    if args.renamed:
        write_csv(selected, removed, args.renamed)

if __name__ == "__main__":
    main()
