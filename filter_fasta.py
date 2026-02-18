#!/usr/bin/env python3

from Bio import SeqIO, Phylo
import argparse
import re
import csv
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
    parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
    parser.add_argument("-f", "--found", type=str, help="Output file for sequences found in search filter")
    parser.add_argument("-n", "--notfound", type=str, help="Output file for input sequences not present in search filter (or removed with --ref_ids option)")

    # Search file options: default keep (-f) or remove (-n) all IDs found in file. More options available with CSV input.
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-l", "--list", type=str, help="File with list of IDs to search for.")
    group.add_argument("-s", "--fasta", type=str, help="Fasta file with IDs to search for.")
    group.add_argument("-t", "--tree", type=str, help="Newick tree with IDs to search for.")
    group.add_argument("-c", "--csv", type=str, help="CSV with IDs to search for (default first column, more options below).")

    # CSV options
    parser.add_argument("-a", "--id_column", type=str, help="Column with fasta IDs")
    parser.add_argument("-r", "--ref_column", type=str, help="Column with reference IDs (eg species name) - keep longest sequence for each reference ID")
    parser.add_argument("-k", "--keep", action='store_true', help="Keep IDs not present in CSV (default remove)")


    return parser.parse_args()


def read_csv(id_file, id_column, ref_column):
    with open(id_file) as file:
        reader = csv.reader(file)
        header = next(reader)
        id_index = 0
        if id_column:
            try:
                id_index = header.index(id_column)
            except ValueError:
                sys.exit(f'Error: {id_column} not found in {id_file} column names')
        if ref_column:
            try:
                ref_index = header.index(ref_column)
            except ValueError:
                sys.exit(f'Error: {ref_column} not found in {id_file} column names')
            id_map = {}
            for row in reader:
                if len(row) >= 2 and row[ref_index] != '' and row[ref_index] != 'NA':
                    id_map.setdefault(row[ref_index], []).append(row[id_index])
            return id_map

        else:
            id_list = [row[id_index].strip() for row in reader]
            return id_list


def read_id_file(id_file, format):
    filter = []
    with open(id_file) as infile:
        if format == 'list':
            lines = infile.readlines()
            for line in lines:
                if line != '\n':
                    filter.append(line.strip())

        elif format == 'fasta':
            recs = SeqIO.parse(id_file, "fasta")
            for rec in recs:
                filter.append(rec.id)

        elif format == 'tree':
            tree = Phylo.read(infile, "newick")
            taxa = tree.get_terminals()
            filter = [taxon.name for taxon in taxa]
    filter = [re.sub(r';frame=\d', '', f) for f in filter]
    print(f'{len(filter)} IDs in {id_file}')
    return filter


def filter_fasta(id_file, filter, input_fasta, found, not_found):
    with open(input_fasta) as file:
        records = list(SeqIO.parse(file, "fasta"))
    print(f'{len(records)} records in {input_fasta}')
    frametag = re.compile(r';frame=\d')

    # For each record, check if it is present in filter file (without frame tag)
    if found:
        found_records = (r for r in SeqIO.parse(input_fasta, "fasta") if frametag.sub('', r.id) in filter)
        count = SeqIO.write(found_records, found, "fasta")
        print(f"{count} of {len(filter)} records from {id_file} found in {input_fasta}, saved to {found}")
    if not_found:
        not_found_records = (r for r in SeqIO.parse(input_fasta, "fasta") if frametag.sub('', r.id) not in filter)
        count = SeqIO.write(not_found_records, not_found, "fasta")
        print(f"{count} records from {input_fasta} not present in {id_file}, saved to {not_found}")


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


def replace_ambiguous(sequence, nuc):
    ambiguous = '[-RYSWKMBDHVN]' if nuc else '[-BZJX]'
    seq = str(sequence.upper())
    seq = re.sub(ambiguous, '', seq)
    return seq


def keep_longest(input_fasta, id_map, id_file, found, not_found, keep):
    cleaned = {}
    longest = []
    checked = []
    with open(input_fasta) as file:
        records = list(SeqIO.parse(file, "fasta"))
    print(f'{len(records)} records in {input_fasta}')
    nuc = is_nucleotide(records)
    frametag = re.compile(r';frame=\d')
    for rec in records:
        cleaned[frametag.sub('', rec.id)] = {'seq': replace_ambiguous(rec.seq, nuc), 'seq_id': rec.id}
    for ref, id_list in id_map.items():
        max_len = 0
        max_ref = None
        for rec_id in id_list:
            try:
                if len(cleaned[rec_id]['seq']) > max_len:
                    max_len = len(cleaned[rec_id])
                    max_ref = rec_id
                checked.append(rec_id)
            except KeyError:
                continue
        if max_ref:
            longest.append(cleaned[max_ref]['seq_id'])

    if keep:
        notfound = [rec.id for rec in records if rec.id not in checked]
        selected = [rec for rec in records if rec.id in longest] + [rec for rec in records if rec.id in notfound]
    else:
        selected = [rec for rec in records if rec.id in longest]
    if found:
        count = SeqIO.write(selected, found, "fasta")
        if count > 0:
            print(f"{count} of {len(id_map)} records from {id_file} found in {input_fasta}, saved to {found}")
        else:
            print(f"{count} of {len(id_map)} records from {id_file} found in {input_fasta}")

    if not_found:
        selected_ids = [rec.id for rec in selected]
        removed = [rec for rec in records if rec.id not in selected_ids]
        if removed != []:
            count = SeqIO.write(removed, not_found, "fasta")
            print(f"{count} records from {input_fasta} removed, saved to {not_found}")

def main():
    args = parse_args()
    if args.csv:
        id_column = args.id_column if args.id_column else None

        if args.ref_column:
            keep = args.keep if args.keep else False
            id_map = read_csv(args.csv, id_column, args.ref_column)
            keep_longest(args.input, id_map, args.csv, args.found, args.notfound, keep)
        else:
            filter = read_csv(args.csv, id_column, None)
            filter_fasta(args.csv, filter, args.input, args.found, args.notfound)

    else:
        if args.list:
            filter = args.list
            format = 'list'
        elif args.fasta:
            filter = args.fasta
            format = 'fasta'
        elif args.tree:
            filter = args.tree
            format = 'tree'
        file = read_id_file(filter, format)

        filter_fasta(filter, file, args.input, found=args.found, not_found=args.notfound)
    

if __name__ == "__main__":
    main()