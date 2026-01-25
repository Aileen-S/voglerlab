#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
from itertools import repeat
import multiprocessing
import time

start_time = time.time()


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


def split_into_chunks(records, num_chunks):
    for i in range(num_chunks):
        yield records[i::num_chunks]


def is_identical(seq_a, seq_b, min_overlap):
    overlap = [
        (a, b) for a, b in zip(seq_a, seq_b)
        if a != '-' and b != '-']
    if len(overlap) > min_overlap:
        return all(a == b for a, b in overlap)


def process_chunk(chunk, seq_dict, min_overlap):
    matches = []
    for rec in chunk:
        hits = []
        loop_dict = seq_dict.copy()
        del loop_dict[rec.id]
        for seq_id, seq in loop_dict.items():
            if is_identical(rec.seq, seq, min_overlap):
                hits.append(seq_id)
        if len(hits) > 0:
            hits.append(rec.id)
            matches.append(hits)
    return matches

def choose_best(matches, seq_dict, nuc):
    ambiguous = '[RYSWKMBDHVN-]' if nuc else '[BZJX-]'
    best = ''
    max = 0
    for match in matches:
        ambi_count = sum(1 for char in seq_dict[match] if char in ambiguous)
        if ambi_count > max:
            best = match
            max = ambi_count
    return best


def keep(records, keep_ids):
    # Add all keep records to output
    return records

# Argument parserremoved and 
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Remove duplicate sequences from fasta.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument('-o', '--output', type=str, help="Output fasta with duplicates removed")
parser.add_argument('-d', '--dups', type=str, help="Optional output file for removed duplicate sequences")
parser.add_argument('-m', '--min_overlap', type=int, help="Minumum overlap length for sequences to be compared and removed if identical (default 400 characters)")
parser.add_argument("-k", "--keep", type=str, help="File with list of IDs to keep, regardless of identical sequences")
parser.add_argument('-t', '--thread', type=str, help="Number of threads (default all available)")


args = parser.parse_args()


def main():

    # Read in and clean sequences
    with open(args.input) as file:
        records = list(SeqIO.parse(file, "fasta"))
    all_ids = set(rec.id for rec in records)
    if len(all_ids) < len(records):
        print('WARNING: Sequences IDs not all unique. This may cause errors in filtering.\n')
        sys.exit()
    nuc = is_nucleotide(records)
    clean = replace_ambiguous(records, nuc)
    seq_dict = {rec.id: rec.seq for rec in clean}

    # Split seqeunces into chunks
    num_threads = multiprocessing.cpu_count()
    chunks =  list(split_into_chunks(records, num_threads))


    # Check for identical sequences
    min_overlap = args.min_overlap if args.min_overlap else 400
    chunk_args = zip(chunks, repeat(seq_dict), repeat(min_overlap))
    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.starmap(process_chunk, chunk_args)
    duplicates = [item for sublist in results for item in sublist]


    # Combine match lists
    x = 0
    sorted_dups = set()
    for dup_list in duplicates:
        dup_list.sort()
        sorted_dups.add(tuple(dup_list))
    unique_dups = [list(dup) for dup in sorted_dups]
    dup_ids = [item for sublist in unique_dups for item in sublist]
    for dup_list in unique_dups:
        x += 1
        print(f'Duplicate sequence set {x}: {', '.join(dup_list)}')

    # Choose best (least uncertain characters) for each match set (or sequence on keep list, if present)
    if args.keep:
        keep_ids = [line.strip() for line in open(args.keep)]
    else:
        keep_ids = []

    chosen_ids = []
    for dup_list in unique_dups:
        if args.keep:
            if not any(seq_id in keep_ids for seq_id in dup_list):
                best = choose_best(dup_list, seq_dict, nuc)
                chosen_ids.append(best)
        else:
            best = choose_best(dup_list, seq_dict, nuc)
            chosen_ids.append(best)

    # Write to file
    final_ids = list(set([i for i in all_ids if i not in dup_ids] + chosen_ids + keep_ids))
    final_set = [rec for rec in records if rec.id in final_ids]

    with open(args.output, 'w') as output:
        SeqIO.write(final_set, output, "fasta")
    print(f"Saved {len(final_set)} records of {len(records)} from {args.input} to {args.output}")

    if args.dups:
        final_ids = [rec.id for rec in final_set]
        removed = [rec for rec in records if rec.id not in final_ids]
        with open(args.dups, 'w') as output:
            SeqIO.write(removed, output, "fasta")
        print(f"Saved {len(removed)} duplicate records from {args.input} to {args.dups}")

    end_time = time.time()
    elapsed = (end_time - start_time)
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60

    if hours > 0:
        print(f'Time elapsed: {hours} hours, {minutes} minutes, {seconds:.4f} seconds')
    elif minutes > 0:
        print(f'Time elapsed: {minutes} minutes, {seconds:.4f} seconds')
    else:
        print(f'Time elapsed: {seconds:.4f} seconds')


if __name__ == "__main__":
    main()

    

