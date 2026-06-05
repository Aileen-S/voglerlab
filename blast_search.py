#!/usr/bin/env python3

import argparse
from Bio import Entrez, SeqIO
import subprocess
import time
import os


# Paths to the files and directories
parser = argparse.ArgumentParser(description="Search BLAST database with profile fasta and extract longest match for each record")
parser.add_argument("-db", "--database", type=str, help="Path to BLAST database")
parser.add_argument("-p", "--profile", type=str, help="BLAST search profile")
parser.add_argument("-o", "--output", type=str, help="Output fasta")
parser.add_argument("-s", "--hits", type=str, help="Save BLAST hits file")

args = parser.parse_args()

def blast(args):
    db = args.database
    query = args.profile
    form = '6 sacc sstart send'
    threads = str(os.cpu_count())
    command = ['blastn', '-db', db, '-query', query, '-max_target_seqs', '100000000', '-outfmt', form, '-num_threads', threads, '-evalue', '1e-5']
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        print("BLAST failed:")
        print("stderr:", result.stderr)
        return []
    if args.hits:
        with open(args.hits, 'w') as file:
            file.write(result.stdout)
            print(f'{len(result.stdout.splitlines())} blast hits written to {args.hits}')
    else:
        print(f'{len(result.stdout.splitlines())} blast hits')

    return result.stdout.splitlines()

def get_hits(blast_result):
    blast_hits = {}
    for line in blast_result:
        line = line.strip()
        record = line.split("\t") # Record format 'GBID start end'
        gbid = record[0]
        start = int(record[1])
        end = int(record[2])

        # Swap the values if the alignment is on the reverse strand
        if start > end:
            start, end = end, start

        # Check for too long matches (errors from genomes)
        if end - start > 4000:
            break

        if gbid in blast_hits:  
            # Check if start is lower than current start
            if start < blast_hits[gbid]['start']:
                blast_hits[gbid]['start'] = start
            # Check if end is higher than current end
            if end > blast_hits[gbid]['end']:
                blast_hits[gbid]['end'] = end
        else:
            blast_hits[gbid] = {'start': start, 'end': end}
    print(f'{len(blast_hits)} blast records with hits > 200bps')
    return blast_hits

def extract_sequences(blast_hits, args):
    with open(args.output, 'w') as file:
        print(f'Writing sequences to {args.output}')
        for gbid, rec in blast_hits.items():
            rec_len = int(rec['end']) - int(rec['start'])
            if rec_len > 200 and rec_len < 5000:
                command = ['blastdbcmd', '-db', args.database, '-entry', gbid, '-range', f'{rec["start"]}-{rec["end"]}', '-outfmt', '%f']
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode != 0 or not result.stdout.strip():
                    print(f"Failed to retrieve sequence: {gbid}")
                    print("stderr:", result.stderr)
                    return []
                file.write(result.stdout)


def main():
    blast_result = blast(args)
    blast_hits = get_hits(blast_result)
    extract_sequences(blast_hits, args)

if __name__ == "__main__":
    main()
