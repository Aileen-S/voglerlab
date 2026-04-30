#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from io import StringIO
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import time



def print_message(*message):
    print(time.strftime("%d-%m-%Y %H:%M:%S", time.gmtime()) + "\t" + " ".join(map(str, message)))

# Add option to provide input run list, rather than all in directory
def parse_args():
    parser = argparse.ArgumentParser(description="Align and clean BUSCO v6 output sequences. Requires MAFFT, CIAlign and trimAl")
    parser.add_argument("-i", "--input", type=str, help="Input directory containing completed BUSCO runs", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output directory to store results", required=True)
    parser.add_argument("-it", "--inclusion_threshold", type=float, help="Minimum fraction of BUSCO genes present for inclusion in supermatrix (default=0.5)", default=0.5)
    parser.add_argument("-tt", "--trimal_threshold", type=float, help="Minimum non-gap fraction per column for inclusion in supermatrix (default=0.5)", default=0.5)
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use")
    args = parser.parse_args()

    return args

def find_BUSCOs(input_directory):
    # Run directory paths list
    busco_runs = []
    # Run directory names list
    busco_run_names = []
    for d in os.listdir(input_directory):
        d = os.path.join(input_directory, d)
        print(d)
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.startswith("run_"):
                    busco_runs.append(os.path.join(d, f))
                    busco_run_names.append(os.path.basename(d))
                    break

    if len(busco_runs) == 0:
        print_message("ERROR. Didn't find any BUSCO directories")
        sys.exit()
    else:
        print_message("Found", len(busco_runs), "BUSCO directories:")
    
    return busco_runs, busco_run_names

def find_sequences(busco_runs, busco_run_names):
    # BUSCO genes list
    all_buscos = set()
    # List of run paths for each BUSCO gene
    buscos_complete_per_gene = defaultdict(list)
    buscos_fragmented_per_gene = defaultdict(list)
    # List of BUSCO genes for each run
    buscos_per_run = defaultdict(list)

    # Find single-copy BUSCOs
    for busco_run_name, busco_run in zip(busco_run_names, busco_runs):
        buscos_per_run[busco_run_name] = []
        path = Path(busco_run) / "busco_sequences" / "single_copy_busco_sequences"

        for f in os.listdir(path):
            if f.endswith(".faa"):
                busco_gene = f.replace(".faa", "")
                all_buscos.add(busco_gene)
                buscos_per_run[path].append(busco_gene)
                buscos_complete_per_gene[busco_gene].append(path)

    # Find fragmented BUSCOs
    for busco_run_name, busco_run in zip(busco_run_names, busco_runs):
        buscos_per_run[busco_run_name] = []
        path = Path(busco_run) / "busco_sequences" / "fragmented_busco_sequences"

        for f in os.listdir(path):
            if f.endswith(".faa"):
                busco_gene = f.replace(".faa", "")
                all_buscos.add(busco_gene)
                buscos_per_run[path].append(busco_gene)
                buscos_fragmented_per_gene[busco_gene].append(path)
    
    print_message(f"Found {len(all_buscos)} unique BUSCO genes")

    return all_buscos, buscos_per_run, buscos_complete_per_gene, buscos_fragmented_per_gene


def align_busco(gene, buscos_complete_per_gene, buscos_fragmented_per_gene, threads):
    try:
        complete_records = []
        for path in buscos_complete_per_gene[gene]:
            record = SeqIO.read(f"{path}/{gene}.faa", "fasta")
            taxon_dir = path.parent.parent.parent
            taxon_name = taxon_dir.name
            record.id = taxon_name
            complete_records.append(record)
        fragmented_records = []
        for path in buscos_fragmented_per_gene[gene]:
            record = SeqIO.read(f"{path}/{gene}.faa", "fasta")
            taxon_dir = path.parent.parent.parent
            taxon_name = taxon_dir.name
            record.id = taxon_name
            fragmented_records.append(record)

        fragments_only = True if len(complete_records) == 0 else False
        if fragments_only:
            complete_records = fragmented_records
            print(f'Aligning BUSCO gene {gene} with MAFFT: 0 complete sequences, {len(complete_records)} fragmented sequences')
        else:
            print(f'Aligning BUSCO gene {gene} with MAFFT: {len(complete_records)} complete sequences, {len(fragmented_records)} fragmented sequences')

        # Temp directory for mafft
        os.makedirs("tmp", exist_ok=True)

        with tempfile.TemporaryDirectory(prefix="mafft_temp_", dir="tmp") as temp_dir:
            if complete_records != []:
                complete_path = os.path.join(temp_dir, "complete.fasta")
                with open(complete_path, "w") as temp_input:
                    SeqIO.write(complete_records, temp_input, "fasta")

                # Align complete genes with MAFFT
                command = ['mafft', '--auto', '--anysymbol', '--thread', str(threads), complete_path]
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode != 0 or not result.stdout.strip():
                    print("MAFFT failed:")
                    print("stderr:", result.stderr)
                    return []
                aligned_io = StringIO(result.stdout)
                results = list(SeqIO.parse(aligned_io, "fasta"))
            
                if len(fragmented_records) == 0 or fragments_only:
                    return results
                
                aligned_path = os.path.join(temp_dir, "aligned.fasta")
                with open(aligned_path, "w") as temp_input:
                    SeqIO.write(results, temp_input, "fasta")

                fragments_path = os.path.join(temp_dir, "fragments.fasta")
                with open(fragments_path, "w") as temp_input:
                    SeqIO.write(fragmented_records, temp_input, "fasta")

                # Add fragmented genes
                command = ['mafft', '--addfragments', fragments_path, '--retree', '1', '--maxiterate', '0', '--anysymbol', '--thread', str(threads), aligned_path]
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode != 0 or not result.stdout.strip():
                    print("MAFFT failed:")
                    print("stderr:", result.stderr)
                    return []
                aligned_io = StringIO(result.stdout)
                results = list(SeqIO.parse(aligned_io, "fasta"))
                return results
    
    except Exception as e:
        print(f"Error: {e}")
        return []

def trimal(gene, alignment, trimal_threshold):
    with tempfile.TemporaryDirectory(prefix="trimal_temp_", dir="tmp") as temp_dir:
        path = os.path.join(temp_dir, "input.fasta")
        with open(path, "w") as temp_input:
            SeqIO.write(alignment, temp_input, "fasta")

        command = ["trimal", "-in", path, "-gt", str(trimal_threshold)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout.strip():
            print("trimAl failed:")
            print("stderr:", result.stderr)
            return []
        aligned_io = StringIO(result.stdout)
        results = list(SeqIO.parse(aligned_io, "fasta"))

    # Write to file here?
    return results


def process_buscos(gene, buscos_complete_per_gene, buscos_fragmented_per_gene, threads, trimal_threshold, inclusion_threshold, busco_runs):
    try:
        runs_for_gene = len(buscos_complete_per_gene[gene]) + len(buscos_fragmented_per_gene[gene])
        fraction = runs_for_gene / busco_runs
        if fraction >= inclusion_threshold:
            aligned = align_busco(gene, buscos_complete_per_gene, buscos_fragmented_per_gene, threads)
            trimmed = trimal(gene, aligned, trimal_threshold)
            return gene, trimmed
        else:
            print(f'Skipping {gene} because it is present in less than {int(inclusion_threshold * 100)} of taxa')
            return gene, None

    except Exception as e:
        print(f"Error: {e}")
        return gene, None


def process_buscos_parallel(busco_runs, busco_genes, buscos_complete_per_gene, buscos_fragmented_per_gene,
                            threads, trimal_threshold=0.5, inclusion_threshold=0.5):

    parallel_jobs = int(threads / 2) if threads > 1 else 1
    threads_per_job = int(threads / parallel_jobs)
    print(f"Starting BUSCO alignments: {parallel_jobs} parallel jobs with {threads_per_job} threads each")

    tasks = [
        (gene, buscos_complete_per_gene, buscos_fragmented_per_gene,
         threads_per_job, trimal_threshold, inclusion_threshold, len(busco_runs))
        for gene in busco_genes
    ]

    completed, failed = 0, 0
    alignments = {}

    with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
        futures = {executor.submit(process_buscos, *t): t[0] for t in tasks}
        for future in as_completed(futures):
            gene, trimmed = future.result()
            if trimmed is None:
                failed += 1
            else:
                completed += 1
                alignments[gene] = trimmed

    print(f"Done: {completed} succeeded, {failed} failed")
    return alignments


def concatenate(busco_run_names, alignments):
    supermatrix = []
    partitions = {}
    busco_genes = list(alignments.keys())
    first = True
    p = 0

    for run in busco_run_names:
        record = SeqRecord(Seq(""), id=run)
        for gene in busco_genes:
            found = False
            for seq in alignments[gene]:
                if seq.id == run:
                    record.seq += str(seq.seq)
                    found = True
                    break
            
            # If taxon is missing from this gene, add gaps
            if not found:
                length = len(alignments[gene][0].seq)
                record.seq  += "-" * length

            if first:
                length = len(record.seq)
                start = p + 1
                end = p + length
                partitions[gene] = {'start' : start,
                                'end' : end}
                p += length
                print(p)

        first = False
                
        supermatrix.append(record)

    return supermatrix, partitions


def main():

    # Get input options
    args = parse_args()
    threads = str(os.cpu_count())
    threads = 1 if threads is None else threads
    input_directory = os.path.abspath(args.input)
    output_directory = os.path.abspath(args.output)
    print_message(f"Starting BUSCO output processing with {threads} threads")
    print_message(f"Looking for BUSCO runs in {input_directory}")
    print_message(f"Keeping BUSCO genes present in at least {args.inclusion_threshold} of taxa")
    print_message(f"Keeping columns with at least {args.trimal_threshold} non-gap fraction in supermatrix")

    # Find BUSCO runs
    busco_runs, busco_run_names = find_BUSCOs(input_directory)
    total_runs = len(busco_runs)
    all_buscos, buscos_per_run, buscos_complete_per_gene, buscos_fragmented_per_gene = find_sequences(busco_runs, busco_run_names)

    # Get thread count
    threads = args.threads if args.threads else os.cpu_count()
    threads = 1 if threads is None else threads

    alignments = process_buscos_parallel(busco_runs, all_buscos, buscos_complete_per_gene, buscos_fragmented_per_gene, threads, args.trimal_threshold, args.inclusion_threshold)

    # Write to file
    output_path = os.makedirs(output_directory, exist_ok=True)
    for aln in alignments:
        with open(f"{output_directory}/{aln}.fasta", "w", ) as output:
            SeqIO.write(alignments[aln], output, "fasta")

    supermatrix, partitions = concatenate(busco_run_names, alignments)

    with open(f"{output_directory}/supermatrix.fasta", "w", ) as output:
        SeqIO.write(supermatrix, output, "fasta")

    # RAxML format, eg: "AA, 84030at7041=1-984"
    with open (f"{output_directory}/partitions.txt", "w") as output:
        for busco_gene, coordinates in partitions.items():
            output.write(f"AA, {busco_gene} = {coordinates['start']}-{coordinates['end']};\n")


if __name__ == "__main__":
    main()

