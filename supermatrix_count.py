#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import csv
import string

def parse_args():
    parser = argparse.ArgumentParser(description="Count number of genes and number of bases per gene in supermatrix")
    parser.add_argument("-f", "--fasta", type=str, help="Input supermatrix fasta")
    parser.add_argument("-p", "--partitions", type=str, help="optional catfasta2phyml output partition file")
    parser.add_argument("-o", "--output", type=str, help="Output csv file")

    return parser.parse_args()


# Read partition file
def read_partitions(partition_file):
    with open(partition_file) as file:
        partitions = {}
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            part, crds = line.split(' = ')
            # Get partition location and length
            crds = crds.split('-')
            start = int(crds[0]) - 1
            end = int(crds[1])
            length = end - start
            # Get partition name
            part = part.split('/')[-1]
            partitions[part] = [start, end, length]

    print(f'{len(list(partitions.keys()))} partitions in partition file')
    return partitions

def process_fasta(input_fasta, partitions):
    records = SeqIO.parse(input_fasta, "fasta")
    count = {}
    # With partition file
    for rec in records:
        gaps = rec.seq.count('-') + rec.seq.count('N') + rec.seq.count('X')
        count[rec.id] = {'total_nucleotides': len(rec.seq) - gaps}
        row = [rec.id, '', '']
        g = 0
        prec = [] # gene prescence/absence annotation for tree
        for gene in partitions.keys():
            gene_seq = rec.seq[partitions[gene][0]:partitions[gene][1]]
            gaps = gene_seq.count('-') + gene_seq.count('N') + gene_seq.count('X') + gene_seq.count('*')
            length = partitions[gene][2] - gaps
            count[rec.id][gene] = length
            if length > 0:
                g += 1
        count[rec.id]['total_partitions'] = g
    return count

def write_csv(count, output_csv, partitions):
    # Write CSV gene length file
    with open(output_csv, "w") as file:
        writer = csv.writer(file)
        partition_names = list(partitions.keys()) if partitions else []
        partition_codes = {name: string.ascii_uppercase[i % 26] for i, name in enumerate(partition_names)}
        print('Writing partition representation string for each sequence:')
        for name in partition_names:
            print(f'{partition_codes[name]} = {name}')
        writer.writerow(["taxon", "partition_representation", "total_partitions", "total_nucleotides"] + [n.replace('.fasta', '') for n in partition_names])
        x = 0
        for rec in count.keys():
            part_rep = []
            x += 1
            row = [rec, count[rec]['total_partitions'], count[rec]['total_nucleotides']]
            for name in partition_names:
                length = count[rec][name]
                row.append(length)
                # Write partition_rep string
                part_rep.append(partition_codes[name]) if length > 0 else part_rep.append('-')
            row[1] = ''.join(part_rep)
            writer.writerow(row)
    print(f'{x} records written to {output_csv}')

def main():
    args = parse_args()
    if args.partitions:
        partitions = read_partitions(args.partitions)
        partitions_genes = partitions.keys()    
    gene_counts = process_fasta(args.fasta, partitions)
    write_csv(gene_counts, args.output, partitions)

if __name__ == "__main__":
    main()