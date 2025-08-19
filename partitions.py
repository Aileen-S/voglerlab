import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Create raxml-compatible partition files from catfasta2phyml.pl output (https://github.com/nylander/catfasta2phyml)")
    parser.add_argument("-i", "--input", type=str, help="Catfasta partition file")
    parser.add_argument('-p', '--prefix', type=str, help="Output file prefix")
    parser.add_argument('-q', '--iqtree', action='store_true', help="Output IQTREE format (standard is RAxML)")
    parser.add_argument('-t', '--sequence_type', choices=['aa', 'nt', 'ry'], help="Specify sequence type: protein, nucleotide or RY-coded")
    # 12S, 16S, 18S and 28S will be treated as non coding genes by default. Any aditional non-coding genes, or differnet naming format, can be added using -n
    parser.add_argument('-n', '--non-coding', type=str, help="Comma delimited list of non-coding genes")
    return parser.parse_args()


def process_partition_file(input_file, sequence_type):
    with open(input_file) as infile:
        partitions = {}
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            line = line.split(' = ')
            partition_name = os.path.basename(line[0])
            coordinates = line[1]
            partitions[partition_name] = coordinates
    return partitions


def write_partitions(partitions, sequence_type, prefix):
    data = data_map[sequence_type]
    # Gene partitions
    with open(f'{prefix}_partitions_{sequence_type}_gene.txt', 'w') as outfile:
        for partition_name, coordinates in partitions.items():
            if not any(g in partition_name for g in non_coding_genes):
                outfile.write(f'{data}, {partition_name} = {coordinates}\n')
            else:
                outfile.write(f'DNA, {partition_name} = {coordinates}\n')

    if sequence_type != 'aa':
        with open(f'{prefix}_partitions_{sequence_type}_codon.txt', 'w') as outfile:
            for partition_name, coordinates in partitions.items():
                if not any(g in partition_name for g in non_coding_genes):
                    coords = coordinates.split('-')
                    start = int(coords[0])
                    end = int(coords[1])
                    for i in range(3):
                        outfile.write(f'{data}, {partition_name}.{i + 1} = {start + i}-{end}\\3\n')
                else:
                    outfile.write(f'DNA, {partition_name} = {coordinates}\n')


non_coding_genes = ['12S', '16S', '18S', '28S']

data_map = {'aa': 'PROT',
        'nt': 'DNA',
        'bin': 'RY'}


def main():
    args = parse_args()
    if args.iqtree:
        data['aa'] = 'AA'
    if args.non_coding:
        non_coding_genes.extend(args.non_coding.split(','))

    partitions = process_partition_file(args.input, args.sequence_type)

    write_partitions(partitions, args.sequence_type, args.prefix)


if __name__ == "__main__":
    main()