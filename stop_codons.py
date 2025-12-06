#!/usr/bin/env python3

from Bio import SeqIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Check fasta alignment for internal stop codons")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--check", type=str, help="Output fasta for sequences with internal stop codons")
    parser.add_argument("-g", "--good", type=str, help="Output fasta for other sequences")
    parser.add_argument("-d", "--data", choices=['nt', 'aa'], default='nt', type=str, help="Nucleotide or amino acid")
    parser.add_argument("-m", "--mito", action='store_true', help="Sequence is invertebrate mitochondrial")
    parser.add_argument("-f", "--frame", type=int, help="Specify reading frame (default 1)")


    return parser.parse_args()


def find_internal_stop_codons(records, data, locus, frame):
    good = []
    check = []
    if data == 'nt':
        stop_codons = ['TAA', 'TAG'] if locus == 'mito' else ['TAA', 'TAG', 'TGA']
        frame = frame - 1
    for rec in records:
        found_stop = False
        seq = str(rec.seq).rstrip('-')
        if data == 'nt':
            codons = [seq[i:i+3] for i in range(frame, len(seq), 3)]
            if len(codons[-1].replace('N', '').replace('-', '')) < 3:
                codons = codons[:-1]
            if any(codon in stop_codons for codon in codons[:-1]):
                check.append(rec)
                found_stop = True
            if found_stop == False:
                good.append(rec)
        elif data == 'aa':
            check.append(rec) if '*' in seq[:-1] else good.append(rec)
    # print(f'{len(good)} sequences without internal stop codons')
    return check, good


def main():
    args = parse_args()
    locus = 'mito' if args.mito else None
    frame = args.frame if args.frame else 1
    records = list(SeqIO.parse(args.input, 'fasta'))
    check, good = find_internal_stop_codons(records=records, data=args.data, locus=locus,frame=frame)

    # Save sequences without internal stop codons
    with open(args.good, 'w') as file:
        SeqIO.write(good, file, 'fasta')
        print(f'{len(good)} sequences without internal stop codons saved to {args.good}')
    # Save sequences containing internal stop codons
    with open(args.check, 'w') as file:
        SeqIO.write(check, file, 'fasta')
        print(f'{len(check)} sequences with internal stop codons written to {args.check}')


if __name__ == "__main__":
    main()
