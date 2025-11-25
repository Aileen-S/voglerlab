#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re



def parse_args():
    parser = argparse.ArgumentParser(description="Check fasta alignment for internal stop codons and partial codons")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-o", "--output", type=str, help="Output fasta")
    parser.add_argument("-s", "--stop", type=str, help="Save seqeunces with stop codons to specified file, instead of cleaning sequences")
    parser.add_argument("-m", "--mito", action='store_true', help="Sequence is invertebrate mitochondrial")
    parser.add_argument("-f", "--frame", type=int, help="Specify reading frame (default 1)")

    return parser.parse_args()

def find_reading_frame(fasta_file, trans_table):
    
    frame_totals = [0, 0, 0]
    stop_codons = ['TAA', 'TAG', 'TGA']
    if trans_table == 5:
        stop_codons.remove('TGA')
    records = SeqIO.parse(fasta_file, 'fasta')
    for rec in records:
        for f in range(3):
            x = 0
            codons = [rec.seq[i:i+3] for i in range(f, len(rec.seq), 3)]
            for codon in codons:
                if codon in stop_codons:
                    x += 1
            frame_totals[f] += x
    frame = frame_totals.index(min(frame_totals)) + 1
    print(f'{frame_totals[0]} stop codons in frame 1')
    print(f'{frame_totals[1]} stop codons in frame 2')
    print(f'{frame_totals[2]} stop codons in frame 3')
    print(f'Best reading frame: {frame}')
    return frame


def replace_partial_and_stop_codons(fasta_file, trans_table, frame, stop=False):
    stop_codons = ['TAA', 'TAG', 'TGA']
    frame = frame - 1
    if trans_table == 5:
        stop_codons.remove('TGA')
        locus = 'invertebrate mitochondrial'
    else:
        locus = 'standard code'
    good = []
    check = []
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    total_records = len(records)
    print(f'Found {total_records} sequences in {fasta_file}')
    print(f'Searching for internal stop codons: {locus} {(", ").join(stop_codons)}')
    recs = []
    stop_seqs = []
    x = 0
    y = 0
    for record in records:
        has_stop_codon = False
        sequence = str(record.seq).upper()
        # Split into codons according to reading frame
        codons = [sequence[i:i+3] for i in range(frame, len(sequence), 3)]
        seq = ''
        # Exclude the last codon
        for codon in codons[0:-1]:
            if codon in stop_codons:
                if stop:
                    has_stop_codon = True
                    stop_seqs.append(record)
                    break
                # else:
                #     codon = '---'
                #     x += 1
            if '-' in codon and any(c.isalpha() for c in codon):
                codon = codon.replace('-', 'N')
                y += 1
            seq = seq + codon
        if has_stop_codon == False:
            record.seq = seq
            recs.append(record)
    if stop:
        print(f'{len(stop_seqs)} seqeunces with internal stop codons')
    else:
        # print(f'{x} internal stop codons and {y} partial internal codons replaced')
        print(f'{y} partial internal codons replaced')

    return recs, stop_seqs



def main():
    args = parse_args()
    trans_table = 5 if args.mito else 1

    if args.frame:
        frame = args.frame
    else:
        frame = find_reading_frame(args.input, trans_table)

    if args.stop:
        recs, stop = replace_partial_and_stop_codons(args.input, trans_table, frame, stop=True)
    else:
        recs, stop = replace_partial_and_stop_codons(args.input, trans_table, frame)

    # Save sequences without internal stop codons
    with open(args.output, 'w') as file:
        SeqIO.write(recs, file, 'fasta')
    print(f'{len(recs)} sequences saved to {args.output}')

    if len(stop) > 0:  
        with open(args.stop, 'w') as file:
            SeqIO.write(stop, file, 'fasta')
        print(f'{len(stop)} sequences with internal stop codons saved to {args.stop}')
    
if __name__ == "__main__":
    main()
