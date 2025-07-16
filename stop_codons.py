from Bio import SeqIO
import argparse, argcomplete


def parse_args():
    parser = argparse.ArgumentParser(description="Check fasta alignment for internal stop codons")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-c", "--check", type=str, help="Output fasta for sequences with internal stop codons")
    parser.add_argument("-g", "--good", type=str, help="Output fasta for other sequences")
    parser.add_argument("-m", "--mito", action='store_true', help="Sequence is invertebrate mitochondrial")
    parser.add_argument("-f", "--frame", type=int, help="Specify reading frame (default 1)")


    argcomplete.autocomplete(parser)
    return parser.parse_args()

def find_internal_stop_codons(fasta_file, locus, frame):
    stop_codons = ['TAA', 'TAG', 'TGA']
    if locus == 'invertebrate mitochondrial':
        stop_codons.remove('TGA')
    else:
        locus = 'standard code'
    start = (frame - 1) if frame else 0
    good = []
    check = []
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    total_records = len(records)
    print(f'Found {total_records} sequences in {fasta_file}')
    print(f'Searching for internal stop codons: {locus} {(", ").join(stop_codons)}')
    for record in records:
        has_stop_codon = False
        sequence = str(record.seq)
        # Split into codons according to reading frame
        codons = [sequence[i:i+3] for i in range(start, len(sequence), 3)]
        # Exclude the last codon
        for codon in codons[0:-1]:
            if codon in stop_codons:
                check.append(record)
                has_stop_codon = True
                break
        if has_stop_codon == False:
            good.append(record)

    return check, good



def main():
    args = parse_args()
    check, good = find_internal_stop_codons(args.input, args.mito and 'invertebrate mitochondrial' or None, args.frame)

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
