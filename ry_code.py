from Bio import SeqIO
import argparse
import re

def main():

    parser = argparse.ArgumentParser(description="RY coding")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-o", "--output", type=str, help="Output fasta")
    parser.add_argument("-n", "--numeric", action='store_true', help="0/1 output instead of RY")

    args = parser.parse_args()

    if args.numeric:
        with open(args.output, "w") as output:
            records = SeqIO.parse(args.input, "fasta")
            for rec in records:
                seq = str(rec.seq.upper())
                seq = re.sub(r'[AGRCTYSWKMDHBVN]', replace_numeric, seq)
                output.write(f'>{rec.id}\n{seq}\n')

    else:
        with open(args.output, "w") as output:
            records = SeqIO.parse(args.input, "fasta")
            for rec in records:
                seq = str(rec.seq.upper())
                seq = re.sub(r'[AGRCTYSWKMDHBVN]', replace_ry, seq)
                output.write(f'>{rec.id}\n{seq}\n')


# 1/0 code
def replace_numeric(match):
    char = match.group(0)
    if char in 'AGR':
        return '1'
    elif char in 'CTY':
        return '0'
    elif char in 'SWKMDHBVN':
        return '-'

# RY code
def replace_ry(match):
    char = match.group(0)
    if char in 'AGR':
        return 'R'
    elif char in 'CTY':
        return 'Y'
    elif char in 'SWKMDHBVN':
        return '-'


if __name__ == "__main__":
    main()