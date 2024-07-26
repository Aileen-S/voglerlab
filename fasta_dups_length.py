from Bio import SeqIO
import argparse, argcomplete
import sys

parser = argparse.ArgumentParser(description="Remove duplicate fasta IDs and/or pad sequences to match longest.")
parser.add_argument("-i", "--input", type=str, help="Input fasta")
parser.add_argument('-o', '--output', type=str, help="Output fasta")
parser.add_argument('-d', '--dups', action='store_true', help='Remove duplicate fasta IDs')
parser.add_argument('-l', '--length', action='store_true', help='Pad sequences to match longest')

argcomplete.autocomplete(parser)
args = parser.parse_args()

if not args.dups and not args.length:
    print("Error: Either the -d flag (remove duplicates) or the -l flag (pad sequences) must be used.")
    sys.exit(1)

seen = []
records = []
dups = []

if args.dups:
# Remove duplicate IDs
    for record in SeqIO.parse(args.input, "fasta"):  
        if str(record.name) not in seen:
            seen.append(str(record.name))
            records.append(record)
        else:
            dups.append(record.name)

    if len(dups) > 0:
        print(f'{len(dups)} duplicate records removed from {args.input}: {dups}')
else:
    for record in SeqIO.parse(args.input, "fasta"):  
        records.append(record)

if args.length:
    # Add gaps to pad all sequences to length of longest
    max_len = max([len(rec.seq) for rec in records])
    print(f'Longest sequence length = {max_len}\n'
        'Adding gaps to pad shorter sequences')
    for rec in records:
        padding = '-'*(max_len - len(rec.seq))
        rec.seq += padding

SeqIO.write(records, args.output, 'fasta')
print(f'Output written to {args.output}')