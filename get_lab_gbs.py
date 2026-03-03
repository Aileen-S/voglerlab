import glob
import argparse

parser = argparse.ArgumentParser(description="Get lab mitogenome genbanks from ID list")
parser.add_argument("-i", "--input", type=str, help="Input ID list file")
parser.add_argument("-o", "--output", type=str, help="Name output file: default is <inputfile>.gb")
parser.add_argument('-v', '--version', type=str, help='Database version, eg: gbmaster_2024-11-15')

args = parser.parse_args()

file = open(args.input)
if args.output:
    output = open(args.output, "w")
else:
    output = open(f'{args.input}.gb', 'w')

x = 0
lines = file.readlines()
print(f'Searching /mbl/share/workspaces/groups/voglerlab/MMGdatabase/{args.version}/ for IDs in {args.input}')
for line in lines:
    line = line.strip()
    try:
        name = f'/mbl/share/workspaces/groups/voglerlab/MMGdatabase/{args.version}/{line}.gb*'
        name = glob.glob(name)
        x += 1
        record = open(name[0])
        record = record.read()
        output.write(record)
    except FileNotFoundError:
        print(f'No record found for {line}')

print(f'{x} records saved to {args.output}' if args.output else
      f'{x} records saved to {args.input}.gb')
