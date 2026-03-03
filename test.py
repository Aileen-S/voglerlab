#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re


filter = []
with open('test.ids') as infile:
    lines = infile.readlines()
    for line in lines:
        if line != '\n':
            filter.append(line.strip())


with open('/home/aileen/onedrive/treebuilding/251119carabidae/lab/mmg.gb') as file:
    records = list(SeqIO.parse(file, "genbank"))
print(f'{len(records)} records found')

found = []
for rec in records:
    print(rec)
    if rec.name in filter:
        found.append(rec)


with open('test.gb', 'w') as outfile:
    SeqIO.write(found, outfile, "genbank")
