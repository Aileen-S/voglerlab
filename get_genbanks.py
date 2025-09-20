#!/usr/bin/env python3

import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
import time
import re
from Bio.Seq import Seq, UndefinedSequenceError

# Function definitions

def get_feat_name(feat):
    featnames = []
    nametags = ['gene', 'product', 'label', 'standard_name']  # Search these four keys for gene name
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                featnames.append(featname)
    return featnames


def set_feat_name(feat, name):
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys() :
                feat.qualifiers[t][0] = name
    return feat


# Search GenBank and retrive list of GBIDs.
# Args:
#     term (str): Search term for nucleotides.
#     chunk (optional): Number of results per request (default: 10000).
#     retries (optional): Number of retries on HTTP errors (default: 3).
#     delay (optional): Delay in seconds between retries (default: 30).
# Yields:
#     list: List of GenBank IDs from the search.
def get_gbids(query, chunk=10000, retries=10, delay=30):
    gbids = set()
    # Define search term 
    terms = [f"txid{txid}" for txid in query] if isinstance(query, list) else [query]
    if isinstance(query, list):
        print('Retriving GenBank record IDs for input taxon ID list')
    else:
        print(f'Retriving GenBank record IDs for {query}')
    for term in terms:
        for attempt in range(retries):
            try:
                # Get record count for search term
                searchhand = Entrez.esearch(db="nucleotide", term=term, retmax=0)
                searchrec = Entrez.read(searchhand)
                count = int(searchrec["Count"])
                # Get GBIDs
                for start in range(0, count, chunk):
                    searchhand = Entrez.esearch(db="nucleotide", term=term, retstart=start, retmax=chunk)
                    searchrec = Entrez.read(searchhand)
                    gbids.update(searchrec['IdList'])
            # If HTTP error, pause and try again
            except Entrez.HTTPError:
                print(f"HTTP Error: retrying in {delay} seconds")
                time.sleep(delay)
    if gbids == []:
        print(f"Failed to retrieve records")
        return None
    print(f'Found {len(gbids)} GenBank IDs')
    gbids = list(gbids)
    return gbids

# Accessions to GBIDs

def accs_to_gbids(ids, chunk_size=500, retries=10, delay=30):
    total = len(ids)
    all_ids = []
    for i in range(0, total, chunk_size):
        chunk = ids[i:i + chunk_size]
        for attempt in range(retries):
            try:
                term = " OR ".join(chunk)
                handle = Entrez.esearch(db="nucleotide", term=term)
                record = Entrez.read(handle)
                all_ids.extend(record["IdList"])
                break  # successful attempt
            except Entrez.HTTPError:
                print(f"HTTP error for chunk {i}-{i+chunk_size}; retrying in {delay}s")
                time.sleep(delay)
        else:
            print(f"Failed to retrieve GBIDs {i}-{i+chunk_size}")

    return all_ids



# Search GenBank with ID list
def search_genbank(ids, chunk_size=500, retries=10, delay=30, save=False, output="records.gb", exclude=[]):
    total = len(ids)
    processed = 0
    print(f'Downloading records for {total} GenBank IDs')

    if save:
        outfile = open(output, "w")

    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]

        for attempt in range(retries):
            try:
                with Entrez.efetch(db="nucleotide", id=','.join(chunk), rettype="gb", retmode="text") as handle:
                    results = SeqIO.parse(handle, "gb")
                    for record in results:
                        try:
                            processed += 1
                            if processed % 500 == 0:
                                print(f"Downloaded {processed} of {total} records")
                            if processed == total:
                                print(f"Downloaded {processed} of {total} records")
                            yield record

                            if save:
                                SeqIO.write(record, outfile, "genbank")
                        except HTTPException:
                            print(f"Incompete read error with record {record.name}")
                    break # Stop addition retries if successful
            except Entrez.HTTPError:
                print("HTTP error fetching records; retrying in 20 seconds")
                time.sleep(delay)
        else:
            print(f"Failed to retrieve records for chunk {i}-{i+chunk_size}")



def genbank_metadata(rec):
    # NCBI taxon ID
    db_xref = rec.features[0].qualifiers.get("db_xref", [])
    txid = ""
    for ref in db_xref:
        if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
            txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value

    # Taxonomy
    # Replace the following characters: > < . ( ) ; : ' ,
    spec = re.sub(r"[><.();:'\"]", "", rec.annotations["organism"]).replace(",", "")
    spec_parts = [part for part in spec.split(" ") if not re.search(r'\d', part) and not part.isupper()]
    spec = " ".join(spec_parts)
    #specfasta = spec.replace(" ", "_")

    taxonomy = ['', '', '', '', '', '']
    for tax in rec.annotations["taxonomy"]:
        if tax in suborders: taxonomy[0] = tax
        if tax.endswith('formia'): taxonomy[1] = tax
        if tax.endswith('oidea'): taxonomy[2] = tax
        if tax.endswith('idae'): taxonomy[3] = tax
        if tax.endswith('inae'): taxonomy[4] = tax
        if tax.endswith('ini'): taxonomy[5] = tax
    #taxonomy.append(spec.split(' ')[0])
    #fastatax = f"{txid}_{taxonomy[2]}_{taxonomy[3]}_{taxonomy[4]}_{specfasta}"

    # Location
    if "country" in rec.features[0].qualifiers:
        location = rec.features[0].qualifiers["country"][0]
        if ":" in location:
            country, region = location.split(":", 1)
        else:
            country = location
    else:
        country = ""
        region = ""
    if "lat_lon" in rec.features[0].qualifiers:
        latlon = rec.features[0].qualifiers["lat_lon"][0]
        try:
            ll_list = latlon.split(" ")
            if ll_list[1] == "N":
                lat = ll_list[0]
            else:
                lat = "-" + ll_list[0]
            if ll_list[3] == "E":
                long = ll_list[2]
            else:
                long = "-" + ll_list[2]
        except IndexError:
            lat = ""
            long = ""
    else:
        latlon = ""
        lat = ""
        long = ""

    # References
    refs = []
    if "references" in rec.annotations:
        first = rec.annotations['references'][0]
        refs.append(first.authors)
        refs.append(first.title)
        refs.append(first.journal)
    output = {"gbid": rec.name,
              "txid": txid,
              "description": rec.description,
              "spec_id": rec.annotations["organism"],
              "spec": spec,
              "date": rec.annotations["date"],
              "taxonomy": taxonomy,
              "country": country,
              "lat": lat,
              "long": long,
              "refs": refs,
              "row": [txid, rec.name, '', '', ''] + taxonomy + [spec, country, lat, long] + refs}
    return output

# Argument parser
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta.")

# Input options
# Choose either taxon or file (file of IDs with one per line)
# If using file, specify accessions or taxon IDs
# Add your email address for NCBI
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-f', '--file', type=str, help="Input file with accession or taxon ID list")
parser.add_argument('-r', '--ref', choices=['txid', 'gbid'], help="If using --file option, specify accessions or taxon IDs.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")

# Filtering options
parser.add_argument("-x", "--exclude", type=str, help="Input file with list of accession to skip (will not search or save these records)")
parser.add_argument('-m', '--mito', action='store_true', help='Save only mitochondrial protein-coding genes')
parser.add_argument("-g", "--gene", type=str, help="Save specified gene only (use gene name format from genes dist below)")
parser.add_argument("-l", "--longest", action='store_true', help="Save only longest sequence for each NCBI taxid")

# Output options
# Optional output of genbank format file as well as fasta
parser.add_argument("-s", "--save", type=str, help="Output genbank file with initial search results")

args = parser.parse_args()         # Process input args from command line


genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "RRNS", "SSU", "RRN12", "S-RRNA", "12S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "RRNL", "LSU", "RRN16", "L-RRNA", "16S LARGE SUBUNIT RIBOSOMAL RNA", "LARGE SUBUNIT RIBOSOMAL RNA"],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA", "SMALL SUBUNIT RIBOSOMAL RNA"],
         "28S": ["28S", "28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA", "LARGE SUBUNIT RIBOSOMAL RNA", "28S LARGE SUBUNIT"],

         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],

         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I',   'CYTOCHROME C OXIDASE SUBUNIT I',   'COXI',   'CO1', 'COI',   'CYTOCHROME COXIDASE SUBUNIT I',   'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE I', 'CYTOCHROME OXYDASE SUBUNIT 1', 'CYTOCHROME OXIDASE C SUBUNIT I', 'COX 1', 'COX1', 'CYTCHROME OXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE 1'],
         "COX2": ['CYTOCHROME C OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE SUBUNIT II',  'CYTOCHROME C OXIDASE SUBUNIT II',  'COXII',  'CO2', 'COII',  'CYTOCHROME COXIDASE SUBUNIT II',  'CYTOCHROME OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE II', 'CYTOCHROME C OXIDASE II', 'CYTOCHROME OXYDASE C SUBUNIT 2', 'CYTOCHROME OXIDASE C SUBUNIT 2', 'COX2', 'CYTOCHROME OXIDASE (CO) II'],
         "COX3": ['CYTOCHROME C OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE SUBUNIT III', 'CYTOCHROME C OXIDASE SUBUNIT III', 'COXIII', 'CO3', 'COIII', 'CYTOCHROME COXIDASE SUBUNIT III', 'CYTOCHROME OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE III', 'CYTOCHROME OXIDASE C SUBUNIT 3',  'COX3', 'CYTOMCHROME C OXIDASE SUBUNIT 1'],
         
         "CYTB": ['CYTOCHROME B', 'CYB', 'COB', 'COB / CYTB', 'CYTB', "COB/CYTB"],
         "ND1": ['NAD1', 'NSD1', 'NADH1', 'NADH DEHYDROGENASE SUBUNIT I', 'NADH DEHYDROGENASE SUBUNIT 1', 'NADH DESHYDROGENASE SUBUNIT 1', 'NAD1-0', 'ND1'],
         "ND2": ['NAD2', 'NSD2', 'NADH2', 'NADH DEHYDROGENASE SUBUNIT II', 'NADH DEHYDROGENASE SUBUNIT 2', 'NADH DESHYDROGENASE SUBUNIT 2', 'NAD2-0', 'ND2'],
         "ND3": ['NAD3', 'NSD3', 'NADH3', 'NADH DEHYDROGENASE SUBUNIT III', 'NADH DEHYDROGENASE SUBUNIT 3', 'NADH DESHYDROGENASE SUBUNIT 3', 'NAD3-0', 'ND3'],
         "ND4": ['NAD4', 'NSD4', 'NADH4', 'NADH DEHYDROGENASE SUBUNIT IV', 'NADH DEHYDROGENASE SUBUNIT 4', 'NADH DESHYDROGENASE SUBUNIT 4', 'NAD4-0', 'ND4'],
         "ND4L": ['NAD4L', 'NSD4L', 'NADH4L', 'NADH DEHYDROGENASE SUBUNIT IVL', 'NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DESHYDROGENASE SUBUNIT 4L', 'NAD4L-0', 'ND4L'],
         "ND5": ['NAD5', 'NSD5', 'NADH5', 'NADH DEHYDROGENASE SUBUNIT V', 'NADH DEHYDROGENASE SUBUNIT 5', 'NADH DESHYDROGENASE SUBUNIT 5', 'NAD5-0', 'ND5'],
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6'],
         
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "CAD": ["CAD", "CAD FRAGMENT 1", "CARBAMOYLPHOSPHATE SYNTHETASE"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA", "EF-1A"],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "RNApol": ["RNA POL II", "RNA POL2", "RNA POLYMERASE II LARGE SUBUNIT"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"]}

mito = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
nuc = ['AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']
rna = ['12S', '16S', '18S', '28S']
cds = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']


if args.email:
    Entrez.email = args.email

if args.ref == 'gbid':
    gbids = []
    file = open(args.file)
    lines = file.readlines()
    for line in lines:
        acc = line.strip()
        gbids.append(acc)
    print(f'{len(gbids)} IDs found in {args.file}')



else:
    if args.ref  == 'txid':
        taxids = []
        file = open(args.file)
        lines = file.readlines()
        for line in lines:
            taxid = line.strip()
            taxids.append(taxid)
        print(f'{len(taxids)} IDs found in {args.file}')
        gbids = get_gbids(taxids)
        #print(f"Found {len(gbids)} records for input ID list")
    if args.taxon:
        basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"

        # Retrieve GBIDs for search term
        gbids = get_gbids(basesearch)
        #print(f"Found {len(gbids)} records for {args.taxon}")

# Remove 'exlude' GBIDs from search list
exc_accs = []
exc = ''
if args.exclude:
    exclude = open(args.exclude)
    lines = exclude.readlines()
    for line in lines:
        acc = line.strip()
        exc_accs.append(acc)
    print(f'{len(exc_accs)} IDs found in {args.exclude} with be excluded from search results')

    # Convert to from accessions to GBIDs
    exc = accs_to_gbids(exc_accs)
    filtered = [id for id in gbids if id not in exc]
    gbids = filtered
    print(f'{len(gbids)} IDs remaining after exclusions removed')

# Search through GBIDs
meta = {}
seqs = {}
nohits = []
other_type = set()
misc_feature = set()
unrec_genes = {}
unrec_species = []
x = 0  # Count taxids

if args.save:
    results = search_genbank(gbids, save=True, output=args.save, exclude=exc)
else:
    results = search_genbank(gbids, exclude=exc)
for rec in results:
    if args.taxon:
        if args.taxon not in rec.annotations["taxonomy"]:
            unrec_species.append(rec.name)
            continue
    output = genbank_metadata(rec)
    meta[rec.name] = output
    g = 0
    for feature in rec.features:
        type = feature.type
        if type not in ('CDS', 'rRNA'):
            if type == 'misc_feature':
                if 'note' in feature.qualifiers:
                    misc_feature.add(tuple(feature.qualifiers['note']))
            else:
                other_type.add(type)
            continue
        names = get_feat_name(feature)                       # Find gene name
        stdname = ""
        for k, v in genes.items():
            for name in names:
                if name in v:
                    stdname = k
                    g += 1
        if stdname == '':
            for name in names:
                if name in unrec_genes:
                    unrec_genes[name].append(rec.name)
                else:
                    unrec_genes[name] = [rec.name]
            continue

        if args.mito:
            if stdname not in mito:
                continue
        if args.gene:
            if stdname != args.gene:
                continue
        frame = ''
        if stdname in cds:
            if 'codon_start' in feature.qualifiers:
                frame = feature.qualifiers["codon_start"][0]
            else:
                print(f"Reading frame missing from record {rec.name}, {stdname}.")
        # try:
        seq = feature.extract(rec.seq)
        # except UndefinedSequenceError:
        #     print(f"Error extracting sequence for record '{rec.name}', feature '{names}')")
        #     continue
        gene_output = {"gbid": rec.name,
                       "gene": stdname,
                       "length": len(seq),
                       "seq": seq,
                       "frame": frame}
        if stdname in seqs:                              # If taxon ID in dict
            if output['txid'] in seqs[stdname]:                 # If gene in dict for that taxon ID
                seqs[stdname][output['txid']].append(gene_output)    # Add gene info list to dict
                x += 1
            else:
                seqs[stdname][output['txid']] = [gene_output]      # Otherwise add to dict with new key
                x += 1
        else:
            seqs[stdname] = {output['txid']: [gene_output]}      # Otherwise add to dict with new key
            x += 1
    if g == 0:
        nohits.append(rec.name)



if args.longest:
    print("Saving longest sequence for each gene for each NCBI taxonomy ID")
    # Set record length as 0, iterate through records and replace whenever another sequence is longer.
    def findmax(x):
        max = x[0]["length"]
        maxrec = x[0]
        for record in x:
            if record["length"] > max:
                max = record["length"]
                maxrec = record
        return maxrec


    # Dict for longest sequences, key is gene stdname, value is list of records
    saved_gbids = []
    saved_recs = {}
    for stdname, tax in seqs.items():
        for tax, records in tax.items():
            chosen = findmax(records)
            saved_gbids.append(chosen['gbid'])
            if gene in saved_recs:
                if tax in saved_recs[gene]:
                    saved_recs[gene][tax].append(chosen)
                else:
                    saved_recs[gene][tax] = [chosen]
            else:
                saved_recs[gene] = {tax: [chosen]}

else:
    saved_recs = seqs
    saved_gbids = gbids

# for gene, records in saved_recs.items():
#     print(f"{len(records)} records found for {gene}")

# Write fastas
noframe = {}
for gene, tax in saved_recs.items():
    file = open(f"{gene}.fasta", "w")
    x = 0
    y = 0
    for tax, records in tax.items():
        for rec in records:
            try:
                seq = rec['seq']
                if gene in rna:
                    file.write(f">{rec['gbid']}\n{seq}\n")
                    x += 1
                else:
                    file.write(f">{rec['gbid']};frame={rec['frame']}\n{seq}\n")
                    x += 1
                    if rec['frame'] == '': 
                        if gene in noframe: noframe[gene].append(rec['gbid'])
                        else: noframe[gene] = [rec['gbid']]
            except UndefinedSequenceError:
                print(f"Error extracting sequence for record '{rec['gbid']}', {gene})")
    print(f'{x} records written to {gene}.fasta')


# Write CSV metadata file
added = []
with open("metadata.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow(
        ["ncbi_taxid", "genbank_accession", "bold_id", "bold_bin", "lab_id", "suborder", "infraorder", "superfamily", "family", 
        "subfamily", "tribe", "species", "country", "latitude", "longitude", "ref_authoer", "ref_title", "ref_journal"])
    for gbid, rec in meta.items():
        if gbid in gbids:
            if gbid not in added:
                added.append(rec['gbid'])
                row = rec['row']
            writer.writerow(row)
    print("Metadata saved to metadata.csv")


if args.file:
    if len(nohits) > 0:
        print(f'\nNo requested genes found in the following records: {nohits}')

if unrec_genes != {}:
    print("\nUnrecognised genes printed below, and to 'other_genes.csv")
    with open("other_genes", "w") as file:
        writer = csv.writer(file)
        writer.writerow(['gene', 'count', 'records'])
        for gene, recs in unrec_genes.items():
            print(f'{gene}: {len(recs)} record' if len(recs) == 1 else f'{gene}: {len(recs)} records')
            #recs = ', '.join(recs)
            writer.writerow([gene, len(recs), recs])

print('\nMisc Features:')
print(misc_feature)
print("\nOther Feature Types:")
print(other_type)
print("\nUnrecognised Species:")
print(unrec_species)
print("\nMissing Reading Frames:")
for gene, gbids in noframe.items():
    print(f"{gene}: {gbids}")