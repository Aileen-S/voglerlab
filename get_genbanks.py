import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
import time
import re

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



def get_gbids(term, chunk=10000, retries=10, delay=30):
    # Search GenBank and retrive list of GBIDs.
    # Args:
    #     term (str): Search term for nucleotides.
    #     chunk (optional): Number of results per request (default: 10000).
    #     retries (optional): Number of retries on HTTP errors (default: 3).
    #     delay (optional): Delay in seconds between retries (default: 30).
    # Yields:
    #     list: List of GenBank IDs from the search.
    gbids = []
    print(f"Searching Genbank for {term}")
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
                gbids.extend(searchrec['IdList'])
            return gbids
        # If HTTP error, pause and try again
        except Entrez.HTTPError:
            print(f"HTTP Error: retrying in {delay} seconds")
            time.sleep(delay)

    print(f"Failed to retrieve records after {retries} attempts.")
    return None


# Search GenBank with ID list
def search_genbank(ids, chunk_size=500, retries=10, delay=30):
    count = ids.split(',')
    total = len(count)
    processed = 0
    print(f'Downloading {total} records')
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]

        for attempt in range(retries):
            try:
                handle = Entrez.efetch(db="nucleotide", id=','.join(chunk), rettype="gb", retmode="text")
                results = SeqIO.parse(handle, "gb")
                for record in results:
                    processed += 1
                    if processed % 500 == 0:
                        print(f"Downloaded {processed} of {total} records")
                    if processed == total:
                        print(f"Downloaded {processed} of {total} records")
                    yield record
                break
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
        ll_list = latlon.split(" ")
        if ll_list[1] == "N":
            lat = ll_list[0]
        else:
            lat = "-" + ll_list[0]
        if ll_list[3] == "E":
            long = ll_list[2]
        else:
            long = "-" + ll_list[2]
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
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-f', '--file', type=str, help="Input file with accession or taxon ID list")
parser.add_argument('-r', '--ref', choices=['txid', 'gbid'], help="If using --file option, specify accessions or taxon IDs.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")
parser.add_argument('-m', '--mito', action='store_true', help='Save only mitochondrial protein-coding genes')
parser.add_argument('-c', '--coi', action='store_true', help='Save only COX1')
parser.add_argument("-g", "--gene", type=str, help="Save specified gene only")

args = parser.parse_args()         # Process input args from command line
#args = argparse.Namespace(taxon='Amphizoidae', mpc=True, email='aileen.scott@nhm.ac.uk', nuclear=False) # This is how I step through the script interactively

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "SSU"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "LSU"],
         "ATP6": ['ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6', 'ATP6'],
         "ATP8": ['ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8', 'ATP8'],
         "COX1": ['CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I', 'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1', 'COX1', 'CYTOCHROME OXIDASE I', 'CYTOCHROME OXIDASE C SUBUNIT I'],
         "COX2": ['CYTOCHROME C OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE SUBUNIT II', 'CYTOCHROME C OXIDASE SUBUNIT II', 'COXII', 'CO2', 'COII', 'CYTOCHROME COXIDASE SUBUNIT II', 'CYTOCHROME OXIDASE SUBUNIT 2', 'CYTOCHROME OXYDASE SUBUNIT 2', 'COX2'],
         "COX3": ['CYTOCHROME C OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE SUBUNIT III', 'CYTOCHROME C OXIDASE SUBUNIT III', 'COXII', 'CO3', 'COIII', 'CYTOCHROME COXIDASE SUBUNIT III', 'CYTOCHROME OXIDASE SUBUNIT 3', 'CYTOCHROME OXYDASE SUBUNIT 3', 'COX3'],
         "CYTB": ['CYTOCHROME B', 'CYB', 'COB', 'COB / CYTB', 'CYTB', "COB/CYTB"],
         "ND1": ['NAD1', 'NSD1', 'NADH1', 'NADH DEHYDROGENASE SUBUNIT I', 'NADH DEHYDROGENASE SUBUNIT 1', 'NADH DESHYDROGENASE SUBUNIT 1', 'NAD1-0', 'ND1'],
         "ND2": ['NAD2', 'NSD2', 'NADH2', 'NADH DEHYDROGENASE SUBUNIT II', 'NADH DEHYDROGENASE SUBUNIT 2', 'NADH DESHYDROGENASE SUBUNIT 2', 'NAD2-0', 'ND2'],
         "ND3": ['NAD3', 'NSD3', 'NADH3', 'NADH DEHYDROGENASE SUBUNIT III', 'NADH DEHYDROGENASE SUBUNIT 3', 'NADH DESHYDROGENASE SUBUNIT 3', 'NAD3-0', 'ND3'],
         "ND4": ['NAD4', 'NSD4', 'NADH4', 'NADH DEHYDROGENASE SUBUNIT IV', 'NADH DEHYDROGENASE SUBUNIT 4', 'NADH DESHYDROGENASE SUBUNIT 4', 'NAD4-0', 'ND4'],
         "ND4L": ['NAD4L', 'NSD4L', 'NADH4L', 'NADH DEHYDROGENASE SUBUNIT IVL', 'NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DESHYDROGENASE SUBUNIT 4L', 'NAD4L-0', 'ND4L'],
         "ND5": ['NAD5', 'NSD5', 'NADH5', 'NADH DEHYDROGENASE SUBUNIT V', 'NADH DEHYDROGENASE SUBUNIT 5', 'NADH DESHYDROGENASE SUBUNIT 5', 'NAD5-0', 'ND5'],
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6'],
         "18S": ["18S", "18S RIBOSOMAL RNA", "18S RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA"],
         "28S": ["28S RIBOSOMAL RNA", "28S RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA"],
         "AK": ["AK", "ARGININE KINASE", "ARGK", "ARGKIN", "ARGS", "ARK"],
         "CAD": ["CAD", "CAD FRAGMENT 1", "CARBAMOYLPHOSPHATE SYNTHETASE"],
         "EF1A": ["EF1-ALPHA", "EF1A", "ELONGATION FACTOR 1 ALPHA", "ELONGATION FACTOR 1-ALPHA"],
         "H3": ["H3", "HISTONE 3", "HISTONE H3", "HIS3"],
         "RNApol": ["RNA POL II", "RNA POL2", "RNA POLYMERASE II LARGE SUBUNIT"],
         "Wg": ["WG", "WINGLESS", "WNG", "WNT", "WNT1", "WNT-4"]}

mito = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
nuc = ['AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']
rna = ['12S', '16S', '18S', '28S']
cds = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

unrec_genes = []
unrec_species = []
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

        # Generate search term to get all sequences in the search taxonomy
    if args.taxon:
        basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"

        # Retrieve GBIDs for search term
        gbids = get_gbids(term=basesearch)
        print(f"Found {len(gbids)} records for {args.taxon}")

# Search through GBIDs
meta = {}
seqs = {}
nohits = []
other_type = set()
misc_feature = set()
x = 0  # Count taxids

gbid_str = ",".join(gbids)   
results = search_genbank(gbid_str)
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
            unrec_genes.append(name)
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
        seq = feature.extract(rec.seq)
        gene_output = {"gbid": rec.name,
                       "gene": stdname,
                       "length": len(seq),
                       "seq": seq,
                       "frame": frame}
        if output['txid'] in seqs:                              # If taxon ID in dict
            if stdname in seqs[output['txid']]:                 # If gene in dict for that taxon ID
                seqs[output['txid']][stdname].append(gene_output)    # Add gene info list to dict
                x += 1
            else:
                seqs[output['txid']][stdname] = [gene_output]      # Otherwise add to dict with new key
                x += 1
        else:
            seqs[output['txid']] = {stdname: [gene_output]}      # Otherwise add to dict with new key
            x += 1
    if g == 0:
        nohits.append(rec.name)

# for rec in species.values():
#     for r in rec.values():
#         print(r)

print(f"\n{x} sequences found for requested genes\n"
      f"Saving longest sequence for each gene for each NCBI taxonomy ID")

#print("\nUnrecognised Species")
#print(unrec_species)


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
gbids = []
longest = {}
for tax, stdname in seqs.items():
    for gene, records in stdname.items():
        chosen = findmax(records)
        gbids.append(chosen['gbid'])
        if gene in longest:
            longest[gene].append(chosen)
        else:
            longest[gene] = [chosen]

            
# Write fastas
for gene, records in longest.items():
    file = open(f"{gene}.fasta", "w")
    x = 0
    y = 0
    for rec in records:
        if gene in rna:
            file.write(f">{rec['gbid']}\n{rec['seq']}\n")
            x += 1
        else:
            if rec['frame'] == '':
                if y == 0:
                    rf = open(f"{gene}.rf", "w")
                rf.write(f">{rec['gbid']}\n{rec['seq']}\n")
                y += 1
            else:
                file.write(f">{rec['gbid']};frame={rec['frame']}\n{rec['seq']}\n")
                x += 1

    print(f'{x} records written to {gene}.fasta')
    if y > 0:
        print(f'{y} records without reading frame written to {gene}.rf')

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

print("\nUnrecognised Genes")
counter = Counter(unrec_genes)
print(counter)
print('Misc Features')
print(misc_feature)
print("Other Feature Types")
print(other_type)
