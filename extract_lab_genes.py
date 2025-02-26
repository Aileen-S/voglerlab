import argparse, argcomplete
import csv
from Bio import Entrez
from Bio import SeqIO
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
    specfasta = spec.replace(" ", "_")

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
              "row": [txid, rec.name, '', '', ''] + taxonomy + [rec.annotations["organism"], country, lat, long] + refs}
    return output


# Argument parser
# Add option to find only mito genes, or only selected genes.
parser = argparse.ArgumentParser(description="Search GenBank file, retrieve gene sequences and save as fasta.")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-g', '--gb_file', type=str, help="Input genbank format file")
#parser.add_argument('-m', '--metadata', action='store_true', help='Save metadata')
parser.add_argument('-i', '--fasta_id', choices=['gbid', 'txid', 'both'], help="Choose identifiers for output fastas. Default is gbid.")
parser.add_argument('-l', '--list', type=str, help="Limit to list of db_ids in file")

argcomplete.autocomplete(parser)
args = parser.parse_args()         # Process input args from command line

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

mito = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
rna = ['12S', '16S', '18S', '28S']

genes = {"12S": ["12S", "12S RIBOSOMAL RNA", "12S RRNA", "SSU", "SMALL SUBUNIT RIBOSOMAL RNA", "RRNS"],
         "16S": ["16S", "16S RIBOSOMAL RNA", "16S RRNA", "LARGESUBUNITRIBOSOMALRNA", "LSU", "LARGE SUBUNIT RIBOSOMAL RNA", "RRNL"],
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
         "ND6": ['NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI', 'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0', 'ND6']}

cds = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'H3', 'RNApol', 'Wg']

ids = []
if args.list:
    file = open(args.list)
    lines = file.readlines()
    for line in lines:
        ids.append(line.strip())


# Search through GBIDs
unrec_genes = []
species = {}
#count = {}
with open(args.gb_file) as file:
    record = SeqIO.parse(file, "gb")
    for rec in record:
        if args.list:
            if rec.name not in ids:
                continue
        output = {}
        # if args.metadata:
        #     output = genbank_metadata(rec)
        g = 0
        for feature in rec.features:
            #type = feature.type
            if feature.type in ('CDS', 'rRNA'):
                names = get_feat_name(feature)                       # Find gene name
                #print(names)
                stdname = "none"
                for k, v in genes.items():
                    for name in names:
                        if name in v:
                            stdname = k
                            g += 1

                if stdname == 'none':
                    unrec_genes.append(name)
                    print(f'unrec: {name}')
                if stdname in cds:
                    if 'codon_start' in feature.qualifiers:
                        frame = feature.qualifiers["codon_start"][0]
                    else:
                        print(f"Reading frame missing from record {rec.name}, {stdname}.")
                else:
                    frame = ''
                seq = feature.extract(rec.seq)
                output = {"gbid": rec.name,
                          "gene": stdname,
                          "length": len(seq),
                          "seq": seq,
                          "frame": frame}
                if stdname not in species:
                    species[stdname] = []
                species[stdname].append(output)
        #count[rec.name] = g



# Write CSV metadata file
# if args.metadata:
#     gbids = []
#     with open("metadata.csv", "w") as file:
#         writer = csv.writer(file)
#         writer.writerow(
#             ["ncbi_taxid", "genbank_accession", "bold_id", "bold_bin", "lab_id", "suborder", "infraorder", "superfamily", "family", 
#             "subfamily", "tribe", "species", "country", "latitude", "longitude", "ref_authoer", "ref_title", "ref_journal"])
#         for gene, records in species.items():
#             for rec in records:
#                 if rec['gbid'] not in gbids:
#                     gbids.append(rec['gbid'])
#                     row = [rec['txid'], '', '', '', rec['gbid']] + rec['taxonomy'] + [rec["spec"], rec['country'], rec['lat'], rec['long']] + rec['refs']
#                 writer.writerow(row)
#         print("Metadata saved to metadata.csv")

for gene, records in species.items():
    x = 0
    y = 0
    if gene in rna:
        file = open(f"{gene}.fasta", "w")
        for rec in records:
            fasta_id = f">{rec['gbid']}\n{rec['seq']}\n"
            file.write(fasta_id)
            x += 1

    else:
        file = open(f"{gene}.fasta", "w")
        for rec in records:
            if rec['frame'] == '':
                fasta_id = f">{rec['gbid']}\n{rec['seq']}\n"
            else:
                fasta_id = f">{rec['gbid']};frame={rec['frame']}\n{rec['seq']}\n"

            if 'frame' in fasta_id:
                file.write(fasta_id)
                x += 1
            else:
                if y == 0:
                    rf = open(f"{gene}.rf", "w")
                rf.write(fasta_id)
                y += 1
    print(f'{x} records written to {gene}.fasta')
    if y > 0:
        print(f'{y} records without reading frame written to {gene}.rf')


