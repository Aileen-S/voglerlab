#!/usr/bin/env python3

import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
import time
import re
from Bio.Seq import Seq, UndefinedSequenceError
from collections import defaultdict

suborders = ['Adephaga', 'Polyphaga', 'Myxophaga', 'Archostemata']

# Argument parser
parser = argparse.ArgumentParser(description="Search GenBank, retrieve gene sequences and save as fasta.")
parser.add_argument('-i', '--input', type=str, help="Input genbank file, rather than searching NCBI")
parser.add_argument('-o', '--output', type=str, help="Output CSV file")
parser.add_argument("-s", "--save", type=str, help="Output genbank file with initial search results - provide file path")
parser.add_argument("-t", "--taxon", type=str, help="Taxon of interest")
parser.add_argument('-id', '--id_list', type=str, help="Input file with accession or taxon ID list")
parser.add_argument('-r', '--ref', choices=['txid', 'gbid'], help="If using --id_list option, specify accessions or taxon IDs.")
parser.add_argument("-e", "--email", type=str, help="Your email registered with NCBI")


args = parser.parse_args()         # Process input args from command line



def get_gbids(query, chunk=10000, retries=10, delay=30):
    gbids = set()
    terms = [f"txid{txid}" for txid in query] if isinstance(query, list) else [query]
    
    if isinstance(query, list):
        print('Retrieving GenBank record IDs for input taxon ID list')
    else:
        print(f'Retrieving GenBank record IDs for {query}')

    for term in terms:
        for attempt in range(retries):
            try:
                # Search for records
                with Entrez.esearch(db="nucleotide", term=term, retmax=0, usehistory="y") as search_handle:
                    search_record = Entrez.read(search_handle)
                count = int(search_record["Count"])
                web_env = search_record["WebEnv"]
                query_key = search_record["QueryKey"]

                # Fetch record IDs in chunks
                for start in range(0, count, chunk):
                    with Entrez.esearch(db="nucleotide", term=term, retstart=start, retmax=chunk, idtype="acc", webenv=web_env, query_key=query_key) as fetch_handle:
                        fetch_record = Entrez.read(fetch_handle)
                    gbids.update(fetch_record['IdList'])
                    time.sleep(1)
                break
            
            except (Entrez.HTTPError, RuntimeError) as e:
                print(f"Error for {term} on attempt {attempt + 1}: {e} Retrying in {delay} seconds...")
                time.sleep(delay)
        else:
            print(f"Failed to retrieve records for {term} after {retries} attempts.")
            return None

    if not gbids:
        print("Failed to retrieve any records.")
        return None

    print(f'Found {len(gbids)} total GenBank IDs.')
    gbids = {id.split('.')[0] for id in gbids}

    return gbids


def set_search_parameters(args):
    if args.ref == 'gbid':
        gbids = []
        file = open(args.id_list)
        lines = file.readlines()
        for line in lines:
            acc = line.strip()
            gbids.append(acc)
        print(f'{len(gbids)} IDs found in {args.id_list}')

    elif args.ref  == 'txid':
        taxids = []
        file = open(args.id_list)
        lines = file.readlines()
        for line in lines:
            taxid = line.strip()
            taxids.append(taxid)
        print(f'{len(taxids)} IDs found in {args.id_list}')
        gbids = get_gbids(taxids)

    if args.taxon:
        basesearch = f"(\"{args.taxon}\"[Organism] OR \"{args.taxon}\"[All Fields])"
        gbids = get_gbids(basesearch)
        #print(f"Found {len(gbids)} records for {args.taxon}")
       
    return list(gbids)


# Search GenBank with ID list
def search_genbank(ids, chunk_size=500, retries=10, delay=30, save=False, output="records.gb"):
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
            except (Entrez.HTTPError, RuntimeError) as e:
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
    spec = rec.annotations["organism"]
    taxonomy = ['', '', '', '', '', '']
    for tax in rec.annotations["taxonomy"]:
        if tax in suborders: taxonomy[0] = tax
        if tax.endswith('formia'): taxonomy[1] = tax
        if tax.endswith('oidea'): taxonomy[2] = tax
        if tax.endswith('idae'): taxonomy[3] = tax
        if tax.endswith('inae'): taxonomy[4] = tax
        if tax.endswith('ini'): taxonomy[5] = tax
    taxonomy_string = '|'.join(rec.annotations["taxonomy"])
    date = rec.annotations["date"]

    # Location
    if "geo_loc_name" in rec.features[0].qualifiers:
        location = rec.features[0].qualifiers["geo_loc_name"][0]
        if ":" in location:
            country, region = location.split(":", 1)
            country = country.strip()
            region = region.strip()
        else:
            country = location
            region = ""
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
              "date": date,
              "taxonomy": taxonomy,
              "country": country,
              "region": region,
              "lat": lat,
              "long": long,
              "refs": refs,
              "row": [txid, rec.name, date, '', '', ''] + taxonomy + [spec, taxonomy_string, country, region, lat, long] + refs}
    return output


def find_genes(results, args):

    meta = {}
    seqs = {}
    nohits = []
    other_type = set()
    misc_feature = set()
    unrec_genes = {}
    unrec_species = []
    x = 0  # Count taxids

    gene_names = defaultdict(int)
    for rec in results:
        if args.taxon:
            if args.taxon not in rec.annotations["taxonomy"]:
                unrec_species.append(rec.name)
                continue
        output = genbank_metadata(rec)
        meta[rec.name] = output

    return meta



def main():

    # Get records
    if args.input:
        with open(args.input) as file:
            results = list(SeqIO.parse(file, "genbank"))
            print(f'{len(list(results))} records in {args.input}')            
  
    if not args.input:
        Entrez.email = args.email if args.email else None
        gbids = set_search_parameters(args)
        if args.save:
            results = search_genbank(gbids, save=True, output=args.save)
        else:
            results = search_genbank(gbids)

    meta = {}
    for rec in results:
        if args.taxon:
            if args.taxon not in rec.annotations["taxonomy"]:
                unrec_species.append(rec.name)
                continue
        output = genbank_metadata(rec)
        meta[rec.name] = output


    # Write CSV metadata file
    x = 0
    with open(args.output, "w") as file:
        writer = csv.writer(file)
        writer.writerow(
            ["ncbi_taxid", "genbank_accession", "date", "bold_id", "bold_bin", "lab_id", "suborder", "infraorder", "superfamily", "family", 
            "subfamily", "tribe", "species", "taxonomy", "country", "region", "latitude", "longitude", "ref_author", "ref_title", "ref_journal"])
        for gbid, rec in meta.items():
            writer.writerow(rec['row'])
            x += 1
        print(f"Metadata for {x} records saved as CSV file to {args.output}")

if __name__ == "__main__":
    main()