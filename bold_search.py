#!/usr/bin/env python3

# Script to download sequences and associated metadata from BOLD v5
from Bio import Entrez
import argparse
import requests
import sys
import time
from datetime import datetime
import json

def get_search_taxa(taxon):
    ncbi_ranks = {"superkingdom": 1, "kingdom": 2, "subkingdom": 3, "superphylum": 4, "phylum": 5, "subphylum": 6, "infraphylum": 7, 
                  "superclass": 8, "class": 9, "subclass": 10, "infraclass": 11, "cohort": 12, "subcohort": 13, "superorder": 14, 
                  "order": 15, "suborder": 16, "infraorder": 17, "parvorder": 18, "superfamily": 19, "family": 20, "subfamily": 21, 
                  "tribe": 22, "subtribe": 23, "genus": 24, "subgenus": 25, "section": 26, "subsection": 27, "series": 28, "subseries": 29, 
                  "species group": 30, "species subgroup": 31, "species": 32, "forma specialis": 33, "subspecies": 34, "varietas": 35, 
                  "subvariety": 36, "forma": 37, "serogroup": 38, "serotype": 39, "strain": 40, "isolate": 41}
    bold_ranks = [2, 5, 9, 15, 20, 21, 22, 24, 32]
    try:
        with Entrez.esearch(db="Taxonomy", term=taxon) as handle:
            record = Entrez.read(handle)
        txid = record["IdList"][0]
    except IndexError:
        print(f"{taxon} not found on NCBI. Exiting.")
        sys.exit()
    handle = Entrez.efetch(db="Taxonomy", id=txid, retmode="xml")
    record = Entrez.read(handle)
    rank = record[0]["Rank"]
    taxon_rank = ncbi_ranks[rank]
    print(f"Taxon rank: {rank}")
    if taxon_rank in bold_ranks:
        return [taxon]
    else:
        for rank_no in bold_ranks:
            if rank_no > taxon_rank:
                index = rank_no
                break
        sub_rank = list(ncbi_ranks.keys())[list(ncbi_ranks.values()).index(index)]
        search_term = f"txid{txid}[Subtree] AND {sub_rank}[RANK]"
        with Entrez.esearch(db="Taxonomy", term=search_term, retmax=500) as handle:
            record = Entrez.read(handle)
        ids = ",".join(record["IdList"])
        with Entrez.efetch(db="Taxonomy", id=ids, retmode="xml") as handle:
            tax_records = Entrez.read(handle)
        sub_rank_list = [r["ScientificName"] for r in tax_records]
        print(f"Searching for {sub_rank} records in NCBI subtree: {sub_rank_list}")
        return sub_rank_list

def find_bold_taxon(tax):
    # Find taxon
    search_1 = f"https://portal.boldsystems.org/api/query/preprocessor?query=tax:{tax}"
    response_1 = requests.get(search_1, timeout=60).json()
    term_1  = response_1["successful_terms"][0]["matched"]
    if term_1.startswith("tax:"):
        print(f"{tax} found in BOLD taxonomy: {term_1}")
        return True
    else:
        print(f"{tax} not found in BOLD taxonomy. Checking NCBI.")
        return False

def retrieve_bold_taxon(tax, retries=20, delay=5):
    # Find taxon
    search_1 = f"https://portal.boldsystems.org/api/query/preprocessor?query=tax:{tax}"
    response_1 = requests.get(search_1, timeout=60).json()
    term_1  = response_1["successful_terms"][0]["matched"]
    # Get search ID
    search_2 = f"https://portal.boldsystems.org/api/query?query={term_1}&extent=full"
    response_2 = requests.get(search_2, timeout=60).json()
    term_2 = response_2['query_id']
    # print(f'Recieved search ID: {response_2["query_id"]}')
    
    print(f"Downloading BOLD data for taxon: {tax}")
    for attempt in range(retries):
        try:
            # Download metadata
            search_3 = f"https://portal.boldsystems.org/api/documents/{term_2}==/download?format=tsv"
            response_3 = requests.get(search_3)
            if response_3.status_code == 200:
                print(f"Records found: {len(response_3.text.splitlines()) - 1}")
                return response_3.text
                break
            else:
                print(f"HTTP error fetching records on attempt {attempt + 1} of {retries}; retry in {delay} seconds")
                time.sleep(delay)

        except (requests.exceptions.RequestException, json.JSONDecodeError, KeyError, IndexError) as e:
            print(f'Error on attempt {attempt + 1}: {e}')
            sleep(delay)
    print(f"Failed to retrieve records after {retries} attempts.")
    sys.exit()

def search_bold_ids(ids, chunk_size=100):
    # Split the list of IDs into chunks
    chunks = [ids[i:i + chunk_size] for i in range(0, len(ids), chunk_size)]

    # Initialize an empty list to store the results
    results = []

    # Make API calls for each chunk
    for i, chunk in enumerate(chunks):
        # API URL
        combined_url = f"http://www.boldsystems.org/index.php/API_Public/combined?ids={'|'.join(chunk)}&format=tsv"
        # print(combined_url)
        # Download metadata
        print(f"Downloading data for IDs in chunk {i}")
        response = requests.get(combined_url)
        tsv_data = response.text

        # If this is not the first chunk, skip the header line
        if i > 0:
            tsv_data = '\n'.join(tsv_data.split('\n')[1:])

        results.append(tsv_data)

    # Join the results together
    combined_results = '\n'.join(results)

    return combined_results


parser = argparse.ArgumentParser(description="Download barcodes and metadata from BOLD v5. Input either taxon or list of BOLD IDs.")
parser.add_argument("-e", "--email", help="email address for NCBI query, if using --taxon")
parser.add_argument("-t", "--taxon", help="Specify taxa (single taxon or comma delimited list)")
parser.add_argument("-i", "--ids", help="File with list of BOLD IDs, one per line")
parser.add_argument("-o", "--output", help="Output TSV file path")
args = parser.parse_args()


def main():
    if args.taxon:
        Entrez.email = args.email if args.email else ''
        taxa = args.taxon.split(',')
        t = 0
        data = ''
        search_taxa = []
        for taxon in taxa:
            if find_bold_taxon(taxon):
                search_taxa.append(taxon)
            else:
                search_taxa += get_search_taxa(taxon)
        for taxon in search_taxa:
            d = retrieve_bold_taxon(taxon)
            if t > 0: 
                d = '\n'.join(d.split('\n')[1:])
            t += 1
            data += d

    elif args.ids:
        with open(args.ids, 'r') as file: 
            ids = file.read().splitlines()
            #ids = '|'.join(ids)
        data = search_bold_ids(ids)

    # Write to file
    if args.output:
        output = args.output
    else:
        date = datetime.today().strftime('%y%m%d')
        output = f"bold_metadata_{date}.tsv"

    with open(output, "w") as file:
        file.write(data)
        print(f"Saved to {output}")


if __name__ == "__main__":
    main()