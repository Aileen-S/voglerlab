### Script to download sequences and associated metadata from BOLD v5 ###

import argparse
import requests
import sys
import time
from datetime import datetime
import json


def search_bold_taxon(tax, retries=20, delay=5):
    # API URL
    #combined_url = f"http://www.boldsystems.org/index.php/API_Public/combined?taxon={tax}&format=tsv"

    # Find taxon
    search_1 = f"https://portal.boldsystems.org/api/query/preprocessor?query=tax:{tax}"
    response_1 = requests.get(search_1, timeout=60).json()
    term_1  = response_1["successful_terms"][0]["matched"]
    if term_1.startswith("tax:"):
        print(f"{tax} found in BOLD taxonomy: {term_1}")
    else:
        print(f"{tax} not found in BOLD taxonomy. Exiting.")
        sys.exit()

    # Get search ID
    search_2 = f"https://portal.boldsystems.org/api/query?query={term_1}&extent=full"
    response_2 = requests.get(search_2, timeout=60).json()
    term_2 = response_2['query_id']
    print(f'Recieved search ID: {response_2["query_id"]}')
    
    print(f"Downloading data for taxon: {tax}")
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


def main():
    parser = argparse.ArgumentParser(description="Download barcodes and metadata from BOLD Systems.")
    parser.add_argument("-t", "--taxon", help="Specify taxa (single taxon or comma delimited list)")
    parser.add_argument("-i", "--ids", help="File with list of BOLD IDs, one per line")
    parser.add_argument("-o", "--output", help="Output TSV file path")
    args = parser.parse_args()

    if args.taxon:
        taxa = args.taxon.split(',')
        t = 0
        data = ''
        for taxon in taxa:
            d = search_bold_taxon(taxon)
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