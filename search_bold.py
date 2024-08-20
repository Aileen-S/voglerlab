### Script to download sequences and associated metadata from BOLD ###

import argparse
import requests

def search_bold(taxon):
    # API URL
    combined_url = f"http://www.boldsystems.org/index.php/API_Public/combined?taxon={taxon}&format=tsv"

    # Download metadata
    print(f"Downloading data for taxon: {taxon}")
    response = requests.get(combined_url)
    tsv_data = response.text
    return tsv_data

def main():
    parser = argparse.ArgumentParser(description="Download barcodes and metadata from BOLD Systems.")
    parser.add_argument("-t", "--taxon", required=True, help="Specify taxon")
    parser.add_argument("-o", "--output", help="Output TSV file path")

    args = parser.parse_args()

    data = search_bold(args.taxon)

    # Write to file
    if args.output:
        output = args.output
    else:
        output = f"{args.taxon}_metadata.tsv"

    with open(output, "w") as file:
        file.write(data)

    if args.output:
        print(f"Metadata saved to {args.output}")
    else:
        print(f"Metadata saved to {args.taxon}_metadata.tsv")

if __name__ == "__main__":
    main()
