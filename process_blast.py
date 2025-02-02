import argparse, argcomplete
from Bio import Entrez, SeqIO
import time

# Search GenBank, wait and try again if 'too many requests' error
def search_genbank(ids, retries=60, delay=40):
    for attempt in range(0, retries):
        try:
            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            results = SeqIO.parse(handle, "gb")
            return results
        except Entrez.HTTPError:
            print(f"HTTP error fetching records; retry in {delay} seconds")
            time.sleep(delay)
    print(f"Failed to retrieve records after {retries} attempts.")
    return None

# Paths to the files and directories
parser = argparse.ArgumentParser(description="Filter BLAST file to get max and min sequence position from hits. Gets longest for each TXID" 
                                "BLAST outfmt '6 staxids sacc sstart send'. Result can be used to extract full sequences from BLAST database.")
parser.add_argument("-i", "--input", type=str, help="BLAST output file")
parser.add_argument("-o", "--output", type=str, help="Output file with match coordinates")
parser.add_argument("-l", "--longest", action='store_true', help="Only keep longest sequence for each txid")
parser.add_argument("-e", "--email", type=str, help="Your email address for NCBI")

argcomplete.autocomplete(parser)
args = parser.parse_args()

Entrez.email = args.email

input = args.input
output = args.output

blast_hits = {}

with open(input, "r") as file:
    data = file.readlines()
    for line in data:
        line = line.strip()
        record = line.split("\t") # Record format 'TXID GBID start end'
        txid = record[0]
        gbid = record[1]
        start = int(record[2])
        end = int(record[3])
        
        if args.longest:
            # Check for multiple TXIDs
            if ';' in txid:
                results = search_genbank(gbid)
                try:
                    for r in results:
                        if "db_xref" in r.features[0].qualifiers:
                            db_xref = r.features[0].qualifiers["db_xref"]
                            for ref in db_xref:
                                if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                                    txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
                        else:
                            print(f'Error getting TXID for {gbid}')
                except TypeError:
                    print(f'Error getting TXID for {gbid}')
                    continue

        # Swap the values if the alignment is on the reverse strand
        if start > end:
            start, end = end, start

        # Get match coordinates. Add to other matches for same record.    
        if txid in blast_hits:
            if gbid in blast_hits[txid]:
                # Check for too long matches (errors from genomes)
                if max([end, blast_hits[txid][gbid]['end']]) - min([start, blast_hits[txid][gbid]['start']]) < 4000:
                    # Check if start is lower than current start
                    if start < int(blast_hits[txid][gbid]['start']):
                        blast_hits[txid][gbid]['start'] = int(start)
                    # Check if end is higher than current end
                    if end > int(blast_hits[txid][gbid]['end']):
                        blast_hits[txid][gbid]['end'] = int(end)
                else:
                    # If max coordinated of current and previous match are > 4000bp, save as separate sequence
                    # Merge with current _2 sequence if possible
                    if f'{gbid}_2' in blast_hits[txid]:
                        if max([end, blast_hits[txid][f'{gbid}_2']['end']]) - min([start, blast_hits[txid][gbid]['start']]) < 4000:
                            # Find first start and last end match for each gbid
                            if start < int(blast_hits[txid][f'{gbid}_2']['start']):
                                blast_hits[txid][f'{gbid}_2']['start'] = int(start)
                            if end > int(blast_hits[txid][f'{gbid}_2']['end']):
                                blast_hits[txid][f'{gbid}_2']['end'] = int(end)
                    else:
                        # Save as seperate sequence under TXID. Longest will be chosen in next loop, and suffix removed if necessary.
                        blast_hits[txid][f'{gbid}_2'] = {'start': int(start), 'end': int(end)}
            else:
                blast_hits[txid][gbid] = {'start': int(start), 'end': int(end)}
        else:
            blast_hits[txid] = {gbid: {'start': int(start), 'end': int(end)}}

write = []

if args.longest:
    for txid, gbid in blast_hits.items():
        max_len = 200
        max_acc = ''
        for gb, rec in gbid.items():
            # Get longest match for each TXID
            rec_len = int(rec['end']) - int(rec['start'])
            if rec_len > max_len:
                max_len = rec_len
                max_acc = gb
        #print(f'max_acc = {max_acc}, max_len = {max_len}')
        max_acc = max_acc.replace('_2', '')
        if max_acc != '':
            write.append([max_acc, (blast_hits[txid][max_acc]['start']), blast_hits[txid][max_acc]['end']])
else:
    for txid, gbid in blast_hits.items():
        print(txid)
        for gb, rec in gbid.items():
            gb = gb.replace('_2', '')
            rec_len = int(rec['end']) - int(rec['start'])
            if rec_len > 200:
                write.append([gb, (blast_hits[txid][gb]['start']), blast_hits[txid][gb]['end']])

with open (output,'w') as outfile:
    for rec in write:
        outfile.write(f'{rec[0]}\t{rec[1]}\t{rec[2]}\n')

print(f'Output written to {output}')

# Then:

# Retrieve from BLAST database and save to fasta
# for file in $gene.out
# do
#   input=$file
#   output=$gene.fasta
#   > "$output"
#   while IFS= read -r line
#   do
#     accession=$(echo "$line" | awk '{print $1}')
#     start=$(echo "$line" | awk '{print $2}')
#     end=$(echo "$line" | awk '{print $3}')
 
#     # Run blastdbcmd and append the output to the FASTA file
#     blastdbcmd -db nt -entry "$accession" -range "$start-$end" -target_only -outfmt %f >> "$output"
#   done < "$input"
# done 
