import argparse, argcomplete
import os
import tempfile
import subprocess
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def parse_args():

    parser = argparse.ArgumentParser(description="Clean and align fasta sequences")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-o", "--output", type=str, help="Output fasta for good sequences")
    parser.add_argument("-f", "--filtered", type=str, help="Optional output fasta for sequences that did not meet threshold")

    parser.add_argument("-a", "--aligned", action='store_true', help="Input fasta is already aligned and contains filter seqeunce labelled 'reference")
    parser.add_argument("-g", "--gene", choices=['COX1a', '18S', '28S'], default='COX1a', help="Specify gene to run alignment and filter from reference filter sequence")
    parser.add_argument("-p", "--profile", action='store_true', help="Optional profile fasta to improve alignment")

    argcomplete.autocomplete(parser)
    return parser.parse_args()


reference_seqs = {"COX1a": "AATAAATAATATAAGATTTTGACTNCTACCCCCNTCNTTAACNTTNCTANTAATAAGAAGAATAGTAGAAAGAGGAGCAGGAACAGGATGAACAGTTTACCCNCCNCTATCATCNAATATTGCNCATGGAGGANCNTCNGTNGATTTAGCNATTTTTAGNCTNCATTTAGCAGGAATTTCNTCAATTTTAGGAGCNGTAAATTTTATTACNACAGTNATTAATATACGACCNANAGGAATAACNTTTGATCGAATACCTTTATTTGTTTGAGCAGTAGNAATTACAGCNNTTCTACTACTNCTATCNTTACCNGTTTTAGCAGGAGCNATTACTATACTATTAACAGATCGAAATTTAAATACNTCATTTTTTGACCCNGCAGGAGGAGGAGACCCAATTCTATACCAACATTTATTT"}


def align_sequences(records, profile=False):
    # Write temporary input file for mafft
    print('Aligning sequences with MAFFT')
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        SeqIO.write(records, temp_input, "fasta")
        temp_input.flush()
        input = temp_input.name
    # Get thread count
    threads = str(os.cpu_count())
    threads = 1 if threads is None else threads
    # Call MAFFT
    if profile:
        command = ['mafft', '--addfragments', input, '--retree', '1', '--maxiterate', '0', '--adjustdirection', '--anysymbol', '--thread', str(threads), profile]
    else:
        command = ['mafft', '--retree', '1', '--maxiterate', '0', '--adjustdirection', '--anysymbol', '--thread', str(threads), input]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        print("MAFFT failed:")
        print("stderr:", result.stderr)
        return []
    aligned_io = StringIO(result.stdout)
    results = list(SeqIO.parse(aligned_io, "fasta"))
    return results


def filter_to_reference(records):
    print('Removing sequences without reference sequence')
    # Find reference start and end position
    reference = next(rec for rec in records if 'reference' in rec.id)
    start = next((i for i, c in enumerate(reference.seq) if c != '-'), None)
    end = len(reference.seq) - next((i for i, c in enumerate(reversed(reference.seq)) if c != '-'), None)
    ref_length = end - start

    # Filter seqeunces
    output = []
    removed = []
    for rec in records:
        start = next((i for i, c in enumerate(reference.seq) if c != '-'), None)
        end = len(reference.seq) - next((i for i, c in enumerate(reversed(reference.seq)) if c != '-'), None)
        region = rec.seq[start:end].replace('N', '').replace('-', '')

        # Keep those with at least 95% of reference region
        if len(region) >= ref_length * 0.95:
            output.append(rec)
        else:
            removed.append(rec)
    print(f'{len(output)} sequence remain after filtering')
    return removed, output


def main():
    args = parse_args()
    records = list(SeqIO.parse(args.input, 'fasta'))
    print(f'{len(records)} sequences read from {args.input}')
    if args.gene:
        records.append(SeqRecord(Seq(reference_seqs[args.gene]), id='reference', description='reference'))
        records = align_sequences(records, profile=args.profile)
    removed,records = filter_to_reference(records)

    SeqIO.write(records, args.output, 'fasta')
    print(f'{len(records)} sequences written to {args.output}')
    if args.filtered:
        SeqIO.write(removed, args.filtered, 'fasta')
        print(f'{len(removed)} sequences written to {args.filtered}')


if __name__ == "__main__":
    main()
