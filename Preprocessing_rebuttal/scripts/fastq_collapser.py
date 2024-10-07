import argparse
import gzip
from collections import Counter
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Collapse identical sequences in a FASTA/FASTQ file.")
    parser.add_argument('-i', type=str, default=None, help="Input FASTA/Q file. Defaults to STDIN.")
    parser.add_argument('-o', type=str, default=None, help="Output FASTA/Q file. Defaults to STDOUT.")
    parser.add_argument('-v', action='store_true', help="Verbose mode - print summary of input/output counts.")
    return parser.parse_args()


def collapse_sequences(input_file, output_file, verbose):
    sequence_counter = Counter()

    # Open input file, either gzipped or plain
    if input_file:
        if input_file.endswith('.gz'):
            handle = gzip.open(input_file, 'rt')
        else:
            handle = open(input_file, 'r')
    else:
        handle = sys.stdin

    # Parse the input sequences
    fastx_iterator = SeqIO.parse(handle, "fasta" if input_file.endswith(".fasta") else "fastq")

    # Count identical sequences
    for record in fastx_iterator:
        sequence_str = str(record.seq)
        sequence_counter[sequence_str] += 1

    handle.close()

    # Open output file, either gzipped or plain
    if output_file:
        if output_file.endswith('.gz'):
            out_handle = gzip.open(output_file, 'wt')
        else:
            out_handle = open(output_file, 'w')
    else:
        out_handle = sys.stdout

    # Sort the sequences by their abundance
    sorted_sequences = sorted(sequence_counter.items(), key=lambda x: x[1], reverse=True)

    # Write the collapsed sequences
    counter = 0
    total_reads = 0
    for seq, count in sorted_sequences:
        counter += 1
        total_reads += count
        out_handle.write(f">{counter}-{count}\n{seq}\n")

    if verbose:
        print(f"Input: {len(sequence_counter)} unique sequences (representing {total_reads} reads).")
        print(f"Output: {counter} sequences.")

    if output_file:
        out_handle.close()


def main():
    args = parse_args()
    collapse_sequences(args.i, args.o, args.v)


if __name__ == "__main__":
    main()
