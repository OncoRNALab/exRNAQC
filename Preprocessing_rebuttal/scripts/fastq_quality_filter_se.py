import argparse
import gzip
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Filter FASTQ reads by quality.")
    parser.add_argument('-q', type=int, required=True, help="Minimum quality score to keep.")
    parser.add_argument('-p', type=int, required=True, help="Minimum percent of bases that must have [-q] quality.")
    parser.add_argument('-z', action='store_true', help="Compress output with GZIP.")
    parser.add_argument('-i', type=str, default=None, help="Input FASTA/Q file. Defaults to STDIN.")
    parser.add_argument('-o', type=str, default=None, help="Output FASTA/Q file. Defaults to STDOUT.")
    parser.add_argument('-v', action='store_true', help="Verbose mode - report the number of sequences processed.")
    return parser.parse_args()


def calculate_percentile_quality(qualities, min_quality, min_percent):
    """Calculate if a sequence meets the required quality score and percentage."""
    total_bases = len(qualities)
    bases_above_quality = sum(q >= min_quality for q in qualities)
    percent_above_quality = (bases_above_quality / total_bases) * 100
    return percent_above_quality >= min_percent


def filter_fastq_records(input_file, output_file, min_quality, min_percent, compress_output, verbose):
    num_input_reads = 0
    num_output_reads = 0

    if input_file:
        if input_file.endswith('.gz'):
            handle = gzip.open(input_file, 'rt')
        else:
            handle = open(input_file, 'r')
    else:
        handle = sys.stdin

    if output_file:
        if compress_output:
            out_handle = gzip.open(output_file, 'wt')
        else:
            out_handle = open(output_file, 'w')
    else:
        out_handle = sys.stdout

    fastq_iterator = SeqIO.parse(handle, "fastq")

    for record in fastq_iterator:
        num_input_reads += 1
        qualities = record.letter_annotations["phred_quality"]

        if calculate_percentile_quality(qualities, min_quality, min_percent):
            SeqIO.write(record, out_handle, "fastq")
            num_output_reads += 1

    if verbose:
        discarded = num_input_reads - num_output_reads
        report = (
            f"Quality cut-off: {min_quality}\n"
            f"Minimum percentage: {min_percent}\n"
            f"Input: {num_input_reads} reads.\n"
            f"Output: {num_output_reads} reads.\n"
            f"Discarded {discarded} ({(discarded * 100) / num_input_reads:.2f}%) low-quality reads.\n"
        )
        if output_file:
            print(report)
        else:
            print(report, file=sys.stderr)

    handle.close()
    out_handle.close()


def main():
    args = parse_args()
    filter_fastq_records(args.i, args.o, args.q, args.p, args.z, args.v)


if __name__ == "__main__":
    main()
