
import os
import re
import argparse


parser = argparse.ArgumentParser(description='Submit RNA Seq preprocessing jobs to the UGent HPC (slurm)')

parser.add_argument('-i', required=True, help='Folder containing all the fastq files', metavar='base_dir')
parser.add_argument('-o', required=True, help='Directory where output should be created', metavar='output_dir')
parser.add_argument('-m', required=False, default='RNA', help='String to match in filenames', metavar='match_string')

args = parser.parse_args()
base_dir = args.i
out_dir = args.o
output_file = os.path.join(out_dir, "file_list.txt")
match_string = args.m

def find_paired_files(base_dir, match_string):
    r1_files = []
    r2_files = []

    print(f"Searching for files in: {base_dir}")

    # List all files and match those containing 'RNA'
    for root, _, files in os.walk(base_dir):
        # print(f"Inspecting directory: {root}")  # Debug print
        for file in files:
            pattern = rf"^({match_string}\d+)_"  # Match the sample ID
            match = re.match(pattern, file)
            if match:
                sampleID = match.group(1)
                # print(f"Sample ID: {sampleID}")
            if match_string in file and file.endswith('.fastq.gz'):
                if '_R1' in file:
                    r1_files.append((sampleID, os.path.join(root, file)))
                elif '_R2' in file:
                    r2_files.append((sampleID, os.path.join(root, file)))

    # Ensure that for every R1 there is an R2
    r1_files.sort(key=lambda x: x[0])
    r2_files.sort(key=lambda x: x[0])

    paired_files = []
    
    for (sampleID, r1) in r1_files:
        if len(r2_files)>0:
            r2_file_name = os.path.basename(r1).replace('_R1', '_R2')
            r2 = os.path.join(os.path.dirname(r1), r2_file_name)
            if any(r2 == r2_file for _, r2_file in r2_files):
                paired_files.append((sampleID, r1, r2))
        else:
            paired_files.append((sampleID, r1, None))
    return paired_files


def main():
    paired_files = find_paired_files(base_dir, match_string=match_string)
    with open(output_file, 'w') as f:
        for sampleID, r1, r2 in paired_files:
            f.write(f"{sampleID}\t{r1}\t{r2}\n")

    print(f"Paired files written to {output_file}")

if __name__ == "__main__":
    main()
