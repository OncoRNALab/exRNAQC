---

# exRNAQC - Data Preprocessing

This repository contains all the scripts and supplementary files required for preprocessing full length RNA-seq and small RNA-seq data as part of the **exRNAQC** project. The repository includes Singularity files to construct container images, main pipelines for data preprocessing, annotation files, and small test datasets to validate the pipelines.

## Table of Contents
1. [Project Structure](#project-structure)
2. [Features](#features)
3. [Installation](#installation)
4. [Usage](#usage)
    - [1. Build the Singularity images](#1-build-the-singularity-images)
    - [2. Run the containers](#2-run-the-containers)
    - [3. Update configuration files](#3-update-configuration-files)
    - [4. Prepare input data](#4-prepare-input-data)
    - [5. Execute the preprocessing pipelines](#5-execute-the-preprocessing-pipelines)
5. [License](#license)
6. [Contact](#contact)

## Project Structure

The repository is structured as follows:

```
.
├── SingSmallRNA.def             # Singularity definition for small RNA-seq preprocessing
├── SingFullRNA.def              # Singularity definition for full RNA-seq preprocessing
├── scripts/
│   ├── FullRNA_preprocessing.py # Main pipeline for full RNA-seq preprocessing
│   ├── SmallRNA_preprocessing.py# Main pipeline for small RNA-seq preprocessing
│   ├── get_sample_list.py       # Helper script to generate sample list
│   └── ...                      # Additional scripts needed for the pipelines
├── resources/
│   └── ...                      # Annotation files required for preprocessing
├── FullRNA_DataTest/            # Test dataset for full RNA-seq pipeline
├── SmallRNA_DataTest/           # Test dataset for small RNA-seq pipeline
└── README.md                    # This readme file
```

## Features

- **Full RNA-seq and Small RNA-seq Data Preprocessing**: Comprehensive pipelines for preprocessing RNA sequencing data.
- **Singularity Container Support**: Containerized environments for managing software dependencies.
- **Customizable**: Configurable YAML files for paths to reference indexes, annotation files, and other pipeline parameters.
- **Sample Testing**: Includes small RNA-seq and full RNA-seq test datasets to validate the pipelines.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/OncoRNALab/exRNAQC.git
   cd exRNAQC
   ```

2. Ensure you have [Singularity](https://sylabs.io/singularity/) installed on your system. This is necessary to build and run the containers used for preprocessing.

## Usage

### 1. Build the Singularity images

To preprocess RNA-seq data, first build the Singularity images for small RNA and full RNA using the provided definition files.

For **small RNA-seq**:

```bash
singularity --debug build small.sif SingSmallRNA.def
```

For **full RNA-seq**:

```bash
singularity --debug build full.sif SingFullRNA.def
```

### 2. Run the containers

Run the containers with sufficient memory (recommended: 60 GB).

For **full RNA-seq** preprocessing:

```bash
singularity run full.sif -mem 60000M
```

For **small RNA-seq** preprocessing:

```bash
singularity run small.sif -mem 60000M
```

### 3. Update configuration files

Navigate to the `scripts` directory and review the YAML configuration files (`config_full.yaml` and `config_small.yaml`). You will need to:

1. Set the correct paths to the reference index and annotation files (e.g., `.gtf` files).
2. Set `repeat_analysis` to `no` unless repeat analysis is required.

```bash
cd scripts/
```

### 4. Prepare input data

Before running the main pipelines, create a text file containing the sample names and paths to the FASTQ files (paired-end or single-end). You can use the `get_sample_list.py` script to generate this file.

For example, to generate a sample list for small RNA data:

```bash
python3 get_sample_list.py -i ../SmallRNA_DataTest/ -o ../SmallRNA_DataTest/
```

Ensure that your sample list format matches the example in `file_list.txt`.

### 5. Execute the preprocessing pipelines

Once everything is set up, you can run the pipelines as follows:

For **full RNA-seq**:

```bash
python3 FullRNA_preprocessing.py --config config_full.yaml
```

For **small RNA-seq**:

```bash
python3 SmallRNA_preprocessing.py --config config_small.yaml
```

## Output Folder Structure

The pipeline generates an output folder for each sample (e.g., `RNA022999`). Here's an overview of the structure and contents of the output directory for each sample:

```
Out_FullRNA/
└── RNA022999/                            # Folder named after the sample
    ├── cutadapt/                         # Contains fastq files after adaptor trimming
    │   ├── RNA022999_1_trim.fastq.gz
    │   ├── RNA022999_2_trim.fastq.gz
    ├── trimming/                         # Contains fastq files after quality filtering and length trimming
    │   ├── RNA022999_1_trim_len_qc.fastq.gz
    │   ├── RNA022999_2_trim_len_qc.fastq.gz
    ├── FASTQC/                           # Quality control reports from FASTQC
    │   ├── RNA022999_1_trim_fastqc.html
    │   ├── RNA022999_2_trim_fastqc.html
    ├── dedup_clumpify-subs_atstart/    # Main output directory after deduplication
    │   ├── RNA022999_1_trim_len_qc_subs.fastq.gz # Quality filtered, subsampled, and deduplicated fastq
    │   ├── RSeQC_output.txt              # Output from RSeQC for strandedness analysis
    │   ├── RNA022999_htout/              # Gene counts output
    │   │   ├── RNA022999.Aligned.sortedByName.out_htseq_counts.txt
    │   ├── RNA022999_idxout/             # Counts per chromosome/spike
    │   │   ├── RNA022999.Aligned.sortedByCoord.out_idxstat.txt
    │   ├── RNA022999_klout/           # Kallisto abundance estimation files
    │   │   ├── abundance.tsv
    │   │   ├── abundance.h5
    │   ├── RNA022999_srout/              # STAR alignment outputs
    │       ├── RNA022999_Aligned.sortedByCoord.out.bam  # Sorted by coordinate
    │       ├── RNA022999_Aligned.sortedByName.out.bam   # Sorted by name
    │       ├── RNA022999.Aligned.sortedByCoord.out.bam.bai # Indexed BAM file
```
```
├── RNA004105						# Main sample directory
│   ├── dedup_none-subs_atstart						# Subdirectory (dependent on subsampling choices)
│   │   ├── logs							# folder with error and output files for each script
│   │   │   ├── 00_combinecopy_RNA004105L1.err
│   │   │   ├── ...
│   │   │   └── 08_zipfastq_RNA004105L1.out
│   │   ├── redundant_notann.txt					# unique read ids that are both in miR/contam and in notann
│   │   ├── RNA004105_allspikes.txt				# which unique read sequence maps to which spike and how mnay times does the sequence occur
│   │   ├── RNA004105_collapsed_nospikes.fa			# same as collapse.fa, but without spike reads
│   │   ├── RNA004105_collapse.fa				# FASTA file where identical sequences of FASTQ are collapsed: each unique sequence gets a unique id (first number behind ">") followed by the number of reads with this sequence
│   │   ├── RNA004105_contam.txt				# contaminant reads: rRNA, snRNA, misc_RNA, tRNA, piRNA counts (sample_id, read_id, counts, type, ensembl_id) 
│   │   ├── RNA004105_isomiRs.txt				# isomiR reads 
(sample_id, read_id, counts, MIMAT_id, 5'offset, 3'offset, sequence)
│   │   ├── RNA004105_mapped.sam				# bowtie output with mapped reads (no spikes)
│   │   ├── RNA004105_miRs.txt				# counts per miRNA (MIMAT) = sum of isomiR counts
│   │   ├── RNA004105_not_annotated_nomiR_nocontam.txt	# not annotated (but mapped) reads without redundant reads (sse redundant_notann.txt)
│   │   ├── RNA004105_not_annotated.txt			# not annotated (but mapped) reads (read_id, counts, sequence, chrommosome, start, end, strand)
│   │   ├── RNA004105.pdf					# length distributions per type
│   │   ├── RNA004105_qc_subs.fastq.gz			# subsampled (quality filtered) FASTQ
│   │   ├── RNA004105_read_count_new.txt			# nr of mapped (without spikes!), miR, contam, and notann reads
│   │   ├── RNA004105_read_length_new.txt			# read length per unique read
│   │   ├── RNA004105_spikes.txt				# sum of all reads per spike
│   │   ├── RNA004105_technical_artifacts.txt		# technical artifacts
│   │   └── scripts							# folder with all the scripts (to repeat/check analyses)
│   │       ├── 00_combinecopy_RNA004105L1.sh
│   │       ├── ...
│   │       └── 08_zipfastq_RNA004105L1.sh
│   ├── FASTQC_original							# FASTQC of original FASTQ
│   │   ├── RNA004105_fastqc.html
│   │   └── RNA004105_fastqc.zip
│   ├── FASTQC_qc							# FASTQC of trimmed + quality filtered FASTQ
│   │   ├── RNA004105_qc_fastqc.html
│   │   └── RNA004105_qc_fastqc.zip
│   ├── RNA004105_clipped.fastq.gz				# FASTQ after adapter trimming
│   ├── RNA004105.fastq.gz					# original FASTQ
│   ├── RNA004105_qc.fastq.gz				# FASTQ after adapter trimming & quality control
│   └── RNA004105_line_count.txt				# number of lines in the different FASTQ files (should be divided by 4 to get the number of reads)
```
## License


## Contact

---

### Notes:
