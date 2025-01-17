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
    - [3. Build Indexes](#3-Build-Indexes)
    - [4. Prepare input data](#4-prepare-input-data)
    - [5. Update configuration files](#5-update-configuration-files)
    - [6. Execute the preprocessing pipelines](#6-execute-the-preprocessing-pipelines)
    - [7. Subsampling and Repeat Analysis](#7-Subsampling-and-Repeat-Analysis)
5. [Pipeline Workflow](#Pipeline-Workflow)
6. [Output Folder Structure](#output-folder-structure)
7. [License](#license)
8. [Contact](#contact)
   
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
│   └── supl_files.tar.gz        # fasta files required for building the indexes
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

### 1. Clone the repository:
   ```bash
   git clone https://github.com/OncoRNALab/exRNAQC.git
   cd exRNAQC/Preprocessing_rebuttal
   ```

### 2. Ensure you have [Singularity](https://sylabs.io/singularity/) installed on your system

This is necessary to build and run the containers used for preprocessing.

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

Run the containers with sufficient memory (60 GB).

For **full RNA-seq** preprocessing:

```bash
singularity run full.sif -mem 60000M
```

For **small RNA-seq** preprocessing:

```bash
singularity run small.sif -mem 60000M
```
### 3. Build Indexes:

1. Extract files in the resource subdirectory
```bash
cd resources/
tar -xzvf resources_files.tar.gz 
```
3. STAR index (v2.7.11b)
   
The UCSC hg38.fa file was concatenated with ribosomal, spike-ins and chrIS fasta files to make the final fasta and gtf files for building the index. 

```bash
cat hg38.fa rDNA_2.fa ERCC92.fa chrIS.fa > GRCh38_ucsc_lift_ERCC_chrIS_r45S.fa
cat Homo_sapiens.GRCh38.100_ucsc-named.gtf ERCC92.gtf RNAsequins.v2.2.gtf > GRCh38_ucsc_lift_ERCC_chrIS_r45S.gtf
```
Build BED file
```bash
awk -F"\t" 'BEGIN{OFS="\t"} $3 == "exon" { print $1,$4-1,$5,$6,$8,$7 }' GRCh38_ucsc_lift_ERCC_chrIS_r45S.gtf > GRCh38.109_exons_unsorted.bed
sortBed -i GRCh38.109_exons_unsorted.bed > GRCh38.109_exons_sorted.bed
mergeBed -i GRCh38.109_exons_sorted.bed -s -c 4,5,6 -o distinct > Homo_sapiens.GRCh38.109.chrIS_spikes_45S_exons_sorted_merged.bed
```
Build STAR index

```bash
STAR   --runMode genomeGenerate   \
    --runThreadN 12 \
    --genomeDir path/to/star_index  \
    --genomeFastaFiles path/to/GRCh38_ucsc_lift_ERCC_chrIS_r45S.fa \
    --sjdbGTFfile path/to/GRCh38_ucsc_lift_ERCC_chrIS_r45S.gtf \
    --sjdbOverhang 74   \
    --genomeChrBinNbits 14 
```

2. Kallisto index

Similarly, the fasta file for building the index was done concatenating Ensembl version 91 of the human genome, ERCC and Sequin spike sequences.

```bash
kallisto index -i path/to/index path/to/fasta
```
3. Bowtie index

Similarly, the UCSC hg38.fa with added spike-ins sequences was used to built the index. 

```bash
bowtie-build -f path/to/fasta path/to/index
```

### 4. Prepare input data

Before running the main pipelines, create a text file containing the sample names and paths to the FASTQ files (paired-end or single-end). You can use the `get_sample_list.py` script to generate this file.

For example, to generate a sample list for small RNA data:

```bash
python3 get_sample_list.py -i ../SmallRNA_DataTest/ -o ../SmallRNA_DataTest/
```

### 5. Update configuration files

Navigate to the `scripts` directory and review the YAML configuration files (`config_full.yaml` and `config_small.yaml`). You will need to:

1. Set the correct paths to the reference index and annotation files (e.g., `.gtf` files).
2. Set `repeat_analysis` to `no` unless repeat analysis is required.
3. Set the path to the text file containing sample names and paths

```bash
cd scripts/
```

Ensure that your sample list format matches the example in `file_list.txt`.

### 6. Execute the preprocessing pipelines

Once everything is set up, you can run the pipelines as follows:

For **full RNA-seq**:

```bash
python3 FullRNA_preprocessing.py --config config_full.yaml
```

For **small RNA-seq**:

```bash
python3 SmallRNA_preprocessing.py --config config_small.yaml
```
### 7. Subsampling and Repeat Analysis

1. Calculate the level of subsampling desired. Navigate to the output directory and run:

```bash
#fullRNAseq
for sample in $(ls | grep "RNA0"); do lines=4; qcfil_lines=`zcat ${sample}/trimming/${sample}_1_trim_len_qc.fastq.gz | wc -l`; echo $((qcfil_lines/lines)) $(echo $sample); done >> lines_fastq.txt
sort -nk1 lines_fastq.txt | head
#smallRNAseq
for sample in $(ls | grep "RNA0"); do lines=4; qcfil_lines=`zcat ${sample}/${sample}_qc.fastq.gz | wc -l`; echo $((qcfil_lines/lines)) $(echo $sample); done >> lines_fastq.txt
sort -nk1 lines_fastq.txt | head
```
2. Add the number of reads for subsampling to the yaml config file.
3. Change the repeat_analysis to yes in the config yaml file
```
subsample_to_nr: "700000"
repeat_analysis: "yes"
```
5. Run the pipeline again.
- preprocess of FulRNA test dataset (with 100K reads) takes aprox. 16min with 64gb of RAM
- preprocess of smallRNA test dataset (with 10K reads) takes aprox. 1min with 64gb of RAM

## Pipeline Workflow
### Full length RNA 
#### 1.	Load Parameters:
- Parse command-line arguments to get the path to a YAML configuration file.
- Load the configuration parameters from the YAML file.
#### 2.	Initialize Variables:
- Extract various parameters from the configuration file, such as output directory, data type, deduplication method, subsampling timing, and others.
#### 3.	Load STAR Index:
- Load the STAR index into memory for alignment purposes.
#### 4.	Read Sample List:
- Read the sample list file to get pairs of sample IDs and their corresponding read files (R1 and R2).
#### 5.	Process Each Sample:
- For each sample, perform the following steps:
  - Print the sample ID and read files being processed.
  - Determine the deduplication mode and subsampling timing.
  - Perform quality control (QC) and trimming:
  - Run FastQC on the original reads.
  - If `adapter contamination` is specified, trim the adapters.
  - Trim the reads to a specified length.
  - Optionally, remove reads with bad quality.
  - Run FastQC on the trimmed reads.
  - If `repeat analysis` is enabled, skip the QC and trimming steps.
  - If `single-end data` and `fragment length mean` or `standard deviation` are not provided, exit with an error.
  - Perform `subsampling` if specified.
  - Perform `deduplication`:
    - If using Clumpify, align reads with STAR, infer strandedness, remove duplicates, and run Kallisto quantification.
    - If using Picard, align reads with STAR, deduplicate with Picard, convert BAM to FASTQ, and run Kallisto quantification.
  - Sort BAM files by name.
  - Run HTSeq for quantification and generate index statistics.
#### 6.	Unload STAR Index:
- Unload the STAR index from memory after processing all samples.
#### 7.	Handle Subsampling Level Determination:
  - If `subsampling` is set to start but no subsampling level is provided, exit with instructions to determine the subsampling level after quality filtering.

### SmallRNA
#### 1. Directory Setup:
- For each sample in the `file_list.txt`, the script will:
  - Create the necessary directories for output, scripts, and logs.
  - Determine the subsampling mode based on the `subsampling` parameter.

#### 2. Initial Processing:
- If `repeat_analysis` is set to `'no'`:
  - Combine and copy files using **combinecopy**.
  - Perform initial quality control using **fastqc**.
  - Trim adapters using **cutadapt_smallRNA**.
  - Filter out low-quality reads using **removebadqc**.
  - Perform another round of quality control using **fastqc**.
  
- If `repeat_analysis` is set to `'yes'`:
  - Prepare for repeat analysis using **preprepeat**.

#### 3. Subsampling:
- If `subsampling` is `'yes'` and `subsample_to_nr` is not `'0'`:
  - Subsample reads using **subsample**.
  - Perform quality control on subsampled reads using **fastqc**.

#### 4. Further Processing:
- Collapse identical sequences using **prepcounts**.
- Map sequences using **mapbowtie**.
- Annotate mapped reads using **annotate**.
- Gzip all FASTQ files using **zipfastq**.

#### 5. Completion:
- Print a message indicating the completion of processing for each sample.
- If `subsampling` is `'yes'` and `subsample_to_nr` is `'0'`, exit the program with a message to determine the subsampling level after quality filtering.

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
Out_SmallRNA/
└── RNA014937/                            # Folder named after the sample
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
