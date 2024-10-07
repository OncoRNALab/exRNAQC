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
   git clone https://github.com/your-username/exRNAQC.git
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

## License


## Contact

---

### Notes:
