# exRNAQC

Quick access to HTML reports with results for each study and the list of figures: [here](https://oncornalab.github.io/exRNAQC/))

## Full Length RNAseq
### Rationale
Used for RNA access libraries, but with adaptations one can also use it for other library preps. Before starting the analysis you need to know if the data is stranded or unstranded and single or paired end sequenced. The latter, you will notice in the fastq files. If the strandedness is unsure, you can always load a STAR output file in IGV.

### Setup
- To run it in the HPC use latest script: RNASeq_preprocessing.py (see Preprocessing/FullLengthRNASeq). For more info on how to run it, look at README.md in same folder
- To run it locally see instructions in Preprocessing_rebuttal

### MultiQC
For quality control, you can run MultiQC in the parent directory when the pipeline is finished.
```bash
multiqc -f .
```

### Further processing
See Rmarkdown files (exRNAQC004 and exRNAQC005).

You can download the data locally on your computer with RSync (PICARD_output.txt, .rds files for coverage analysis, abundance.tsv files for kallisto counts):
```
Rsync -zvar --exclude="*gz" --exclude="*fa" --exclude="*sam" --exclude="*bam" --exclude="*/*srout*/_STAR*" --exclude="*/*srout*/Unmapped*" --exclude="*/fastqc*" hpc:/PATH .
```

## Small RNAseq
### Setup
- To run it in the HPC use latest script: smallRNASeq_preprocessing.py (see Preprocessing/SmallRNAseq). For more info on how to run it, look at README.md in same folder
- To run it locally see instructions in Preprocessing_rebuttal

### Further processing
See Rmarkdown files (exRNAQC011 and exRNAQC013).


## Preprocessing exceptions
### In case of manual demultiplexing
The pipeline only runs if the folder layout has a specific structure: `NSQ_RunXXX/RNA002019-123456/RNA002019-123456_S1_L001_R1_001.fastq.gz`. If the run is not downloaded from Basespace, but all fastq files are in 1 folder (so the folder should only contain files with `*fastq.gz`, try running:
	
```bash
rename _R1_001.fastq.gz _L001_R1_001.fastq.gz *gz && # Optional, if the lane is not in the filename (e.g. demultiplexing without splitting lanes).
for sample in $(ls *fastq.gz | cut -d'_' -f 1 | uniq); 
 do
   current_time=$(date "+%H%M%S") && 
   mkdir ${sample}-${current_time} &&
   mv $sample*R1_001.fastq.gz ${sample}-${current_time}; 
 done
 ```
 
### No RNA numbers in file name
It can occur that the filename does not contain an RNA number, but a sample name (e.g. `exRNAQC-005-D3-Citrate-t16_S1_L001_R1_001.fastq.gz` instead of `RNA004898S3_S1_L001_R1_001.fastq.gz`
There is a different `BaseSpaceDownloader.py` located under `./Resources` in this repo that changes the folder names to this RNA number. The files can then be renamed with `renameFq.py`
