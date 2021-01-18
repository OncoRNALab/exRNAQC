#### Running the small RNA seq pipeline

You need to load Python3

```
python smallRNASeq_preprocessing.py -h
usage: smallRNASeq_preprocessing.py [-h] -t {se,pe} -b base_dir -o output_dir
                                    -m string_match [--mismatch {0,1,2,3}]
                                    [-l L] [--subs {no,yes}] [-n N]
                                    [--repeat {yes,no}] -u user_email

Submit smallRNA Seq preprocessing jobs to the UGent HPC (slurm)

optional arguments:
  -h, --help            show this help message and exit
  -t {se,pe}            Single end (se) or paired end (pe) data
  -b base_dir           Base directory where the sample subdirectories are
                        located
  -o output_dir         Directory where output should be created
  -m string_match       String match to select sample folders in base
                        directory e.g. RNA0
  --mismatch {0,1,2,3}  Number of mismatches [default: 0] allowed in first L
                        (option -l) bases
  -l L                  Mismatches (--mismatch) allowed in first L bases
                        [default: 25, should be >= 5]
  --subs {no,yes}       Subsampling needed? [default: no]
  -n N                  Nr of reads to subsample to. If you need the analysis
                        to stop in order to determine the subsampling level
                        first, do not enter a number here (or -n 0)
  --repeat {yes,no}     Repeated analysis? [default: no] If yes: general fastq
                        copying and filtering will not be repeated
  -u user_email         Submitter email address (used for error reporting)
```

If you need to do subsampling, but you do not know upfront to which number: run the pipeline with --subs {yes} but without the -n command (or with -n 0). This will only execute the commands up until the subsampling step and allow you to determine the min number of lines in the fastq e.g. via the following command:
```bash
cat */*_read_count.txt | grep "RNA0" | sort -nk4
```
The number in the last column of the first line will equal the lowest nr of lines in all fastqs -> if you divide this by 4 (fastq has 4 lines for every read) and round down to e.g. the nearest (half a) million, you get the number to subsample to (-n option)
You can now continue the pipeline by running the same command line with the extra -n <number to subsample> option as well as --repeat yes (the latter will avoid copying and filtering the fastqs again)

Multiple jobs/sample:
1) clipping adaptors, quality filtering and collapsing
2) fastqc
3) obtain spikes (spikes.txt, allspikes.txt)
4) annotation of other reads
miRs.txt --> miRNAs
isomiRs.txt --> isomiRs
contam.txt --> contaminants such as tRNAs, piwiRNAs, ...
technical_artifacts.txt
5) some quality control (.pdf, read_count_new.txt, read_length_new.txt)

#### Example for exRNAQC011
```bash
module load Python/3.6.6-intel-2018b 
python smallRNASeq_preprocessing.py -m RNA0 -t se --mismatch 0 -l 25 --subs yes --repeat no -b /data/20181005_exRNAQC011-99504405 -o $VSC_SCRATCH_VO/exRNAQC011 -u example@ugent.be

# subsampling level = 1.5M, restart pipeline:
python smallRNASeq_preprocessing.py -m RNA0 -t se --mismatch 0 -l 25 --subs yes -n 1500000 --repeat yes -b /data/20181005_exRNAQC011-99504405 -o $VSC_SCRATCH_VO/exRNAQC011 -u example@ugent.be
```

#### Output
```bash
├── RNA004105L1-193593406						# Main sample directory
│   ├── dedup_none-subs_atstart						# Subdirectory (dependent on subsampling choices)
│   │   ├── logs							# folder with error and output files for each script
│   │   │   ├── 00_combinecopy_RNA004105L1.err
│   │   │   ├── ...
│   │   │   └── 08_zipfastq_RNA004105L1.out
│   │   ├── redundant_notann.txt					# unique read ids that are both in miR/contam and in notann
│   │   ├── RNA004105L1-193593406_allspikes.txt				# which unique read sequence maps to which spike and how mnay times does the sequence occur
│   │   ├── RNA004105L1-193593406_collapsed_nospikes.fa			# same as collapse.fa, but without spike reads
│   │   ├── RNA004105L1-193593406_collapse.fa				# FASTA file where identical sequences of FASTQ are collapsed: each unique sequence gets a unique id (first number behind ">") followed by the number of reads with this sequence
│   │   ├── RNA004105L1-193593406_contam.txt				# contaminant reads: rRNA, snRNA, misc_RNA, tRNA, piRNA counts (sample_id, read_id, counts, type, ensembl_id) 
│   │   ├── RNA004105L1-193593406_isomiRs.txt				# isomiR reads 
(sample_id, read_id, counts, MIMAT_id, 5'offset, 3'offset, sequence)
│   │   ├── RNA004105L1-193593406_mapped.sam				# bowtie output with mapped reads (no spikes)
│   │   ├── RNA004105L1-193593406_miRs.txt				# counts per miRNA (MIMAT) = sum of isomiR counts
│   │   ├── RNA004105L1-193593406_not_annotated_nomiR_nocontam.txt	# not annotated (but mapped) reads without redundant reads (sse redundant_notann.txt)
│   │   ├── RNA004105L1-193593406_not_annotated.txt			# not annotated (but mapped) reads (read_id, counts, sequence, chrommosome, start, end, strand)
│   │   ├── RNA004105L1-193593406.pdf					# length distributions per type
│   │   ├── RNA004105L1-193593406_qc_subs.fastq.gz			# subsampled (quality filtered) FASTQ
│   │   ├── RNA004105L1-193593406_read_count_new.txt			# nr of mapped (without spikes!), miR, contam, and notann reads
│   │   ├── RNA004105L1-193593406_read_length_new.txt			# read length per unique read
│   │   ├── RNA004105L1-193593406_spikes.txt				# sum of all reads per spike
│   │   ├── RNA004105L1-193593406_technical_artifacts.txt		# technical artifacts
│   │   └── scripts							# folder with all the scripts (to repeat/check analyses)
│   │       ├── 00_combinecopy_RNA004105L1.sh
│   │       ├── ...
│   │       └── 08_zipfastq_RNA004105L1.sh
│   ├── FASTQC_original							# FASTQC of original FASTQ
│   │   ├── RNA004105L1-193593406_fastqc.html
│   │   └── RNA004105L1-193593406_fastqc.zip
│   ├── FASTQC_qc							# FASTQC of trimmed + quality filtered FASTQ
│   │   ├── RNA004105L1-193593406_qc_fastqc.html
│   │   └── RNA004105L1-193593406_qc_fastqc.zip
│   ├── RNA004105L1-193593406_clipped.fastq.gz				# FASTQ after adapter trimming
│   ├── RNA004105L1-193593406.fastq.gz					# original FASTQ
│   ├── RNA004105L1-193593406_qc.fastq.gz				# FASTQ after adapter trimming & quality control
│   └── RNA004105L1-193593406_line_count.txt				# number of lines in the different FASTQ files (should be divided by 4 to get the number of reads)
```
