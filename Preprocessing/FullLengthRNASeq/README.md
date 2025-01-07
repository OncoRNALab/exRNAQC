#### Running the RNA access pipeline
You need to load Python3.6

```
$ python RNASeq_preprocessing.py -h
usage: RNASeq_preprocessing.py [-h] -t {se,pe} -b base_dir -o output_dir 
                                -m string_match -s {yes,no} 
                                --dedup {no,picard,clumpify} 
                                [--pico {yes,no}]
                                [--subs {no,start,afterclumpify}] [-n N]
                                [--repeat {yes,no}] -u user_email

Submit RNA Seq preprocessing jobs to the UGent HPC (slurm)

optional arguments:
  -h, --help            show this help message and exit
  -t {se,pe}            Single end (se) or paired end (pe) data
  -b base_dir           Base directory where the sample subdirectories are
                        located
  -o output_dir         Directory where output should be created
  -m string_match       String match to select sample folders in base
                        directory e.g. RNA0
  -s {yes,no}           (reverse) stranded sequencing or not (unstranded)?
  -l L                  Normal read length (longer reads will be trimmed)
  --dedup {no,picard,clumpify}
                        Method of duplicate removal [no, clumpify, picard]
  --pico {yes,no}       PicoV2 sequencing? [default: no]
  --polyA {yes,no}      PolyA trimming needed? [default: no]
  --adapter {truseq,nextera,no}
                        Adapter trimming needed? If so, which one [default:
                        no]
  --subs {no,start,afterclumpify}
                        Subsampling needed? [default: no] If so, when? At
                        start or after clumpify duplicate removal
  -n N                  Nr of reads to subsample to. If you need the analysis
                        to stop in order to determine the subsampling level
                        first, do not enter a number here (or -n 0)
  --coverage {yes,no}   Coverage analysis needed? [default: no]
  --fragmentmean FRAGMENTMEAN
                        Mean fragment length (only necessary for kallisto in
                        SE mode)
  --fragmentsd FRAGMENTSD
                        Mean fragment length (only necessary for kallisto in
                        SE mode)
  --repeat {yes,no}     Repeated analysis? [default: no] If yes: general fastq
                        copying and filtering will not be repeated
  -u user_email         Submitter email address (used for error reporting)
```

If you need to do subsampling, but you do not know upfront to which number: run the pipeline with --subs {start, afterclumpify} but without the -n command (or with -n 0). This will only execute the commands up until the subsampling step and allow you to determine the min number of lines in the fastq e.g. via the following command (executed from the output directory)
```bash
$ for sample in $(ls | grep "RNA0"); do qcfil_lines=`wc -l ${sample}/${sample}_1_qc.fastq`; echo $qcfil_lines $sample; done >> lines_fastq.txt
$ sort -nk1 lines_fastq.txt | head
```
OR (in case of dup removal after clumpify) 
```bash
$ for sample in $(ls | grep "RNA0"); do qcfil_lines=`wc -l ${sample}/dedup_clumpify-subs_afterdedup/${sample}_1_qc_clumped.fastq`; echo $qcfil_lines $sample; done >> lines_fastq.txt
$ sort -nk1 lines_fastq.txt | head
```
The number on the first line will equal the lowest nr of lines in all fastqs -> if you divide this by 4 (fastq has 4 lines for every read) and round down to e.g. the nearest million, you get the number to subsample to (-n option)
You can now continue the pipeline by running the same command line with the extra -n <number to subsample> option as well as --repeat yes (the latter will avoid copying and filtering the fastqs again)
  
#### Example for exRNAQC004
all samples start with "RNA0", paired end, stranded, duplicate removal with Clumpify/Picard after subsampling
```
### Start Clumpify pipeline with --subs start option, but no -n option yet (-> results in subfolder dedup_clumpify-subs_atstart)
$ python RNASeq_preprocessing.py -m RNA0 -t pe -s yes --dedup clumpify --subs start --repeat no -b $VSC_DATA_VO/NSQ_Run479_Run481_Run482 -o $VSC_SCRATCH_VO/exRNAQC004/preprocessing -u example@ugent.be

### Once all jobs are finished: determine subsampling level 
$ for sample in $(ls | grep "RNA0"); do qcfil_lines=`wc -l ${sample}/${sample}_1_qc.fastq`; echo $qcfil_lines $sample; done >> lines_fastq.txt
sort -nk1 lines_fastq.txt | head

### Restart Clumpify pipeline with -n <nr to subsample> and --repeat option
$ python RNASeq_preprocessing.py -m RNA0 -t pe -s yes --dedup clumpify --subs start -n 21000000 --repeat yes -b $VSC_DATA_VO/exRNAQC/data/exRNAQC004/NSQ_Run479_Run481_Run482 -o $VSC_SCRATCH_VO/projects/exRNAQC004/preprocessing -u annelien.morlion@ugent.be

### Same options can be used to run Picard deduplication in parallel (-> new subfolder dedup_picard-subs_atstart), only --dedup option changed to picard
$ python RNASeq_preprocessing.py -m RNA0 -t pe -s yes --dedup picard --subs start -n 21000000 --repeat yes -b $VSC_DATA_VO/exRNAQC/data/exRNAQC004/NSQ_Run479_Run481_Run482 -o $VSC_SCRATCH_VO/projects/exRNAQC004/preprocessing -u annelien.morlion@ugent.be
```

#### Count number of reads
After the pipeline has finished, the number of initial reads, after subsampling and after duplicate removal can be calculated with: 
```bash
for sample in $(ls | grep "RNA0"); do lines=4; qcfil_lines=`zcat ${sample}/${sample}_1.fastq.gz | wc -l`; echo $((qcfil_lines/lines)) $(echo $sample | cut -d'-' -f1); done >> lines_fastq_pre.txt # initial
for sample in $(ls | grep "RNA0"); do lines=4; qcfil_lines=`zcat ${sample}/dedup_clumpify-subs_atstart/${sample}_1_qc_subs.fastq.gz | wc -l`; echo $((qcfil_lines/lines)) $(echo $sample | cut -d'-' -f1); done >> lines_fastq_subs.txt # after SS
for sample in $(ls | grep "RNA0"); do lines=4; qcfil_lines=`zcat ${sample}/dedup_clumpify-subs_atstart/${sample}_1_qc_subs_clumped.fastq.gz | wc -l`; echo $((qcfil_lines/lines)) $(echo $sample | cut -d'-' -f1); done >> lines_fastq_post.txt # after dedup
```

#### MultiQC
For quality control, you can run MultiQC in the parent directory when the pipeline is finished.
```bash
ml MultiQC/1.7-intel-2018b-Python-3.6.6
ml swap matplotlib/2.2.3-intel-2018b-Python-3.6.6
multiqc -f .
```

#### Output
In the output_dir, there will be a number of files created. Depending on the chosen deduplication and subsampling options, output for each sample will be placed in a different subfolder
1. with clumpify deduplication: dedup_clumpify-subs_none, dedup_clumpify-subs_atstart, dedup_clumpify-subs_afterdedup
2. with picard deduplication: dedup_picard-subs_none, dedup_picard-subs_atstart
3. without deduplication: dedup_none-subs_none, dedup_none-subs_atstart

The directory structure will look like this (similar for no dedup or no/later subsampling):

```bash
├── RNA003584L1-173868704                                         # Main sample directory
│   ├── dedup_clumpify-subs_atstart                               # Deduplication + subsampling options directory
│   │   ├── 25_klout                                              # Folder containing Kallisto abundance files 
│   │   │   ├── abundance.h5                                                # (based on deduplicated fastqs)
│   │   │   ├── abundance.tsv
│   │   │   └── run_info.json
│   │   ├── 26_srout                                              # Folder containing deduplicated bamfiles
│   │   │   ├── Aligned.sortedByCoord.out.bam                     # Sorted by coordinate
│   │   │   ├── Aligned.sortedByCoord.out.bam.bai                 # Indexed bam
│   │   │   ├── Aligned.sortedByName.out.bam                      # Sorted by name
│   │   │   ├── Log.final.out
│   │   │   ├── Log.out
│   │   │   ├── Log.progress.out
│   │   │   ├── SJ.out.tab
│   │   │   ├── _STARgenome
│   │   │   ├── _STARpass1
│   │   │   ├── Unmapped.out.mate1.gz                             # Unmapped reads of first fastq
│   │   │   └── Unmapped.out.mate2.gz                             # Unmapped reads of second fastq
│   │   ├── 27_htout
│   │   │   └── RNA003584L1_htseq_counts.txt                      # HTSeq counts (alternative to Kallisto)
│   │   ├── 28_idxout
│   │   │   └── RNA003584L1_idxstat.txt                           # Counts per chromosome/spike
│   │   ├── 30_coverage                                           # Directory with files for coverage analyses
│   │   │   ├── RNA003584L1_coverage_per_nt_all_strand.txt
│   │   │   ├── RNA003584L1_coverage_per_nt_minus_strand.txt
│   │   │   └── RNA003584L1_coverage_per_nt_plus_strand.txt
│   │   ├── FASTQC_qc_subs                                        # FASTQC of quality filtered + subsampled fastqs
│   │   │   ├── RNA003584L1-173868704_1_qc_subs_fastqc.html
│   │   │   ├── ...
│   │   ├── FASTQC_qc_subs_clumped                                # FASTQC of quality filtered + subs + dedup fastqs
│   │   │   ├── RNA003584L1-173868704_1_qc_subs_clumped_fastqc.html
│   │   │   ├── ...
│   │   ├── logs                                                  # Log directory: output and error file for each job
│   │   │   ├── 00_combinecopy_RNA003584L1.err
│   │   │   ├── 00_combinecopy_RNA003584L1.out
│   │   │   ├── ...
│   │   ├── RNA003584L1-173868704_1_qc_subs_clumped.fastq.gz      # Quality filtered + subsampled + clumpify duplicate removed fastq 1
│   │   ├── RNA003584L1-173868704_1_qc_subs.fastq.gz              # Quality filtered + subsampled fastq 1
│   │   ├── RNA003584L1-173868704_2_qc_subs_clumped.fastq.gz      # Quality filtered + subsampled + clumpify duplicate removed fastq 2
│   │   ├── RNA003584L1-173868704_2_qc_subs.fastq.gz              # Quality filtered + subsampled fastq 2
│   │   ├── RSeQC_output_all.txt                                  # RSeQC output for strandedness
│   │   ├── RSeQC_output.txt
│   │   └── scripts                                               # Scripts directory, contains slurm script for each step
│   │       ├── 00_combinecopy_RNA003584L1.sh                               # (can be useful if steps need to be rerun)
│   │       ├── 01_fastqc_RNA003584L1_orig.sh
│   │       ├── ...
│   ├── dedup-picard_subs-atstart                                 # Deduplication + subsampling options directory
│   │   ├── 23_srout
│   │   │   ├── Aligned.sortedByCoord.out.bam                     # bam file after mapping qc+subs fastq
│   │   │   ├── Aligned.sortedByCoord.picard.bam                  # Picard deduplicated bam file (sorted by coordinate)
│   │   │   ├── Aligned.sortedByCoord.picard.bam.bai              # Indexed bam
│   │   │   ├── Aligned.sortedByName.picard.bam                   # Picard deduplicated bam file (sorted by name)
│   │   │   ├── Log.final.out
│   │   │   ├── Log.out
│   │   │   ├── Log.progress.out
│   │   │   ├── SJ.out.tab
│   │   │   ├── _STARgenome
│   │   │   ├── _STARpass1
│   │   │   ├── Unmapped.out.mate1.gz                             # Unmapped reads of first fastq
│   │   │   └── Unmapped.out.mate2.gz                             # Unmapped reads of second fastq
│   │   ├── 24_picardout
│   │   │   └── picard_dup.metrics                                # Picard output file (contains the % of duplication)
│   │   ├── 25_klout                                              # Folder containing Kallisto abundance files 
│   │   │   ├── abundance.h5                                      # (based on deduplicated fastqs)
│   │   │   ├── abundance.tsv
│   │   │   └── run_info.json
│   │   ├── 26_htout                                                
│   │   │   └── RNA003584L1_htseq_counts.txt                      # HTSeq counts (alternative to Kallisto)
│   │   ├── 27_idxout
│   │   │   └── RNA003584L1_idxstat.txt                           # Counts per chromosome/spike
│   │   ├── 29_coverage                                           # Directory with files for coverage analyses
│   │   │   ├── RNA003584L1_coverage_per_nt_all_strand.txt
│   │   │   ├── RRNA003584L1_coverage_per_nt_minus_strand.txt
│   │   │   └── RNA003584L1_coverage_per_nt_plus_strand.txt
│   │   ├── FASTQC_qc_subs                                        # FASTQC of quality filtered + subsampled fastqs
│   │   │   ├── RNA003584L1-173868704_1_qc_subs_fastqc.html
│   │   │   ├── ...
│   │   ├── logs                                                  # Log directory: output and error file for each job
│   │   │   ├── 20_preprepeat_RNA003584L1.err
│   │   │   ├── 20_preprepeat_RNA003584L1.out
│   │   │   ├── ...
│   │   ├── RNA003584L1-173868704_1_picard.fastq.gz               # Fastq 1 after Picard dup removal (bam to fastq)
│   │   ├── RNA003584L1-173868704_1_qc_subs.fastq.gz              # Quality filtered + subsampled fastq 1
│   │   ├── RNA003584L1-173868704_2_picard.fastq.gz               # Fastq 2 after Picard dup removal (bam to fastq)
│   │   ├── RNA003584L1-173868704_1_qc_subs.fastq.gz              # Quality filtered + subsampled fastq 2
│   │   ├── RSeQC_output_all.txt
│   │   ├── RSeQC_output.txt
│   │   ├── scripts                                               # Scripts directory, contains slurm script for each step
│   │   │   ├── 20_preprepeat_RNA003584L1.sh                                # (can be useful if steps need to be rerun)
│   │   │   ├── ...
│   ├── FASTQC_original                                           # FASTQC of original fastqs
│   │   ├── RNA003584L1-173868704_1_fastqc.html
│   │   ├── ...
│   ├── FASTQC_qc                                                 # FASTQC of quality filtered fastqs
│   │   ├── RNA003584L1-173868704_1_qc_fastqc.html
│   │   ├── ...
│   ├── RNA003584L1-173868704_1.fastq.gz                          # Original fastq 1 (concatenated fastq of different lanes)
│   ├── RNA003584L1-173868704_1_qc.fastq.gz                       # Quality filtered fastq 1
│   ├── RNA003584L1-173868704_2.fastq.gz                          # Original fastq 2 (concatenated fastq of different lanes)
│   └── RNA003584L1-173868704_2_qc.fastq.gz                       # Quality filtered fastq 2

```
