# make sure you use python 3 (module load Python/3.6.6-intel-2018b)
# python RNASeq_preprocessing.py <other args>
# e.g. python RNASeq_preprocessing.py -m RNA0 -s yes -t pe -l 75 --dedup clumpify --subs start -n 20000000 --repeat no -b $VSC_DATA_VO/NSQ_run453-84636552_Run458-86026942 -o $VSC_SCRATCH_VO/test_pipeline -u example@ugent.be

import subprocess
import os
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Submit RNA Seq preprocessing jobs to the UGent HPC (slurm)')

# Read arguments
parser.add_argument('-t', nargs=1, choices=['se','pe'], required=True, help='Single end (se) or paired end (pe) data')
parser.add_argument('-b', nargs=1, required=True, help='Base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-o', nargs=1, required=True, help='Directory where output should be created', metavar='output_dir')
parser.add_argument('-m', nargs=1, required=True, help='String match to select sample folders in base directory e.g. RNA0', metavar='string_match')
parser.add_argument('-s', nargs=1, choices=['yes','no'], required=True, help='(reverse) stranded sequencing or not (unstranded)?')
parser.add_argument('-l', nargs=1, required=True, help='Normal read length (longer reads will be trimmed)')
parser.add_argument('--dedup', nargs=1, choices=['no', 'picard', 'clumpify'], required=True, help='Method of duplicate removal [no, clumpify, picard]')
parser.add_argument('--pico', nargs=1, choices=['yes','no'], default=['no'], help='PicoV2 sequencing? [default: no]')
parser.add_argument('--polyA', nargs=1, choices=['yes','no'], default=['no'], help='PolyA trimming needed? [default: no]')
parser.add_argument('--adapter', nargs=1, choices=['truseq','nextera','no'], default=['no'], help='Adapter trimming needed? If so, which one [default: no]')
parser.add_argument('--subs', nargs=1, choices=['no','start','afterclumpify'], default=['no'], help='Subsampling needed? [default: no] If so, when? At start or after clumpify duplicate removal')
#parser.add_argument('--seq', nargs=1, choices=['novaseq','hiseq','nextseq'], default=['nextseq'], help='Sequencing machine? [default: nextseq] Is important for duplicate removal', metavar=seq_machine)
parser.add_argument('-n', nargs=1, default=['0'], help='Nr of reads to subsample to. If you need the analysis to stop in order to determine the subsampling level first, do not enter a number here (or -n 0)')
parser.add_argument('--coverage', nargs=1, choices=['yes','no'], default=['no'], help='Coverage analysis needed? [default: no]')
parser.add_argument('--fragmentmean', nargs=1, default=['-1'], help='Mean fragment length (only necessary for kallisto in SE mode)')
parser.add_argument('--fragmentsd', nargs=1, default=['-1'], help='Mean fragment length (only necessary for kallisto in SE mode)')
parser.add_argument('--repeat', nargs=1, choices=['yes','no'], default=['no'], help='Repeated analysis? [default: no] If yes: general fastq copying and filtering will not be repeated')
parser.add_argument('-u', nargs=1, required=True, help='Submitter email address (used for error reporting)', metavar='user_email')

# Parse arguments
args = parser.parse_args()
print(args)
data_type = args.t[0]
base_dir = args.b[0].rstrip("/")
string_match = args.m[0]
email = args.u[0]
deduplication = args.dedup[0]
picoV2 = args.pico[0]
stranded = args.s[0]
readlength = args.l[0]
clumpifylength = str(int(0.8*int(readlength))) #use 80% of read (do not take into account last 20%, more chance on sequencing errors)
adaptercontam = args.adapter[0]
polyAtrim = args.polyA[0]
#print(picoV2)
subs_timing = args.subs[0]
coverage = args.coverage[0]
repeat_analysis = args.repeat[0]
fragmentmean = args.fragmentmean[0]
fragmentsd = args.fragmentsd[0]
subsample_to_nr = args.n[0] # nr of reads you want to subsample to (based on floor of wc -l divided by 4), will be ignored if no subsampling asked for in command line
#seq_machine = args.seq[0]

kal_index = "Kallisto_index_hg38_91_withspikes_rDNA45S"
#star_index = "Genome_hg38_spikes_MTr45S_star_index"
star_index = "Genome_hg38_spikes_chrIS_MTr45S_star_index"
#gtf = "Homo_sapiens.GRCh38.91.gtf"
gtf = "Homo_sapiens.GRCh38.91_withspikes.gtf"
chrom = "Homo_sapiens.GRCh38.chr_spikes.txt"
intergenic_bed="Homo_sapiens.GRCh38.91_intergenic_stranded.bed"
intron_bed="Homo_sapiens.GRCh38.91_introns.bed"
exon_bed = "Homo_sapiens.GRCh38.91_exons_sorted_merged2.bed"
exonminus_bed = "Ensembl_GRCh38_91_exons_minus.bed"
exonplus_bed = "Ensembl_GRCh38_91_exons_plus.bed"
tasks = 1
mem = 40


def sbatch(job_name, command, index, mem = mem, tasks = tasks, workdir = '.', time='4:00:00', dep=''): 
	if dep != '':
		dep = ' --dependency=afterok:{} --kill-on-invalid-dep=yes'.format(dep) 
	printlines = [
		"#!/bin/bash",
		"",
		"#SBATCH -J {}".format(job_name+'_'+dedup_mode),
		"#SBATCH -D {}".format(workdir),
		"#SBATCH --mem={}G".format(mem),
		"#SBATCH --cpus-per-task={}".format(tasks), #nr of processors per task needed
		"#SBATCH -t {}".format(time),
		"#SBATCH --mail-user={}".format(email),
		"#SBATCH --mail-type=FAIL",
		"#SBATCH -o {1}/{2}/logs/{3:02d}_{0}.out".format(job_name, workdir, dedup_mode, index),
		"#SBATCH -e {1}/{2}/logs/{3:02d}_{0}.err".format(job_name, workdir, dedup_mode, index),
		#"#SBATCH --test-only", #Uncomment this to test instead of directly submitting job
		""
	]
	
	printlines.extend(command.split("; "))
	#print(printlines)
	printjobfile('{0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, dedup_mode, index, job_name),printlines)
	#print(job_name)
	sbatch_command = 'sbatch{4} {0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, dedup_mode, index, job_name, dep)
	#print(sbatch_command)
	sbatch_response = subprocess.getoutput(sbatch_command) 
	#print(sbatch_response)
	job_id = sbatch_response.split(' ')[3].strip() 
	#job_id = job_name
	#print(job_id)
	return job_id

def printjobfile(filename, printlines):
        with open(filename, 'w') as the_file:
                for line in printlines:
                        the_file.write(line+"\n")

### Copy from DATA to SCRATCH
def combinecopy(sampleID,dep='',jobsuffix=''):
	""" Copy from original directory to working directory + unzip
	Check if the fastq files per lane are already combined ([[ -s file ]]: True if file has a Size greater than zero),
	if not concatenate them
	"""
	if data_type == 'se':
		command = '[[ -s {0}/{2}.fastq.gz ]] && cp {0}/{2}.fastq.gz {1}/{2}.fastq.gz || cat {0}/*R1_*.fastq.gz > {1}/{2}.fastq.gz; '.format(input_dir, output_dir, sampleID)
	elif data_type == 'pe':
		command = '[[ -s {0}/{2}_1.fastq.gz ]] && cp {0}/{2}_1.fastq.gz {1}/{2}_1.fastq.gz || cat {0}/*R1_*.fastq.gz > {1}/{2}_1.fastq.gz; '.format(input_dir, output_dir, sampleID)
		command = command + '[[ -s {0}/{2}_2.fastq.gz ]] && cp {0}/{2}_2.fastq.gz {1}/{2}_2.fastq.gz || cat {0}/*R2_*.fastq.gz > {1}/{2}_2.fastq.gz; '.format(input_dir, output_dir, sampleID)
	command = command + 'gunzip {0}/{1}*.fastq.gz; '.format(output_dir, sampleID)
	job_id = sbatch('combinecopy_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)
	return job_id

### If not first analysis: unzip the fastq files for further use (if they are gzipped)
def preprepeat(sampleID,dep='',jobsuffix=''):
	""" If the quality filtered files are already present, unzip them (if still needed) for further analysis
	Test if the file is gzipped [[ -s file.gz ]]: True if file has a Size greater than zero -> gunzip
	Else: test whether the unzipped file exists (if it doesn't, the function will throw an error)
	"""
	if data_type == 'se':
		command = '[[ -s {0}/{1}_qc.fastq.gz ]] && gunzip {0}/{1}_qc.fastq.gz || head -1 {0}/{1}_qc.fastq; '.format(output_dir, sampleID)
	elif data_type == 'pe':
		command = '[[ -s {0}/{1}_1_qc.fastq.gz ]] && gunzip {0}/{1}_1_qc.fastq.gz || head -1 {0}/{1}_1_qc.fastq; '.format(output_dir, sampleID)
		command = command + '[[ -s {0}/{1}_2_qc.fastq.gz ]] && gunzip {0}/{1}_2_qc.fastq.gz || head -1 {0}/{1}_2_qc.fastq; '.format(output_dir, sampleID)
	job_id = sbatch('preprepeat_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep)
	return job_id	
	
### FastQC
def fastqc(sampleID,suffix='',outDIR='.',subDIR='',dep='',jobsuffix=''): 
	# Build the command for fastqc 
	command = 'module purge; module load FastQC/0.11.8-Java-1.8; mkdir {}/{}{}; '.format(output_dir, subDIR, outDIR)
	if data_type == 'se':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}{4}.fastq; '.format(output_dir, subDIR, outDIR, sampleID, suffix)
	elif data_type == 'pe':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_1{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_2{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
	job_id = sbatch('fastqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)
	return job_id

def trimlength(sampleID,suffix='',dep='',jobsuffix=''):
	""" Remove trailing bases from reads
	--minimum-length 20: discard reads shorter than 20
	-l: shorten each read down to a certain length, use the --length option or the short version -l (we want to trim last nt off 76nt reads)
	-q x,y: quality-trim read with a threshold of x from 3' and with threshold of y from 5'
	-o & -p: output file for first read (-o) and second read (-p) of pair
	--pair-filter=any: (default) read pair is discarded if one of the reads (R1 or R2) fulfills the filtering criterion (e.g. one read is shorter than 20)
	"""
	command = 'module purge; module load cutadapt/1.18-intel-2018b-Python-3.6.6; '
	#command = command + 'mkdir {0}/{1}/{2:02d}_cutadaptout; '.format(output_dir, dedup_mode, index)
	if data_type == 'se':
		command = command + 'cutadapt -l {3} --minimum-length=20 -o {0}/{1}_len.fastq {0}/{1}{2}.fastq; '.format(output_dir, sampleID, suffix, readlength)
		command = command + 'gzip {0}/{1}{2}.fastq; '.format(output_dir, sampleID, suffix)
	elif data_type == 'pe':
		command = command + 'cutadapt -l {3} --minimum-length=20 --pair-filter=any -o {0}/{1}_1{2}_len.fastq -p {0}/{1}_2{2}_len.fastq {0}/{1}_1{2}.fastq {0}/{1}_2{2}.fastq; '.format(output_dir, sampleID, suffix, readlength)
		command = command + 'gzip {0}/{1}_1{2}.fastq; gzip {0}/{1}_2{2}.fastq; '.format(output_dir, sampleID, suffix)
	job_id = sbatch('trimlength_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep,tasks=4)
	return job_id

def removebadqc(sampleID,suffix='',percentage=str(80),quality=str(19),dep='',jobsuffix=''):
	"""
	Keep only read pairs where each read has at least p% of bases with a phred score > q
	-q: minimum quality phred score
	-p: percentage that has the minimum quality phred score
	FASTQC_filter.py: Keep only read pairs where each read has at least p% of bases with a phred score > q
	fastq_quality_filter: remove reads with lower quality, keeping only reads that have at least p% of bases with a quality score of q or more
	-i: FASTA/Q input file
	-o: FASTA/Q output file
	-Q33: If Illumina encoding is >= 1.8, you need to provide this option to fastq_quality_filter. For Illumina encoding <=1.5, option is not required
	
	"""
	command = 'module purge; module load Biopython/1.72-foss-2018b-Python-3.6.6; '
	#command = command + 'mkdir {0}/{1}/{2:02d}_removebadqcout; '.format(output_dir, dedup_mode, index)
	if data_type == 'se':
		command = 'module purge; module load FASTX-Toolkit/0.0.14-intel-2018a; '
		command = command + 'fastq_quality_filter -p {3} -q {4} -Q33 -i {0}/{1}{2}.fastq -o {0}/{1}_qc.fastq; '.format(output_dir, sampleID, suffix, percentage, str(int(quality)+1))
		command = command + 'qc_lines=`wc -l {0}/{1}_qc.fastq`; echo {1} $qc_lines > {0}/qc_lines.txt; '.format(output_dir,sampleID)
	elif data_type == 'pe':
		command = command + 'split -l 10000000 --additional-suffix=.fastq {0}/{1}_1{2}.fastq {0}/temp_1_; '.format(output_dir, sampleID, suffix)
		command = command + 'rm {0}/{1}_1{2}.fastq; '.format(output_dir, sampleID, suffix)
		command = command + 'split -l 10000000 --additional-suffix=.fastq {0}/{1}_2{2}.fastq {0}/temp_2_; '.format(output_dir, sampleID, suffix)
		command = command + 'rm {0}/{1}_2{2}.fastq; '.format(output_dir, sampleID, suffix)
		command = command + 'for suffix in $(find {0}/. -type f -name "*temp_1*" | awk -F \'temp_1_\' \'{{print $2}}\' | cut -d \'.\' -f1); '.format(output_dir, index)
		command = command + 'do; '
		command = command + '   python /data/gent/vo/000/gvo00027/RNA_seq_pipeline/resources_Annelien/FASTQC_filter.py -r1 {0}/temp_1_$suffix.fastq -r2 {0}/temp_2_$suffix.fastq -p {1} -q {2} -o1 {0}/temp_qc_1_$suffix.fastq -o2 {0}/temp_qc_2_$suffix.fastq; '.format(output_dir, percentage, quality)
		command = command + '   rm {0}/temp_1_$suffix.fastq;    rm {0}/temp_2_$suffix.fastq; '.format(output_dir)
		command = command + 'done; '
		command = command + 'cat {0}/temp_qc_1_* > {0}/{1}_1_qc.fastq; rm {0}/temp_qc_1*; '.format(output_dir, sampleID)
		command = command + 'cat {0}/temp_qc_2_* > {0}/{1}_2_qc.fastq; rm {0}/temp_qc_2*; '.format(output_dir, sampleID)
		command = command + 'qc_lines=`wc -l {0}/{1}_1_qc.fastq`; echo {1} $qc_lines > {0}/qc_lines.txt; '.format(output_dir, sampleID)
	job_id = sbatch('removebadqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep, tasks=4)
	return job_id

def trimadapter(sampleID,suffix='',adapterR1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',adapterR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',dep='',jobsuffix=''):
	""" Remove adapters and trim reads
	-a: regular 3' adapter R1 (To reduce the number of falsely trimmed bases, the alignment algorithm requires that, by default, at least three bases match between adapter and read)
	-A: regular 3' adapter R2
	--minimum-length 20: discard reads shorter than 20
	-U: remove x bases from beginning (positive x) or end (negative x) of R2 read (-U 3 needed for picoV2 data)
	-l: shorten each read down to a certain length, use the --length option or the short version -l (we want to trim last nt off 76nt reads)
	-q x,y: quality-trim read with a threshold of x from 3' and with threshold of y from 5'
	-o & -p: output file for first read (-o) and second read (-p) of pair
	--pair-filter=any: (default) read pair is discarded if one of the reads (R1 or R2) fulfills the filtering criterion (e.g. one read is shorter than 20)
	"""
	command = 'module purge; module load cutadapt/1.18-intel-2018b-Python-3.6.6; '
	command = command + 'mkdir {0}/{1:02d}_trimout; '.format(output_dir, index)
	### Normal adapter trimming or additional removal of 3 nt with PicoV2?
	if picoV2=='yes':
		if data_type == 'se':
			command = command + 'cutadapt -a {4} --minimum-length=20 -U 3 -o {0}/{1:02d}_trimout/{2}_adapter.fastq {0}/{2}{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
		elif data_type == 'pe':
			command = command + 'cutadapt --pair-filter=any -a {4} -A {5} --minimum-length=20 -U 3 -o {0}/{1:02d}_trimout/{2}_1_adapter.fastq -p {0}/{1:02d}_trimout/{2}_2_adapter.fastq {0}/{2}_1{3}.fastq {0}/{2}_2{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
	else: #no picoV2 trimming necessary
		if data_type == 'se':
			command = command + 'cutadapt -a {4} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_adapter.fastq {0}/{2}{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
		elif data_type == 'pe':
			command = command + 'cutadapt --pair-filter=any -a {4} -A {5} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_1_adapter.fastq -p {0}/{1:02d}_trimout/{2}_2_adapter.fastq {0}/{2}_1{3}.fastq {0}/{2}_2{3}.fastq; '.format(output_dir,index,sampleID,suffix,adapterR1,adapterR2)
	### PolyA adapter trimming necessary?
	if polyAtrim=='yes':
		if data_type == 'se':
			command = command + 'cutadapt -a T{{100}} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_adapterpA.fastq {0}/{1:02d}_trimout/{2}_adapter.fastq; '.format(output_dir,index,sampleID) #remove polyA tails (adapter sequence is considered to be a sequence of 100A for the R2 read and a sequence of 100T for the R1 read)
			command = command + 'mv {0}/{1:02d}_trimout/{2}_adapterpA.fastq {0}/{2}{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
		elif data_type == 'pe':
			command = command + 'cutadapt --pair-filter=any -a T{{100}} -A A{{100}} --minimum-length=20 -o {0}/{1:02d}_trimout/{2}_1_adapterpA.fastq -p {0}/{1:02d}_trimout/{2}_2_adapterpA.fastq {0}/{1:02d}_trimout/{2}_1_adapter.fastq {0}/{1:02d}_trimout/{2}_2_adapter.fastq; '.format(output_dir,index,sampleID) #remove polyA tails (adapter sequence is considered to be a sequence of 100A for the R2 read and a sequence of 100T for the R1 read)
			command = command + 'mv {0}/{1:02d}_trimout/{2}_1_adapterpA.fastq {0}/{2}_1{3}_trim.fastq; mv {0}/{1:02d}_trimout/{2}_2_adapterpA.fastq {0}/{2}_2{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
	else: #no polyA trimmming necessary
		if data_type == 'se':
			command = command + 'mv {0}/{1:02d}_trimout/{2}_adapter.fastq {0}/{2}{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
		elif data_type == 'pe':
			command = command + 'mv {0}/{1:02d}_trimout/{2}_1_adapter.fastq {0}/{2}_1{3}_trim.fastq; mv {0}/{1:02d}_trimout/{2}_2_adapter.fastq {0}/{2}_2{3}_trim.fastq; '.format(output_dir,index,sampleID,suffix) #give them back original name, to be able to run sequel of pipeline just like for capture data
	command = command + 'rm -r {0}/{1:02d}_trimout; '.format(output_dir, index) #rm temp dir
	job_id = sbatch('trimadapter_'+shortname+jobsuffix, command, index,workdir=output_dir,dep=dep,tasks=4)
	return job_id

### Subsampling
def subsample(sampleID,nreads,suffix='',subDIR='',dep='',jobsuffix=''):
	""" Seqtk subsampling
	Downsample fastq files to x number of reads (take min of all samples)
	-s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of code
	"""
	command = 'module purge; module load seqtk/1.3-foss-2018b; '
	#if subsample_to_nr == '0': #if no subsampling nr given
	#	command = command + 'echo \"Determine the subsampling level first! (i.e. floored minimum of wc-l fastq divided by 4)\"; exit 1; '
	if data_type == 'se':
		command = command + 'seqtk sample -s100 {0}/{1}{2}{3}.fastq {4} > {0}/{5}/{2}{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
	elif data_type == 'pe':
		command = command + 'seqtk sample -s100 {0}/{1}{2}_1{3}.fastq {4} > {0}/{5}/{2}_1{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
		command = command + 'seqtk sample -s100 {0}/{1}{2}_2{3}.fastq {4} > {0}/{5}/{2}_2{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,dedup_mode)
	job_id = sbatch('subsample_'+shortname+jobsuffix,command,index,workdir=output_dir,mem=40, time='8:00:00', dep=dep,tasks=2)
	return job_id

def clumpifydupremoval_trimmed(sampleID, suffix='', subDIR='',subst='2',kmersize='31', passes='20',dep='',jobsuffix=''):
	""" Duplicate removal with clumpify
	Cutadapt: trim every read to 80% of its length (-l) and gzip file for processing with Clumpify
	Clumpify:
	Needs gzipped input files and can be run in PE and SE mode (for SE: no in2=, out2=)
	dedupe: Remove duplicate reads. For pairs, both must match. By default, deduplication does not occur (dedupe=f)
	(default allowns=t: No-called bases will not be considered substitutions)
	Duplicates are removed within clumps (that are based on reads sharing a specific kmer)
	subs: number of substitutions allowed (if subs=0, only exact matches are considered as duplicates, default 2, note that this allows 2 substitutions in both ends of the pair)
	subrate=0.0: If set, the number of substitutions allowed will be max(subs, subrate*min(length1, length2)) for 2 sequences.
	k: kmer size, reads that share a kmer of given length end in same clump (default 31), lower k will result in more reads per clump and thus more chance to eliminate duplicates
	passes: number of times a different kmer is selected for seeding clumps (default 1), eventually, any pair of duplicates will land in same clump given enough passes if they share a single kmer
	For plasma data (with a lot of duplicates), plateau is reached after about 10 passes but additional passes do not take much longer -> to be safe, run with 15 to 20 passes
	dedupedist=40: (dist) Max distance to consider for optical duplicates. Higher removes more duplicates but is more likely to remove PCR rather than optical duplicates. This is platform-specific;
		recommendations: NextSeq = 40 (and spany=t); HiSeq 1T = 40; HiSeq 2500 = 40; HiSeq 3k/4k = 2500; Novaseq = 12000
	spany=f: Allow reads to be considered optical duplicates if they are on different tiles, but are within dupedist in the y-axis.  Should only be enabled when looking for tile-edge duplicates (as in NextSeq).
	"""
	command = 'module purge; module load BBMap/38.26-foss-2018b; module load cutadapt/1.18-foss-2018b-Python-3.6.6; '
	command = command + 'mkdir {0}/{1}/{2:02d}_clumpout; '.format(output_dir, dedup_mode, index)
	if data_type == 'se':
		command = command + 'cutadapt -l {6} -o {0}/{1}/{2:02d}_clumpout/{3}_temptrim.fastq {0}/{5}{3}{4}.fastq; '.format(output_dir, dedup_mode, index, sampleID, suffix,subDIR, clumpifylength)
		command = command + 'gzip {0}/{1}/{2:02d}_clumpout/{3}_temptrim.fastq; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'clumpify.sh in={0}/{1}/{2:02d}_clumpout/{3}_temptrim.fastq.gz out={0}/{1}/{2:02d}_clumpout/{3}_tempclumped.fastq.gz dedupe subs={4} k={5} passes={6}; '.format(output_dir, dedup_mode, index, sampleID, subst, kmersize, passes)
		command = command + 'gunzip {0}/{1}/{2:02d}_clumpout/{3}_tempclumped.fastq.gz; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'grep \'^@\' {0}/{1}/{2:02d}_clumpout/{3}_tempclumped.fastq > {0}/{1}/{2:02d}_clumpout/{3}_tempnames.txt; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'grep -Ff {0}/{1}/{2:02d}_clumpout/{3}_tempnames.txt -A3 {0}/{5}{3}{4}.fastq | grep -v -- \"^--$\" > {0}/{1}/{3}{4}_clumped.fastq; '.format(output_dir, dedup_mode, index, sampleID, suffix,subDIR)
		command = command + 'rm -r {0}/{1}/{2:02d}_clumpout; '.format(output_dir, dedup_mode, index)
	elif data_type == 'pe':
		command = command + 'cutadapt --pair-filter=both -l {6} -o {0}/{1}/{2:02d}_clumpout/{3}_temptrim_1.fastq -p {0}/{1}/{2:02d}_clumpout/{3}_temptrim_2.fastq {0}/{5}{3}_1{4}.fastq {0}/{5}{3}_2{4}.fastq; '.format(output_dir, dedup_mode, index, sampleID, suffix,subDIR, clumpifylength)
		command = command + 'gzip {0}/{1}/{2:02d}_clumpout/{3}_temptrim*.fastq; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'clumpify.sh in={0}/{1}/{2:02d}_clumpout/{3}_temptrim_1.fastq.gz in2={0}/{1}/{2:02d}_clumpout/{3}_temptrim_2.fastq.gz out={0}/{1}/{2:02d}_clumpout/{3}_tempclumped_1.fastq.gz out2={0}/{1}/{2:02d}_clumpout/{3}_tempclumped_2.fastq.gz dedupe subs={4} k={5} passes={6}; '.format(output_dir, dedup_mode, index, sampleID, subst, kmersize, passes)
		command = command + 'gunzip {0}/{1}/{2:02d}_clumpout/{3}_tempclumped*.fastq.gz; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'grep \'^@\' {0}/{1}/{2:02d}_clumpout/{3}_tempclumped_1.fastq > {0}/{1}/{2:02d}_clumpout/{3}_tempnames_1.txt; grep \'^@\' {0}/{1}/{2:02d}_clumpout/{3}_tempclumped_2.fastq > {0}/{1}/{2:02d}_clumpout/{3}_tempnames_2.txt; '.format(output_dir, dedup_mode, index, sampleID)
		command = command + 'grep -Ff {0}/{1}/{2:02d}_clumpout/{3}_tempnames_1.txt -A3 {0}/{5}{3}_1{4}.fastq | grep -v -- \"^--$\" > {0}/{1}/{3}_1{4}_clumped.fastq; '.format(output_dir, dedup_mode, index, sampleID, suffix,subDIR)
		command = command + 'grep -Ff {0}/{1}/{2:02d}_clumpout/{3}_tempnames_2.txt -A3 {0}/{5}{3}_2{4}.fastq | grep -v -- \"^--$\" > {0}/{1}/{3}_2{4}_clumped.fastq; '.format(output_dir, dedup_mode, index, sampleID, suffix,subDIR)
		command = command + 'rm -r {0}/{1}/{2:02d}_clumpout; '.format(output_dir, dedup_mode, index)
	#command = command + 'gzip {0}/{1}/{3}*_clumped.fastq; '.format(output_dir, dedup_mode, index, sampleID)	
	job_id = sbatch('removedupClumpifytrim_'+shortname+jobsuffix,command,index,workdir=output_dir,mem=40,time='8:00:00',dep=dep, tasks=4)
	return job_id

### Duplicate removal
def picarddupremcoord(sampleID,bamIN='Aligned.sortedByCoord.out.bam',bamOUT='Aligned.sortedByCoord.picard.bam',subdir='.',dep='',jobsuffix=''):
	""" Duplicate removal with Picard
	(When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates)
	ASSUME_SORT_ORDER=coordinate
	REMOVE_DUPLICATES=true: remove all duplicates (optical and sequencing)
	VALIDATION_STRINGENCY=SILENT: (default: STRICT) Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags)
	M: file to write duplication metrics to
	"""
	command = 'module purge; module load picard/2.18.5-Java-1.8.0_162; '
	command = command + 'mkdir {0}/{1}/{2:02d}_picardout; '.format(output_dir, dedup_mode, index)
	command = command + 'java -jar ${{EBROOTPICARD}}/picard.jar MarkDuplicates I={0}/{1}/{2:02d}_srout/{3} O={0}/{1}/{5:02d}_picardout/{4} ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M={0}/{1}/{5:02d}_picardout/{6}_picard_dup.metrics; '.format(output_dir,dedup_mode, alignindex,bamIN,bamOUT,index, sampleID)
	#'rm {}_srout/{}'.format(sampleID,bamOUT) #not needed anymore, just interesting to see how many duplicates picard still finds (picard_dup.metrics)
	command = command + 'mv {0}/{1}/{2:02d}_picardout/{3} {0}/{1}/{4:02d}_srout/{3}; '.format(output_dir, dedup_mode, index,bamOUT,alignindex)
	job_id = sbatch('removedupPicard_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep, time='8:00:00', tasks=4)
	return job_id

def bamtofastqfile2(sampleID,bamIN,dep='',jobsuffix=''):
	""" Convert bam file to fastq
	-fq: output file (SE, or mate 1 for PE); -fq2: output file 2 (mate2 for PE)
	Read pairs should be together for this conversion (if not the case, sortbambyname first)
	Issue: if multimapping is allowed, one read with 10 alignments can be converted to 10 reads!
	With AllBestScore option of STAR, they are all marked as primary -> better to select unique reads (name-based) after conversion (awk)
	"""
	command = 'module purge; module load BEDTools/2.27.1-intel-2018a; '
	if data_type == 'se':
		command = command + 'bedtools bamtofastq -i {0}/{1}/{2:02d}_srout/{3} -fq {0}/{1}/tmp1.fastq; '.format(output_dir, dedup_mode, alignindex, bamIN)
		command = command + 'paste - - - - < {0}/{1}/tmp1.fastq | awk \'!_[$1]++\' | tr \'\t\' \'\n\' > {0}/{1}/{2}_picard.fastq; '.format(output_dir, dedup_mode, sampleID)
		command = command + 'rm {0}/{1}/tmp1.fastq; '.format(output_dir,dedup_mode)
	elif data_type == 'pe':
		command = command + 'bedtools bamtofastq -i {0}/{1}/{2:02d}_srout/{3} -fq {0}/{1}/tmp1.fastq -fq2 {0}/{1}/tmp2.fastq; '.format(output_dir, dedup_mode, alignindex, bamIN)
		command = command + 'paste - - - - < {0}/{1}/tmp1.fastq | awk \'!_[$1]++\' | tr \'\t\' \'\n\' > {0}/{1}/{2}_1_picard.fastq; '.format(output_dir, dedup_mode, sampleID)
		command = command + 'paste - - - - < {0}/{1}/tmp2.fastq | awk \'!_[$1]++\' | tr \'\t\' \'\n\' > {0}/{1}/{2}_2_picard.fastq; '.format(output_dir, dedup_mode, sampleID)
		command = command + 'rm {0}/{1}/tmp1.fastq; rm {0}/{1}/tmp2.fastq; '.format(output_dir,dedup_mode)
	job_id = sbatch('bamtofastqfile2_'+shortname+jobsuffix, command, index, workdir=output_dir, dep=dep)
	return job_id

def bamtofastqfile(sampleID,bamIN,dep='',jobsuffix=''):
	""" Convert bam or sam file to fastq with samtools fastq
	-1: write reads with READ1 FLAG set (and not READ2) to file (SE, or mate 1 for PE); -2: write reads with READ2 FLAG set (and not READ1) to file 2 (mate 2 for PE)
	-0 /dev/null: write reads where both READ1 and READ2 are set; or neither is set to file (/dev/null = remove them) (reads that do not fit in simple paired-end sequencing model)
	-s /dev/null: write singleton reads to file (/dev/null = remove them)
	-n: use read names as they are (instead of adding '/1' and '/2')
	-F 0x900: Do not output alignments with any bits set in INT present in the FLAG field
	0x900 = 0x100 (not primary alignment) & 0x800 (supplementary alignment)
	"""
	command = 'module purge; module load SAMtools/1.8-intel-2018a; '
	if data_type == 'se':
		command = command + 'samtools fastq -1 {0}/{1}/{2}_picard.fastq -0 /dev/null -n -F 0x900 {0}/{1}/{3:02d}_srout/{4}; '.format(output_dir, dedup_mode, sampleID, alignindex, bamIN)
	elif data_type == 'pe':
		command = command + 'samtools fastq -1 {0}/{1}/{2}_1_picard.fastq -2 {0}/{1}/{2}_2_picard.fastq -0 /dev/null -n -F 0x900 {0}/{1}/{3:02d}_srout/{4}; '.format(output_dir, dedup_mode, sampleID, alignindex, bamIN)
	job_id = sbatch('bamtofastqfile_'+shortname+jobsuffix, command, index, workdir=output_dir, dep=dep)
	return job_id

### Mapping and counting (STAR/HTSeq/Kalllisto)
def alignment(sampleID,suffix='',subDIR='',dep='', jobsuffix=''): 
	""" Paired end STAR alignment with output as bam (coordinate sorted)
	--outFileNamePrefix *: output folder
	--outSAMtype BAM SortedByCoordinate: outputs BAM file that is sorted by coordinate ("Aligned.sortedByCoord.out.bam")
	--outReadsUnmapped Fastx: unmapped reads output in separate fastq files (Unmapped.out.mate1/2)
	--twopassMode Basic: run STAR 2-pass mapping for each sample separately
	--outMultimapperOrder Random: outputs multiple alignments for each read in random order, and also also randomizes the choice of the primary alignment from the highest scoring alignments
	--outSAMmultNmax -1: max number of multiple alignments for a read that will be output to the SAM/BAM files (default -1: all alignments (up to outFilterMultimapNmax) will be output)
	--outFilterMultimapNmax 10: (default 10) max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped (counted as "mapped to too many loci")
	--outSAMprimaryFlag AllBestScore: For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in the FLAG (column 2 of the SAM).
	The unmarked alignment is selected from the best ones (i.e. highest scoring). outSAMprimaryFlag AllBestScore option will change default behavior and output all alignments with the best score as primary alignments (i.e. 0x100 bit in the FLAG unset)
	--outFilterScoreMinOverLread 0.66: alignment will be output only if its score normalized over read length (sum of mate lengths for PE reads) is higher than or equal to this value (default 0.66)
	--outFilterMatchNminOverLread 0.66: alignment will be output only if the number of matched bases normalized over read length (sum of mate lengths for PE reads) is higher than or equal to this value (default 0.66)
	--outFilterMatchNmin 20: alignment will be output only if the number of matched bases is higher than or equal to this value (default 0)
	"""
	command = 'module purge; module load STAR/2.6.0c-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_srout; '.format(output_dir, dedup_mode, index)
	#command = command + 'if [ {0} = . ] ; then echo {0}; else mkdir {0}; fi; '.format(subdir) #if subdir is specified, make directory for it (if it does not already exist)
	#'STAR --runThreadN 10 --outFileNamePrefix {}_srout/ --readFilesIn {}.fastq {}.fastq --genomeDir {} --sjdbGTFfile {} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag OneBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20'.format(sampleID,sampleID1,sampleID2,star_index,gtf)
	if data_type == 'se':
		command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/{3:02d}_srout/{4}. --readFilesIn {1}/{7}{4}{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,sampleID,suffix,gtf,subDIR)
	elif data_type == 'pe':
		command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/{3:02d}_srout/{4}. --readFilesIn {1}/{7}{4}_1{5}.fastq {1}/{7}{4}_2{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,sampleID,suffix,gtf,subDIR)
	command = command + 'gzip {0}/{1}/{2:02d}_srout/*Unmapped.*; '.format(output_dir, dedup_mode, index)
	job_id = sbatch('alignSTAR_'+shortname+jobsuffix, command, index, workdir=output_dir, time='8:00:00', mem=40,dep=dep,tasks=8) 
	return job_id

def countskallisto(sampleID,suffix='',subDIR='',dep='', jobsuffix=''):
	""" kallisto quant to run quantification algorithm (paired end mode)
	-t: number of threads to use (default 1)
	-i: index (used for quantification)
	--rf-stranded: strand specific mode: only fragments where the first read in the pair pseudoaligns to the reverse strand of a transcript are processed (vice versa in --fr-stranded)
	-o: output directory
	--single: SE sequencing (fragment length mean and sd must also be supplied for single-end reads using -l and -s)
	Produces 3 output files by default: abundances.tsv (plaintext file of abundance estimates), abundances.h5, run_info.json
	"""
	command = 'module purge; module load kallisto/0.44.0-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_klout; '.format(output_dir, dedup_mode, index)
	if stranded == 'yes':
		if data_type == 'se':
			command = command + 'kallisto quant -t 10 -i {0} --rf-stranded --single -l {7} -s {8} -o {1}/{2}/{3:02d}_klout/ {1}/{6}{4}{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR, fragmentmean, fragmentsd)
		elif data_type == 'pe':
			command = command + 'kallisto quant -t 10 -i {0} --rf-stranded -o {1}/{2}/{3:02d}_klout/ {1}/{6}{4}_1{5}.fastq {1}/{6}{4}_2{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR)
	else: #unstranded
		if data_type == 'se':
			command = command + 'kallisto quant -t 10 -i {0} --single -l {7} -s {8} -o {1}/{2}/{3:02d}_klout/ {1}/{6}{4}{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR, fragmentmean, fragmentsd)
		elif data_type == 'pe':
			command = command + 'kallisto quant -t 10 -i {0} -o {1}/{2}/{3:02d}_klout/ {1}/{6}{4}_1{5}.fastq {1}/{6}{4}_2{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR)
	job_id = sbatch('countskallisto_'+shortname+jobsuffix,command, index, workdir=output_dir, dep=dep,tasks=2)
	return job_id

def countskallisto_withbam(sampleID,suffix='',subDIR='',dep='', jobsuffix=''):
	""" kallisto quant to run quantification algorithm (paired end mode)
	-t: number of threads to use (default 1)
	-i: index (used for quantification)
	--rf-stranded: strand specific mode: only fragments where the first read in the pair pseudoaligns to the reverse strand of a transcript are processed (vice versa in --fr-stranded)
	-o: output directory
	--genomebam: visualize alignments by outputting this in bam file (this option needs a gtf and chromosomes file): output is pseudoalignments.bam and pseudoalignments.bam.bai
	Produces 3 other output files by default: abundances.tsv (plaintext file of abundance estimates), abundances.h5, run_info.json
	"""
	command = 'module purge; module load kallisto/0.44.0-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_klout; '.format(output_dir, dedup_mode, index)
	if stranded == 'yes':
		if data_type == 'se':
			command = command + 'kallisto quant -t 10 -i {0} --rf-stranded --single -l {8} -s {9} -o {1}/{2}/{3:02d}_klout/ --genomebam --gtf {7} --chromosomes {8} {1}/{6}{4}{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR,gtf,chrom, fragmentmean, fragmentsd)
		elif data_type == 'pe':
			command = command + 'kallisto quant -t 10 -i {0} --rf-stranded -o {1}/{2}/{3:02d}_klout/ --genomebam --gtf {7} --chromosomes {8} {1}/{6}{4}_1{5}.fastq {1}/{6}{4}_2{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR,gtf,chrom)
	else: #unstranded
		if data_type == 'se':
			command = command + 'kallisto quant -t 10 -i {0} --single -l {8} -s {9} -o {1}/{2}/{3:02d}_klout/ --genomebam --gtf {7} --chromosomes {8} {1}/{6}{4}{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR,gtf,chrom, fragmentmean, fragmentsd)
		elif data_type == 'pe':
			command = command + 'kallisto quant -t 10 -i {0} -o {1}/{2}/{3:02d}_klout/ --genomebam --gtf {7} --chromosomes {8} {1}/{6}{4}_1{5}.fastq {1}/{6}{4}_2{5}.fastq; '.format(kal_index,output_dir,dedup_mode,index,sampleID,suffix,subDIR,gtf,chrom)
	job_id = sbatch('countskallistobam_'+shortname+jobsuffix,command, index, workdir=output_dir, dep=dep,tasks=2)
	return job_id

def sortbambyname(bamIN, bamOUT, dep='', jobsuffix=''):
	""" Some algorithms only run on name sorted bam files instead of coordinate sorted ones, this function makes the conversion from coo to name sorted.
	-n: sort by read names instead of chromosomal coordinates
	-o: output (sam, bam, cram format is deduced from filename extension -> make sure it ends on .bam for bam file
	"""
	command = 'module purge; module load SAMtools/1.8-intel-2018a; '
	command = command + 'samtools sort -o {0}/{1}/{2:02d}_srout/{3} -n {0}/{1}/{2:02d}_srout/{4}; '.format(output_dir,dedup_mode,alignindex,bamOUT, bamIN)
	job_id = sbatch('sortbambyname_'+shortname+jobsuffix,command, index, workdir=output_dir, dep=dep, mem=40,tasks=2)
	return job_id

def countshtseq(bamIN, dep='', jobsuffix=''):
	""" HTSeq quantification of STAR mapped bam files.
	--order name: (default) needs name sorted bam files
	--nonunique none: (default) if the union of features for each position in the read is > 1, the read (pair) is counted as ambiguous and not counted for any features +
	if read (pair) aligns to more than one location in reference, it is scored as alignment_not_unique (for each location)
	--stranded reverse: for PE reads, the second read has to be on same strand and the first read has to be on opposite strand (for stranded=yes it is vice versa)
	"""
	command = 'module purge; module load HTSeq/0.11.0-foss-2018b-Python-2.7.15; '
	command = command + 'mkdir {0}/{1}/{2:02d}_htout; '.format(output_dir, dedup_mode, index)
	if stranded == 'yes':
		command = command + 'htseq-count --format bam --order name --nonunique none --stranded reverse {0}/{1}/{2:02d}_srout/{3} {4} > {0}/{1}/{6:02d}_htout/{5}_htseq_counts.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,gtf,shortname,index)
	else:
		command = command + 'htseq-count --format bam --order name --nonunique none --stranded no {0}/{1}/{2:02d}_srout/{3} {4} > {0}/{1}/{6:02d}_htout/{5}_htseq_counts.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,gtf,shortname,index)
	#'rm {0}_srout/{1}'.format(sampleID,bamIN) #remove name sorted bam for memory reasons (not needed anymore in pipeline)
	job_id = sbatch('countshtseq_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep, tasks=2)
	return job_id

def idxstat(bamIN,dep='', jobsuffix=''):
	""" Perform samtools idxstats
	Retrieve and print stats in the index file corresponding to the input file. The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
	Before calling idxstats, the input BAM file should be indexed by samtools index.
	"""
	command = 'module purge; module load SAMtools/1.8-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_idxout; '.format(output_dir, dedup_mode, index)
	command = command + 'samtools index {0}/{1}/{2:02d}_srout/{3}; '.format(output_dir,dedup_mode,alignindex,bamIN) #index bam file
	command = command + 'samtools idxstats {0}/{1}/{2:02d}_srout/{3} > {0}/{1}/{4:02d}_idxout/{5}_idxstat.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,index,shortname)
	job_id = sbatch('idxstat_'+shortname+jobsuffix, command, index,workdir=output_dir, dep=dep)
	return job_id

def strandedness(sampleID, bamIN,dep='',jobsuffix=''):
	""" RSeQC to retrieve % correct strandedness.
	Grep the line that shows what fraction of reads is explained by fr-firststrand (reverse in htseq)
	(1+-,1-+,2++,2-- category, e.g. 1+- read 1 '+' mapped to + strand while gene is on '-' strand is what we expect for reverse stranded)
	"""
	command = 'module purge; module load Python/2.7.14-intel-2018a; module load bx-python/0.8.1-intel-2018a-Python-2.7.14; module load RSeQC/2.6.4-intel-2018a-Python-2.7.14; '
	#command = command + 'mkdir {0}/{1}/{2:02d}_strandedness; '.format(output_dir, dedup_mode, index)
	if deduplication == 'clumpify':
		# run STAR first on FASTQ without deduplication
		command = command + 'module load STAR/2.6.0c-intel-2018a; '
		command = command + 'mkdir {0}/{1}/nodedup_srout; '.format(output_dir, dedup_mode, index)
		#command = command + 'if [ {0} = . ] ; then echo {0}; else mkdir {0}; fi; '.format(subdir) #if subdir is specified, make directory for it (if it does not already exist)
		#'STAR --runThreadN 10 --outFileNamePrefix {}_srout/ --readFilesIn {}.fastq {}.fastq --genomeDir {} --sjdbGTFfile {} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag OneBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20'.format(sampleID,sampleID1,sampleID2,star_index,gtf)
		#samplename, suffix=fastq_suffix, subDIR=sub_dir, dep=subsclump_jobid
		#star_index,output_dir,dedup_mode,index,sampleID,suffix,gtf,subDIR)
		if data_type == 'se':
			command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/nodedup_srout/ --readFilesIn {1}/{7}{4}{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,samplename,fastq_suffix,gtf,sub_dir)
		elif data_type == 'pe':
			command = command + 'STAR --runThreadN 10 --outFileNamePrefix {1}/{2}/nodedup_srout/ --readFilesIn {1}/{7}{4}_1{5}.fastq {1}/{7}{4}_2{5}.fastq --genomeDir {0} --sjdbGTFfile {6} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outMultimapperOrder Random --outSAMmultNmax -1 --outFilterMultimapNmax 10 --outSAMprimaryFlag AllBestScore --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20; '.format(star_index,output_dir,dedup_mode,index,samplename,fastq_suffix,gtf,sub_dir)
		command = command + 'rm {1}/{2}/nodedup_srout/Unmapped.*; '.format(output_dir, dedup_mode, index)
		command = command + 'infer_experiment.py -r {3} -i {0}/{1}/nodedup_srout/{2} > {0}/{1}/{4}_RSeQC_output_all.txt; '.format(output_dir,dedup_mode,bamIN,exon_bed, sampleID) #output all metrics
		command = command + 'out=`cat {0}/{1}/{2}_RSeQC_output_all.txt | grep "1+-" | cut -d":" -f2`; '.format(output_dir,dedup_mode,sampleID) #grep the percentage we are interested it (how many are on correct strand?)
		command = command + 'echo {2} $out > {0}/{1}/{2}_RSeQC_output.txt; '.format(output_dir,dedup_mode,sampleID)
		command = command + 'rm -r {0}/{1}/nodedup_srout; '.format(output_dir,dedup_mode)
		job_id = sbatch('strandedness_'+shortname+jobsuffix, command, index, workdir=output_dir, time='8:00:00', mem=40,dep=dep,tasks=8) 
		return job_id
	else: # strandedness on BAM file (without duplicate removal) that was already created with STAR
		command = command + 'infer_experiment.py -r {4} -i {0}/{1}/{2:02d}_srout/{3} > {0}/{1}/{5}_RSeQC_output_all.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exon_bed,sampleID) #output all metrics
		command = command + 'out=`cat {0}/{1}/{2}_RSeQC_output_all.txt | grep "1+-" | cut -d":" -f2`; '.format(output_dir,dedup_mode,sampleID) #grep the percentage we are interested it (how many are on correct strand?)
		command = command + 'echo {2} $out > {0}/{1}/{2}_RSeQC_output.txt; '.format(output_dir,dedup_mode,sampleID)
		job_id = sbatch('strandedness_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep)
		return job_id
	
def coverage(bamIN,dep='',jobsuffix=''):
	""" BEDTools: for downstream coverage analyses
	"""
	command = 'module purge; module load BEDTools/2.27.1-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_coverage; '.format(output_dir, dedup_mode, index)
	command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bga -split -strand + | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_plus_strand.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonplus_bed,index,shortname)
	command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bga -split -strand - | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_minus_strand.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonminus_bed,index,shortname)
	#command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bga -split | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_all_strand.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonplus_bed,index,shortname)
	job_id = sbatch('coverage_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep, tasks=2)
	

def coverageNoA(bamIN,dep='',jobsuffix=''):
	""" BEDTools: for downstream coverage analyses
	"""
	command = 'module purge; module load BEDTools/2.27.1-intel-2018a; '
	command = command + 'mkdir {0}/{1}/{2:02d}_coverage; '.format(output_dir, dedup_mode, index)
	command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bg -split -strand + | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_plus_strand_noA.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonplus_bed,index,shortname)
	command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bg -split -strand - | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_minus_strand_noA.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonminus_bed,index,shortname)
	#command = command + 'genomeCoverageBed -ibam {0}/{1}/{2:02d}_srout/{3} -bg -split | intersectBed -a {4} -wo -b - > {0}/{1}/{5:02d}_coverage/{6}_coverage_per_nt_all_strand_noA.txt; '.format(output_dir,dedup_mode,alignindex,bamIN,exonplus_bed,index,shortname)
	job_id = sbatch('coverageNoA_'+shortname+jobsuffix,command,index,workdir=output_dir, dep=dep, tasks=2)

def zipfastq(dep='',jobsuffix=''):
	""" Gzip fastq files and change permissions """
	command = 'gzip {0}/{1}/*.fastq; gzip {0}/*.fastq; '.format(output_dir, dedup_mode) #gzip all fastq files
	command = command + 'chmod -R 774 {0}/{1}; chgrp -R gvandesompele_lab {0}/{1}; '.format(output_dir, dedup_mode)
	job_id = sbatch('zipfastq_'+shortname+jobsuffix,command,index, workdir=output_dir, tasks=2, time='8:00:00', dep=dep)
	return job_id

### Select the functions you want to use + make sure the dependencies match!

## Make sure you have a file with the samples you want (e.g. by command line: ls | grep "RNA" > listsamples.txt)
#samples = [line.rstrip() for line in open(origdir+"/"+sys.argv[2])]
#for samplename in samples:
for samplename in os.listdir(base_dir): ##alternative approach if you want all samples with RNA in the name immediately
	if os.path.isdir(os.path.join(base_dir,samplename)):
		if re.search(string_match,samplename):
			print(samplename)
			output_dir = args.o[0].rstrip("/")
			input_dir = base_dir+"/"+samplename
			output_dir = output_dir+"/"+samplename
			shortname=str.split(samplename,"-")[0] #retrieve only the first part of the name to use in job name
			
			if deduplication == 'clumpify':
				dedup_mode = 'dedup_clumpify'
			elif deduplication == 'picard':
				dedup_mode = 'dedup_picard'
			else:
				dedup_mode = 'dedup_none'
			
			# name of subdirectory to put files in
			if subs_timing == 'start':
				dedup_mode += '-subs_atstart'
			elif subs_timing == 'afterclumpify':
				dedup_mode += '-subs_afterdedup'
			else:
				dedup_mode += '-subs_none'
				
			
			#make subdir for each sample in working directory (if it does not exist yet)
			os.makedirs(output_dir, exist_ok=True)
			os.makedirs(output_dir+"/"+dedup_mode+"/scripts", exist_ok=True)
			#os.makedirs(output_dir+"/"+dedup_mode+"/tmp", exist_ok=True)
			os.makedirs(output_dir+"/"+dedup_mode+"/logs", exist_ok=True)
			
			index = 0
			sub_dir = ''
			
			#### General fastq generation (if needed)
			# copy and unzip everything to this working directory
			if repeat_analysis == 'no':
				# copy and unzip everything to this working directory
				copycombine_jobid = combinecopy(samplename)
				index += 1
				fastqc_jobid = fastqc(samplename,'','FASTQC_original', dep=copycombine_jobid, jobsuffix='_orig')
				index += 1
				fastq_suffix = ''
				
				#lentrim_jobid = cutadapt_l75(samplename, suffix=fastq_suffix,dep=copycombine_jobid)
		
				#index += 1
				if adaptercontam != 'no':
					if adaptercontam == 'truseq':
						adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
						adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
					elif adaptercontam == 'nextera':
						adapter1 = 'CTGTCTCTTATACACATCT' # For Nextera/TruSight, the same sequence is used for both reads (cf Illumina website)
						adapter2 = 'CTGTCTCTTATACACATCT'
					trimadapt_jobid = trimadapter(samplename, suffix=fastq_suffix, adapterR1=adapter1, adapterR2=adapter2, dep=copycombine_jobid, jobsuffix='')
					fastq_suffix += '_trim'
				else:
					trimadapt_jobid = copycombine_jobid
				trimlength_jobid = trimlength(samplename, suffix=fastq_suffix, dep=trimadapt_jobid)
				index += 1
				fastq_suffix += '_len'
				
				removebadqc_jobid = removebadqc(samplename, suffix=fastq_suffix, dep=trimlength_jobid)
				index += 1
				fastq_suffix = '_qc'
				fastqcfilt_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, dep=removebadqc_jobid, jobsuffix='_filt') #fastq file not in subdir yet
				index +=1
				### end of general fastq generation
				
				# if deduplication == 'clumpify': #already determine strandedness on original fastq (for others, it will be determined after running STAR)
				# 	# Determine strandedness (first STAR will map fastq without duplicate removal, then strandedness determined)
				# 	strandcoordsorted_jobid = strandedness('Aligned.sortedByCoord.out.bam',dep=removebadqc_jobid,jobsuffix='')
				# 	index += 1
			
			else: #not first analysis (repeat_analysis = 'yes')
				#skip the steps above and immediately start with the other tasks (by putting the dependency as '')
				index = 20 #start from a new index (to make clear it is a repeat)
				fastq_suffix = '_qc'
				preprepeat_jobid = preprepeat(samplename)
				#preprepeat_jobid = ''
				index += 1
				removebadqc_jobid = preprepeat_jobid
			####
			
			if (data_type == 'se') & ((fragmentmean == '-1') | (fragmentsd =='-1')): #if SE and options for fragment length mean and sd were not changed:
				sys.exit("ERROR: For single end sequencing, kallisto requires a fragment length mean and sd")
				
			if subs_timing == 'start':
				if subsample_to_nr == '0':
					continue #go to next samplename
					#sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/${sample}_1_qc.fastq; done > lines_qcfilt_fastq.txt")
					# subsample to minimum based on wc -l divided by 4 and floored to a million
				else:
					#sub_dir = ''
					subsstart_jobid = subsample(samplename, str(subsample_to_nr), suffix=fastq_suffix, subDIR=sub_dir, dep=removebadqc_jobid, jobsuffix='')
					fastq_suffix += '_subs'
					sub_dir = dedup_mode+'/'
					index += 1
					#fastqctrim_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, jobsuffix='_trim')
					#index += 1
			else:
				#sub_dir = ''
				subsstart_jobid = removebadqc_jobid
			
			#names of bam files
			bam_coord = samplename+'.Aligned.sortedByCoord.out.bam'
			bam_name = samplename+'.Aligned.sortedByName.out.bam'
			
			if deduplication == 'clumpify':
				# Determine strandedness (first STAR will map fastq without duplicate removal, then strandedness determined)
				strandcoordsorted_jobid = strandedness(samplename, bam_coord,dep=subsstart_jobid,jobsuffix='')
				index += 1
				clumpify_jobid = clumpifydupremoval_trimmed(samplename,suffix=fastq_suffix,subDIR=sub_dir,dep=strandcoordsorted_jobid)
				#fastqc_clumped_jobid = fastqc(samplename+'_1_trim_qcfil80_subs_clumped',samplename+'_2_trim_qcfil80_subs_clumped','FASTQC_clumped', subdir='subs_dedup', dep=clumpify_jobid, jobsuffix='_clumped')
				fastq_suffix += '_clumped'
				sub_dir = dedup_mode+'/'
				index += 1
				
				if subs_timing == 'afterclumpify':
					if subsample_to_nr == '0':
						continue #go to next samplename
						#sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/dedup_clumpify-subs_afterdedup/{sample}_1_qc_clumped.fastq; done > lines_qcfilt_fastq.txt")
						# subsample to minimum based on wc -l divided by 4 and floored to a million
					else:
						subsclump_jobid = subsample(samplename, str(subsample_to_nr), suffix=fastq_suffix, subDIR=sub_dir, dep=clumpify_jobid, jobsuffix='')
						fastq_suffix += '_subs'
						index += 1
				else:
					subsclump_jobid = clumpify_jobid
				
				fastqcclump_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, subDIR=sub_dir, dep=subsclump_jobid, jobsuffix='_clump')
				index += 1
				
				# Kallisto quantification on Clumpify duplicate removed fastq files
				kallistoclump_jobid = countskallisto(samplename,suffix=fastq_suffix,subDIR=sub_dir,dep=subsclump_jobid, jobsuffix='_clump')
				index += 1
				
				# Mapping with STAR of Clumpify duplicate removed FASTQ
				star_jobid = alignment(samplename, suffix=fastq_suffix, subDIR=sub_dir, dep=subsclump_jobid, jobsuffix='_clump')
				alignindex = index
				endfastq_jobid=star_jobid
				namesortedbam_jobid = sortbambyname(bamIN=bam_coord, bamOUT=bam_name, dep=star_jobid, jobsuffix='_clump')
				index += 1
			
			else: #no duplicate removal or picard dup removal
				# Mapping with STAR of low quality reads filtered FASTQ
				star_jobid = alignment(samplename, suffix=fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, jobsuffix='_noClump')
				alignindex = index
				index +=1
				
				# RSeQC: to retrieve % correct strandedness (based on original BAM - no duplicate removal)i
				strandcoordsorted_jobid = strandedness(samplename, bam_coord,dep=star_jobid,jobsuffix='')
				index += 1
			
				if deduplication == 'picard':
					# update bam file names for further use
					bam_coord_nodedup = bam_coord
					bam_coord = samplename+'.Aligned.sortedByCoord.picard.bam'
					bam_name = samplename+'.Aligned.sortedByName.picard.bam'
					# duplicate removal based on coordinate sorted bam (Picard)
					picardcoord_jobid = picarddupremcoord(samplename, bamIN=bam_coord_nodedup, bamOUT=bam_coord, dep=star_jobid, jobsuffix='_picard') #dep=starsubs_jobid
					namesortedbam_jobid = sortbambyname(bamIN=bam_coord, bamOUT=bam_name, dep=picardcoord_jobid, jobsuffix='_picard')
					index += 1
					# Kallisto quantification (based on fastq files after duplicate removal -> bamtofastq conversion needed)
					bamtofastq_jobid = bamtofastqfile(samplename, bamIN=bam_name, dep=namesortedbam_jobid, jobsuffix='_picard')
					fastq_suffix = '_picard'
					sub_dir = dedup_mode + '/'
					kallisto_jobid = countskallisto(samplename, suffix=fastq_suffix, subDIR=sub_dir, dep=bamtofastq_jobid, jobsuffix='_picard')
					index += 1
					# fastqc of picard dup removed + converted files
					fastqc_picard_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, subDIR=sub_dir, dep=bamtofastq_jobid, jobsuffix='_picard')
					endfastq_jobid=fastqc_picard_jobid
					index += 1
					
				else: #no duplicate removal
					namesortedbam_jobid = sortbambyname(bamIN=bam_coord, bamOUT=bam_name, dep=star_jobid, jobsuffix='_nodedup')
					index += 1
					# Kallisto quantification of original (qc filtered) fastq files
					kallisto_jobid = countskallisto(samplename,suffix=fastq_suffix,subDIR=sub_dir,dep=star_jobid, jobsuffix='_nodedup')			
					index += 1
					endfastq_jobid=kallisto_jobid
			
			# HTSeq quantification (bam needs to be sorted by name)
			htseq_jobid = countshtseq(bamIN=bam_name, dep=namesortedbam_jobid, jobsuffix='')
			index += 1
			# Run idxstats (samtools)
			idxstat_jobid = idxstat(bamIN=bam_coord, dep=namesortedbam_jobid, jobsuffix='')
			index += 1
			
			if coverage == 'yes':
				# BEDTools: for downstream coverage analyses
				covcoordsorted_jobid = coverage(bam_coord,dep=namesortedbam_jobid, jobsuffix='')
				index += 1
				#covcoordsortednoAclump_jobid = coverageNoA('Aligned.sortedByCoord.out.bam',subdir='subs_dedup',dep=starsubs_jobid,jobsuffix='_clump')
			
			# Gzip all fastq files
			zipfastq_jobid = zipfastq(dep=endfastq_jobid)

if (subs_timing == 'start') & (subsample_to_nr == '0'):
	sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/${sample}_1_qc.fastq; done > lines_qcfilt_fastq.txt")

elif (subs_timing == 'afterclumpify') & (subsample_to_nr == '0'):
	sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/dedup_clumpify-subs_afterdedup/{sample}_1_qc_clumped.fastq; done > lines_qcfilt_fastq.txt")
