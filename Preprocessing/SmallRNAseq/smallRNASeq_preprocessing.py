# make sure you use python 3
# python smallRNASeq_preprocessing.py <other args>
# e.g. python smallRNASeq_preprocessing.py -m RNA0 -t se --subs start -n 20000000 --repeat no -b $VSC_DATA_VO/NSQ_run453-84636552_Run458-86026942 -o $VSC_SCRATCH_VO/test_pipeline -u example@ugent.be

import subprocess
import os
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Submit smallRNA Seq preprocessing jobs to the UGent HPC (slurm)')

# Read arguments
parser.add_argument('-t', nargs=1, choices=['se','pe'], required=True, help='Single end (se) or paired end (pe) data')
parser.add_argument('-b', nargs=1, required=True, help='Base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-o', nargs=1, required=True, help='Directory where output should be created', metavar='output_dir')
parser.add_argument('-m', nargs=1, required=True, help='String match to select sample folders in base directory e.g. RNA0', metavar='string_match')
parser.add_argument('--mismatch', choices=['0','1','2','3'], nargs=1, default=['0'], help='Number of mismatches [default: 0] allowed in first L (option -l) bases')
parser.add_argument('-l', nargs=1, default=['25'], help='Mismatches (--mismatch) allowed in first L bases [default: 25, should be >= 5]')
parser.add_argument('--subs', nargs=1, choices=['no','yes'], default=['no'], help='Subsampling needed? [default: no]')
parser.add_argument('-n', nargs=1, default=['0'], help='Nr of reads to subsample to. If you need the analysis to stop in order to determine the subsampling level first, do not enter a number here (or -n 0)')
parser.add_argument('--repeat', nargs=1, choices=['yes','no'], default=['no'], help='Repeated analysis? [default: no] If yes: general fastq copying and filtering will not be repeated')
parser.add_argument('-u', nargs=1, required=True, help='Submitter email address (used for error reporting)', metavar='user_email')

# Parse arguments
args = parser.parse_args()
print(args)
data_type = args.t[0]
base_dir = args.b[0].rstrip("/")
string_match = args.m[0]
email = args.u[0]
subsampling = args.subs[0]
repeat_analysis = args.repeat[0]
subsample_to_nr = args.n[0] # nr of reads you want to subsample to (based on floor of wc -l divided by 4), will be ignored if no subsampling asked for in command line
mismatch = args.mismatch[0]
mismatchregion = args.l[0]

tasks=1
adapter = "TGGAATTCTCGGGTGCCAAGG"
#refgenome = "GRCh37"
refgenome = "bt1_hg38"

def sbatch(job_name, command, index, mem='', tasks = tasks, workdir = '.', time='4:00:00', dep=''): 
	if dep != '':
		dep = ' --dependency=afterok:{} --kill-on-invalid-dep=yes'.format(dep) 
	printlines = [
		"#!/bin/bash",
		"",
		"#SBATCH -J {}".format(job_name+'_'+subs_mode),
		"#SBATCH -D {}".format(workdir),
		#"#SBATCH --mem={}G".format(mem),
		"#SBATCH --cpus-per-task={}".format(tasks), #nr of processors per task needed
		"#SBATCH -t {}".format(time),
		"#SBATCH --mail-user={}".format(email),
		"#SBATCH --mail-type=FAIL",
		"#SBATCH -o {1}/{2}/logs/{3:02d}_{0}.out".format(job_name, workdir, subs_mode, index),
		"#SBATCH -e {1}/{2}/logs/{3:02d}_{0}.err".format(job_name, workdir, subs_mode, index),
		#"#SBATCH --test-only", #Uncomment this to test instead of directly submitting job
		""
	]
	
	printlines.extend(command.split("; "))
	#print(printlines)
	printjobfile('{0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, subs_mode, index, job_name),printlines)
	#print(job_name)
	sbatch_command = 'sbatch{4} {0}/{1}/scripts/{2:02d}_{3}.sh'.format(output_dir, subs_mode, index, job_name, dep)
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
		command = '[[ -s {0}/{2}.fastq.gz ]] && cp {0}/{2}.fastq.gz {1}/{2}.fastq.gz || cat {0}/*R1*.fastq.gz > {1}/{2}.fastq.gz; '.format(input_dir, output_dir, sampleID)
	elif data_type == 'pe':
		command = '[[ -s {0}/{2}_1.fastq.gz ]] && cp {0}/{2}_1.fastq.gz {1}/{2}_1.fastq.gz || cat {0}/*R1*.fastq.gz > {1}/{2}_1.fastq.gz; '.format(input_dir, output_dir, sampleID)
		command = command + '[[ -s {0}/{2}_2.fastq.gz ]] && cp {0}/{2}_2.fastq.gz {1}/{2}_2.fastq.gz || cat {0}/*R2*.fastq.gz > {1}/{2}_2.fastq.gz; '.format(input_dir, output_dir, sampleID)
	command = command + 'gunzip {0}/{1}*.fastq.gz; '.format(output_dir, sampleID)
	job_id = sbatch('combinecopy_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)
	return job_id

### FastQC
def fastqc(sampleID,suffix='',outDIR='.',subDIR='',dep='',jobsuffix=''): 
	# Build the command for fastqc 
	command = 'module load FastQC/0.11.8-Java-1.8; mkdir {}/{}{}; '.format(output_dir, subDIR, outDIR)
	if data_type == 'se':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}{4}.fastq; '.format(output_dir, subDIR, outDIR, sampleID, suffix)
	elif data_type == 'pe':
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_1{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
		command = command + 'fastqc -o {0}/{1}{2} {0}/{1}{3}_2{4}.fastq; '.format(output_dir, subDIR, outDIR,sampleID,suffix)
	job_id = sbatch('fastqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep)
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

### Cutadapt 
def cutadapt_smallRNA(sampleID,suffix='',subDIR='',dep='',jobsuffix=''):
	""" Remove adapters and trim reads
	--discard-untrimmed: if this option is used with a linked adapter, read is considered to be trimmed only if all required adapters were found
	-a: regular 3' adapter R1 (To reduce the number of falsely trimmed bases, the alignment algorithm requires that, by default, at least three bases match between adapter and read)
	-A: regular 3' adapter R2
	-q x,y: quality-trim read with a threshold of x from 3' and with threshold of y from 5'
	-o & -p: output file for first read (-o) and second read (-p) of pair
	-m: minimum length
	-e: maximum error rate in adapter
	--pair-filter=any: (default) read pair is discarded if one of the reads (R1 or R2) fulfills the filtering criterion (e.g. one read is shorter than 20)
	"""
	command = 'module purge; module load cutadapt/1.16-intel-2018a-Python-3.6.4; module load FASTX-Toolkit/0.0.14-intel-2018a; '
	command = command + 'cutadapt --discard-untrimmed -m 15 -e 0.15 -q 20 -a {2} -o {0}/{1}_clipped.fastq {0}/{1}.fastq; '.format(output_dir,sampleID,adapter)
	job_id = sbatch('cutadaptsmall_'+shortname+jobsuffix, command, index,workdir=output_dir, dep=dep,tasks=4)
	return job_id

def removebadqc(sampleID,suffix='',percentage=str(80),quality=str(19),dep='',jobsuffix=''):
	"""
	fastq_quality_filter: remove reads with lower quality, keeping only reads that have at least p% of bases with a quality score of q or more
	FASTQC_filter.py: Keep only read pairs where each read has at least p% of bases with a phred score > q
	-q: minimum quality phred score
	-p: percentage that has the minimum quality phred score
	-i: FASTA/Q input file
	-o: FASTA/Q output file
	-Q33: If Illumina encoding is >= 1.8, you need to provide this option to fastq_quality_filter. For Illumina encoding <=1.5, option is not required
	"""
	#command = command + 'mkdir {0}/{1}/{2:02d}_removebadqcout; '.format(output_dir, subs_mode, index)
	if data_type == 'se':
		command = 'module purge; module load FASTX-Toolkit/0.0.14-intel-2018a; '
		command = command + 'fastq_quality_filter -p {3} -q {4} -Q33 -i {0}/{1}{2}.fastq -o {0}/{1}_qc.fastq; '.format(output_dir, sampleID, suffix, percentage, str(int(quality)+1))
		command = command + 'echo -e "sampleID\\torig_lines\\tclipped_lines\\tqc_lines" > {0}/{1}_line_count.txt; sampleID={1}; orig_lines=`wc -l {0}/{1}.fastq | cut -f1 -d\' \'`; clipped_lines=`wc -l {0}/{1}_clipped.fastq | cut -f1 -d\' \'`; qc_lines=`wc -l {0}/{1}_qc.fastq | cut -f1 -d\' \'`; echo -e "$sampleID\\t$orig_lines\\t$clipped_lines\\t$qc_lines" >> {0}/{1}_line_count.txt; '.format(output_dir,sampleID)
	elif data_type == 'pe':
		command = 'module purge; module load Biopython/1.72-foss-2018b-Python-3.6.6; '
		command = command + 'python FASTQC_filter.py -r1 {0}/{1}_1{2}.fastq -r2 {0}/{1}_2{2}.fastq -p {3} -q {4} -o1 {0}/{1}_1_qc.fastq -o2 {0}/{1}_2_qc.fastq; '.format(output_dir, sampleID, suffix, percentage, quality)
		command = command + 'rm {0}/{1}_1{2}.fastq; rm {0}/{1}_2{2}.fastq; '.format(output_dir, sampleID, suffix)
	job_id = sbatch('removebadqc_'+shortname+jobsuffix,command,index,workdir=output_dir,dep=dep, tasks=4)
	return job_id

def prepcounts(sampleID, suffix='',subDIR='',dep='',jobsuffix=''):
	"""
	fastq_collapser: collapses identical sequences in a FASTQ or FASTA file into a single sequence with the nr of occurences in line above sequence (behind the unique read number)
	-i: FASTA/Q input file. default is STDIN
	-o: FASTA/Q output file. default is STDOUT
	-v: verbose (short summary of input/output counts)
	-Q33: If Illumina encoding is >= 1.8, you need to provide this option. For Illumina encoding <=1.5, option is not required
	Getspikes.py: retrieves the number of sequences of each spike (using the collapsed fastq information) and outputs a fasta without spikes
	"""
	if data_type == 'se':
		command = 'module purge; module load FASTX-Toolkit/0.0.14-intel-2018a; '
		command = command + 'fastx_collapser -i {0}/{1}{2}{3}.fastq -o {0}/{4}/{2}_collapse.fa -v -Q33; '.format(output_dir,subDIR,sampleID,suffix,subs_mode)
		command = command + 'ml Biopython/1.71-intel-2018a-Python-3.6.4; '
		command = command + 'python Getspikes.py spike_sequences.csv {0}/{1}/{2}_collapse.fa {0}/{1}/{2}; '.format(output_dir,subs_mode,sampleID)
	elif data_type == 'pe':
		command = ''.format() ####
	job_id = sbatch('prepcounts_'+shortname+jobsuffix, command, index, workdir=output_dir,dep=dep,tasks=4)
	return job_id

def mapbowtie(sampleID,dep='',jobsuffix=''):
	"""
	Calculate counts with bowtie + general statistics
	Bowtie:
	-f: FASTA as input format (-q would be FASTQ)
	-k: Report up to <int> valid alignments per read or pair (default: 1)
	-n: When the -n option is specified (which is the default), bowtie determines which alignments are valid according to the following policy:
	   * alignments may have no more than N mismatches (where N is a number 0-3, set with -n) in the first L bases (where L is a number 5 or greater, set with -l) on the high-quality (left) end of the read. first L bases = seed.
	   * sum of the Phred quality values at all mismatched positions (not just in the seed) may not exceed E (set with -e). Where qualities are unavailable (e.g. if the reads are from a FASTA file), the Phred quality defaults to 40.
	--best: Running Bowtie in --best mode eliminates strand bias by forcing Bowtie to select one strand or the other with a probability that is proportional to the number of best sites on the strand.
	"""
	if data_type == 'se':
		command = 'module purge; module load Bowtie/1.2.2-intel-2018a; '
		command = command + 'bowtie -f -k 10 -n {4} -l {5} --best {3} {0}/{1}/{2}_collapsed_nospikes.fa {0}/{1}/{2}_mapped.sam; '.format(output_dir,subs_mode,sampleID,refgenome, mismatch, mismatchregion)
		command = command + 'module load BioPerl/1.7.2-intel-2018a-Perl-5.26.1; '
	elif data_type == 'pe':
		command = ''.format() ####
	job_id = sbatch('mapbowtie_'+shortname+jobsuffix, command, index, workdir=output_dir, dep=dep, tasks=4)
	return job_id

def annotate(sampleID,dep='',jobsuffix=''):
	"""
	Based on bowtie output file: annotate miRs, isomiRs, contaminants
	When a read is assigned both to contam/miR AND to not_annotated, remove these reads from not_annotated.txt file (based on unique read ID)
	In the end, calculate number of mapped reads
	"""
	if data_type == 'se':
		command = 'module purge; module load BioPerl/1.7.2-intel-2018a-Perl-5.26.1; '
		command = command + 'cd {0}/{1}; '.format(output_dir, subs_mode)
		command = command + 'perl match_bowtie_outputv5.pl {0}/{1} {2}; '.format(output_dir,subs_mode,sampleID)
		command = command + 'cat {0}/{1}/{2}_not_annotated.txt | sort -nk1,1 -u | awk \'{{print $1}}\' | sort -nk1 > {0}/{1}/uniq_reads_notann.txt; '.format(output_dir, subs_mode, sampleID)
		command = command + 'cat {0}/{1}/{2}_contam.txt | sort -nk2,2 -u | awk \'{{print $2}}\' | sort -nk1 > {0}/{1}/uniq_reads_contam.txt; '.format(output_dir, subs_mode, sampleID)
		command = command + 'cat {0}/{1}/{2}_isomiRs.txt | sort -nk2,2 -u | awk \'{{print $2}}\' | sort -nk1 > {0}/{1}/uniq_reads_miRs.txt; '.format(output_dir, subs_mode, sampleID)
		command = command + 'sort -nk1 {0}/{1}/uniq_* | uniq -d > {0}/{1}/redundant_notann.txt; '.format(output_dir, subs_mode, sampleID)
		command = command + 'awk \'NR==FNR {{a[$0]=$0}} NR>FNR {{ if($1!=a[$1]) print $0}}\' {0}/{1}/redundant_notann.txt {0}/{1}/{2}_not_annotated.txt > {0}/{1}/{2}_not_annotated_nomiR_nocontam.txt; '.format(output_dir, subs_mode, sampleID)
		command = command + 'rm {0}/{1}/uniq_*; '.format(output_dir, subs_mode, sampleID)
		command = command + 'perl get_read_length.pl {0}/{1} {2}; '.format(output_dir,subs_mode,sampleID)
		command = command + 'perl get_num_mapped_reads.pl {0}/{1} {2}; '.format(output_dir,subs_mode,sampleID)
		command = command + 'notann_counts=`cat {0}/{1}/{2}_not_annotated_nomiR_nocontam.txt | sort -k1,1 -u | awk \'{{sum+=$2}} END {{print sum}}\'`; '.format(output_dir, subs_mode, sampleID)
		command = command + 'echo -e "mapped_notann\\t$notann_counts" >> {0}/{1}/{2}_read_count_new.txt; '.format(output_dir,subs_mode, sampleID)
		command = command + 'module load R/3.4.4-intel-2018a-X11-20180131; '
		command = command + 'Rscript plot_read_length_distribution.R {0}/{1} {2}; '.format(output_dir,subs_mode,sampleID)
		#command = command + 'Rscript plot_read_count.R'.format(output_dir,subs_dir,sampleID)
	elif data_type == 'pe':
		command = ''.format() ####
	job_id = sbatch('annotate_'+shortname+jobsuffix, command, index, workdir=output_dir, dep=dep, tasks=4)
	return job_id

### Subsampling
def subsample(sampleID,nreads,suffix='',subDIR='',dep='',jobsuffix=''):
	""" Seqtk subsampling
	Downsample fastq files to x number of reads (take min of all samples)
	-s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of code
	"""
	command = 'module purge; module load seqtk/1.3-foss-2018a; '
	#if subsample_to_nr == '0': #if no subsampling nr given
	#	command = command + 'echo \"Determine the subsampling level first! (i.e. floored minimum of wc-l fastq divided by 4)\"; exit 1; '
	if data_type == 'se':
		command = command + 'seqtk sample -s100 {0}/{1}{2}{3}.fastq {4} > {0}/{5}/{2}{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,subs_mode)
	elif data_type == 'pe':
		command = command + 'seqtk sample -s100 {0}/{1}{2}_1{3}.fastq {4} > {0}/{5}/{2}_1{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,subs_mode)
		command = command + 'seqtk sample -s100 {0}/{1}{2}_2{3}.fastq {4} > {0}/{5}/{2}_2{3}_subs.fastq; '.format(output_dir,subDIR,sampleID,suffix,nreads,subs_mode)
	job_id = sbatch('subsample_'+shortname+jobsuffix,command,index,workdir=output_dir,mem=40, time='4:00:00', dep=dep,tasks=2)
	return job_id

def zipfastq(dep='',jobsuffix=''):
	""" Gzip fastq files and change permissions """
	command = 'gzip {0}/{1}/*.fastq; gzip {0}/*.fastq; '.format(output_dir, subs_mode) #gzip all fastq files
	command = command + 'chmod -R 774 {0}/{1}; chgrp -R gvandesompele_lab {0}/{1}; '.format(output_dir, subs_mode)
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
					
			if subsampling == 'no':
				subs_mode = 'dedup_none-subs_none'
			elif subsampling == 'yes':
				subs_mode = 'dedup_none-subs_atstart'			#make subdir for each sample in working directory (if it does not exist yet)
		
			os.makedirs(output_dir, exist_ok=True)
			os.makedirs(output_dir+"/"+subs_mode+"/scripts", exist_ok=True)
			#os.makedirs(output_dir+"/"+subs_mode+"/tmp", exist_ok=True)
			os.makedirs(output_dir+"/"+subs_mode+"/logs", exist_ok=True)
				
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
				
				cutadapt_jobid = cutadapt_smallRNA(samplename, suffix=fastq_suffix,dep=copycombine_jobid)
				fastq_suffix += '_clipped'
				index += 1
				removebadqc_jobid = removebadqc(samplename, suffix=fastq_suffix, dep=cutadapt_jobid)
				index += 1
				fastq_suffix = '_qc'
				fastqcfilt_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, dep=removebadqc_jobid, jobsuffix='_filt') #fastq file not in subdir yet
				### end of general fastq generation
			
			else: #not first analysis (repeat_analysis = 'yes')
				#skip the steps above and immediately start with the other tasks (by putting the dependency as '')
				index = 20 #start from a new index (to make clear it is a repeat)
				fastq_suffix = '_qc'
				preprepeat_jobid = preprepeat(samplename)
				index += 1
				removebadqc_jobid = preprepeat_jobid
			####
			
			# Subsample if needed
			if subsampling == 'yes':
				if subsample_to_nr == '0':
					continue #go to next samplename
					#sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/${sample}_1_qc.fastq; done > lines_qcfilt_fastq.txt")
					# subsample to minimum based on wc -l divided by 4 and floored to a million
				else:
					#sub_dir = ''
					subsstart_jobid = subsample(samplename, str(subsample_to_nr), suffix=fastq_suffix, subDIR=sub_dir, dep=removebadqc_jobid, jobsuffix='')
					fastq_suffix += '_subs'
					sub_dir = subs_mode+'/'
					index += 1
					#fastqctrim_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, jobsuffix='_trim')
					#index += 1
			else:
				#sub_dir = ''
				subsstart_jobid = removebadqc_jobid
			
			# Process the files
			prepcounts_jobid = prepcounts(samplename,fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid)
			index += 1
			mapbowtie_jobid = mapbowtie(samplename, dep=prepcounts_jobid)
			index += 1
			annotate_jobid = annotate(samplename, dep=mapbowtie_jobid)
			index += 1
			# Gzip all fastq files
			zipfastq_jobid = zipfastq(dep=mapbowtie_jobid)

if (subsampling == 'yes') & (subsample_to_nr == '0'):
	sys.exit("Determine the subsampling level after quality filtering")
