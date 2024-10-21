
import os
import subprocess
import re
import sys
import yaml
import argparse


def read_file_to_tuples(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [tuple(line.strip().split()) for line in lines]

def sbatch(job_name, command, index, mem='', tasks=1, workdir='.', time='4:00:00', dep=''):
	print(f"Running job: {job_name}")
	subprocess.run(command, shell=True, check=True)

def combinecopy(sampleID, f1, f2, jobsuffix='', shortname="", index="", output_dir="", data_type='se'):
	if data_type == 'se':
		command = f'[[ -s {f1}]] && cp {f1} {output_dir}/{sampleID}.fastq.gz || cat {f1} > {output_dir}/{sampleID}.fastq.gz; gunzip {output_dir}/{sampleID}*.fastq.gz'
	elif data_type == 'pe':
		command = f'[[ -s {f1}]] && cp {f1} {output_dir}/{sampleID}_1.fastq.gz || cat {f1} > {output_dir}/{sampleID}_1.fastq.gz; [[ -s {f2} ]] && cp {f2} {output_dir}/{sampleID}_2.fastq.gz || cat {f2} > {output_dir}/{sampleID}_2.fastq.gz; gunzip {output_dir}/{sampleID}*.fastq.gz'
	sbatch(f'combinecopy_{shortname}{jobsuffix}', command, index)
	

def fastqc(sampleID, suffix='', outDIR='.', subDIR='', dep='', jobsuffix='', data_type='se'):
    command = f'mkdir -p {output_dir}/{subDIR}{outDIR};'
    if data_type == 'se':
        command = command + f'fastqc -o {output_dir}/{subDIR}{outDIR} {output_dir}/{subDIR}{sampleID}{suffix}.fastq'
    elif data_type == 'pe':
        command = command + f'fastqc -o {output_dir}/{subDIR}{outDIR} {output_dir}/{subDIR}{sampleID}{suffix}_1.fastq;'
        command = command + f'fastqc -o {output_dir}/{subDIR}{outDIR} {output_dir}/{subDIR}{sampleID}{suffix}_2.fastq'
    sbatch(f'fastqc_{shortname}{jobsuffix}', command, index)


def preprepeat(sampleID, dep='', jobsuffix='', data_type='se'):
	if data_type == 'se':
		command = f'[[ -s {output_dir}/{sampleID}_qc.fastq.gz ]] && gunzip {output_dir}/{sampleID}_qc.fastq.gz || head -1 {output_dir}/{sampleID}_qc.fastq'
	elif data_type == 'pe':
		command = f'[[ -s {output_dir}/{sampleID}_qc_1.fastq.gz ]] && gunzip {output_dir}/{sampleID}_qc_1.fastq.gz || head -1 {output_dir}/{sampleID}_qc_1.fastq; [[ -s {output_dir}/{sampleID}_qc_2.fastq.gz ]] && gunzip {output_dir}/{sampleID}_qc_2.fastq.gz || head -1 {output_dir}/{sampleID}_qc_2.fastq'
	sbatch(f'preprepeat_{shortname}{jobsuffix}', command, index)

def cutadapt_smallRNA(sampleID, suffix='', subDIR='', dep='', jobsuffix='', data_type='se'):
    if data_type == 'se':
        command = f'cutadapt --discard-untrimmed -m 15 -e 0.15 -q 20 -a {adapter} -o {output_dir}/{sampleID}_clipped.fastq {output_dir}/{sampleID}.fastq'
    elif data_type == 'pe':
        command = f'cutadapt --discard-untrimmed -m 15 -e 0.15 -q 20 -a {adapter1} -A {adapter2} -o {output_dir}/{sampleID}_1_clipped.fastq -p {output_dir}/{sampleID}_2_clipped.fastq {output_dir}/{sampleID}_1.fastq {output_dir}/{sampleID}_2.fastq'
    sbatch(f'cutadaptsmall_{shortname}{jobsuffix}', command, index)
	
def removebadqc(sampleID, suffix='', percentage='80', quality='19', dep='', jobsuffix='', data_type='se'):
	if data_type == 'se':
		command = f'python3 fastq_quality_filter_se.py -p {percentage} -q {int(quality)+1} -i {output_dir}/{sampleID}{suffix}.fastq -o {output_dir}/{sampleID}_qc.fastq; echo -e "sampleID\\torig_lines\\tclipped_lines\\tqc_lines" > {output_dir}/{sampleID}_line_count.txt; sampleID={sampleID}; orig_lines=$(wc -l {output_dir}/{sampleID}.fastq | cut -f1 -d\' \'); clipped_lines=$(wc -l {output_dir}/{sampleID}_clipped.fastq | cut -f1 -d\' \'); qc_lines=$(wc -l {output_dir}/{sampleID}_qc.fastq | cut -f1 -d\' \'); echo -e "$sampleID\\t$orig_lines\\t$clipped_lines\\t$qc_lines" >> {output_dir}/{sampleID}_line_count.txt'
	elif data_type == 'pe':
		command = f'python3 FASTQC_filter.py -r1 {output_dir}/{sampleID}_1{suffix}.fastq -r2 {output_dir}/{sampleID}_2{suffix}.fastq -p {percentage} -q {quality} -o1 {output_dir}/{sampleID}_qc_1.fastq -o2 {output_dir}/{sampleID}_qc_2.fastq; rm {output_dir}/{sampleID}_1{suffix}.fastq; rm {output_dir}/{sampleID}_2{suffix}.fastq'
	sbatch(f'removebadqc_{shortname}{jobsuffix}', command, index)

def prepcounts(sampleID, suffix='', subDIR='', dep='', jobsuffix='', spike_seqs=".", data_type='se'):
	if data_type == 'se':
		command = f'python3 fastq_collapser.py -i {output_dir}/{subDIR}{sampleID}{suffix}.fastq -o {output_dir}/{subs_mode}/{sampleID}_collapse.fa -v; python Getspikes.py {spike_seqs} {output_dir}/{subs_mode}/{sampleID}_collapse.fa {output_dir}/{subs_mode}/{sampleID}'
	sbatch(f'prepcounts_{shortname}{jobsuffix}', command, index)

def prepcounts_test(sampleID, suffix='', subDIR='', dep='', jobsuffix='', data_type='se'):
	if data_type == 'se':
		command = f'python3 fastq_collapser.py -i {output_dir}/{subDIR}{sampleID}{suffix}.fastq -o {output_dir}/{subs_mode}/{sampleID}_collapse.fa -v'
	sbatch(f'prepcounts_{shortname}{jobsuffix}', command, index)

def mapbowtie(sampleID, dep='', jobsuffix='', refgenome="", data_type='se'):
	if data_type == 'se':
		command = f'bowtie -f -k 10 -n {mismatch} -l {mismatchregion} --best {refgenome} {output_dir}/{subs_mode}/{sampleID}_collapsed_nospikes.fa {output_dir}/{subs_mode}/{sampleID}_mapped.sam'
	sbatch(f'mapbowtie_{shortname}{jobsuffix}', command, index)

def annotate(sampleID, dep='', jobsuffix='', data_type='se'):
	if data_type == 'se':
		command = f'cd {output_dir}/{subs_mode}; perl ../../../scripts/match_bowtie_outputv5.pl . {sampleID}; cat {sampleID}_not_annotated.txt | sort -nk1,1 -u | awk \'{{print $1}}\' | sort -nk1 > uniq_reads_notann.txt; cat {sampleID}_contam.txt | sort -nk2,2 -u | awk \'{{print $2}}\' | sort -nk1 > uniq_reads_contam.txt; cat {sampleID}_isomiRs.txt | sort -nk2,2 -u | awk \'{{print $2}}\' | sort -nk1 > uniq_reads_miRs.txt; sort -nk1 uniq_* | uniq -d > redundant_notann.txt; awk \'NR==FNR {{a[$0]=$0}} NR>FNR {{ if($1!=a[$1]) print $0}}\' redundant_notann.txt {sampleID}_not_annotated.txt > {sampleID}_not_annotated_nomiR_nocontam.txt; rm uniq_*; perl ../../../scripts/get_read_length.pl . {sampleID}; perl ../../../scripts/get_num_mapped_reads.pl . {sampleID}; notann_counts=$(cat {sampleID}_not_annotated_nomiR_nocontam.txt | sort -k1,1 -u | awk \'{{sum+=$2}} END {{print sum}}\'); echo -e "mapped_notann\\t$notann_counts" >> {sampleID}_read_count_new.txt'
	sbatch(f'annotate_{shortname}{jobsuffix}', command, index)
	
def makeplots(sampleID):
	f'Rscript ../../../scripts/plot_read_length_distribution.R . {sampleID}'

def subsample(sampleID, nreads, suffix='', subDIR='', dep='', jobsuffix='', data_type='se'):
	if data_type == 'se':
		command = f'seqtk sample -s100 {output_dir}/{subDIR}{sampleID}{suffix}.fastq {nreads} > {output_dir}/{subs_mode}/{sampleID}{suffix}_subs.fastq'
	elif data_type == 'pe':
		command = f'seqtk sample -s100 {output_dir}/{subDIR}{sampleID}{suffix}_1.fastq {nreads} > {output_dir}/{subs_mode}/{sampleID}{suffix}_subs_1.fastq; seqtk sample -s100 {output_dir}/{subDIR}{sampleID}{suffix}_2.fastq {nreads} > {output_dir}/{subs_mode}/{sampleID}{suffix}_subs_2.fastq'
	sbatch(f'subsample_{shortname}{jobsuffix}', command, index)

def zipfastq(dep='', jobsuffix=''):
	command = f'gzip {output_dir}/{subs_mode}/*.fastq; gzip {output_dir}/*.fastq'
	sbatch(f'zipfastq_{shortname}{jobsuffix}', command, index)

#parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('--config',required=True, help='path to yaml config file')
yaml_file = parser.parse_args()
with open(yaml_file.config, "r") as fp:
    args = yaml.safe_load(fp)
output_dir_abs = args["output_dir_abs"]
subsampling = args["subsampling"]
subsample_to_nr = args["subsample_to_nr"]
repeat_analysis = args["repeat_analysis"]
adapter = args["adapter"]
mismatch = args["mismatch"]
mismatchregion = args["mismatchregion"]
spike_seqs = args["spike_seqs"]
bowtie_index = args["bowtie_index"]
data_type = args["data_type"]
SampleList = args["SampleList"]
paired_files = read_file_to_tuples(SampleList) 

for samplename, f1, f2 in paired_files:
    output_dir = output_dir_abs + "/" + samplename
    shortname = samplename
    if subsampling == 'no':
        subs_mode = 'dedup_none-subs_none'
    elif subsampling == 'yes':
        subs_mode = 'dedup_none-subs_atstart'
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir + "/" + subs_mode + "/scripts", exist_ok=True)
    os.makedirs(output_dir + "/" + subs_mode + "/logs", exist_ok=True)
    index = 0
    sub_dir = ''
    if repeat_analysis == 'no':
        copycombine_jobid = combinecopy(samplename, f1, f2, jobsuffix='', shortname=shortname, index=index, output_dir=output_dir, data_type=data_type)
        index += 1
        fastqc(samplename, '', 'FASTQC_original', dep=copycombine_jobid, jobsuffix='_orig')
        index += 1
        fastq_suffix = ''
		# adapter trimming
        cutadapt_jobid = cutadapt_smallRNA(samplename, suffix=fastq_suffix, dep=None, data_type=data_type)
        fastq_suffix += '_clipped'
        index += 1
        # Filter fastq reads per quality score
        removebadqc_jobid = removebadqc(samplename, suffix=fastq_suffix, percentage='80', quality='19', dep='', jobsuffix='', data_type=data_type)
        index += 1
        fastq_suffix = '_qc'
        fastqcfilt_jobid = fastqc(samplename,fastq_suffix,'FASTQC'+fastq_suffix, dep=removebadqc_jobid, jobsuffix='_filt') #fastq file not in subdir yet
        continue		
    else:  # not first analysis (repeat_analysis = 'yes')
        index = 20  # start from a new index (to make clear it is a repeat)
        fastq_suffix = '_qc'
        preprepeat_jobid = preprepeat(samplename)
        index += 1
        removebadqc_jobid = preprepeat_jobid
    # Subsample if needed
    if subsampling == 'yes':
        if subsample_to_nr == '0':
            continue 
        else:
			#subsampling
            subsstart_jobid = subsample(samplename, str(subsample_to_nr), suffix=fastq_suffix, subDIR=sub_dir, dep=removebadqc_jobid, jobsuffix='')
            fastq_suffix += '_subs'
            sub_dir = subs_mode+'/'
            index += 1
            fastqctrim_jobid = fastqc(samplename, fastq_suffix, 'FASTQC' + fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, jobsuffix='_trim')
            index += 1
    else:
        subsstart_jobid = removebadqc_jobid
        # Process the files
	# Collapse identical sequences in a FASTA/FASTQ file.
    prepcounts_jobid = prepcounts(samplename,fastq_suffix, subDIR=sub_dir, dep=subsstart_jobid, spike_seqs=spike_seqs)
    index += 1
	# mapping with bowtie
    mapbowtie_jobid = mapbowtie(samplename, refgenome=bowtie_index, dep=prepcounts_jobid, data_type=data_type)
    index += 1
	# Annotated mapped reads as isomiRs, mature miRNA, contaminants, tRNA, piRNA contaminants, read that map to stop oligo, not annotated
	# for details check the script match_bowtie_outputv5.pl
    annotate_jobid = annotate(samplename, dep=mapbowtie_jobid)
    index += 1
    # Gzip all fastq files
    zipfastq_jobid = zipfastq(dep=mapbowtie_jobid)
    print(f"Finished processing {samplename}")
if (subsampling == 'yes') & (subsample_to_nr == '0'):
    sys.exit("Determine the subsampling level after quality filtering")

