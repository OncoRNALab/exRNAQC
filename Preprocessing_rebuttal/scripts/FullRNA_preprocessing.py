import os
import subprocess
import re
import sys
import yaml
import argparse

#defining functions

def read_file_to_tuples(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [tuple(line.strip().split()) for line in lines]

def fastqc(fastq1, fastq2, sufix, outDIR='.'):
    # Create the output directory if it doesn't exist
    output_path = os.path.join(outDIR, "FASTQC", F"FASTQC{sufix}")
    os.makedirs(output_path, exist_ok=True)
    # Build the command for FastQC
    command = ['fastqc', '-o', output_path, fastq1, fastq2]
    # Run the command
    print(command)
    subprocess.run(command, check=True)
    print(f"FastQC analysis complete. Results are in {output_path}")

def trimlength_local(sampleID, fastq1, fastq2, outDIR='.', suffix='', readlength=75):
    """ Remove trailing bases from reads and save the output in a 'trimming' subdirectory """
    # Create the output directory if it doesn't exist
    output_path = os.path.join(outDIR, 'trimming')
    os.makedirs(output_path, exist_ok=True)
    
    # Build the command for cutadapt
    trimmed_r1 = os.path.join(output_path, f'{sampleID}_1{suffix}_len.fastq')
    command = ['cutadapt', '-l', str(readlength), '--minimum-length=20', '-o', trimmed_r1]
    
    if fastq2:
        trimmed_r2 = os.path.join(output_path, f'{sampleID}_2{suffix}_len.fastq')
        command.extend(['-p', trimmed_r2, fastq1, fastq2])
    else:
        trimmed_r2 = None
        command.append(fastq1)
    
    # Run the command
    subprocess.run(command, check=True)
    
    # Compress the original fastq files
    subprocess.run(['gzip', fastq1], check=True)
    if fastq2:
        subprocess.run(['gzip', fastq2], check=True)
    
    print(f"Trimming complete. Results are in {output_path}")
    print(trimmed_r1, trimmed_r2)
    return trimmed_r1, trimmed_r2

def trimadapter_local(sampleID, fastq1, fastq2=None, outDIR='.', adapterR1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', adapterR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', picoV2='no', polyAtrim='no'):
    """ Remove adapters and trim reads """
    output_path = os.path.join(outDIR, 'cutadapt')
    os.makedirs(output_path, exist_ok=True)
    
    # Build the base command for cutadapt
    command = ['cutadapt', '-a', adapterR1, '--minimum-length=20']
    
    if picoV2 == 'yes':
        command.extend(['-U', '3'])
    
    if fastq2:
        command.extend(['-A', adapterR2, '-o', os.path.join(output_path, f'{sampleID}_1_adapter.fastq'), '-p', os.path.join(output_path, f'{sampleID}_2_adapter.fastq'), fastq1, fastq2])
    else:
        command.extend(['-o', os.path.join(output_path, f'{sampleID}_adapter.fastq'), fastq1])
    
    # Run the command
    subprocess.run(command, check=True)
    
    if polyAtrim == 'yes':
        polyA_command = ['cutadapt', '-a', 'T{100}', '--minimum-length=20']
        if fastq2:
            polyA_command.extend(['-A', 'A{100}', '-o', os.path.join(output_path, f'{sampleID}_1_adapterpA.fastq'), '-p', os.path.join(output_path, f'{sampleID}_2_adapterpA.fastq'), os.path.join(output_path, f'{sampleID}_1_adapter.fastq'), os.path.join(output_path, f'{sampleID}_2_adapter.fastq')])
        else:
            polyA_command.extend(['-o', os.path.join(output_path, f'{sampleID}_adapterpA.fastq'), os.path.join(output_path, f'{sampleID}_adapter.fastq')])
        
        # Run the polyA trimming command
        subprocess.run(polyA_command, check=True)
        
        # Rename the files back to original names and compress them
        if fastq2:
            os.rename(os.path.join(output_path, f'{sampleID}_1_adapterpA.fastq'), os.path.join(output_path, f'{sampleID}_1_trim.fastq'))
            os.rename(os.path.join(output_path, f'{sampleID}_2_adapterpA.fastq'), os.path.join(output_path, f'{sampleID}_2_trim.fastq'))
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_1_trim.fastq')], check=True)
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_2_trim.fastq')], check=True)
            r1 = os.path.join(output_path, f'{sampleID}_1_trim.fastq.gz')
            r2 = os.path.join(output_path, f'{sampleID}_2_trim.fastq.gz')
            return r1, r2
        else:
            os.rename(os.path.join(output_path, f'{sampleID}_adapterpA.fastq'), os.path.join(output_path, f'{sampleID}_trim.fastq'))
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_trim.fastq')], check=True)
            r1 = os.path.join(output_path, f'{sampleID}_trim.fastq.gz')
            return r1, None
    else:
        # Rename the files back to original names and compress them
        if fastq2:
            os.rename(os.path.join(output_path, f'{sampleID}_1_adapter.fastq'), os.path.join(output_path, f'{sampleID}_1_trim.fastq'))
            os.rename(os.path.join(output_path, f'{sampleID}_2_adapter.fastq'), os.path.join(output_path, f'{sampleID}_2_trim.fastq'))
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_1_trim.fastq')], check=True)
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_2_trim.fastq')], check=True)
            r1 = os.path.join(output_path, f'{sampleID}_1_trim.fastq.gz')
            r2 = os.path.join(output_path, f'{sampleID}_2_trim.fastq.gz')
            return r1, r2
        else:
            os.rename(os.path.join(output_path, f'{sampleID}_adapter.fastq'), os.path.join(output_path, f'{sampleID}_trim.fastq'))
            subprocess.run(['gzip', os.path.join(output_path, f'{sampleID}_trim.fastq')], check=True)
            r1 = os.path.join(output_path, f'{sampleID}_trim.fastq.gz')
            return r1, None

def removebadqc(sampleID, fastq1, fastq2, outDIR='.', suffix="", percentage=80, quality=19):
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
    output_path = os.path.join(outDIR, 'trimming')
    os.makedirs(output_path, exist_ok=True)
    print(f'Quality Filtering: {fastq1}{fastq2}')

    if fastq2:
        # Paired-end processing
        temp_fastq1 = os.path.join(output_path, f'{sampleID}_temp_1.fastq')
        temp_fastq2 = os.path.join(output_path, f'{sampleID}_temp_2.fastq')
        qc_fastq1 = os.path.join(output_path, f'{sampleID}_1{suffix}_qc.fastq')
        qc_fastq2 = os.path.join(output_path, f'{sampleID}_2{suffix}_qc.fastq')

        # Split the fastq files into smaller chunks
        subprocess.run(['split', '-l', '10000000', '--additional-suffix=.fastq', fastq1, os.path.join(output_path, 'temp_1_')], check=True)
        subprocess.run(['split', '-l', '10000000', '--additional-suffix=.fastq', fastq2, os.path.join(output_path, 'temp_2_')], check=True)

        # Remove the original fastq files
        os.remove(fastq1)
        os.remove(fastq2)

        # Process each chunk
        for suffix in [f.split('_')[-1].split('.')[0] for f in os.listdir(output_path) if f.startswith('temp_1_')]:
            temp_chunk1 = os.path.join(output_path, f'temp_1_{suffix}.fastq')
            temp_chunk2 = os.path.join(output_path, f'temp_2_{suffix}.fastq')
            qc_chunk1 = os.path.join(output_path, f'temp_qc_1_{suffix}.fastq')
            qc_chunk2 = os.path.join(output_path, f'temp_qc_2_{suffix}.fastq')

            # Run FASTQC_filter.py on each chunk
            try:
                subprocess.run(['python3', 'FASTQC_filter.py', '-r1', temp_chunk1, '-r2', temp_chunk2, '-p', str(percentage), '-q', str(quality), '-o1', qc_chunk1, '-o2', qc_chunk2], check=True)
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while running FASTQC_filter.py on chunk {suffix}: {e}")

            # Remove the temporary chunks
            os.remove(temp_chunk1)
            os.remove(temp_chunk2)

        # Concatenate the processed chunks
        with open(qc_fastq1, 'wb') as f_out:
            for chunk in [f for f in os.listdir(output_path) if f.startswith('temp_qc_1_')]:
                with open(os.path.join(output_path, chunk), 'rb') as f_in:
                    f_out.write(f_in.read())
                os.remove(os.path.join(output_path, chunk))

        with open(qc_fastq2, 'wb') as f_out:
            for chunk in [f for f in os.listdir(output_path) if f.startswith('temp_qc_2_')]:
                with open(os.path.join(output_path, chunk), 'rb') as f_in:
                    f_out.write(f_in.read())
                os.remove(os.path.join(output_path, chunk))

        # Compress the final processed files
        subprocess.run(['gzip', qc_fastq1], check=True)
        subprocess.run(['gzip', qc_fastq2], check=True)

        return f'{qc_fastq1}.gz', f'{qc_fastq2}.gz'
    
def subsample(sampleID, fastq1, fastq2=None, nreads=1000, suffix='', output_dir=".", dedup_mode='clumpify', data_type='pe'):
    """ Seqtk subsampling
    Downsample fastq files to x number of reads (take min of all samples)
    -s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of code
    """
    output_path = os.path.join(output_dir, dedup_mode)
    os.makedirs(output_path, exist_ok=True)
    
    if data_type == 'se':
        command = f'seqtk sample -s 100 {fastq1} {nreads} | gzip > {output_path}/{sampleID}{suffix}_subs.fastq.gz'
        subprocess.run(command, shell=True, check=True)
        r1 = os.path.join(output_path, f'{sampleID}{suffix}_subs.fastq.gz')
        return os.path.abspath(r1), None
    elif data_type == 'pe':
        command1 = f'seqtk sample -s 100 {fastq1} {nreads} | gzip > {output_path}/{sampleID}_1{suffix}_subs.fastq.gz'
        command2 = f'seqtk sample -s 100 {fastq2} {nreads} | gzip > {output_path}/{sampleID}_2{suffix}_subs.fastq.gz'
        subprocess.run(command1, shell=True, check=True)
        subprocess.run(command2, shell=True, check=True)
        r1 = os.path.join(output_path, f'{sampleID}_1{suffix}_subs.fastq.gz')
        r2 = os.path.join(output_path, f'{sampleID}_2{suffix}_subs.fastq.gz')
        return os.path.abspath(r1), os.path.abspath(r2)

def subsample_v2(sampleID, fastq1, fastq2=None, nreads=1000, suffix='', output_dir=".", dedup_mode='clumpify', data_type='pe'):
    """ Seqtk subsampling
    Downsample fastq files to x number of reads (take min of all samples)
    -s100: sets the seed to a fixed number (here 100, but may be any number) in order to have the same result if you rerun this part of code
    """
    output_path = os.path.join(output_dir, dedup_mode)
    os.makedirs(output_path, exist_ok=True)
    
    if data_type == 'se':
        command = f'seqtk sample -s 100 {fastq1} {nreads} > {output_path}/{sampleID}{suffix}_subs.fastq'
        subprocess.run(command, shell=True, check=True)
        r1 = os.path.join(output_path, f'{sampleID}{suffix}_subs.fastq')
        return os.path.abspath(r1), None
    elif data_type == 'pe':
        command1 = f'seqtk sample -s 100 {fastq1} {nreads} > {output_path}/{sampleID}_1{suffix}_subs.fastq'
        command2 = f'seqtk sample -s 100 {fastq2} {nreads} > {output_path}/{sampleID}_2{suffix}_subs.fastq'
        subprocess.run(command1, shell=True, check=True)
        subprocess.run(command2, shell=True, check=True)
        r1 = os.path.join(output_path, f'{sampleID}_1{suffix}_subs.fastq')
        r2 = os.path.join(output_path, f'{sampleID}_2{suffix}_subs.fastq')
        return os.path.abspath(r1), os.path.abspath(r2)
    


def clumpifydupremoval_trimmed(sampleID, fastq1, fastq2=None, outDIR='.', subst='2', kmersize='31', passes='20', dedup_mode='clumpify', data_type='pe'):
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
    output_path = os.path.join(outDIR, dedup_mode, f"{sampleID}_clumpout")
    output_path_gzip_files = os.path.join(outDIR, dedup_mode)
    os.makedirs(output_path, exist_ok=True)
    
    if data_type == 'se':
        # Single-end processing
        trimmed_fastq = os.path.join(output_path, f'{sampleID}_temptrim.fastq')
        clumped_fastq = os.path.join(output_path, f'{sampleID}_clumped.fastq')
        
        # Trim reads
        subprocess.run(['cutadapt', '-l', '80', '-o', trimmed_fastq, fastq1], check=True)
        subprocess.run(['gzip', trimmed_fastq], check=True)
        
        # Clumpify
        subprocess.run(['clumpify.sh', f'in={trimmed_fastq}.gz', f'out={clumped_fastq}.gz', 'dedupe', f'subs={subst}', f'k={kmersize}', f'passes={passes}'], check=True)
        subprocess.run(['gunzip', f'{clumped_fastq}.gz'], check=True)
        
        # Extract original reads
        tempnames = os.path.join(output_path, f'{sampleID}_tempnames.txt')
        with open(tempnames, 'w') as f:
            subprocess.run(['grep', '^@', clumped_fastq], stdout=f, check=True)
        
        final_fastq = os.path.join(output_path, f'{sampleID}_final.fastq')
        with open(final_fastq, 'w') as f:
            subprocess.run(['grep', '-Ff', tempnames, '-A3', fastq1], stdout=f, check=True)
        
        subprocess.run(['gzip', final_fastq], check=True)
        return f'{final_fastq}.gz', None
    
    elif data_type == 'pe':
        # Paired-end processing
        trimmed_fastq1 = os.path.join(output_path, f'{sampleID}_temptrim_1.fastq')
        trimmed_fastq2 = os.path.join(output_path, f'{sampleID}_temptrim_2.fastq')
        tempclumped_fastq1 = os.path.join(output_path, f'{sampleID}_tempclumped_1.fastq')
        tempclumped_fastq2 = os.path.join(output_path, f'{sampleID}_tempclumped_2.fastq')
        clumped_fastq1 = os.path.join(output_path_gzip_files, f'{sampleID}_clumped_1.fastq')
        clumped_fastq2 = os.path.join(output_path_gzip_files, f'{sampleID}_clumped_2.fastq')
        
        # Trim reads
        subprocess.run(['cutadapt', '--pair-filter=both', '-l', '80', '-o', trimmed_fastq1, '-p', trimmed_fastq2, fastq1, fastq2], check=True)
        subprocess.run(['gzip', trimmed_fastq1], check=True)
        subprocess.run(['gzip', trimmed_fastq2], check=True)
        
        # Clumpify
        subprocess.run(['clumpify.sh', f'in={trimmed_fastq1}.gz', f'in2={trimmed_fastq2}.gz', f'out={tempclumped_fastq1}.gz', f'out2={tempclumped_fastq2}.gz', 'dedupe', f'subs={subst}', f'k={kmersize}', f'passes={passes}'], check=True)
        subprocess.run(['gunzip', f'{tempclumped_fastq1}.gz'], check=True)
        subprocess.run(['gunzip', f'{tempclumped_fastq2}.gz'], check=True)
        
        # Extract original reads
        tempnames1 = os.path.join(output_path, f'{sampleID}_tempnames_1.txt')
        tempnames2 = os.path.join(output_path, f'{sampleID}_tempnames_2.txt')
        with open(tempnames1, 'w') as f:
            subprocess.run(f'grep ^@ {os.path.abspath(tempclumped_fastq1)}', shell=True, stdout=f, check=True)
        with open(tempnames2, 'w') as f:
            subprocess.run(f'grep ^@ {os.path.abspath(tempclumped_fastq2)}', shell=True, stdout=f, check=True)
        
        subprocess.run(['gunzip', f'{trimmed_fastq1}.gz'], check=True)
        subprocess.run(['gunzip', f'{trimmed_fastq2}.gz'], check=True)

        with open(clumped_fastq1, 'w') as f:
            subprocess.run(f'grep -Ff {os.path.abspath(tempnames1)} -A3 {os.path.abspath(trimmed_fastq1)} | grep -v -- "^--$"', shell=True, stdout=f, check=True)
        with open(clumped_fastq2, 'w') as f:
            subprocess.run(f'grep -Ff {os.path.abspath(tempnames2)} -A3 {os.path.abspath(trimmed_fastq2)} | grep -v -- "^--$"', shell=True, stdout=f, check=True)
        
        temp_files = f'rm -r {output_path}'
        subprocess.run(temp_files, shell=True, check=True)

        subprocess.run(['gzip', clumped_fastq1], check=True)
        subprocess.run(['gzip', clumped_fastq2], check=True)

        return f'{clumped_fastq1}.gz', f'{clumped_fastq2}.gz'

def bamtofastqfile(sampleID, bamIN, output_dir, dedup_mode, data_type):
    """ Convert bam or sam file to fastq with samtools fastq
    -1: write reads with READ1 FLAG set (and not READ2) to file (SE, or mate 1 for PE); -2: write reads with READ2 FLAG set (and not READ1) to file 2 (mate 2 for PE)
    -0 /dev/null: write reads where both READ1 and READ2 are set; or neither is set to file (/dev/null = remove them) (reads that do not fit in simple paired-end sequencing model)
    -s /dev/null: write singleton reads to file (/dev/null = remove them)
    -n: use read names as they are (instead of adding '/1' and '/2')
    -F 0x900: Do not output alignments with any bits set in INT present in the FLAG field
    0x900 = 0x100 (not primary alignment) & 0x800 (supplementary alignment)
    """
    output_path = os.path.join(output_dir, dedup_mode)
    os.makedirs(output_path, exist_ok=True)
    
    if data_type == 'se':
        fastq_out = os.path.join(output_path, f'{sampleID}_picard.fastq')
        command = [
            'samtools', 'fastq', '-1', fastq_out, '-0', '/dev/null', '-n', '-F', '0x900', bamIN
        ]
    elif data_type == 'pe':
        fastq_out1 = os.path.join(output_path, f'{sampleID}_1_picard.fastq')
        fastq_out2 = os.path.join(output_path, f'{sampleID}_2_picard.fastq')
        command = [
            'samtools', 'fastq', '-1', fastq_out1, '-2', fastq_out2, '-0', '/dev/null', '-n', '-F', '0x900', bamIN
        ]
    
    subprocess.run(command, check=True)
    return fastq_out if data_type == 'se' else (fastq_out1, fastq_out2)

def countskallisto_local(sampleID, fastq1, fastq2, output_dir, dedup_mode, data_type, kal_index, fragmentmean=None, fragmentsd=None, stranded='no'):
    """ kallisto quant to run quantification algorithm (paired end mode)
    -t: number of threads to use (default 1)
    -i: index (used for quantification)
    --rf-stranded: strand specific mode: only fragments where the first read in the pair pseudoaligns to the reverse strand of a transcript are processed (vice versa in --fr-stranded)
    -o: output directory
    --single: SE sequencing (fragment length mean and sd must also be supplied for single-end reads using -l and -s)
    Produces 3 output files by default: abundances.tsv (plaintext file of abundance estimates), abundances.h5, run_info.json
    """
    output_path = os.path.join(output_dir, dedup_mode, f"{sampleID}_klout")
    os.makedirs(output_path, exist_ok=True)

    if stranded == 'yes':
        if data_type == 'se':
            command = [
                'kallisto', 'quant', '-i', kal_index, '--rf-stranded', '--single',
                '-l', str(fragmentmean), '-s', str(fragmentsd), '-o', output_path, fastq1
            ]
        elif data_type == 'pe':
            command = [
                'kallisto', 'quant', '-i', kal_index, '--rf-stranded', '-o', output_path, fastq1, fastq2
            ]
    else:  # unstranded
        if data_type == 'se':
            command = [
                'kallisto', 'quant', '-i', kal_index, '--single', '-l', str(fragmentmean),
                '-s', str(fragmentsd), '-o', output_path, fastq1
            ]
        elif data_type == 'pe':
            command = [
                'kallisto', 'quant',  '-i', kal_index, '-o', output_path, fastq1, fastq2
            ]

    subprocess.run(command, check=True)
    return output_path

def alignment(sampleID, fastq1, fastq2, output_dir, dedup_mode, star_index, gtf, data_type):
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
    output_path = os.path.join(output_dir, dedup_mode, f"{sampleID}_srout")
    os.makedirs(output_path, exist_ok=True)

    # Unzip the input fastq files if they are gzipped
    if fastq1.endswith('.gz'):
        subprocess.run(['gunzip', fastq1], check=True)
        fastq1 = fastq1[:-3]
    if fastq2 and fastq2.endswith('.gz'):
        subprocess.run(['gunzip', fastq2], check=True)
        fastq2 = fastq2[:-3]

    if data_type == 'se':
        command = [
            'STAR', '--runThreadN 10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', fastq1,
            '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outReadsUnmapped', 'Fastx', '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random',
            '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10', '--outSAMprimaryFlag', 'AllBestScore',
            '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
        ]
    elif data_type == 'pe':
        command = [
            'STAR', '--runThreadN 10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', fastq1, fastq2,
            '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outReadsUnmapped', 'Fastx', '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random',
            '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10', '--outSAMprimaryFlag', 'AllBestScore',
            '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
        ]

    subprocess.run(command, check=True)
    subprocess.run(['gzip', f'{output_path}/*Unmapped.*'], shell=True)
    return output_path
def load_star_index(star_index):
    """ Load STAR index into memory """
    command = ['STAR', '--genomeLoad', 'LoadAndExit', '--genomeDir', star_index]
    subprocess.run(command, check=True)

def unload_star_index(star_index):
    """ Unload STAR index from memory """
    command = ['STAR', '--genomeLoad', 'Remove', '--genomeDir', star_index]
    subprocess.run(command, check=True)

def strandedness(sampleID, bamIN, output_dir, dedup_mode, exon_bed, data_type, star_index, gtf):
    """ RSeQC to retrieve % correct strandedness.
    Grep the line that shows what fraction of reads is explained by fr-firststrand (reverse in htseq)
    (1+-,1-+,2++,2-- category, e.g. 1+- read 1 '+' mapped to + strand while gene is on '-' strand is what we expect for reverse stranded)
    """
    output_path = os.path.join(output_dir, dedup_mode)
    os.makedirs(output_path, exist_ok=True)

    if dedup_mode == 'clumpify':
        # run STAR first on FASTQ without deduplication
        if data_type == 'se':
            command = [
                'STAR', '--runThreadN', '10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', f'{output_path}/{sampleID}.fastq',
                '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outReadsUnmapped', 'Fastx',
                '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random', '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10',
                '--outSAMprimaryFlag', 'AllBestScore', '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
            ]
        elif data_type == 'pe':
            command = [
                'STAR', '--runThreadN', '10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', f'{output_path}/{sampleID}_1.fastq', f'{output_path}/{sampleID}_2.fastq',
                '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outReadsUnmapped', 'Fastx',
                '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random', '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10',
                '--outSAMprimaryFlag', 'AllBestScore', '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
            ]
        subprocess.run(command, check=True)
        subprocess.run(['rm', f'{output_path}/Unmapped.*'], check=True)

    command = [
        'infer_experiment.py', '-r', exon_bed, '-i', bamIN
    ]
    output_file = os.path.join(output_path, f'{sampleID}_RSeQC_output_all.txt')
    with open(output_file, 'w') as f:
        subprocess.run(command, stdout=f, check=True)
    percentage = "0.0"
    with open(output_file, 'r') as f:
        for line in f:
            if "1+-" in line:
                percentage = line.split(":")[1].strip()
                break

    final_output = os.path.join(output_path, f'{sampleID}_RSeQC_output.txt')
    with open(final_output, 'w') as f:
        f.write(f'{sampleID} {percentage}\n')

    return final_output

def star_alignment_strandness(sampleID, fastq1, fastq2, output_dir, dedup_mode, data_type, star_index, gtf):
        """ Run STAR alignment and return the path to the resulting BAM file """
        output_path = os.path.join(output_dir, dedup_mode)
        os.makedirs(output_path, exist_ok=True)

        if data_type == 'se':
            command = [
                'STAR', '--runThreadN', '10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', fastq1,
                '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outReadsUnmapped', 'Fastx',
                '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random', '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10',
                '--outSAMprimaryFlag', 'AllBestScore', '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
            ]
        elif data_type == 'pe':
            command = [
                'STAR', '--runThreadN', '10', '--outFileNamePrefix', f'{output_path}/{sampleID}.', '--readFilesIn', fastq1, fastq2,
                '--genomeDir', star_index, '--sjdbGTFfile', gtf, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outReadsUnmapped', 'Fastx',
                '--twopassMode', 'Basic', '--outMultimapperOrder', 'Random', '--outSAMmultNmax', '-1', '--outFilterMultimapNmax', '10',
                '--outSAMprimaryFlag', 'AllBestScore', '--outFilterScoreMinOverLread', '0.66', '--outFilterMatchNminOverLread', '0.66', '--outFilterMatchNmin', '20'
            ]
        subprocess.run(command, check=True)
        #removing unmmaped and Log files
        unmmaped = f'rm {os.path.join(output_path, f"{sampleID}.Unmapped*")}'
        subprocess.run(unmmaped, shell=True, check=True)
        log = f'rm {os.path.join(output_path, f"{sampleID}.Log*")}'
        subprocess.run(log, shell=True, check=True)
        extra = f'rm {os.path.join(output_path, f"{sampleID}.SJ*")}'
        subprocess.run(extra, shell=True, check=True)
        bam_file = os.path.join(output_path, f'{sampleID}.Aligned.sortedByCoord.out.bam')
        return bam_file

def run_infer_experiment(sampleID, bamIN, output_dir, dedup_mode, exon_bed):
    """ Run infer_experiment.py to get strandedness and return the path to the output file."""
    output_path = os.path.join(output_dir, dedup_mode)
    os.makedirs(output_path, exist_ok=True)

    command = [
        'infer_experiment.py', '-r', exon_bed, '-i', bamIN
    ]
    output_file = os.path.join(output_path, f'{sampleID}_RSeQC_output_all.txt')
    with open(output_file, 'w') as f:
        subprocess.run(command, stdout=f, check=True)
    percentage = "0.0"
    with open(output_file, 'r') as f:
        for line in f:
            if "1+-" in line:
                percentage = line.split(":")[1].strip()
                break

    final_output = os.path.join(output_path, f'{sampleID}_RSeQC_output.txt')
    with open(final_output, 'w') as f:
        f.write(f'{sampleID} {percentage}\n')

    return final_output

def sortbambyname(sampleID, bamIN, bamOUT, output_dir, dedup_mode):
    """ Some algorithms only run on name sorted bam files instead of coordinate sorted ones, this function makes the conversion from coo to name sorted.
    -n: sort by read names instead of chromosomal coordinates
    -o: output (sam, bam, cram format is deduced from filename extension -> make sure it ends on .bam for bam file
    """
    output_path = os.path.join(output_dir, dedup_mode, f"{sampleID}_srout")
    os.makedirs(output_path, exist_ok=True)
    
    command = [
        'samtools', 'sort', '-n', '-o', os.path.join(output_path, bamOUT), bamIN
    ]
    subprocess.run(command, check=True)
    return os.path.join(output_path, bamOUT)

def picarddupremcoord(sampleID, bamIN='Aligned.sortedByCoord.out.bam', bamOUT='Aligned.sortedByCoord.picard.bam', output_dir='.', dedup_mode='clumpify'):
    """ Duplicate removal with Picard
    (When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates)
    ASSUME_SORT_ORDER=coordinate
    REMOVE_DUPLICATES=true: remove all duplicates (optical and sequencing)
    VALIDATION_STRINGENCY=SILENT: (default: STRICT) Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags)
    M: file to write duplication metrics to
    """

    output_path = os.path.join(output_dir, dedup_mode, f"{sampleID}_srout")
    output_path2 = os.path.join(output_dir, dedup_mode, f"{sampleID}_picardout")
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path2, exist_ok=True)
    
    bamOUT_path = os.path.join(output_path, bamOUT)
    metrics_file = os.path.join(output_path2, f'{sampleID}_picard_dup.metrics')
    
    command = [
        'java', '-jar', '/opt/picard/picard.jar', 'MarkDuplicates',
        f'I={bamIN}', f'O={bamOUT_path}',
        'ASSUME_SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=SILENT',
        'REMOVE_DUPLICATES=true', f'M={metrics_file}'
    ]
    subprocess.run(command, check=True)
    return bamOUT_path

def countshtseq(sampleID, bamIN, gtf, output_dir, dedup_mode, stranded='no'):
    """ HTSeq quantification of STAR mapped bam files.
    --order name: (default) needs name sorted bam files
    --nonunique none: (default) if the union of features for each position in the read is > 1, the read (pair) is counted as ambiguous and not counted for any features +
    if read (pair) aligns to more than one location in reference, it is scored as alignment_not_unique (for each location)
    --stranded reverse: for PE reads, the second read has to be on same strand and the first read has to be on opposite strand (for stranded=yes it is vice versa)
    """
    output_path = os.path.join(output_dir, dedup_mode, f"{sampleID}_htout")
    os.makedirs(output_path, exist_ok=True)
    
    htseq_output = os.path.join(output_path, f'{os.path.basename(bamIN).replace(".bam", "")}_htseq_counts.txt')
    
    if stranded == 'yes':
        command = [
            'htseq-count', '--format', 'bam', '--order', 'name', '--nonunique', 'none', '--stranded', 'reverse',
            bamIN, gtf
        ]
    else:
        command = [
            'htseq-count', '--format', 'bam', '--order', 'name', '--nonunique', 'none', '--stranded', 'no',
            bamIN, gtf
        ]
    
    with open(htseq_output, 'w') as f:
        subprocess.run(command, stdout=f, check=True)
    
    return htseq_output

def idxstat(sampleID, bamIN, output_dir, dedup_mode):
    """ Perform samtools idxstats
    Retrieve and print stats in the index file corresponding to the input file. The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
    Before calling idxstats, the input BAM file should be indexed by samtools index.
    """
    output_path = os.path.join(output_dir, dedup_mode, f'{sampleID}_idxout')
    os.makedirs(output_path, exist_ok=True)
    
    # Index the BAM file
    subprocess.run(['samtools', 'index', bamIN], check=True)
    
    # Run idxstats and save the output
    idxstat_output = os.path.join(output_path, f'{os.path.basename(bamIN).replace(".bam", "")}_idxstat.txt')
    with open(idxstat_output, 'w') as f:
        subprocess.run(['samtools', 'idxstats', bamIN], stdout=f, check=True)
    
    return idxstat_output

def main():
    #Load Parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',required=True, help='path to yaml config file')
    yaml_file = parser.parse_args()
    with open(yaml_file.config, "r") as fp:
        args = yaml.safe_load(fp)
    output_dir = args["output_dir"]
    data_type = args["data_type"]
    deduplication=args['deduplication']
    subs_timing=args['subs_timing']
    subsample_to_nr = args['subsample_to_nr']
    kal_index = args["kal_index"]
    fragmean=args["fragmean"]
    fragsd=args["fragsd"]
    star_index = args["star_index"]
    gtf = args["gtf"]
    adaptercontam=args["adaptercontam"]
    exon_bed=args["exon_bed"]
    repeat_analysis=args["repeat_analysis"]
    SampleList = args["SampleListFile"]
    # Load the STAR index into memory before processing
    load_star_index(star_index)
    paired_files = read_file_to_tuples(SampleList)    

    for sampleID, r1, r2 in paired_files:        
        print(f"Processing {sampleID}: {r1} and {r2}")
        output_dir_sample = os.path.join(output_dir, sampleID)
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
        # fastq_suffix = 'FASTQC'
        # First step QC and trimming
        if repeat_analysis == 'no':
            fastqc(r1, r2, '_original', output_dir_sample)
            if adaptercontam != 'no':
                if adaptercontam == 'truseq':
                    adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
                    adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
                elif adaptercontam == 'nextera':
                    adapter1 = 'CTGTCTCTTATACACATCT' # For Nextera/TruSight, the same sequence is used for both reads (cf Illumina website)
                    adapter2 = 'CTGTCTCTTATACACATCT'
                r1, r2 = trimadapter_local(sampleID, r1, r2, output_dir_sample, adapterR1=adapter1, adapterR2=adapter2, picoV2='no', polyAtrim='no')
                fastq_suffix = '_trim'

            trimL1, trimL2 = trimlength_local(sampleID, r1, r2, output_dir_sample, suffix=fastq_suffix, readlength=75)
            fastq_suffix += '_len'
            if (data_type == 'pe'):
                trim1, trim2 = removebadqc(sampleID, trimL1, trimL2, outDIR=output_dir_sample, suffix=fastq_suffix, percentage=80, quality=19)
                fastq_suffix += "_qc"
            else:
                trim1, trim2 = trimL1, trimL2
            fastqc(trim1, trim2, fastq_suffix, output_dir_sample)
            continue
        else: 
            if (data_type == 'pe'):
                trim1 = os.path.join(output_dir_sample, "trimming", f'{sampleID}_1_trim_len_qc.fastq.gz') 
                trim2 = os.path.join(output_dir_sample, "trimming", f'{sampleID}_2_trim_len_qc.fastq.gz') 
                fastq_suffix = '_trim_len_qc'
            else:
                trim1 = os.path.join(output_dir_sample, "trimming", f'{sampleID}_1_trim_len.fastq.gz') 
                fastq_suffix = '_trim_len'

        if (data_type == 'se') & ((fragmean == '-1') | (fragsd =='-1')): #if SE and options for fragment length mean and sd were not changed:
            sys.exit("ERROR: For single end sequencing, kallisto requires a fragment length mean and sd")     
        # # # Step 3 subsample
        if subs_timing == 'start':
            if subsample_to_nr == '0':
                continue # Skip subsampling
            else:
                trim1, trim2 = subsample_v2(sampleID, trim1, trim2, int(subsample_to_nr), fastq_suffix, output_dir_sample, dedup_mode, data_type)
                print(f"Subsampled to {subsample_to_nr} reads, files: {trim1} and {trim2}")
                fastq_suffix += '_subs'
        
        bam_coord = sampleID+'.Aligned.sortedByCoord.out.bam'
        bam_name = sampleID+'.Aligned.sortedByName.out.bam'
        # Step 4 deduplication
        if deduplication == 'clumpify':
                # get strandedness
                bamIN = star_alignment_strandness(sampleID, trim1, trim2, output_dir_sample, dedup_mode, data_type, star_index, gtf)
                run_infer_experiment(sampleID, bamIN, output_dir_sample, dedup_mode, exon_bed)
                #clumpify duplicate removal
                r1_dep, r2_dep = clumpifydupremoval_trimmed(sampleID, trim1, trim2, outDIR=output_dir_sample, dedup_mode=dedup_mode, data_type=data_type)
                print(f"Clumpify deduplication complete. Files: {r1_dep} and {r2_dep}")
                fastq_suffix += '_clumped'
                # Run kallisto quantification
                if subs_timing == 'afterclumpify':
                    if subsample_to_nr == '0':
                        continue #go to next samplename
						#sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/dedup_clumpify-subs_afterdedup/{sample}_1_qc_clumped.fastq; done > lines_qcfilt_fastq.txt")
						# subsample to minimum based on wc -l divided by 4 and floored to a million
                    else:
                        r1_dep, r2_dep = subsample(sampleID, r1_dep, r1_dep, int(subsample_to_nr), fastq_suffix, output_dir_sample, dedup_mode, data_type)
                        fastq_suffix += '_subs'
                fastqc(r1_dep, r2_dep, fastq_suffix, output_dir_sample)
                kal_out = countskallisto_local(sampleID, r1_dep, r2_dep, output_dir_sample, dedup_mode, data_type, kal_index, fragmentmean=fragmean, fragmentsd=fragsd, stranded='no')
                print(f'kallisto: {kal_out}')
                # STAR alignment
                star = alignment(sampleID, r1_dep, r2_dep, output_dir_sample, dedup_mode, star_index, gtf, data_type)
                bam_coord = os.path.join(star, bam_coord)
                print(f'STAR alignment: {star}')
                # sort bam file by name
                bam_name = sortbambyname(sampleID, bamIN=bam_coord, bamOUT=bam_name, output_dir=output_dir_sample, dedup_mode=dedup_mode)
                print(f'Sorted BAM file: {bam_name}')
        else:
            # star alignment
            star = alignment(sampleID, trim1, trim2, output_dir_sample, dedup_mode, star_index, gtf, data_type)
            bam_IN = os.path.join(star, bam_coord)
            # get stradedness
            strandedness(sampleID, bam_IN, output_dir_sample, dedup_mode, exon_bed, data_type, star_index, gtf)
            if deduplication == 'picard':
                bam_coord_nodedup = os.path.join(star, bam_coord)
                # Deduplication with Picard
                bam_coord = sampleID+'.Aligned.sortedByCoord.picard.bam'
                star_picard = picarddupremcoord(sampleID, bamIN=bam_coord_nodedup, bamOUT=bam_coord, output_dir=output_dir_sample, dedup_mode=dedup_mode)
                bam_coord = star_picard
                bam_name = sampleID+'.Aligned.sortedByName.picard.bam'
                bam_name = sortbambyname(sampleID, bamIN=star_picard, bamOUT=bam_name, output_dir=output_dir_sample, dedup_mode=dedup_mode)
                star_fq1, star_fq2 = bamtofastqfile(sampleID, bam_name, output_dir_sample, dedup_mode, data_type)
                fastq_suffix = '_picard'
                # Kallisto quantification
                kal_out = countskallisto_local(sampleID, star_fq1, star_fq2, output_dir_sample, dedup_mode, data_type, kal_index, fragmentmean=fragmean, fragmentsd=fragsd, stranded='no')
            else:
                star_sort = sortbambyname(bamIN=star_picard, bamOUT=bam_name, output_dir=output_dir_sample, dedup_mode=dedup_mode)
    
    #     # HTSeq quantification (bam needs to be sorted by name)
        countshtseq(sampleID, bam_name, gtf, output_dir_sample, dedup_mode, stranded='no')
        idxstat(sampleID, bam_coord, output_dir_sample, dedup_mode)
        print(f"Finished processing {r1} and {r2}")
    # # Unload the STAR index from memory after processing
    unload_star_index(star_index)
    

    if (subs_timing == 'start') & (subsample_to_nr == '0'):
        sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/${sample}_1_qc.fastq; done > lines_qcfilt_fastq.txt")

    elif (subs_timing == 'afterclumpify') & (subsample_to_nr == '0'):
        sys.exit("Determine the subsampling level after quality filtering \n(i.e. floored minimum of wc-l fastq divided by 4)\n   for sample in $(ls); do wc -l ${sample}/dedup_clumpify-subs_afterdedup/{sample}_1_qc_clumped.fastq; done > lines_qcfilt_fastq.txt")


if __name__ == "__main__":
    main()
