#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from sys import argv
import re

infile = open(argv[1], 'r',encoding='utf-8-sig')
spikes = pd.read_csv(infile)
#print(spikes)

fasta_sequences = SeqIO.parse(open(argv[2]),'fasta')

all_spikes=pd.DataFrame()
fasta_records=[]
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    match_spikes= spikes[spikes["sequence"].str.contains(sequence)].copy()
    if not match_spikes.empty:
        match_spikes['count'] = int(fasta.id.split("-")[1])
        print(sequence)
        match_spikes['sequence'] = str(sequence)
        all_spikes = all_spikes.append(match_spikes)
    else:
        fasta_records.append(fasta)

output_fasta = str(argv[3]) + "_collapsed_nospikes.fa"
SeqIO.write(fasta_records, output_fasta, "fasta-2line") #"fasta" cuts the read sequence after 80(?) characters, "fasta-2line" keeps fasta record on 2 lines (1 for name, 1 for sequence)

if not all_spikes.empty:
    sum_all_spikes = all_spikes.groupby(['spikeID']).sum()
    output_spikes = str(argv[3]) + "_spikes.txt"
    sum_all_spikes.to_csv(output_spikes,sep='\t')
    output_spikes = str(argv[3]) + "_allspikes.txt"
    all_spikes.to_csv(output_spikes,sep='\t')
else:
    output_spikes = str(argv[3]) + "_spikes.txt"
    open(output_spikes,"w").close()
    output_spikes = str(argv[3]) + "_allspikes.txt"
    open(output_spikes,"w").close()

#	open(output_spikes, 'wp' as fp:
#	fp.write("No spikes found")
#    output_spikes = str(argv[3]) + "_allspikes.txt"
