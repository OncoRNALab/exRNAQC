import sys
import os
import re
import argparse
import random
import subprocess

parser = argparse.ArgumentParser(description='Subsample miRNA vector')

# Read arguments
parser.add_argument('-b', nargs=1, required=True, help='Base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-s', nargs=1, required=True, help='Subdirectory where you can find miR count files', metavar='sub_dir')
parser.add_argument('-o', nargs=1, required=True, help='Directory where output should be created', metavar='output_dir')
parser.add_argument('-m', nargs=1, required=True, help='String match to select sample folders in base directory e.g. RNA0', metavar='string_match')
parser.add_argument('-n', nargs=1, required=True, help='Subsampling level', metavar='subs_level')
args = parser.parse_args()

base_dir = args.b[0].rstrip("/")
subs_level = args.n[0]
string_match = args.m[0]
sub_dir = args.s[0].rstrip("/")

def sample(iterable,n):
    """
    Returns @param n random items from @param iterable.
    Reservoir sampling cf. Knuth, 1981 
    """
    reservoir = []
    for t,item in enumerate(iterable):
        #print(t)
        #print(item)
        if t<n:
            reservoir.append(item)
        else:
            m=random.randint(0,t)
            if m<n:
                reservoir[m] = item
    return reservoir
# https://gregable.com/2007/10/reservoir-sampling.html

for samplename in os.listdir(base_dir): ##alternative approach if you want all samples with RNA in the name immediately
    if os.path.isdir(os.path.join(base_dir,samplename)):
        if re.search(string_match,samplename):
            print(samplename)
            output_dir = args.o[0].rstrip("/")
            input_dir = base_dir+"/"+samplename+"/"+sub_dir
            shortname=str.split(samplename,"-")[0]
            #open file
            input_file = input_dir+"/"+samplename+"_miRs.txt" #file with miR counts
            #duplicate MIMATIDs for the nr of times this MIMAT is detected (counts in miR counts file), this file will be removed afterwards
            tmp_file = output_dir+"/"+samplename+"_miRs_tmp.txt" 
            COMMAND = "awk '{{for(i=1;i<=$3;i++) print $2}}' OFS='\t' {0} > {1}".format(input_file, tmp_file) #new file (_miRs_tmp.txt), result of duplicating MIMATID times counts
            subprocess.call(COMMAND, shell=True)
            #open file with duplicated MIMATIDs and save to list
            duplMIMAT = open(tmp_file, 'r').readlines()
            
            tmp_file2 = output_dir+"/"+samplename+"_miRs_tmp2.txt"
            subsMIMAT = open(tmp_file2, 'w')
            subsMIMAT.writelines( '%s\n' % item.rstrip() for item in sample(duplMIMAT,int(subs_level))) #do reservoir sampling with the tmp list and write the outcome to new file (_miRs_tmp2.txt), which is again a file with duplicated MIMATIDs on each line
            subsMIMAT.close()
            
            output_file = output_dir+"/"+samplename+"_miRs_subs.txt"
            COMMAND = "rm {2}; awk 'NR==FNR {{dups[$1]++}} END {{for (mimat in dups) {{print mimat,dups[mimat]}} }}' {0} > {1} ; rm {0}".format(tmp_file2, output_file, tmp_file) #reduce the number of lines by counting how much lines you have of each MIMATID = new subs counts
            #COMMAND = "rm {2}; sort {0} | uniq -c > {1}".format(tmp_file2, output_file, tmp_file)
            subprocess.call(COMMAND, shell=True)
exit()
