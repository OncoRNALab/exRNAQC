
# Advanced data analyses phase 1/2
exRNAQC020 encompasses all extra analysis done in the context of phase 1 and phase 2 of the exRNAQC study. 

## Fraction of circRNAs in plasma
The unique structure of circRNAs makes them resistant agains exonucleases and are therefore hypothesized to remain relatively stable in a RNAse-rich environment, such as plasma. We investigated the differences between tubes in how the fraction of circRNAs changes over time. 

### circRNA detection from mRNA capture data
First, we need to re-analyse our sequencing data and look for reads mapping to circRNAs. We used our in-house circRNA detection pipeline (which can be found [here](https://github.ugent.be/vandesompelelab/circRNA)) to detect circRNA reads in the mRNA capture sequencing data for the tested tubes and incubation times. We used the following command to run the pipeline:

```
submit_circRNA_slurm.py -p CIRCexplorer2_TopHat -t pe -d clumpify -b base_dir -g hg38 -o output_dir -s string_match -u user_email
```

This command institutes the pipeline to use TopHat for mapping, clumpify for deduplication and CIRCexplorer2 for circRNA detection. 

### Determination of unique circRNA and linear RNA reads
Next, we need to determine which reads are mapping to the backsplice junctions of circRNAs and to linear-unique splice junctions. Using these reads, we can determine the total fraction of circRNA molecules relative to linear RNAs. To facilitate this, we use CiLiQuant (which is described in a recently uploaded [preprint](https://github.com/OncoRNALab/CiLiQuant). For each sample, we used the following command to run CiLiQuant: 

```
tail -n +2 CIRCexplorer2_TopHat_clumpify/04_thfusionout/junctions.bed > junctions_rmline.bed
python CiLiQuant.py -b CIRCexplorer2_TopHat_clumpify/CIRCexplorer2_circRNAs.txt -bc 8 -fc 5 -j junctions_rmline.bed -v 11 -e $exon_bed -o cq_output -n ${sample}
```

The ```tail``` command removes the first line of the junctions.bed file which can be found in one of the circRNA output folders. The second command runs CiLiQuant using the generated circRNAs and junctions, as well as a BED file containing all the exons. 

