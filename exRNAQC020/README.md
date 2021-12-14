
# Advanced data analyses phase 1/2
exRNAQC020 encompasses all extra analysis done in the context of phase 1 and phase 2 of the exRNAQC study. 

## Fraction of circRNAs in plasma
The unique structure of circRNAs makes them resistant agains exonucleases and are therefore hypothesized to remain relatively stable in a RNAse-rich environment, such as plasma. We investigated the differences between tubes in how the fraction of circRNAs changes over time. 

### circRNA detection f
We used our in-house circRNA detection pipeline, which can be found [here](https://github.ugent.be/vandesompelelab/circRNA). We used the following command to run the pipeline:

```submit_circRNA_slurm.py -p CIRCexplorer2_TopHat -t pe -d clumpify -b base_dir -g hg38 -o output_dir -s string_match -u user_email```


