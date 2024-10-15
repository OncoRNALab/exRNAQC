# Phase 2

## exRNAQC017

For the evaluation of different tube-kit combinations in the second phase of exRNAQC, blood was 
drawn from 5 healthy volunteers. We tested 2 different kits in combination with 3 different 
tubes that delivered the best results in phase 1: 

Kits:

* miRNeasy Serum/Plasma Kit (abbreviated to MIR; Qiagen, 217184)
* QIAamp ccfDNA/RNA Kit (abbreviated to CCF; Qiagen, 55184)

Tubes:

* EDTA
* serum
* citrate

## repeated measures analysis

For every metric used to compare kit-tube combinations, we tested the possibility of 
interaction. We build a linear mixed-effects model with tube, RNA isolation kit and timelapse 
as fixed effects and donorID as random effect.

- Run the Rmarkdown repeated_measures_smallRNA.Rmd to reproduce the subplot (b) in Figure 6.
- Run the Rmarkdown repeated_measures.Rmd to reproduce the subplot (a) in Figure 6.
