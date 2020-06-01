+++
title="Usage"
description="Describes the way to run the different modes"
weight=3
sort_by="weight"
+++

The `fs-tool` is a command-line program that can essentially be used in two modes:

1. **Allele Mode**: input is a file(s)  of peptide predictions for the HLA class I alleles and proteome of interest. It performs the fraction shared calculations pairwise for each allele in the input file and generates an output file which lists the fraction shared by the different metrics. 
2. **Cohort Mode**: The same input as in `allele mode` is provided but in addition the cohort to be studied and the index allele(s) of interest are provided. It performs the fraction shared calculations for each metric for each individual in the cohort.

Use `allele mode` if you are interested in studying the attributes of an allele independent of a particular cohort, `cohort mode` if you want to do regression analysis to investigate cell types underlying an association between the index allele and outcome in a particular cohort. 

`fs-tool` automatically toggles between modes depending on the input files you provide (i.e. if you provide a cohort and index allele as well as the binding predictions it will automatically perform `cohort mode`, if you only provide binding predictions it will execute `allele mode`).

