+++
title="Additional Arguments"
description="Describes all the additional parameters that can be provided"
weight=0
sort_by="weight"
+++

In addition, to the arguments described in previous sections ([Binding Predictions](@/input/binding_predictions.md) and [Cohort](@/input/cohort.md)), detailed descriptions of the remaining arguments are:

---

#### `-o` or `--output`

Specifies the output directory where all the results will be stored. This is a required argument.  If the directory does not exist, it will be created. 

---

#### `-i` or `--index`
Required only when `-c` or `--cohort` is provided.  

Each genotype in the cohort will be compared to this allele.  Multiple index alleles can be provides as space delimited values (e.g. `-i A02:01 C08:01`).

---

#### `-m` or `--measure`

In addition, to the default calculations based on TCR/KIR interaction positions it is possible to provide a custom motif and name.  All calculations will also be performed for the supplied measures.

One or multiple measures can be provided as `Name:pos,pos,pos...` 

For example `-m Example1:1,2,3 Example2:4,5,6`, will produce additional entries in the output for measures named `Example1` and `Example2` based on the specified peptide positions/motifs.

---

#### `--drop-default`

If this flag is supplied, no calculations will be made on default motifs representing the TCR and KIR "interaction space".

---
#### `-p` or `--peptide-length`

By default only 9-mer peptides are used for calculations.  However, it is also possible to perform separate calculation for each k-mer as long as it is present in the binding predictions.  

Multiple peptide lengths can be supplied at the same time as `-p 8 9 10 11`, which will perform independent calculation for each of the 4 k-mer sizes.

---

#### `-u` or `--unique`

If the same motif occurs across diferent bound peptides for a given allele, they will be counted/considered only once.  

---

#### `--prefix`

This will prepend a identifier of choice to all the output files.  This is useful if you want to name the outputs after the specific cohort and proteome combination (e.g. `--prefix HTLV_RandomCohort`).

---

#### `--threads`

By default 4 threads are used for calculations but more can be provided, for very large proteomes this will make significant speed improvements if increased.

---

#### `--update`

This flag will ensure that new EBI kir ligand motif data is downloaded before any calculations.  

---
