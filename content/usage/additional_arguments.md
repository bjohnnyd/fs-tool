+++
title="Additional Arguments"
description="Describes all the additional parameters that can be provided"
weight=0
sort_by="weight"
+++

The full list of arguments are shown below:

```
fstool 0.2.5
Calculates fraction of shared bound motifs between HLA alleles while incorporating KIR ligand and LILRB binding
information.

USAGE:
    fs-tool [FLAGS] [OPTIONS] --binding-predictions <binding-predictions>... --output <output>

FLAGS:
        --drop-default    Drop default measures based on TCR and KIR motifs
    -h, --help            Prints help information
    -q, --quiet           Disables any information being printed to terminal (except errors)
        --settings        Lists default measure names and motif positions as well as the default location updated kir
                          ligand will be stored
    -u, --unique          Whether only unique peptide/motif sequences should be considered in the calculations
        --update          Updates the current kir ligand group data
    -V, --version         Prints version information
    -v, --verbose         Determines verbosity of the processing, can be specified multiple times -vvv

OPTIONS:
    -b, --binding-predictions <binding-predictions>...
            Path to file containing predicted Class I affinity data (NetMHCpan results)

    -c, --cohort <cohort>
            Cohort of individuals for which all measures will be calculated

    -i, --index <index>...
            Index allele used for cohort calculations only, all individuals will be compared to these alleles

    -m, --measure <measure>...
            Custom motif positions to use for calculations (format `Name:index,index..` e.g. KIR:2,7,8,9)

    -o, --output <output>                                 Directory to store outputs
    -p, --peptide-length <peptide-length>...
            Which length of input peptide sequence to consider [default: 9]  [possible values: 8, 9, 10, 11]

        --prefix <prefix>                                 Prefix to assign to all outputs
    -t, --threads <threads>                               Number of threads [default: 4]
```

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

One or multiple measures can be provided, `-m Example1:1,2,3 Example2:4,5,6`. The output will then contain additional results for `Example1` and `Example2`.

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
