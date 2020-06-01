+++
title="Cohort"
description="Describes the cohort data structure"
weight=2
+++

The path to the cohort data is specified using the `-c` or `--cohort` flag and specifying this input requires -`i` or `--index` to be supplied.  Every individual in the cohort will be compared to each allele passed under the `--index`flag.

```
    -c, --cohort <cohort>
        Cohort of individuals for which all measures will be calculated
    -i, --index <index>...
	Index allele used for cohort calculations only, all individuals will be compared to these alleles
```

The required column is `id` or `ID`. The HLA allele columns can be any of `A1`, `A2`, `B1`, `B2`, `C1`, `C2` where `1` and `2` represent the two alleles an individual has for that specific gene/locus. In addition to the HLA allele the KIR columns should be specified of the format `KIR2DL1` or `2DL1`, `KIR2DL2` or `2DL1` etc., see the [Nomenclature](@/input/nomenclature.md) section for the KIR columns names that are allowed and the HLA allele names that can be specified in the `A1` to `C2` columns.

Alternative names are allowed for some of the columns and are listed below:

- `ID` = `id`, `sample`, `rowid`
- `A1` = `A.1`
- `A2` = `A.2`
- `B1` = `B.1`
- `B2` = `B.2`
- `C1` = `C.1`
- `C2` = `C.2`
- `KIR2DL1` = `2DL1`
- ..
- `KIR3DL2` = `3DL2`

The KIR columns can have the following values represent the presence or absence of a gene:

- `Presence` = `1`, `TRUE`, `Y`, `T`
- `Absence` = `0`, `FALSE`, `N`, `F`


An example cohort file with the column names can be seen [HERE][cohort].

[cohort]: https://github.com/bjohnnyd/fs-tool/blob/master/tests/input/cohorts/example_cohort.csv




