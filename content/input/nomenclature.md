+++
title="Nomenclature"
description="Describes the Nomenclature standard throughtout"
weight=0
+++


`fs-tool` in all inputs expects a certain level on nomenclature standards being maintained.

## HLA Class I

HLA alllele nomenclature can vary substantially across tools, laboratories or even agencies.  
There are certain times where it will take a considerable time to standardize the names where the benefit is not obvious.
Example is the `HLA Prefix` or the `*` separator, however the field separator `:` and maintaining the Gene (`A0201` vs  `0201`) as part of the allele name ensure certain safeguards.

An obvious one, is that using different programs and dealing with values as `0201` or `A0201` could lead to unwanted removal of the `0` digit as well as cause problems during parsing when such inputs are supplied to various programs.

For this reason in all inpus to `fs-tool`, the `:` (colon) field separator must be present in both the `--index` allele supplied and in the cohort HLA genotype columns. This means that the following is true:

| Example | Valid |
| :-----: | :----:|
| A0201   |   No   |
| 0201    |   No   |
| A02:01  |  Yes  |
| A*02:01 |  Yes  |
| HLA-A*02:01  |  Yes  |
| HLA-A*24:09N  |  Yes  |
| A30:14L  |  Yes  |


For further information on HLA nomenclature, please visit [http://hla.alleles.org/nomenclature/naming.html](http://hla.alleles.org/nomenclature/naming.html).

## KIR 

KIR alleles supplied are valid with or without acronym and with or without  allele specifed. However, alleles speciefed have to have at least the three letter series specified (`KIR2DL1*003` vs `KIR2DL1*3`).  
The `KIR acronym` is not required and can be ommited in the column names of cohort data.

| Example | Valid |
| :-----: | :----:|
| 2DL1*03   |   No   |
| KIR2DL1   |   Yes   |
| 2DL1   |   Yes   |
| KIR2DL5A   |   Yes   |
| 2DL5B   |   Yes   |
| 2DL1*0030202   |   Yes   |


The following KIR genes are currently supported:

- 2DL1
- 2DL2
- 2DL3
- 2DL4
- 2DL5A/2DL5B
- 2DS1
- 2DS2
- 2DS3
- 2DS4
- 2DS5
- 3DL1/3DS1
- 3DL2
- 3DL3
- 2DP1
- 3DP1

For more information, on the KIR gene/alelle nomenclature,  please visit  [https://www.ebi.ac.uk/ipd/kir/alleles.html](https://www.ebi.ac.uk/ipd/kir/alleles.html)
