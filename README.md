# fs-tool 

[![Release][ico-version]][link-version]
[![Build Status][ico-travis]][link-travis]
[![Software License][ico-license]](LICENSE.md)

Command-line tool to calculate fraction of shared bound peptides between HLA alleles from NetMHCpan binding predictions.  The tool currently reports fraction of shared peptides based on default motifs in the peptide but additional positions can be supplied by the user.

Current HLA allele ligand group assignments that is included with this tool was obtained from `https://www.ebi.ac.uk/` on `2019-12-29`.
The kir ligand motifs can be updated using the tool.

Detailed instructions and descriptions are available in the [Documentation](https://bjohnnyd.github.io/fs-tool/public) (documentation will be available soon).



## Install

### Binary

The simplest way to install is using the precompiled binaries provided below:

| ![picture](static/64px-Tux.png) | ![picture](static/64px-MacOS_logo.png)  | ![picture](static/64px-Windows_logo_2012.png) |
| :-----------------------------: | :-------------------------------------: |:--------------------------------------------: |
| [TAR](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-unknown-linux-gnu.tar.gz) | (COMING SOON)  | [TAR](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.tar.gz) |
| [ZIP](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip) | (COMING SOON)  | [ZIP](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip) |

Using the command line you can obtain the binary using (Linux):

``` bash
$ curl -sL https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-unknown-linux-gnu.tar.gz | tar xvz && chmod +x fstool && ./fstool -h
```

### Build 

To compile from source rustup is required and can be obtained [HERE](https://rustup.rs/).  After installing rustup download the release archive file and build:

```bash
$ git clone https://github.com/bjohnnyd/fs-tool.git && cd fs-tool && cargo build --release --bin fs-tool
```

The compiled binary can then be ran using:

``` bash
$ ./target/release/fstool -h
```

All releases and associated binaries and archives are accessible here [Releases](https://github.com/bjohnnyd/fs-tool/releases/latest).

## Usage

Running `fstool -h` will list all possible arguments:

``` bash
$ ./fstool -h
```

```

fstool 0.2.0
Calculates fraction of shared bound motifs between HLA alleles while incorporating KIR ligand and LILRB binding
information.

USAGE:
    fstool [FLAGS] [OPTIONS] --binding-predictions <binding-predictions> --output <output>

FLAGS:
        --drop-default    Drop default measures based on TCR and KIR motifs.
    -h, --help            Prints help information
    -q, --quiet           Disables any information being printed to terminal (except errors)
        --settings        Lists default measure names and motif positions as well as the default location updated kir
                          ligand will be stored
    -u, --unique          Whether only unique peptide/motif sequences should be considered in the calculations
        --update          Updates the current kir ligand group data
    -V, --version         Prints version information
    -v, --verbose         Determines verbosity of the processing, can be specified multiple times -vvv

OPTIONS:
    -b, --binding-predictions <binding-predictions>
            Path to file containing predicted Class I affinity data (NetMHCpan results)

    -c, --cohort <cohort>                              Cohort of individuals for which all measures will be calculated
    -i, --index <index>...
            Index allele used for cohort calculations only, all individuals will be compared to these alleles

    -m, --measure <measure>...
            Custom motif positions to use for calculations (format `Name:index,index..` e.g. KIR:2,7,8,9)

    -o, --output <output>                              Directory to store outputs
    -p, --peptide-length <peptide-length>...
            Which length of input peptide sequence to consider [default: 9]  [possible values: 8, 9, 10, 11]

        --prefix <prefix>                              Prefix to assign to all outputs
    -t, --threads <threads>                            Number of threads [default: 4]


```

## Example

### Allele calculations only 

To run comparison on positions `1,3,7`, and to name the output measure `Example` while updating data from EBI:

``` bash
$ ./fstool cargo run --bin fstool --  -b tests/netmhcpan/netmhcpan_wBA.txt  --prefix "example_cohort_Gag_180_209" -o example_result
```


to drop the default measures `TCR` and `KIR` the flag `--drop-default-measures` can be used.

### Cohort calculations

To perform calculations, using the default measures, for `A02:01` and `C08:02` the following command can be ran:

``` bash
$ ./fstool cargo run --bin fstool --  -b tests/netmhcpan/netmhcpan_wBA.txt  --prefix "example_cohort_Gag_180_209" -o example_result -i A03:01 C08:02 -c tests/example_cohort.csv
```

### Output

The created directory `example_result` will contain the following output: 

| File name | Description  | 
| :-----------------------------: | :-------------------------------------: |
| **example_cohort_Gag_180_209_allele_binding_summary.csv** |  summary of allele peptide binding counts per protein  |  |
| **example_cohort_Gag_180_209_allele_fs_result.csv** | fraction shared calculation results for all combinations of alleles in the binding predictions |  |
| **example_cohort_Gag_180_209_allele_metadata.csv** | lists netmhcpan nearest neighbour information and what allele was the ligand motif assignment based on  |  |
| **example_cohort_Gag_180_209_cohort_result.csv** |  lists per cohort subject calculations for each index allele, peptide length, measure |  |


Additional information can be found in the [Output](https://bjohnnyd.github.io/fs-tool/public/output) section of the [Documentation](https://bjohnnyd.github.io/fs-tool/public) (documentation will be available soon).


## Authors and Citation

Please cite [eLife 2020;9:e54558](https://doi.org/10.7554/eLife.54558).

- [Bisrat Johnathan Debebe][link-author]
- [Becca Asquith][link-author1]
- [Lies Boelen][link-author2]

## Citation

## License

The MIT License (MIT). Please see [License File](LICENSE.md) for more information.

[ico-version]: https://img.shields.io/github/v/release/bjohnnyd/fs-tool?include_prereleases
[ico-license]: https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square
[ico-travis]: https://img.shields.io/travis/com/bjohnnyd/fs-tool/master?style=flat-square
[ico-downloads]: https://img.shields.io/packagist/dt/:vendor/fs-tool.svg?style=flat-square

[link-version]: https://github.com/bjohnnyd/fs-tool/releases
[link-travis]: https://travis-ci.com/bjohnnyd/fs-tool
[link-downloads]: https://packagist.org/packages/bjohnnyd/fs-tool
[link-author]: https://github.com/bjohnnyd
[link-author1]: https://github.com/becca-asquith
[link-author2]: https://github.com/liesb
