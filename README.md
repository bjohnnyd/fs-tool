# fs-tool 

[![Software License][ico-license]](LICENSE.md)
[![Build Status](https://travis-ci.com/bjohnnyd/fs-tool.svg?branch=dev_fs)](https://travis-ci.com/bjohnnyd/fs-tool)

Command-line tool to calculate fraction of shared bound peptides between HLA alleles from NetMHCpan binding predictions.  The tool currently reports fraction of shared peptides based on default motifs in the peptide but additional positions can be supplied by the user.

Current HLA allele ligand group assignments that is included with this tool was obtained from `https://www.ebi.ac.uk/` on `2019-12-29`.
The kir ligand motifs can be updated using the tool.

## Install

The simplest way to install is using the precompiled binaries provided below:

| ![picture](static/64px-Tux.png) | ![picture](static/64px-MacOS_logo.png)  | ![picture](static/64px-Windows_logo_2012.png) |
| :-----------------------------: | :-------------------------------------: |:--------------------------------------------: |
| [TAR](path to release archive) | [TAR](path to release archive)  | [TAR](path to release archive) |
| [ZIP](path to release archive) | [ZIP](path to release archive)  | [ZIP](path to release archive) |

Using the command line you can obtain the binary using:

``` bash
$ curl -LO https://github.com/bjohnnyd/fs-tool/releases/latest/fstool && chmod +x fstool && ./fs-tool -h
```

To compile from source rustup is required and can be obtained [HERE](https://rustup.rs/).  After installing rustup download the release archive file and build:

``` bash
$ curl -sL https://github.com/bjohnnyd/fs-tool/archive/v0.2.0.tar.gz |  tar xvz && cd fs-tool-0.2.0 && cargo build --release --bin fs-tool
```

The compiled binary can then be ran using:

``` bash
$ ./target/release/fs-tool -u
```

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
        --location        Returns the directory and files where the kir ligand data is stored
    -q, --quiet           Disables any information being printed to terminal (except errors)
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
$ ./fstool -n resources/netmhcpan/example_netmhcpan_wBA.txt -m "Random:1,3,7" -u -o example
```


to drop the default measures `TCR` and `KIR` the flag `--drop-default-measures` can be used.

The created directory `example` will contain the following results:



## Usage

Detailed usage instructions are available [HERE](https://bjohnnyd.github.io/fs-tool/public/)

## Authors and Citation

Please cite [eLife 2020;9:e54558](https://doi.org/10.7554/eLife.54558).

- [Bisrat Johnathan Debebe][link-author]
- [Becca Asquith][link-author1]
- [Lies Boelen][link-author2]

## Citation

## License

The MIT License (MIT). Please see [License File](LICENSE.md) for more information.

[ico-version]: https://img.shields.io/packagist/v/:vendor/fs-tool.svg?style=flat-square
[ico-license]: https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square
[ico-travis]: https://img.shields.io/travis/:vendor/fs-tool/master.svg?style=flat-square
[ico-scrutinizer]: https://img.shields.io/scrutinizer/coverage/g/:vendor/fs-tool.svg?style=flat-square
[ico-code-quality]: https://img.shields.io/scrutinizer/g/:vendor/fs-tool.svg?style=flat-square
[ico-downloads]: https://img.shields.io/packagist/dt/:vendor/fs-tool.svg?style=flat-square

[link-packagist]: https://packagist.org/packages/:vendor/fs-tool
[link-travis]: https://travis-ci.org/:vendor/fs-tool
[link-scrutinizer]: https://scrutinizer-ci.com/g/:vendor/fs-tool/code-structure
[link-code-quality]: https://scrutinizer-ci.com/g/:vendor/fs-tool
[link-downloads]: https://packagist.org/packages/:vendor/fs-tool
[link-author]: https://github.com/bjohnnyd
[link-author1]: https://github.com/becca-asquith
[link-author2]: https://github.com/liesb
