# fs-tool 

[![Software License][ico-license]](LICENSE.md)
[![Build Status](https://travis-ci.com/bjohnnyd/fs-tool.svg?branch=dev_fs)](https://travis-ci.com/bjohnnyd/fs-tool)

Command-line tool to calculate fraction of shared bound peptides between HLA alleles from NetMHCpan binding predictions.
The tool currently reports fraction of shared peptides based on default motifs in the peptide but additional positions can be supplied by the user.

Current HLA allele ligand group assignments that is included with this tool was obtained from `https://www.ebi.ac.uk/` on `2019-12-29`.
The kir ligand motifs can be updated using the tool.

## Install

The simplest way to install is using the precompiled binaries provided below:

| ![picture](static/64px-Tux.png) | ![picture](static/64px-MacOS_logo.png)  | ![picture](static/64px-Windows_logo_2012.png) |
| :-----------------------------: | :-------------------------------------: |:--------------------------------------------: |
| [tar.gz](path to release archive) | [tar.gz](path to release archive)  | [tar.gz](path to release archive) |
| [zip](path to release archive) | [zip](path to release archive)  | [zip](path to release archive) |

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
    fs-tool 0.2.0
    Calculates fraction of shared peptides between HLA alleles based on NetMHCpan predictions

    USAGE:
        fs-tool [FLAGS] [OPTIONS]

    FLAGS:
        -d, --debug
            --drop-default-measures
        -h, --help                     Prints help information
            --update-ligand-groups
        -V, --version                  Prints version information

    OPTIONS:
        -m, --measures <measures>...
        -n, --netmhcpan <netmhcpan>
        -o, --output <output>
        -p, --peptide-length <peptide-length>...     [default: 9 10 11]
        -t, --threads <threads>                      [default: 4]
```

## Example

### Allele calculations only 
To run comparison on positions `1,3,7`, and to name the output measure `Example` while updating data from EBI:

``` bash
$ ./fstool -n resources/netmhcpan/example_netmhcpan_wBA.txt -m "Random:1,3,7" -u -o example
```


to drop the default measures `TCR` and `KIR` the flag `--drop-default-measures` can be used.

The created directory `example` will contain the following results:




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
