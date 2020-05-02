# fs-tool 

[![Software License][ico-license]](LICENSE.md)

Command-line tool to calculate fraction of shared bound peptides between HLA alleles from NetMHCpan binding predictions.
The tool currently reports fraction of shared peptides based on default motifs in the peptide but additional positions can be supplied by the user.

Upcoming releases will add the ability to compare fraction shared between an HLA allele and an individual taking into account the HLA and KIR genotypes.

Current HLA allele ligand group assignments that is included with this tool was obtained from `https://www.ebi.ac.uk/` on `2019-12-29`.
When the tool is ran, the first time after installation, it will download current ligand group information and will report where locally it has been stored (OS dependent location).

## Install

Precompiled binary that will run on Linux 64 systems is available for each release. You can download the latest binary from the Releases tab.

To get the latest binary through CLI and download the latest ligand group assignments from EBI run the following command:

``` bash
$ curl -LO https://github.com/bjohnnyd/fs-tool/releases/download/v0.1.3/fs-tool && chmod +x fs-tool && ./fs-tool
```

To compile from source rustup is required and can be obtained [HERE](https://rustup.rs/).  

In addition, [OpenSSL](https://www.openssl.org) needs to be present:

 * macOS with Homebrew:
            ```
            $ brew install openssl pkg-config
            ``` 
 * Linux:
            ```
            $ sudo apt-get install openssl libssl-dev pkg-config
            ```


After installing rustup download the release archive file and build:

``` bash
$ curl -sL https://github.com/bjohnnyd/fs-tool/archive/v0.1.3.tar.gz |  tar xvz && cd fs-tool-0.1.3 && RUSTFLAGS="-Awarnings" cargo build --release --bin fs-tool
```

The resulting binary can then be ran to download the updated ligand data with:

``` bash
$ ./target/release/fs-tool
```

## Usage

Running `fs-tool -h` will list all possible arguments:

``` bash
$ ./fs-tool -h
```

```
    fs-tool 0.1.1
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

To run comparison on positions `1,3,7`, and to name the output measure `Random` while updating data from EBI:

``` bash
$ ./fs-tool -n resources/netmhcpan/example_netmhcpan_wBA.txt -m "Random:1,3,7" --update-ligand-groups -o random_result.tsv
```

will result in:

```
    Measure Index   NonIndex        FS      PeptideLength   IndexBound      NonIndexBound
    TCR     B*2705  A*0301  1       9       1       Bw4-80T 1       A3
    TCR     A*0301  B*2705  1       9       1       A3      1       Bw4-80T
    TCR     B*2705  A*0301  0       10      0       Bw4-80T 1       A3
    TCR     A*0301  B*2705  0       10      1       A3      0       Bw4-80T
    TCR     B*2705  A*0301  0       11      0       Bw4-80T 1       A3
    TCR     A*0301  B*2705  0       11      1       A3      0       Bw4-80T
    Random  B*2705  A*0301  1       9       1       Bw4-80T 1       A3
    Random  A*0301  B*2705  1       9       1       A3      1       Bw4-80T
    Random  B*2705  A*0301  0       10      0       Bw4-80T 1       A3
    Random  A*0301  B*2705  0       10      1       A3      0       Bw4-80T
    Random  B*2705  A*0301  0       11      0       Bw4-80T 1       A3
    Random  A*0301  B*2705  0       11      1       A3      0       Bw4-80T
    KIR     B*2705  A*0301  1       9       1       Bw4-80T 1       A3
    KIR     A*0301  B*2705  1       9       1       A3      1       Bw4-80T
    KIR     B*2705  A*0301  0       10      0       Bw4-80T 1       A3
    KIR     A*0301  B*2705  0       10      1       A3      0       Bw4-80T
    KIR     B*2705  A*0301  0       11      0       Bw4-80T 1       A3
    KIR     A*0301  B*2705  0       11      1       A3      0       Bw4-80T
```

 to drop the default measures `TCR` and `KIR` the flag `--drop-default-measures` can be used.

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
