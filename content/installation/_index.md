+++
title="Installation"
description="Detailed instructions on installing on different OS"
weight=1
+++


The simplest way to use `fs-tool` is to download the provided binary for your OS.

## Binary

- Windows
    - [Zip](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip)
    - [Tar](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.tar.gz)

- Linux
    - [Zip](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip)
    - [Tar](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.tar.gz)

- MacOS
    - [Zip](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-apple-darwin.zip)
    - [Tar](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-apple-darwin.tar.gz)
    
After downloading and uncompressing the archive to test the download open Terminal/Command Prompt and navigate to the download location:

```sh
cd <path/to/uncompressed archive/>
```

and run the command:

```sh
./fs-tool --settings
```

or on Windows:

```ps1
fs-tool.exe --settings
```

the output should look similar to this:

```
   Default Measures:
      - TCR:2,3,4,5,6,9
      - KIR,2,7,8,9
   
   Updated data will be storead at: '<path/to/your/os/data/dir>'
```

where `path/to/your/os/data/dir` will be an OS specific directory where global data will be stored (not the output just metadata).


## Install Manually


The following sections describe how to install from source for each OS separately.

To install manually (build from source) you will first need to download and install [Rustup](https://rustup.rs/).

Building from source should be needed in rare cases that the above binary does not work, the following installation instructions on building and installing is for those rare fail cases (exception is CentOS, where some releases will not work with the above Linux binary and will be needed to build from source).