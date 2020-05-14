+++
title="Installation"
description="Detailed instructions on installing on different OS"
+++

The following sections describe how to install from but the in most cases the available binaries should work as is.

Therefore, the simplest installation steps are to download the appropriate version of the binary for you OS:

- Windows
    - [Zip](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip)
    - [Tar](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.tar.gz)
- Linux
    - [Zip](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.zip)
    - [Tar](https://github.com/bjohnnyd/fs-tool/releases/latest/download/x86_64-pc-windows-gnu.tar.gz)
    
- MacOS (Not available yet)
    
After downloading and uncompressing the archive to test the download open a Terminal or Command Prompt on Windows and navigate to the download location:

```sh
cd <path/to/uncompressed archive>
```

and run the command:

```sh
./fstool --settings
```

or on Windows:

```sh
fstool --settings
```

the output should look similar to this:

```
   Default Measures:
      - TCR:2,3,4,5,6,9
      - KIR,2,7,8,9
   
   Updated data will be storead at: 'path/to/your/os/data/dir'
```

where `path/to/your/os/data/dir` will be an OS specific directory where global data will be stored.


The following are detailed instructions on building and installing on different OS if the above fails (Some releases of CentOS will not work with the above Linux binary).

