+++
title="Installation"
description="Detailed instructions on installing on different OS"
weight=1
+++

The following sections describe how to install from source but in most cases the available binaries should work as is.  Therefore, the simplest installation steps are to download the appropriate version for your OS:

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


The following installation instructions on building and installing is if the above fails (some releases of CentOS will not work with the above Linux binary).

