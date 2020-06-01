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

To install manually (build from source) you will need to:

1. Download and install [Rustup](https://rustup.rs/).
2. Download the latest archive of the repository [HERE](https://github.com/bjohnnyd/fs-tool/archive/master.zip)
3. Uncompress the downloaded `master.zip` and change to the created directory (e.g. `cd fs-tool-master`)
4. Inside the directory run the command `cargo build --release --bin fs-tool`
5. Once the build has finished, you can find the binary inside the directory `/target/release` named `fs-tool` on linux/osx and `fs-tool.exe` on windows.
6. You can copy/move the binary to a directory of choice and add it to your PATH so it can be ran from any location. 
