# Command Variables
dcmd := "docker run --rm -it -v $(pwd):/drone/src -u `id -u`:`id -g` -w /drone/src" 
denv := "--env CC=o64-clang --env CXX=o64-clang++ --env LIBZ_SYS_STATIC=1" 
dimg := "joseluisq/rust-linux-darwin-builder:latest"
apple_target := "x86_64-apple-darwin"
linux_target := "x86_64-unknown-linux-musl"
windows_target := "x86_64-pc-windows-gnu"
binary := "fs-tool"

# Aliases
alias macos := build_apple-darwin
alias linux := build_unknown-linux-musl
alias windows := build_pc-windows-gnu

build_apple-darwin:
	@echo "Compiling {{apple_target}} target binary"
	mkdir -p target/binaries
	{{dcmd}} {{dimg}} rustc -vV > target/binaries/latest.darwin.compilation.log
	{{dcmd}} {{denv}} {{dimg}} \
	cargo build --release --target {{apple_target}}
	cd target/{{apple_target}}/release && \
	tar cvzf ../../binaries/{{apple_target}}.tar.gz ./{{binary}} && \
	zip ../../binaries/{{apple_target}}.zip ./{{binary}}
	@echo "Finished compiling {{apple_target}} target binary"

build_unknown-linux-musl:
	@echo "Compiling {{linux_target}} target binary"
	mkdir -p target/binaries
	{{dcmd}} {{dimg}} rustc -vV > target/binaries/latest.linux.compilation.log
	{{dcmd}} {{dimg}} \
	cargo build --release --target {{linux_target}}
	cd target/{{apple_target}}/release && \
	cd target/{{linux_target}}/release  && \
	tar cvzf ../../binaries/{{linux_target}}.tar.gz ./{{binary}} && \
	zip ../../binaries/{{linux_target}}.zip ./{{binary}}
	@echo "Finished compiling {{linux_target}} target binary"

build_pc-windows-gnu:
	@echo "Compiling {{windows_target}} target binary"
	mkdir -p target/binaries
	cargo install --list | grep cross || cargo install cross
	systemctl is-active --quiet docker || (echo "Docker daemon is not running cannot cross" && exit 1)
	cross build --release --target {{windows_target}}
	cd target/{{windows_target}}/release  && \
	tar cvzf ../../binaries/{{windows_target}}.tar.gz ./{{binary}}.exe && \
	zip ../../binaries/{{windows_target}}.zip ./{{binary}}.exe
	@echo "Finished compiling {{windows_target}} target binary"

precommit:
	@echo "Running standard Testing, Linting and MSRV"
	cargo build
	cargo test -q --all
	cargo fmt -q --all -- --check
	cargo clippy --all-targets \
	--all-features -q -- \
	-D warnings \
	-A clippy::cognitive_complexity
