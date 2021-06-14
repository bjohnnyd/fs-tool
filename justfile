# Command Variables
dcmd := "docker run --rm -it -v $(pwd):/drone/src -w /drone/src" 
denv := "--env CC=o64-clang --env CXX=o64-clang++ --env LIBZ_SYS_STATIC=1" 
dimg := "joseluisq/rust-linux-darwin-builder:latest"
apple_target := "x86_64-apple-darwin"
linux_target := "x86_64-unknown-linux-musl"
windows_target := "x86_64-pc-windows-gnu"

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
	tar cvzf ../../binaries/{{apple_target}}.tar.gz ./fs-tool && \
	zip ../../binaries/{{apple_target}}.zip ./fs-tool
	@echo "Finished compiling {{apple_target}} target binary"

build_unknown-linux-musl:
	@echo "Compiling {{linux_target}} target binary"
	mkdir -p target/binaries
	{{dcmd}} {{dimg}} rustc -vV > target/binaries/latest.linux.compilation.log
	{{dcmd}} {{dimg}} \
	cargo build --release --target {{linux_target}}
	cd target/{{linux_target}}/release  && \
	tar cvzf ../../binaries/{{linux_target}}.tar.gz ./fs-tool && \
	zip ../../binaries/{{linux_target}}.zip ./fs-tool
	@echo "Finished compiling {{linux_target}} target binary"

build_pc-windows-gnu:
	@echo "Compiling {{windows_target}} target binary"
	mkdir -p target/binaries
	{{dcmd}} {{dimg}} rustc -vV > target/binaries/latest.linux.compilation.log
	{{dcmd}} {{dimg}} \
	rustup target add {{windows_target}}
	{{dcmd}} {{dimg}} \
	cargo build --release --target {{windows_target}}
	cd target/{{windows_target}}/release  && \
	tar cvzf ../../binaries/{{windows_target}}.tar.gz ./fs-tool && \
	zip ../../binaries/{{windows_target}}.zip ./fs-tool
	@echo "Finished compiling {{windows_target}} target binary"

precommit:
	@echo "Running standard Testing, Linting and MSRV"
	cargo test -q --all
	cargo fmt -q --all -- --check
	cargo clippy --all-targets \
	--all-features -q -- \
	-D warnings \
	-A clippy::cognitive_complexity
