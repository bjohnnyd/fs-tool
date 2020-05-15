FROM ubuntu:xenial

ENV PATH=/root/.cargo/bin:$PATH
ENV RUSTFLAGS="-Awarnings"

# general packages
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    ca-certificates curl file \
    build-essential \
    git \
    autoconf automake autotools-dev libtool xutils-dev && \
    rm -rf /var/lib/apt/lists/*

# necessary for proper linking of openssl
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    openssl \
    libssl-dev \
    pkg-config 
	

# install toolchain
RUN curl https://sh.rustup.rs -sSf | \
    sh -s -- --default-toolchain stable -y

RUN git clone https://github.com/bjohnnyd/fs-tool.git
WORKDIR fs-tool
RUN cargo build --release --bin fstool
ENTRYPOINT ["./target/release/fstool"]
