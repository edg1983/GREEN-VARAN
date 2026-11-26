# Stage 1: Build completely static musl executables for greenvaran and greendb_query
FROM rust:1-alpine AS builder

RUN apk add --no-cache \
    build-base \
    musl-dev \
    clang-dev \
    clang-libclang \
    llvm-dev \
    cmake \
    make \
    perl \
    linux-headers \
    file

WORKDIR /build

# Cache dependencies before copying source files
COPY Cargo.toml Cargo.lock ./
RUN mkdir -p src/bin && \
    echo "pub fn _dummy() {}" > src/lib.rs && \
    echo "fn main() {}" > src/bin/greenvaran.rs && \
    echo "fn main() {}" > src/bin/greendb_query.rs && \
    RUSTFLAGS="-C target-feature=+crt-static" HTS_STATIC=1 \
    cargo build --release --features static || true

# Build actual binaries
COPY src ./src
RUN touch src/lib.rs src/bin/greenvaran.rs src/bin/greendb_query.rs && \
    RUSTFLAGS="-C target-feature=+crt-static" HTS_STATIC=1 \
    cargo build --release --features static && \
    file /build/target/release/greenvaran | grep -q "statically linked" && \
    file /build/target/release/greendb_query | grep -q "statically linked"

# Stage 2: Minimal runtime image (glibc-compatible for workflow/bin helper tools)
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION="1.4"

LABEL author="Edoardo Giacopuzzi"
LABEL contact="edoardo.giacopuzzi@fht.org"
LABEL version="${VERSION}"
LABEL org.opencontainers.image.source="https://github.com/edg1983/GREEN-VARAN"
LABEL org.opencontainers.image.description="GREEN-VARAN: Fast Rust tool suite for non-coding regulatory variant annotation and prioritization with GREEN-DB"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
        zlib1g \
        bzip2 \
        xz-utils \
        tar \
        unzip \
        procps \
    && apt-get -y autoremove \
    && apt-get -y clean all \
    && rm -rf /var/lib/apt/lists/* /var/cache

# Install compiled Rust binaries
WORKDIR /opt/green-varan_bin
COPY --from=builder /build/target/release/greenvaran /opt/green-varan_bin/greenvaran
COPY --from=builder /build/target/release/greendb_query /opt/green-varan_bin/greendb_query
RUN chmod a+x /opt/green-varan_bin/greenvaran /opt/green-varan_bin/greendb_query

# Copy configs, workflow, and accessory scripts
WORKDIR /opt/green-varan
COPY config ./config
COPY resources ./resources
COPY workflow ./workflow
COPY addon_scripts ./addon_scripts
RUN cp /opt/green-varan_bin/greenvaran /opt/green-varan/workflow/bin/greenvaran && \
    chmod a+x /opt/green-varan/workflow/bin/*

# Link executables into /usr/bin
RUN ln -sf /opt/green-varan_bin/greenvaran /usr/bin/greenvaran && \
    ln -sf /opt/green-varan_bin/greendb_query /usr/bin/greendb_query && \
    ln -sf /opt/green-varan/workflow/bin/bgzip /usr/bin/bgzip && \
    ln -sf /opt/green-varan/workflow/bin/tabix /usr/bin/tabix && \
    ln -sf /opt/green-varan/workflow/bin/vcfanno /usr/bin/vcfanno

# Smoke-test binaries during image build
RUN greenvaran version && greendb_query --help > /dev/null

WORKDIR /data
