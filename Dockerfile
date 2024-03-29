FROM ubuntu:22.04
LABEL author="Edoardo Giacopuzzi"
LABEL contact="edoardo.giacopuzzi@fht.org"
LABEL version="1.3.1"

ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION="1.3.1"

# Install software
RUN apt-get update && \
    apt-get install -y \
        curl \
        bzip2 \
        libbz2-dev \
        libsqlite3-dev \
        libcurl4-openssl-dev \
        zlib1g zlib1g-dev \
        liblzma-dev libssl-dev \
        libffi-dev \
        xz-utils \
        libncurses5-dev \
        wget \
        tar \
        unzip \
        procps \
    && apt-get -y autoremove \
    && apt-get -y clean all \
    && rm -rf /var/cache

# Download GREEN-VARAN binaries
WORKDIR /opt/green-varan_bin
RUN wget https://github.com/edg1983/GREEN-VARAN/releases/download/v${VERSION}/greenvaran && \
    wget https://github.com/edg1983/GREEN-VARAN/releases/download/v${VERSION}/greendb_query && \
    chmod a+x greenvaran greendb_query

# Install full GREEN-VARAN suite
WORKDIR /opt
RUN wget https://github.com/edg1983/GREEN-VARAN/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -zxvf v${VERSION}.tar.gz && \
    rm v${VERSION}.tar.gz && \
    mv GREEN-VARAN-${VERSION} green-varan && \
    cd /opt/green-varan/workflow/bin && \
    chmod a+x *

# Link executables to /usr/bin
WORKDIR /usr/bin
RUN ln -s /opt/green-varan_bin/greenvaran ./ && \
    ln -s /opt/green-varan_bin/greendb_query ./ && \
    ln -s /opt/green-varan/workflow/bin/bgzip ./ && \
    ln -s /opt/green-varan/workflow/bin/tabix ./ && \
    ln -s /opt/green-varan/workflow/bin/vcfanno ./