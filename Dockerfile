# Container for this pipeline: Souporcell 2.5, cellsnp-lite 1.2.3, vireoSNP 0.5.8, samtools 1.9, Python 3.10.
# Rebuilds souporcell_v2.5-python_3.10-vireo_v0.5.8-cellsnp-lite_1.2.3-nextflow-compatible.sif (built 2023-11-27
# from a local image that was never pushed). Versions were read from that SIF. In the rebuild, every tool binary
# and Python package is byte-identical to the SIF; cellsnp-lite pileups match it exactly and fetal_score.py gives
# the same calls (scores within 1e-10 relative).

# labels.json: ContinuumIO docker-images 22.11.1, created 2022-12-23 (Debian 11.6)
FROM continuumio/miniconda3:22.11.1

SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive

# Debian 11 has left the live mirrors, so apt reads the snapshot from the SIF's build date.
# dpkg status in the SIF matches it: libc6 2.31-13+deb11u7, linux-libc-dev 5.10.197-1.
ARG DEBIAN_SNAPSHOT="20231127T000000Z"
RUN printf '%s\n' \
      "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/${DEBIAN_SNAPSHOT} bullseye main" \
      "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/${DEBIAN_SNAPSHOT} bullseye-updates main" \
      "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian-security/${DEBIAN_SNAPSHOT} bullseye-security main" \
      > /etc/apt/sources.list && \
    apt-get -o Acquire::Retries=5 update -y && \
    apt-get -o Acquire::Retries=5 install -y wget build-essential \
      libncurses5-dev zlib1g-dev \
      libbz2-dev liblzma-dev procps

ARG MINIMAP2="2.26"
RUN wget -q -O minimap2.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAP2}/minimap2-${MINIMAP2}_x64-linux.tar.bz2 && \
    mkdir -p /opt/minimap2 && \
    tar xjf minimap2.tar.bz2 -C /opt/minimap2 --strip-components 1
ENV PATH="/opt/minimap2:$PATH"

ARG BEDTOOLS2="2.28.0"
RUN wget -q -O bedtools2.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS2}/bedtools-${BEDTOOLS2}.tar.gz && \
    mkdir -p /opt/bedtools2 && \
    tar zxf bedtools2.tar.gz -C /opt/bedtools2 --strip-components 1 && \
    cd /opt/bedtools2 && make
ENV PATH="/opt/bedtools2/bin:$PATH"

# Rust 1.74.0 (multirust-channel-manifest.toml); the SIF used "stable" on 2023-11-16
ENV CARGO_HOME=/opt/rust/.cargo
ENV RUSTUP_HOME=/opt/rust/.cargo
RUN wget -qO- https://sh.rustup.rs | bash -s -- -y --default-toolchain 1.74.0
ENV PATH="/opt/rust/.cargo/bin:$PATH"

# v2.5 ships Cargo.lock for both crates, so --locked reproduces the dependency versions
ARG SOUPORCELL="2.5"
RUN wget -q -O souporcell.tar.gz https://github.com/wheaton5/souporcell/archive/refs/tags/v${SOUPORCELL}.tar.gz && \
    mkdir -p /opt/souporcell && \
    tar zxf souporcell.tar.gz -C /opt/souporcell --strip-components 1 && \
    cd /opt/souporcell/troublet && cargo build --release --locked && \
    cd /opt/souporcell/souporcell && cargo build --release --locked
ENV PATH="/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH"

# htslib/samtools/bcftools 1.9 (version.h)
ARG SAMTOOLS="1.9"
RUN wget -q -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/${SAMTOOLS}/htslib-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/histlib && tar xjf htslib.tar.bz2 -C /opt/histlib --strip-components 1 && \
    cd /opt/histlib && ./configure && make && make install
RUN wget -q -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS}/samtools-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/samtools && tar xjf samtools.tar.bz2 -C /opt/samtools --strip-components 1 && \
    cd /opt/samtools && ./configure && make && make install
RUN wget -q -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/${SAMTOOLS}/bcftools-${SAMTOOLS}.tar.bz2 && \
    mkdir -p /opt/bcftools && tar xjf bcftools.tar.bz2 -C /opt/bcftools --strip-components 1 && \
    cd /opt/bcftools && ./configure && make && make install

ARG FREEBAYES="1.3.6"
RUN wget -q -O /opt/freebayes https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES}/freebayes-${FREEBAYES}-linux-amd64-static.gz && \
    chmod +x /opt/freebayes
ARG VARTRIX="1.1.22"
RUN wget -q -O /opt/vartrix https://github.com/10XGenomics/vartrix/releases/download/v${VARTRIX}/vartrix_linux && \
    chmod +x /opt/vartrix
ENV PATH="/opt:$PATH"

# conda-meta: cellsnp-lite 1.2.3 (bioconda build h6141fd1_2) with htslib 1.18; python moved to 3.10.13
# The SIF took python 3.10.13 from Anaconda's defaults channel (repo.anaconda.com); this uses
# conda-forge only, as defaults is often blocked and falls under Anaconda's commercial terms.
RUN conda install -y --override-channels -c conda-forge -c bioconda \
      "cellsnp-lite=1.2.3=h6141fd1_2" "htslib=1.18" "libdeflate=1.19" "python=3.10.13" && \
    conda clean -afy

# pip dist-info in /opt/conda; setuptools < 58 is needed to build PyVCF
RUN pip install --no-cache-dir "pip==23.3.1" && \
    pip install --no-cache-dir "setuptools==57.5.0" && \
    pip install --no-cache-dir "Cython==3.0.6" && \
    pip install --no-cache-dir "numpy==1.24.3" "scipy==1.11.4" "pysam==0.22.0" "PyVCF==0.6.8" "pyfaidx==0.7.2.2" && \
    pip install --no-cache-dir "pystan==2.19.1.1" && \
    pip install --no-cache-dir "vireoSNP==0.5.8" "pandas==2.1.3" "matplotlib==3.8.2" \
      "contourpy==1.2.0" "cycler==0.12.1" "fonttools==4.45.1" "kiwisolver==1.4.5" "Pillow==10.1.0" \
      "packaging==23.2" "pyparsing==3.1.1" "python-dateutil==2.8.2" "pytz==2023.3.post1" \
      "tzdata==2023.3" "importlib-metadata==6.8.0" "zipp==3.17.0"

# the SIF's runscript: no ENTRYPOINT, CMD /bin/bash (Nextflow needs no entrypoint)
CMD ["/bin/bash"]
