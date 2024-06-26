## Docker image for download databases

FROM rnakato/r_python:2022.12
LABEL maintainer="Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>"

WORKDIR /opt
USER root

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    libboost-all-dev \
    libbz2-dev \
    libcurl4-gnutls-dev \
    libgsl-dev \
    libgtkmm-3.0-dev \
    libgzstream0 \
    libgzstream-dev \
    liblzma-dev \
    libncurses5-dev \
    libz-dev \
    cmake \
    curl \
    pigz \
    && apt-get clean \
    && rm -rf /var/lib/apt/list

RUN git clone --recursive https://github.com/rnakato/ChIPseqTools.git \
    && cd ChIPseqTools \
    && make

COPY samtools-1.15.1.tar.bz2 samtools-1.15.1.tar.bz2
RUN tar xvfj samtools-1.15.1.tar.bz2 \
    && cd samtools-1.15.1 \
    && ./configure \
    && make && make install \
    && rm /opt/samtools-1.15.1.tar.bz2

RUN git clone https://github.com/rnakato/SSP.git \
    && cd SSP \
    && make

COPY scripts scripts
COPY bin bin
COPY OriDB OriDB
COPY T2Tdata T2Tdata

COPY gffread-0.12.7.Linux_x86_64.tar.gz gffread-0.12.7.Linux_x86_64.tar.gz
RUN tar zxvf gffread-0.12.7.Linux_x86_64.tar.gz \
    && mv gffread-0.12.7.Linux_x86_64/gffread /opt/bin/ \
    && rm -rf gffread-0.12.7.Linux_x86_64.tar.gz gffread-0.12.7.Linux_x86_64

COPY sratoolkit.3.0.0-ubuntu64.tar.gz sratoolkit.3.0.0-ubuntu64.tar.gz
RUN tar zxvf sratoolkit.3.0.0-ubuntu64.tar.gz \
    && mv sratoolkit.3.0.0-ubuntu64/ /opt/bin/sratoolkit.3.0.0 \
    && rm sratoolkit.3.0.0-ubuntu64.tar.gz

ENV PATH ${PATH}:/opt/:/opt/scripts:/opt/UCSCbins:/opt/bin:/opt/ChIPseqTools/bin/:/opt/SSP/bin:/opt/SSP/scripts:/opt/bin/sratoolkit.3.0.0/bin/

USER ubuntu
CMD ["download_genomedata.sh"]