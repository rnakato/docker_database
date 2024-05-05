## Docker image for download databases
FROM rnakato/r_python:2024.04 as common

WORKDIR /opt
USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libboost-all-dev \
    libcurl4-gnutls-dev \
    libgsl-dev \
    libgtkmm-3.0-dev \
    libgzstream0 \
    libgzstream-dev \
    liblzma-dev \
    libz-dev \
    cmake \
    curl \
    pigz \
    && apt-get clean \
    && rm -rf /var/lib/apt/list

# ChIPseqTools
RUN git clone --recursive https://github.com/rnakato/ChIPseqTools.git \
    && cd ChIPseqTools \
    && make

# SSP
RUN git clone https://github.com/rnakato/SSP.git \
    && cd SSP \
    && make

RUN pip --no-cache install ffq

COPY scripts/* scripts/
COPY bin bin
COPY OriDB OriDB

COPY gffread-0.12.7.Linux_x86_64.tar.gz gffread-0.12.7.Linux_x86_64.tar.gz
RUN tar zxvf gffread-0.12.7.Linux_x86_64.tar.gz \
    && mv gffread-0.12.7.Linux_x86_64/gffread /opt/bin/ \
    && rm -rf gffread-0.12.7.Linux_x86_64.tar.gz gffread-0.12.7.Linux_x86_64


FROM rnakato/r_python:2024.04 as normal
LABEL maintainer="Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>"
ENV PATH ${PATH}:/opt/:/opt/scripts:/opt/UCSCbins:/opt/bin:/opt/ChIPseqTools/bin/:/opt/SSP/bin:/opt/SSP/scripts:/opt/bin/sratoolkit.3.0.0/bin/

COPY --from=common / /
USER ubuntu
WORKDIR /home/ubuntu
ENTRYPOINT ["/entrypoint.sh"]
CMD ["download_genomedata.sh"]


FROM rnakato/r_python_gpu:2024.04 as gpu
LABEL maintainer="Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>"
ENV PATH ${PATH}:/opt/:/opt/scripts:/opt/UCSCbins:/opt/bin:/opt/ChIPseqTools/bin/:/opt/SSP/bin:/opt/SSP/scripts:/opt/bin/sratoolkit.3.0.0/bin/

COPY --from=common / /
USER ubuntu
WORKDIR /home/ubuntu
ENTRYPOINT ["/entrypoint.sh"]
CMD ["download_genomedata.sh"]