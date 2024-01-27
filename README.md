# docker_database

Docker image to download reference data and make index files.
This is the base image for [Churros](https://github.com/rnakato/Churros) and [RumBall](https://github.com/rnakato/RumBall).

- Ubuntu 22.04

- GPU mode (cuda:11.8.0-cudnn8-runtime)
   - CUDA 11.8
   - cudnn 8

- Perl 5.36.0 (with plenv)
- Python 3.9 (with Miniconda)
    - MACS2-2.2.9.1

- R 4.x
    - BiocManager
    - Rstudio Desktop
    - Rstudio Server

- SAMtools 1.19.2
- SRAtoolkit 3.0.10
- BEDtools 2.31.0
- OpenBLAS 0.3.24

- user:password
    - ubuntu:ubuntu

- DockerHub:
  - https://hub.docker.com/r/rnakato/database
  - https://hub.docker.com/r/rnakato/database_gpu


## Changelog

- 2024.01
  - Updated SAMtools from 1.17 to 1.19.2
  - Updated SRAtoolkit from 3.0.2 to 3.0.10

- 2023.12
  - Update ChIPseqTools (ssp v1.3.1 and drompa+ v1.18.1)

- 2023.11
  - Removed LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/compat/:/usr/local/cuda/lib64

- 2023.10
  - Added LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/compat/:/usr/local/cuda/lib64

- 2023.06
  - Added fish

- 2022.09
    - Added UCSC build as a argument
    - Update bedtools from v2.30.0 to v2.31.0

- Ensembl106
- This image is based on Ensembl version 106


## Usage

Run normal image:

    docker run -it --rm rnakato/database /bin/bash

Run with GPU:

    docker run -it --rm --gpus all rnakato/database_gpu /bin/bash

The default user is `ubuntu`. Add `-u root` if you want to login as root:

    docker run -it --rm --gpus all -u root rnakato/database_gpu /bin/bash

## Build images from Dockerfile

    # normal
    docker build -t youracount/database --target normal .
    # with GPU
    docker build -t youracount/database_gpu --target gpu .
 