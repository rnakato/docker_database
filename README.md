# docker_database
Docker image to download reference data and make index files

The base image for [Churros](https://github.com/rnakato/Churros) and [RumBall](https://github.com/rnakato/RumBall)

DockerHub:
- https://hub.docker.com/r/rnakato/database
- https://hub.docker.com/r/rnakato/database_gpu

## Changelog

- 2023.10
  - Add LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/compat/:/usr/local/cuda/lib64

- 2023.06
  - Add fish
  
- 2022.09
-- Add UCSC build as a argument

- Ensembl106
- This image is based on Ensembl version 106