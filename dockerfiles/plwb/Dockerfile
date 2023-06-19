FROM conda/miniconda3:latest

MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)
ENV DEBIAN_FRONTEND noninteractive

RUN conda install mamba -n base -c conda-forge
ADD https://raw.githubusercontent.com/pierrebarbera/placement-workbench/master/environment.yml /tmp/

RUN mamba env create -f /tmp/environment.yml

SHELL ["conda", "run", "-n", "plwb", "/bin/bash", "-c"]
RUN mamba install 'tabulate=0.8.10'

WORKDIR /plwb/

COPY placement-workbench-master /plwb/

ENTRYPOINT ["conda", "run", "-n", "plwb"]