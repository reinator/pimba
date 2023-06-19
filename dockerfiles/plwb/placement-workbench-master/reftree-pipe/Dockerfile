FROM continuumio/miniconda3:4.10.3-alpine as basecontainer
# FROM condaforge/mambaforge:4.11.0-2 as basecontainer

WORKDIR "/app"

RUN conda install -y git && conda clean -afy

# get nidhoggr and create its environment
# RUN git clone --recursive https://github.com/lczech/nidhoggr.git
# when testing, this obviates the need to commit/pull
COPY . ./nidhoggr

RUN conda env create -f nidhoggr/environment.yaml && \
    conda clean -afy

ENTRYPOINT ["./nidhoggr/search.py"]
