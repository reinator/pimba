FROM mambaorg/micromamba:1.4.0

LABEL maintainer="LEAO Tiago Ferreira <tiago.leao@pq.itv.org>"

WORKDIR /opt

RUN micromamba install -y -n base \
    -c conda-forge -c bioconda prinseq-plus-plus=1.2.4 && \
    micromamba clean --all --yes

# Add Conda bin directory to PATH
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /app

ENTRYPOINT ["prinseq++"]
