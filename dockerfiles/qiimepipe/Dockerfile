FROM ubuntu:latest

MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)
ENV DEBIAN_FRONTEND noninteractive


RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        python3.9-minimal python3-pip\
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN pip3 install setuptools
RUN pip3 install uparse_python
RUN pip3 install die

WORKDIR /qiimepipe/
COPY *.py /qiimepipe/
COPY fastq.py /usr/local/lib/python3.6/dist-packages/uparse_python/
COPY primer.py /usr/local/lib/python3.6/dist-packages/uparse_python/
COPY die.py /usr/local/lib/python3.6/dist-packages/uparse_python/
COPY uc.py /usr/local/lib/python3.6/dist-packages/uparse_python/
COPY fasta.py /usr/local/lib/python3.6/dist-packages/uparse_python/

ENV LANG pt_BR.UTF-8