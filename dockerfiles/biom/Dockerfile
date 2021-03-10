FROM ubuntu:latest

MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)
ENV DEBIAN_FRONTEND noninteractive


RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        python3.6-dev python3-pip gcc\
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN pip3 install setuptools
RUN pip3 install numpy
RUN pip3 install Cython
RUN pip3 install biom-format

WORKDIR /biomformat/

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8