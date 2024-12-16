FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)

RUN apt-get -y update
RUN apt-get -y upgrade
#RUN apt-get -y install --no-install-recommends --no-install-suggests autoconf automake bzip2-devel gcc-c++ ghostscript groff-base make zlib-devel
RUN apt-get -y install --no-install-recommends --no-install-suggests build-essential zlibc libc6 gcc g++ libboost1.67-all-dev cmake
RUN apt-get -y install --no-install-recommends --no-install-suggests wget
RUN apt-get -y install --no-install-recommends --no-install-suggests unzip
RUN apt-get -y autoremove; \
	apt-get clean; \
	rm -rf /var/lib/apt/lists/*


RUN wget --no-check-certificate https://github.com/torognes/vsearch/archive/v2.29.1.tar.gz

RUN tar xzf v2.29.1.tar.gz

RUN cd vsearch-2.29.1 && ./autogen.sh

RUN cd vsearch-2.29.1 && ./configure CFLAGS="-O3" CXXFLAGS="-O3"

RUN cd vsearch-2.29.1 && make ARFLAGS="cr"

RUN cd vsearch-2.29.1 && make install

RUN mv /vsearch-2.29.1/bin/vsearch /usr/local/bin/

RUN rm /v2.29.1.tar.gz
