FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive
MAINTAINER OLIVEIRA Renato (renato.renison@gmail.com)

RUN apt-get -y update
RUN apt-get -y upgrade
#RUN apt-get -y install --no-install-recommends --no-install-suggests build-essential zlibc libc6 gcc g++ libboost1.67-all-dev cmake
RUN apt-get -y install --no-install-recommends --no-install-suggests wget
RUN apt-get -y install --no-install-recommends --no-install-suggests unzip
RUN apt-get -y autoremove; \
	apt-get clean; \
	rm -rf /var/lib/apt/lists/*


RUN wget --no-check-certificate https://github.com/torognes/vsearch/releases/download/v2.15.2/vsearch-2.15.2-linux-x86_64.tar.gz

RUN tar xzf vsearch-2.15.2-linux-x86_64.tar.gz

RUN mv /vsearch-2.15.2-linux-x86_64/bin/vsearch /usr/local/bin/

RUN rm -r /vsearch-2.15.2-linux-x86_64

RUN rm /vsearch-2.15.2-linux-x86_64.tar.gz
