FROM rocker/r-ver:latest

MAINTAINER OLIVEIRA Renato <renato.renison@gmail.com>


RUN set -e \
	&& apt-get -y update \
	&& apt-get -y dist-upgrade \
	&& apt-get -y install --no-install-recommends --no-install-suggests zlib1g-dev libxml2 libglpk40\
	&& apt-get -y autoremove \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*


RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('phyloseq')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('grid')"
RUN Rscript -e "install.packages('gridExtra')"
RUN Rscript -e "install.packages('xfun')"
RUN Rscript -e "install.packages('pvclust')"
RUN Rscript -e "install.packages('dendextend')"

WORKDIR /data/

COPY Phyloseq_pimba.R /data/
COPY pvclust.R /data/