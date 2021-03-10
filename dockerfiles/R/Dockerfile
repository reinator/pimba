FROM rocker/r-ver:latest

MAINTAINER OLIVEIRA Renato <renato.renison@gmail.com>

RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "library(devtools)" -e "install_github(\"tobiasgf/lulu\")"
RUN Rscript -e "install.packages('dplyr')"

WORKDIR /data/

COPY file.R /data/