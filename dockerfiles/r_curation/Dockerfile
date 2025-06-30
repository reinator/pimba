# Base Bioconductor image
FROM bioconductor/bioconductor_docker:latest

LABEL maintainer="ROSA Fabricio dos Anjos Santa <fabricio.rosa@pq.itv.org>"

# From CRAN
RUN Rscript -e "install.packages(c('BiocManager', 'data.table', 'devtools', 'dplyr', 'ggplot2', 'openxlsx', 'readxl', 'stringr', 'tidyr'))"

# From Github
RUN Rscript -e "devtools::install_github(\"ropensci/taxizedb\")"

# From Bioconductor
RUN Rscript -e "BiocManager::install('Biostrings', ask = FALSE)"

# Main working directory
WORKDIR /app

# Copy scripts
COPY main.R /app/
COPY custom_functions.R /app/