# ---------------------------
# COMMAND LINE ARGUMENTS
# ---------------------------
args <- commandArgs(TRUE)  # Get input arguments from command line

# ---------------------------
# PACKAGE DEPENDENCIES
# ---------------------------
required_packages <- c('dplyr', 'tidyr', 'openxlsx', 'readxl', 'data.table', 
                      'stringr', 'taxizedb', 'ggplot2')
sapply(required_packages, require, character.only = TRUE)

# Load Biostrings explicitly
library(Biostrings)

# Clean up
rm(required_packages)
# ---------------------------
# DATA PRE-PROCESSING
# ---------------------------

# Load custom functions
source("/app/custom_functions.R")

# Define input files from arguments
tax_file <- args[1]   
otu_file <- args[2]    
fasta_file <- args[3]  
xlsx <- args[4]        

# Combine taxonomy and OTU data
tax_datasets <- combine_tax_otu(
  tax_file = tax_file,
  otu_file = otu_file,
  fasta_file = fasta_file
)

# Process the combined datasets
parsing_taxa(
  tax_otu_df = tax_datasets$tax_otu_table,
  raw_otus_df = tax_datasets$raw_otus,
  output_xlsx = xlsx
)

