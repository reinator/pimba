#!/usr/bin/env Rscript

# Input parsing/validation 
mult <- TRUE
verbose <- TRUE
if(exists("snakemake")) {
  if( length(snakemake@log) > 0 ) { 
    logfile <- file(snakemake@log[[1]], open = "wt")
    sink(logfile)
    sink(logfile,type="message")
  }
  if( snakemake@threads < 2 ) { 
    mult <- FALSE
  } else {
    setThreadOptions( numThreads = snakemake@threads )
  }
  if( length(snakemake@input) == 2 ) {
    merge <- TRUE
    fnF <- snakemake@input[[1]]
    fnR <- snakemake@input[[2]]
    outdir <- dirname(snakemake@output[[1]])
  } else if ( length(snakemake@input) == 1 ) {
    merge <- FALSE
    fnF <- snakemake@input[[1]]
    outdir <- dirname(snakemake@output[[1]])
  } else {
    cat("
    Must supply either one fastq (merged) or two fastq (unmerged, named forward and reverse)
    \n\n")
    q(save="no")
  }
  out_fasta <- file.path( snakemake@output[['fasta']] )
  out_table <- file.path( snakemake@output[['otu_table']] )
  
} else {
  opts <- commandArgs(trailingOnly = TRUE)
  if( length(opts) == 3 ) {
    merge <- TRUE
    fnF <- opts[1]
    fnR <- opts[2]
    outdir <- opts[3]
  } else if( length(opts) == 2 ) {
    merge <- FALSE
    fnF <- opts[1]
    outdir <- opts[2]
  } else {
    cat("
    Script to run dada2 pipeline on one sample (two fastq files)
    
    Usage:
    (Paired-end reads): ./dada2.R <R1.fastq> <R2.fastq> <outdir>
    (Already merged):   ./dada2.R <seqs.fastq> <outdir>
    \n\n")
    q(save="no")
  }
  out_fasta <- file.path( outdir, "ASVs.fa" )
  out_table <- file.path( outdir, "ASVs_counts.tsv" )
}

library(dada2)
getN <- function(x) sum(getUniques(x))

# First mode: two sequence files, assumed to be unmerged paired end reads
if( merge ) {
  cat("Meging input files\n")
  sample.names <- tools::file_path_sans_ext(basename(fnF))

  filtF <- file.path( outdir, "filtered", paste0(sample.names, "_F_filt.fastq.gz") )
  filtR <- file.path( outdir, "filtered", paste0(sample.names, "_R_filt.fastq.gz") )

  cat("Filtering/trimming...\n")
  out <- filterAndTrim(fnF,
    filtF,
    fnR,
    filtR,
    truncLen=0,
    maxN=0,
    maxEE=c(2,5),
    truncQ=2,
    rm.phix=TRUE,
    minLen=30,
    compress=TRUE,
    multithread=mult)
  cat("Done!\n")

  cat("Learning errors...\n")
  errF <- learnErrors(filtF, multithread=mult)
  errR <- learnErrors(filtR, multithread=mult)
  cat("Done!\n")

  cat("Running dada...\n")
  dadaF <- dada(filtF, err=errF, multithread=mult, verbose=2)
  dadaR <- dada(filtR, err=errR, multithread=mult, verbose=2)
  cat("Done!\n")

  cat("Paired-end merging...\n")
  mergers <- mergePairs(dadaF, filtF, dadaR, filtR,
    minOverlap = 10,
    verbose=verbose
  )
  cat("Done!\n")

  cat("Mangling ASV table...\n")
  seqtab <- makeSequenceTable(mergers)
  print(seqtab)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=verbose)
  track <- cbind(out, getN(dadaF), getN(dadaR), getN(mergers), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  cat("Done!\n")
} else {
# Use with caution, as per https://github.com/benjjneb/dada2/issues/446
# Alternate mode: start straight at merged input sequences
  cat("Input already merged\n")
  # get sample name 
  sample.names <- tools::file_path_sans_ext(basename(fnF))
  print(sample.names)
  filtF <- file.path( outdir, "filtered", paste0(sample.names, "_filt.fastq.gz") )
  print(filtF)
  
  cat("Filtering/trimming...\n")
  out <- filterAndTrim(
    fnF,
    filtF,
    truncLen=0,
    maxN=0,
    maxEE=2,
    truncQ=2,
    rm.phix=TRUE,
    minLen=30,
    compress=TRUE,
    multithread=mult)
  cat("Done!\n")

  cat("Learning errors...\n")
  errF <- learnErrors(filtF, multithread=TRUE, verbose=verbose)
  errF <- inflateErr(getErrors(errF), 3)
  cat("Done!\n")

  cat("Running dada...\n")
  dadaF <- dada(filtF, err=errF, multithread=mult, verbose=2)
  cat("Done!\n")

  cat("Mangling ASV table...\n")
  seqtab <- makeSequenceTable(dadaF)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=verbose)
  track <- cbind(out, getN(dadaF), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
  rownames(track) <- sample.names
  cat("Done!\n")
}

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write( asv_fasta, out_fasta )

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, out_table, sep="\t", quote=F, col.names=FALSE)
