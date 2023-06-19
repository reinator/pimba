Snakemake pipeline for phylogenetic tree inference - specifically geared toward reference trees for phylogenetic placement.


Pipeline Overview
-------------------

**Minimal input:**
  - Reference sequences (fasta or .csv of genbank acessions)

**Process and available tools:**

  - handling of multiple input files for separate runs of the pipeline
  - sequence download
  - sequence quality control
  - multiple sequence alignment
  - (optional) model selection
  - tree inference
  - post-analysis
  - task handling for the cluster (slurm)

**Output:**

  - inferred trees
  - quality metrics
  - bootstrap support values
