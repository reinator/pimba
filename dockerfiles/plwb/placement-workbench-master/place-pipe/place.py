#!/usr/bin/env python3

import multiprocessing
import argparse, sys, os
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(script_dir, '..', 'common'))
import util
from datetime import datetime
import pandas as pd
import snakemake
import platform

# if we're on mac/arm, default back down to x86 compiled packages
# ideally this decision would be per-package availability
if platform.machine() == 'arm64':
  os.environ['CONDA_SUBDIR'] = "osx-64"

parser = argparse.ArgumentParser(description='Wrapper to run the pipeline with some default settings.')
###
#  Input Files
###
input_group = parser.add_argument_group('Input')
input_group.add_argument('--fasta-paths', dest='fasta_files', type=str, nargs='+',
                    help='input fasta files')
input_group.add_argument('--fastq-paths', dest='fastq_files', type=str, nargs='+',
                    help="input fastq files. If there are unmerged (paired-end) samples,"
                    " forward and reverse must be in separate files"
                    " and you MUST specify the naming pattern by which they are differentiated"
                    " using the --merge-pattern argument, otherwise the pipeline will try to place them separately.")
input_group.add_argument('--merge-pattern', dest='merge_pattern', type=str,
                    help="Pattern by which unmerged forward/reverse fastq files differ. Enables paired-end" 
                    " merging in the pipeline. For example"
                    " if the files are called A_f.fq and A_r.fq, the appropriate pattern would be"
                    " '_(f|r).fq'. PLEASE NOTE: the forward file must be first when lexicographically sorted"
                    " and the differentiating part must be at the end of the filename/path, meaning"
                    " separating forward and reverse into folders (/forward/sample1.fq etc.) is not supported.")

input_group.add_argument('--reference-tree', dest='ref_tree', type=str,
                    help='Reference tree, in newick format', required=True)
input_group.add_argument('--reference-msa', dest='ref_msa', type=str,
                    help='Reference MSA, in fasta format', required=True)
input_group.add_argument('--model-file', dest='model_file', type=str,
                    help='Reference tree, in newick format', required=True)
input_group.add_argument('--taxonomy-file', dest='taxon_file', type=str,
                    help='Tab-separated file mapping reference labels to their taxonomic paths', required=False)

###
#   Config Manip
###
pipeline_group = parser.add_argument_group('Pipeline Options')
pipeline_group.add_argument('--out-dir', dest='out_dir', type=str,
                    default=None,
                    help='optional output directory. By default, the script creates a timestamped output directory.')
pipeline_group.add_argument('-d', '--datatype', dest='datatype', type=str, nargs='?',
                    const='nt', default='nt', choices=['nt', 'aa'],
                    help="datatype, 'aa' for protein, 'nt' for DNA data")
pipeline_group.add_argument('-p', '--prefix', dest='prefix', type=str,
                    default=None,
                    help='prefix to fasta paths and output (useful to specify where data was mounted to in docker)')
pipeline_group.add_argument('--no-chunkify', dest='no_chunkify', action='store_true',
                    help='use the chunkify routine to split queries into blocks of more managable size')
pipeline_group.add_argument('--sequence-clustering', dest='sequence_clustering', type=str, nargs='?',
                    action='store', default=None, choices=['swarm', 'dada2'],
                    help="Which tools, if any, to use for query clustering.")

###
#   Cluster options
###
cluster_group = parser.add_argument_group('Cluster Execution')
cluster_group.add_argument('--cluster-exec', dest='on_cluster', action='store_true',
                    help="Starts the pipeline in computing-cluster (slurm/sge etc.) submission mode. "
                    "Highly recommended to do this from a screen/tmux session!")
cluster_group.add_argument('--cluster-env', dest='clust_env', type=str, nargs='?',
                    const='auto', default='auto', choices=['auto','slurm', 'sge'],
                    help="What job submission system we are on. 'auto' attempts to autodetect.")

###
#   Runtime Manip
###
misc_group = parser.add_argument_group('Misc. Options')
misc_group.add_argument('--threads', dest='threads', type=int,
                    default=multiprocessing.cpu_count(),
                    help='number of threads to use')
misc_group.add_argument('-v','--verbose', dest='verbose', action='store_true',
                    help='increase verbosity')
args = parser.parse_args()

num_in_paths = 0 if not args.fasta_files else len(args.fasta_files)
num_in_paths = num_in_paths if not args.fastq_files else num_in_paths + len(args.fastq_files)
if num_in_paths == 0:
  util.fail( "Must supply query fasta/fastq files to be placed!" )


util.expect_file_exists( args.ref_tree )
util.expect_file_exists( args.ref_msa )
util.expect_file_exists( args.model_file )

if args.prefix:
  util.expect_dir_exists( args.prefix )

if args.taxon_file:
  util.expect_file_exists( args.taxon_file )

use_chunkify = not args.no_chunkify


# make a unique output dir, labeled by date and time
out_dir = "run-{}".format(datetime.now().strftime("%Y-%m-%d-%H:%M:%S")) if( not args.out_dir ) else args.out_dir

# build up a samples.tsv file to be used in this execution
# 
# first read in a unique list of all the desired files
file_paths = []
if args.fasta_files:
  file_paths.extend( util.ingest_paths( args.fasta_files,
                                        extensions=['.fa', '.afa', '.fasta'],
                                        allow_gz=True ) )

# also fetch the fastq files in the same way
fastq_file_paths = []
if args.fastq_files:
  fastq_file_paths.extend( util.ingest_paths( args.fastq_files,
                                        extensions=['.fq', '.fastq'],
                                        allow_gz=True ) )
if args.merge_pattern:
  # however these we need to treat differently: we find each F/R pair according to the pattern
  import re
  # clean up the pattern and stuff it in a regex
  # escape unescaped dots
  pattern = re.sub(r"\.", "\.", args.merge_pattern)
  reg = re.compile(rf"(.*){pattern}$")
  unmerged_files  = []
  unmerged_names  = []
  for f in fastq_file_paths:
    match = reg.search( f )
    if match:
      unmerged_files.append( f )
      unmerged_names.append( match.group(1) )
    else:
      file_paths.append( f )

  assert(len(unmerged_files) % 2 == 0)
else:
  file_paths.extend( fastq_file_paths )

sample_names = util.get_unique_names( file_paths ) if file_paths else []

if args.merge_pattern and unmerged_files:
  sample_names.extend( util.get_unique_names( unmerged_names ) )
  file_paths.extend( unmerged_files )

# add file paths to the samples, giving it a sample name corresponding to the file name / directory
samples = pd.DataFrame({
  'sample':sample_names,
  'input_file':file_paths
  }).set_index( 'sample' )

# finally, write the samples.tsv to the output folder
util.make_path( out_dir )
samples_file = os.path.join( out_dir, "samples.tsv" )
samples.to_csv( samples_file, sep='\t' )

# next, we set those config values that we wish to override from the defaults, by creating a dict
# of those values
config_overrrides = {
  'data':
  {
    'samples': samples_file,
    'reference-tree': args.ref_tree,
    'reference-alignment': args.ref_msa
  },
  'settings':
  {
    'datatype': args.datatype,
    'use-chunkify': use_chunkify,
    'outdir': out_dir,
  },
  'params':
  {
    'threads': args.threads
  }
}

if args.taxon_file:
  config_overrrides['data']['taxonomy-file'] = args.taxon_file
if args.sequence_clustering:
  config_overrrides['settings']['clustering-tool'] = args.sequence_clustering

calling_dir = os.path.dirname(os.path.abspath(__file__))

# check if mamba exists
conda_front = 'conda'
if util.is_tool('mamba'):
  conda_front = 'mamba'

# get cluster settings
cluster, cluster_config = (None, None) if not args.on_cluster else util.cluster_settings( args.clust_env, calling_dir )

snakemake.snakemake(
  snakefile=os.path.join( calling_dir, "Snakefile" ),
  use_conda=True,
  conda_frontend=conda_front,
  cores=args.threads,
  config=config_overrrides,
  cluster_config=cluster_config,
  cluster=cluster,
  rerun_triggers=["mtime"]
  )
