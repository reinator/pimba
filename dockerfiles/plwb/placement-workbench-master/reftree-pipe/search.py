#!/usr/bin/env python3

import multiprocessing
import argparse, sys, os
from os.path import join
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, join(script_dir, '..', 'common'))
import util
from datetime import datetime
import pandas as pd
import snakemake
import platform

# if we're on mac/arm, default back down to x86 compiled packages
# ideally this decision would be per-package availability
if( platform.machine() == 'arm64' ):
  os.environ['CONDA_SUBDIR'] = "osx-64"

parser = argparse.ArgumentParser(description='Wrapper to run the pipeline with some default settings.')

input_group = parser.add_argument_group('Input')
input_group.add_argument('--fasta-paths', dest='fasta_files', type=str, nargs='+',
                    help='input fasta files')
input_group.add_argument('--csv-paths', dest='accession_files', type=str, nargs='+',
                    help="""input .csv files specifying accessions to be downloaded. These files must contain 
                    at least two columns: 'accession' for the accessions, and 'label' for the name/label associated
                    with an accession. The label will subsequently be used in the fasta/tree files. Note that the column
                    labels can be customized when using the pipeline directly, using the config.yaml file.
                    """)
input_group.add_argument('--constraint-tree', dest='constraint_tree', type=str,
                    help='Optional constraint tree, in newick format')

pipeline_group = parser.add_argument_group('Pipeline Options')
pipeline_group.add_argument('--out-dir', dest='out_dir', type=str,
                    default=None,
                    help='optional output directory. By default, the script creates a timestamped output directory.')
pipeline_group.add_argument('--phat', dest='do_phat', action='store_true',
                    help="""enable the optional PhAT algorithm, reducing the number of taxa down to a target number, 
                    based on a given taxonomy. The taxonomy file is expected to exist alongside the input fasta/csv file.""")
pipeline_group.add_argument('--phat-target-num', dest='phat_target_num', type=int,
                    default=512,
                    help='target number of taxa that the PhAT algorithm should aim for.')
pipeline_group.add_argument('--taxonomy-file', dest='taxonomy_file', type=str,
                    help='taxonomy file (required for PhAT algorithm)')
pipeline_group.add_argument('--compatible-trees', dest='trees_compatible', action='store_true',
                    help="Indicate that the resulting trees will be compatible (same number of taxa, same labels)"
                    "WARNING: if this turnes out to be false, the run will fail at the 'rf_distances_between_samples' rule")
pipeline_group.add_argument('-a','--align', dest='do_align', action='store_true',
                    help='align the sequences')
pipeline_group.add_argument('--modeltest', dest='do_modeltest', action='store_true',
                    help='automatically find the best model')
pipeline_group.add_argument('-d', '--datatype', dest='datatype', type=str, nargs='?',
                    const='nt', default='nt', choices=['nt', 'aa'],
                    help="datatype, 'aa' for protein, 'nt' for DNA data")
pipeline_group.add_argument('-p', '--prefix', dest='prefix', type=str,
                    default=None,
                    help='prefix to fasta paths and output (useful to specify where data was mounted to in docker)')

cluster_group = parser.add_argument_group('Cluster Execution')
cluster_group.add_argument('--cluster-exec', dest='on_cluster', action='store_true',
                    help="Starts the pipeline in computing-cluster (slurm/sge etc.) submission mode. "
                    "Highly recommended to do this from a screen/tmux session!")
cluster_group.add_argument('--cluster-env', dest='clust_env', type=str, nargs='?',
                    const='auto', default='auto', choices=['auto','slurm', 'sge'],
                    help="What job submission system we are on. 'auto' attempts to autodetect.")
cluster_group.add_argument('--nodes', dest='nodes', type=int,
                    default=1,
                    help='number of nodes to use')

misc_group = parser.add_argument_group('Misc. Options')
misc_group.add_argument('--threads', dest='threads', type=int,
                    default=multiprocessing.cpu_count(),
                    help='number of threads to use')
misc_group.add_argument('-v','--verbose', dest='verbose', action='store_true',
                    help='increase verbosity')
args = parser.parse_args()

# input validation
if( args.prefix ):
  util.expect_dir_exists( args.prefix )

if args.do_phat:
    if not args.taxonomy_file:
        util.fail( "When using PhAT, must also supply a valid taxonomy file" )
    else:
        util.expect_file_exists( args.taxonomy_file )
    
    if args.do_align:
        util.fail( "It does not make sense to combine PhAT with alignment, as PhAT requires aligned sequences already." )

skip_alignment = not args.do_align

# make a unique output dir, labeled by date and time
out_dir = "run-{}".format(datetime.now().strftime("%Y-%m-%d-%H:%M:%S")) if( not args.out_dir ) else args.out_dir

# build up a samples.tsv file to be used in this execution
# 
# first read in a unique list of all the desired files
file_paths = []
if( args.fasta_files ):
  file_paths.extend( util.ingest_paths( args.fasta_files, extensions=['.fa', '.afa', '.fasta', '.gz'] ) )
if( args.accession_files ):
  file_paths.extend( util.ingest_paths( args.accession_files, extensions=['.csv'] ) )

if not file_paths:
  util.fail( "Did not detect any input fasta/csv files. Wrong directory?" )

# add file paths to the samples, giving it a sample name corresponding to the file name / directory
samples = pd.DataFrame({
  'sample':util.get_unique_names(file_paths),
  'input_file':file_paths
  }).set_index( 'sample' )

# finally, write the samples.tsv to the output folder
util.make_path( out_dir )
samples_file = join( out_dir, "samples.tsv" )
samples.to_csv( samples_file, sep='\t' )

# check if we want to use modeltest, and set the model string to auto if thats the case
model_string = "" # empty means defaults are used
if args.do_modeltest:
    model_string = "auto"

# next, we set those config values that we wish to override from the defaults, by creating a dict
# of those values
config_overrrides = {
  'data':
  {
    'samples': samples_file,
    'trees_are_compatible': args.trees_compatible
  },
  'settings':
  {
    'skip_autoref': not args.do_phat,
    'autoref': "phat", # always set phat as an autoref tool. ignored internally on skip
    'skip_alignment': skip_alignment,
    'outdir': out_dir,
    'datatype': args.datatype
  },
  'params':
  {
    'threads': args.threads,
    'model': model_string,
    'gappa':
    {
      'phat':
      {
        'target-size': args.phat_target_num
      }
    }
  }
}

if args.constraint_tree:
  config_overrrides['data']['constraint_tree'] = args.constraint_tree
if args.taxonomy_file:
  config_overrrides['data']['taxonomy'] = args.taxonomy_file

calling_dir = os.path.dirname(os.path.abspath(__file__))

# check if mamba exists
conda_front = 'conda'
if util.is_tool('mamba'):
  conda_front = 'mamba'

# get cluster settings
cluster, cluster_config = (None, None) if not args.on_cluster else util.cluster_settings( args.clust_env, calling_dir )

snakemake.snakemake(
  snakefile=join( calling_dir, "Snakefile" ),
  use_conda=True,
  conda_frontend=conda_front,
  cores=args.threads*args.nodes,
  local_cores=args.threads,
  config=config_overrrides,
  latency_wait=3,
  cluster_config=cluster_config,
  cluster=cluster,
  force_incomplete=True
  )

