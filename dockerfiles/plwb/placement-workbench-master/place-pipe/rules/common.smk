# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, re, sys
import socket, platform
common_dir = os.path.abspath(os.path.join( workflow.current_basedir, "..", "..", "common" ))
sys.path.insert(0, common_dir)
from util import is_fasta, is_fastq, expect_file_exists, config_to_file, fail

# Ensure min Snakemake version
snakemake.utils.min_version("5.7")

# =================================================================================================
#     Basic Configuration
# =================================================================================================

# Load the config. If --directory was provided, this is also loaded from there.
# This is useful to have runs that have different settings, but generally re-use the main setup.
configfile: "config.yaml"
# snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# =================================================================================================
#     Get Samples List
# =================================================================================================

def is_existing_fastq( paths: list ):
    for path in paths:
        expect_file_exists( path )
        if not is_fastq( path ):
            return False
    return True


# Prepare a dictionary of samples, from name (used for our wildcards here, and hence for file
# naming within the pipeline) to absolute file paths.
samples = pd.read_table(config["data"]["samples"], dtype=str).set_index(["sample"], drop=False)

# Get just the sample names, to use as a list of wildcards later
sample_names=list(set(samples.index.get_level_values("sample")))
assert(len(sample_names) > 0)
# next we do some quick sanity check as we allow samples to have multiple files, if those files are 
# unmerged fastq files
pd_df_type = pd.core.series.Series
for n in sample_names:
    entry = samples.loc[n, "input_file"]
    if type(entry) is pd_df_type:
        if len(entry) != 2:
            fail( f"Incorrect number ({len(entry)}) of files specified for sample {n}")
        elif not is_existing_fastq( list(entry) ):
            fail( "Specifying two input files for a sample is only allowed for unmerged fastq files!")
    elif type(entry) is str:
        expect_file_exists( entry )
    else:
        fail( f"Sample {n} has a file of unknown type {type(entry)}!")

# List of used clustering approaches
clusterer_list = config["settings"]["clustering-tool"]

# output prefix
outdir=config["settings"]["outdir"].rstrip("/")

# helper bool for model vs evaluate
use_evaluate = bool( config["params"]["epa-ng"]["model-params"] == "" )

# helper bool for krona
make_krona_plots = ('krona' in config['params']['gappa']['assign'].keys()) and bool(config['params']['gappa']['assign']['krona'])

# persist the used config
config_to_file( config, outdir )

hmmer_datatype_string  = "dna" if config["settings"]["datatype"] == 'nt' else "amino"

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# Get a nicely formatted hostname
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# Some helpful messages
logger.info("===========================================================================")
logger.info("    place-pipe")
logger.info("")
logger.info("    Host:               " + hostname)
logger.info("    Snakefile:          " + (workflow.snakefile))
logger.info("    Base directory:     " + (workflow.basedir))
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Config files:       " + (", ".join(workflow.configfiles)))
logger.info("    Samples:            " + str(len(sample_names)))
logger.info("===========================================================================")
logger.info("")

# =================================================================================================
#     Common File Access Functions
# =================================================================================================

def get_sample_fasta( wildcards ):
    """Get fasta file for a given sample""" 
    entry = samples.loc[wildcards.sample, "input_file"]
    
    if type(entry) is pd_df_type:
        return rules.merge_paired_pear.output
    elif is_fastq( entry ):
        return rules.merged_fastq_to_fasta.output
    elif is_fasta( entry ):
        return entry
    else:
        fail(f"Unknown type for sample file? ({type(entry)})")

def get_sample_fastq( wildcards ):
    """Get fastq file for a given sample"""
    entry = samples.loc[wildcards.sample, "input_file"]

    if type(entry) is pd_df_type:
        return sorted(list(entry))
    else:
        assert( is_fastq( entry ) )
        return entry

def relative_input_path( wildcards, input, output ):
    """Returns the relative path to the input file, from the directory of the output file/directory"""
    return os.path.relpath( str(input), os.path.dirname( str(output) ) )

# def get_all_sample_names():
#     """Get the samples names given to all fasta files"""
#     return list(samples.keys())

# def get_all_sample_paths():
#     """Get the paths to all fasta files"""
#     return list(samples.values())

# =================================================================================================
#     Config Related Functions
# =================================================================================================

from operator import getitem
from functools import reduce

def get_highest_override( tool, key ):
    """From the config, get the value labeled with "key", unless the "tool" overrides that value,
    in which case fetch the override"""

    params = config['params']
    if type(tool) is not list: tool = [ tool ]

    try:
        # try to fetch the key, uder the config path for tool
        return reduce( getitem, tool + [key], params )
    except KeyError:
        # if that fails, try to find the key under the base params
        if key in params:
            return params[key]
        else:
            fail("invalid key for 'config['params']': '{}'".format( key ))

def get_threads( tool ):
    return int( get_highest_override( tool, "threads") )

