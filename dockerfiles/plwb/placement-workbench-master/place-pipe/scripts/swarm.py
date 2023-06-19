# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp
from tempfile import NamedTemporaryFile

shell.executable("bash")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get the output directory
outdir = os.path.dirname( snakemake.output[0] )

# swarm expects the input to be dereplicated
derep_file = NamedTemporaryFile( suffix="dereplicated.fa" )
shell( f"swarm --append-abundance 1 -d 0 -w {derep_file.name} -o /dev/null {snakemake.input[0]} {log}" )

# also, trim out any N's from the input fasta file
stripped_file = NamedTemporaryFile( suffix="stripped.fa" )
shell( f"sed '/^>/ ! s/[^ACGTacgt]//g' \"{derep_file.name}\" > {stripped_file.name}" )


# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "swarm", snakemake )

# Output args
ps.add( snakemake.output.seeds,             "--seeds {}" )
ps.add( snakemake.output.statistics_file,   "--statistics-file {}" )

# General options
ps.add_threads()

# Clustering options
ps.add_opt( "differences",              sp.typ.UINT )
ps.add_opt( "no-otu-breaking",          sp.typ.FLAG )
# Fastidious options
ps.add_opt( "boundary",                 sp.typ.UINT )
ps.add_opt( "ceiling",                  sp.typ.UINT )
ps.add_opt( "fastidious",               sp.typ.FLAG )
ps.add_opt( "bloom-bits",               sp.typ.UINT )

# (other) output options
ps.add_opt( "append-abundance",         sp.typ.UINT )
ps.add_opt( "internal-structure" )
ps.add_opt( "network-file" )
ps.add_opt( "output-file" )
ps.add_opt( "mothur",                   sp.typ.FLAG )
ps.add_opt( "uclust-file" )
ps.add_opt( "usearch-abundance",        sp.typ.FLAG )


# Pairwise alignment advanced options
ps.add_opt( "match-reward",             sp.typ.UINT )
ps.add_opt( "mismatch-penalty",         sp.typ.UINT )
ps.add_opt( "gap-opening-penalty",      sp.typ.UINT )
ps.add_opt( "gap-extension-penalty",    sp.typ.UINT )

# Required fasta input file at the end
ps.add( stripped_file.name, "{}",    sp.typ.FILE )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
