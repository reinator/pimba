# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp

shell.executable("bash")

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "hmmalign", snakemake, ['params','hmmer','hmmalign'] )

# General options
ps.add_opt( "trim",             sp.typ.FLAG )
# special handling of the --amino/--dna/--rna flags
ps.add( snakemake.params.states, "--{}" )
ps.add_opt( "informat" )
ps.add_opt( "outformat" )

# output
ps.add( snakemake.output[0],    "-o {}" )

# Input
ps.add( snakemake.input.msa,    "--mapali {}",  sp.typ.FILE )

# positional arguments at the end
ps.add( snakemake.input.hmmfile,"{}",           sp.typ.FILE )
ps.add( snakemake.input.seqfile,"{}",           sp.typ.FILE )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
