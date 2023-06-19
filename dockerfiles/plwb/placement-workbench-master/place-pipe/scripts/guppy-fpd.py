# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
from util import dirname
import snakeparser as sp

shell.executable("bash")

# Get the output directory
outdir = dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "guppy fpd", snakemake, ['params','guppy','fpd'] )

# General options
ps.add_opt( "pp",               sp.typ.FLAG )
ps.add_opt( "prefix" )
ps.add_opt( "csv",              sp.typ.FLAG )
ps.add_opt( "theta" )
ps.add_opt( "chao-d" )
ps.add_opt( "include-pendant",  sp.typ.FLAG )

# output
ps.add( os.path.basename(snakemake.output[0]), "-o {}" )
ps.add( outdir,                 "--out-dir {}" )

# Input positional argument at the end
ps.add( snakemake.input[0],    "{}",   sp.typ.FILE )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
