# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import snakeparser as sp
from tempfile import NamedTemporaryFile
from ete3 import Tree

def resolve( treefile ):
    t = Tree( treefile )
    t.standardize()
    return str(t.write())

shell.executable("bash")

# Get the output directory
outdir = util.dirname( snakemake.output[0] )

# =================================================================================================
#     File prep / aggregation
# =================================================================================================
# If there is only one file on the input, we assume its a file containing multiple newick trees
if len(snakemake.input) == 1:
    ml_trees = snakemake.input[0]
# If there are multiple, we assume that multiple newick files with single trees are specified
# In this case we combine them into one file first
else:
    tmp = NamedTemporaryFile( suffix="ml_trees.newick", mode='w' )
    ml_trees = tmp.name
    for in_path in snakemake.input:
        with open( in_path, 'r' ) as infile:
            for tree in infile.readlines():
                if len(tree.rstrip()) > 0:
                    tmp.write( f"{resolve( tree )}\n" )
    tmp.flush()                    

# =================================================================================================
#     Parse arguments
# =================================================================================================
if "prefix" in snakemake.params.keys():
    prefix = snakemake.params.prefix
else:
    prefix = os.path.join( outdir, "rf_calc" )

ps = sp.Parser("raxml-ng", snakemake, ['params','raxml-ng','rfdist'])

# select the run mode
ps.add( "--rfdist" )

# Required args
ps.add( ml_trees, "--tree {}", sp.typ.FILE )

# Optional args


# Closing args
ps.add( prefix, "--prefix {}" )

ps.add_threads()


# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

result_file = prefix + ".raxml.rfDistances"
util.expect_file_exists( result_file )

# rename the result file
shell( f"mv {result_file} {snakemake.output[0]}" )
