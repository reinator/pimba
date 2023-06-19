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
ps = sp.Parser( "mafft", snakemake )

# as per https://mafft.cbrc.jp/alignment/software/manual/manual.html
# OPTIONS

# Algorithm
ps.add_opt( "auto",             sp.typ.FLAG )
ps.add_opt( "6merpair",         sp.typ.FLAG )
ps.add_opt( "globalpair",       sp.typ.FLAG )
ps.add_opt( "localpair",        sp.typ.FLAG )
ps.add_opt( "genafpair",        sp.typ.FLAG )
ps.add_opt( "fastapair",        sp.typ.FLAG )
ps.add_opt( "weighti",          sp.typ.UINT )
ps.add_opt( "retree",           sp.typ.UINT )
ps.add_opt( "maxiterate",       sp.typ.UINT )
ps.add_opt( "fft",              sp.typ.FLAG )
ps.add_opt( "nofft",            sp.typ.FLAG )
ps.add_opt( "noscore",          sp.typ.FLAG )
ps.add_opt( "memsave",          sp.typ.FLAG )
ps.add_opt( "parttree",         sp.typ.FLAG )
ps.add_opt( "dpparttree",       sp.typ.FLAG )
ps.add_opt( "fastaparttree",    sp.typ.FLAG )
ps.add_opt( "partsize",         sp.typ.UINT )
ps.add_opt( "groupsize",        sp.typ.UINT )

# Parameter
ps.add_opt( "op",               sp.typ.FLOAT )
ps.add_opt( "ep",               sp.typ.FLOAT )
ps.add_opt( "lop",              sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "lep",              sp.typ.FLOAT )
ps.add_opt( "lexp",             sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "LOP",              sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "LEXP",             sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "bl",               sp.typ.UINT )
ps.add_opt( "jtt",              sp.typ.UINT )
ps.add_opt( "tm",               sp.typ.UINT )
ps.add_opt( "aamatrix",         sp.typ.FILE )
ps.add_opt( "fmodel",           sp.typ.FLAG )

# Output
ps.add_opt( "clustalout",       sp.typ.FLAG )
ps.add_opt( "inputorder",       sp.typ.FLAG )
ps.add_opt( "reorder",          sp.typ.FLAG )
ps.add_opt( "treeout",          sp.typ.FLAG )
ps.add_opt( "quiet",            sp.typ.FLAG )

# Input
ps.add_opt( "nuc",              sp.typ.FLAG )
ps.add_opt( "amino",            sp.typ.FLAG )
ps.add_opt( "seed",             sp.typ.FLAG )

# Runtime / Other (undocumented?)
ps.add_opt( "dash",             sp.typ.FLAG )
ps.add_threads("--thread {}")

# Input/Output (they're positional)
ps.add( snakemake.input[0],     valid_func=sp.typ.FILE )
ps.add( ">" )
ps.add( snakemake.output[0] )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string( log_stdout=False ) )
