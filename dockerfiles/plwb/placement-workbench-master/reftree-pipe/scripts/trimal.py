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
ps = sp.Parser( "trimal", snakemake, arg_ident='-' )

# Input/Output
ps.add( snakemake.input[0],     "-in {}", sp.typ.FILE )
ps.add( snakemake.output[0],    "-out {}" )
ps.add( "-fasta" )

# Common options
ps.add_opt( "backtrans",        sp.typ.FILE )
ps.add_opt( "ignorestopcodon",  sp.typ.FLAG )
ps.add_opt( "splitbystopcodon", sp.typ.FLAG )

ps.add_opt( "matrix",           sp.typ.FILE )

ps.add_opt( "htmlout" )

ps.add_opt( "keepheader",       sp.typ.FLAG )

ps.add_opt( "complementary",    sp.typ.FLAG )
ps.add_opt( "colnumbering",     sp.typ.FLAG )

ps.add_opt( "selectcols" )
ps.add_opt( "selectseqs" )

ps.add_opt( "gapthreshold",     sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "simthreshold",     sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "conthreshold",     sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "cons",             sp.typ.UINT(0,100) )

ps.add_opt( "nogaps",           sp.typ.FLAG )
ps.add_opt( "noallgaps",        sp.typ.FLAG )
ps.add_opt( "keepseqs",         sp.typ.FLAG )

ps.add_opt( "gappyout",         sp.typ.FLAG )
ps.add_opt( "strict",           sp.typ.FLAG )
ps.add_opt( "strictplus",       sp.typ.FLAG )

ps.add_opt( "automated1",       sp.typ.FLAG )

ps.add_opt( "terminalonly",     sp.typ.FLAG )
ps.add_opt( "block",            sp.typ.UINT )

ps.add_opt( "resoverlap",       sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "seqoverlap",       sp.typ.UINT(0,100) )

ps.add_opt( "clusters",         sp.typ.UINT(1) )
ps.add_opt( "maxidentity",      sp.typ.FLOAT(0.0,1.0) )

ps.add_opt( "w",                sp.typ.UINT(1) )
ps.add_opt( "gw",               sp.typ.UINT(1) )
ps.add_opt( "sw",               sp.typ.UINT(1) )
ps.add_opt( "cw",               sp.typ.UINT(1) )

ps.add_opt( "sgc",              sp.typ.FLAG )
ps.add_opt( "sgt",              sp.typ.FLAG )
ps.add_opt( "ssc",              sp.typ.FLAG )
ps.add_opt( "sst",              sp.typ.FLAG )
ps.add_opt( "sfc",              sp.typ.FLAG )
ps.add_opt( "sft",              sp.typ.FLAG )
ps.add_opt( "sident",           sp.typ.FLAG )


# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
