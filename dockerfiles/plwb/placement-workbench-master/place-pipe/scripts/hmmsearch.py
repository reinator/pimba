# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from encodings import search_function
from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp

shell.executable("bash")

# the actual output of hmmsearch is a tsv file, which we parse after
outdir = os.path.dirname( snakemake.output[0] )
hits_file = os.path.join( outdir, "hits.tsv" )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "hmmsearch", snakemake, ['params','hmmer','hmmsearch'] )

# Options directing output:
# ps.add( log_file,    "-o {}" )
# ps.add( msa_out,    "-A {}" )
ps.add( hits_file,   "--tblout {}" )
# ps.add( hits_file, "--domtblout {}" )
# ps.add( hits_file, "--pfamtblout {}" )
ps.add_opt( "acc",          sp.typ.FLAG )
ps.add_opt( "noali",        sp.typ.FLAG )
ps.add_opt( "notextw",      sp.typ.FLAG )
ps.add_opt( "textw",        sp.typ.UINT(120) )

# Options controlling reporting thresholds:
ps.add_opt( "E",            sp.typ.FLOAT, "-E {}" )
ps.add_opt( "T",            sp.typ.FLOAT, "-T {}" )
ps.add_opt( "domE",         sp.typ.FLOAT )
ps.add_opt( "domT",         sp.typ.FLOAT )

# Options controlling inclusion (significance) thresholds:
ps.add_opt( "incE",         sp.typ.FLOAT )
ps.add_opt( "incT",         sp.typ.FLOAT )
ps.add_opt( "incdomE",      sp.typ.FLOAT )
ps.add_opt( "incdomT",      sp.typ.FLOAT )

# Options controlling model-specific thresholding:
ps.add_opt( "cut_ga",       sp.typ.FLAG )
ps.add_opt( "cut_nc",       sp.typ.FLAG )
ps.add_opt( "cut_tc",       sp.typ.FLAG )

# Options controlling acceleration heuristics:4
ps.add_opt( "max",          sp.typ.FLAG )
ps.add_opt( "F1",           sp.typ.FLOAT )
ps.add_opt( "F2",           sp.typ.FLOAT )
ps.add_opt( "F3",           sp.typ.FLOAT )
ps.add_opt( "nobias",       sp.typ.FLAG )

# Other expert options:
ps.add_opt( "nonull2",      sp.typ.FLAG )
ps.add_opt( "Z",            sp.typ.UINT, "-Z {}" )
ps.add_opt( "domZ",         sp.typ.UINT )
ps.add_opt( "seed",         sp.typ.UINT )
ps.add_opt( "tformat" )
ps.add_threads( "--cpu {}" )

# positional arguments
ps.add( snakemake.input.hmmfile, "{}",    sp.typ.FILE )
ps.add( snakemake.input.seqfile, "{}",    sp.typ.FILE )

# =================================================================================================
#     Run the search and make a proper output file
# =================================================================================================

shell( ps.get_shell_string() )

outfile = snakemake.output[0]

shell( f"seqtk subseq {snakemake.input.seqfile} <(grep -v '^#' {hits_file} | awk '{{{{print $1}}}}') > {outfile}" )
