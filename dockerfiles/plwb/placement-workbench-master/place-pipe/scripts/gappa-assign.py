# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp

shell.executable("bash")

# Get the output directory
sample_outdir = os.path.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "gappa examine assign", snakemake, ['params','gappa','assign'] )

# Required args
ps.add( snakemake.input.jplace,     "--jplace-path {}", sp.typ.FILE )
ps.add( snakemake.input.taxon_file, "--taxon-file {}",  sp.typ.FILE )

# Optional args
ps.add_opt( "root-outgroup",            sp.typ.FILE )
ps.add_opt( "taxonomy",                 sp.typ.FILE )
ps.add_opt( "ranks-string" )
ps.add_opt( "sub-taxopath" )
ps.add_opt( "max-level",                sp.typ.UINT )
ps.add_opt( "distribution-ratio",       sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "consensus-thresh",         sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "resolve-missing-paths",    sp.typ.FLAG )
ps.add_opt( "distant-label",            sp.typ.FLAG )

# Output args
ps.add_opt( "file-prefix" )
ps.add_opt( "file-suffix" )
ps.add_opt( "cami",                     sp.typ.FLAG )
ps.add_opt( "krona",                    sp.typ.FLAG )
ps.add_opt( "sativa",                   sp.typ.FLAG )
ps.add_opt( "sample-id" )
ps.add_opt( "per-query-results",        sp.typ.FLAG )
ps.add_opt( "best-hit",                 sp.typ.FLAG )
ps.add_opt( "allow-file-overwriting",   sp.typ.FLAG )
ps.add_opt( "verbose",                  sp.typ.FLAG )

# Closing args
ps.add( sample_outdir, "--out-dir {}" , sp.typ.DIR )
ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
