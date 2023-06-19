# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import snakeparser as sp


shell.executable("bash")

# Get the output directory
sample_outdir = util.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "gappa prepare phat", snakemake, ['params','gappa','phat'] )

# Required args
ps.add( snakemake.input.taxonomy_file,  "--taxonomy-file {}",   sp.typ.FILE )
ps.add( snakemake.input.sequence_file,  "--sequence-file {}",   sp.typ.FILE )
ps.add( snakemake.config['params']['gappa']['phat']['target_size'],
    "--target-size {}",
    sp.typ.UINT )

# Optional args
ps.add_opt( "sub-taxonomy" )
ps.add_opt( "min-subclade-size",        sp.typ.UINT )
ps.add_opt( "max-subclade-size",        sp.typ.UINT )
ps.add_opt( "min-tax-level",            sp.typ.UINT )
ps.add_opt( "allow-approximation",      sp.typ.FLAG )
ps.add_opt( "no-taxa-selection",        sp.typ.FLAG )
ps.add_opt( "consensus-method",         sp.typ.IN(['majorities','cavener','threshold']) )
ps.add_opt( "consensus-threshold",      sp.typ.FLOAT(0.0, 1.0) )

# Closing args
ps.add( sample_outdir, "--out-dir {}", sp.typ.DIR )
ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

result_file = os.path.join( sample_outdir, "consensus_sequences.fasta" )
util.expect_file_exists( result_file )

# also, move the result file to what is expected by the output
shell( "mv {} {}".format( 
    result_file, 
    snakemake.output[0] )
)
