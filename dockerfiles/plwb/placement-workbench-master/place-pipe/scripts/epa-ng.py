# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
from util import expect_file_exists
import snakeparser as sp
from tempfile import TemporaryDirectory

shell.executable("bash")

# =================================================================================================
#     Parse arguments
# =================================================================================================

ps = sp.Parser( "epa-ng", snakemake )

# Options
ps.add_opt( "verbose", sp.typ.FLAG )

# Input
# other inputs delayed until after split, see below
ps.add( snakemake.input.tree,   "--tree {}",    sp.typ.FILE )
ps.add( snakemake.params.model, "--model {}" )

# Output
ps.add_opt( "filter-acc-lwr",   sp.typ.FLOAT(0.0, 1.0) )
ps.add_opt( "filter-min-lwr",   sp.typ.FLOAT(0.0, 1.0) )
ps.add_opt( "filter-min",       sp.typ.UINT(1) )
ps.add_opt( "filter-max",       sp.typ.UINT(1) )
ps.add_opt( "precision",        sp.typ.UINT(1) )
ps.add_opt( "redo",             sp.typ.FLAG )
ps.add_opt( "preserve-rooting", sp.typ.IN(['on','off']) )

# Compute
ps.add_opt( "dyn-heur",         sp.typ.FLOAT(0.0, 1.0) )
ps.add_opt( "fix-heur",         sp.typ.FLOAT(0.0, 1.0) )
ps.add_opt( "baseball-heur",    sp.typ.FLAG )
ps.add_opt( "no-heur",          sp.typ.FLAG )
ps.add_opt( "chunk-size",       sp.typ.UINT(1) )
ps.add_opt( "raxml-blo",        sp.typ.FLAG )
ps.add_opt( "no-pre-mask",      sp.typ.FLAG )
ps.add_opt( "rate-scalers",     sp.typ.IN(['on','off','auto']) )
ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

# epa-ng always uses the same output file names, so for parallel instances of this rule,
# we have to use temp directories where we can do our work.
with TemporaryDirectory() as tempdir:
    # We need to split the combined query+ref msa into individual (temp) files.
    split = sp.Parser( "epa-ng", snakemake )
    split.add( "--split" )
    split.add( snakemake.input.msa,         valid_func=sp.typ.FILE )
    split.add( snakemake.input.query,       valid_func=sp.typ.FILE )
    split.add( tempdir, "--out-dir {}", sp.typ.DIR )
    shell( split.get_shell_string() )

    # The output of this step are two files:
    # query.fasta
    # reference.fasta
    # which we use in the following

    # Now run the actual placement. We here do not use the original reference alignment,
    # but the one resulting from the above split step instead, as this is guaranteed to have
    # the same width as the query alignment. Unfortunately, some aligners (such as version 3 of
    # hmmalign) mess around with gaps, so that the produced alignment does not fit with the
    # original any more. Hence, we cannot use the original here, and use the split one instead.

    ps.add( f"{tempdir}/reference.fasta",  "--msa {}",     sp.typ.FILE )
    ps.add( f"{tempdir}/query.fasta",      "--query {}",   sp.typ.FILE )
    ps.add( tempdir,                       "--out-dir {}", sp.typ.DIR )

    shell( ps.get_shell_string() )

    # Finally, move the result files that we care about to our actual output
    jplace_result = f"{tempdir}/epa_result.jplace"
    expect_file_exists( jplace_result )
    shell( f"mv {jplace_result} {snakemake.output.jplace}" )
