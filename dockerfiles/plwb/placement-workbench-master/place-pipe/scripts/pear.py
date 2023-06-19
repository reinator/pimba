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
ps = sp.Parser( "pear", snakemake )

# input args
ps.add( snakemake.input[0], "--forward-fastq {}", sp.typ.FILE )
ps.add( snakemake.input[1], "--reverse-fastq {}", sp.typ.FILE )

# Output args
# ouput arg to pear actually works like the prefix arg in raxml-ng
# it also outputs a fastq file, so we need to convert it later
name = "merge"
prefix = os.path.join( os.path.dirname( snakemake.output[0] ), name )
ps.add( prefix, "--output {}" )

result_file = f"{prefix}.assembled.fastq"

# General options
ps.add_opt( "p-value",              sp.typ.FLOAT(0.0001, 1.0) )
ps.add_opt( "min-overlap",          sp.typ.UINT )
ps.add_opt( "max-assembly-length",  sp.typ.UINT )
ps.add_opt( "min-assembly-length",  sp.typ.UINT )
ps.add_opt( "min-trim-length",      sp.typ.UINT )
ps.add_opt( "quality-threshold",    sp.typ.UINT )
ps.add_opt( "max-uncalled-base",    sp.typ.FLOAT(0.0, 1.0) )
ps.add_opt( "test-method",          sp.typ.IN([1, 2]) )
ps.add_opt( "empirical-freqs",      sp.typ.FLAG )
ps.add_opt( "score-method",         sp.typ.IN([1, 2, 3]) )
ps.add_opt( "phred-base",           sp.typ.UINT )
ps.add_opt( "memory",               sp.typ.UINT )
ps.add_opt( "cap",                  sp.typ.UINT )
ps.add_opt( "nbase",                sp.typ.FLAG )

ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

# convert the result file to fasta
shell( f"seqtk seq -A {result_file} > {snakemake.output[0]}" )
