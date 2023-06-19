from snakemake.shell import shell
import sys, os
import util

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell.executable("bash")

"""
Symlinks the input to the output, with correct relative paths
"""
if len(sys.argv) > 1:
    in_path = sys.argv[1]
else:
    in_path = str(snakemake.input[0])

if len(sys.argv) > 2:
    out_path = sys.argv[2]
else:
    out_path = str(snakemake.output[0])

rel_input = os.path.relpath( in_path, os.path.dirname( out_path ) )

shell( f'ln -s "{rel_input}" "{out_path}" {log}' )
