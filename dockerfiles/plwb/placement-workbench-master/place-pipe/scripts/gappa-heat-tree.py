# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
from tempfile import TemporaryDirectory
import os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp
from util import fail, dirname, extension, filename

shell.executable("bash")

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "gappa examine heat-tree", snakemake, ['params','gappa','heat-tree'] )

# Get all files types of trees that we want to produce by their extensions in the output files.
# We also get the outdir alongside, and check if all outdirs are identical
exts = []
outdir = None
for file in snakemake.output:
    exts.append( extension( file ).split('.')[-1] )
    if outdir:
        if outdir != dirname(file):
            fail(
                "Output file paths for different file types in gappa heat-tree rule have to "
                "be identical except for their file extension."
            )
    else:
        outdir = dirname(file)

# Create command line arguments and check that the extensions are valid.
for ext in exts:
    if ext not in [ "newick", "nexus", "phyloxml", "svg" ]:
        fail(f"Unknown gappa heat-tree file extension: {ext}")
    ps.add( f"--write-{ext}-tree" )

# Normal options
# Input
ps.add( snakemake.input, "--jplace-path {}", sp.typ.FILES )

# Settings
ps.add_opt( "point-mass",               sp.typ.FLAG )
ps.add_opt( "ignore-multiplicities",    sp.typ.FLAG )

# Color
ps.add_opt( "color-list" )
ps.add_opt( "reverse-color",    sp.typ.FLAG )
ps.add_opt( "under-color",      sp.typ.STRING )
ps.add_opt( "clip-under",       sp.typ.FLAG )
ps.add_opt( "over-color",       sp.typ.STRING )
ps.add_opt( "clip-over",        sp.typ.FLAG )
ps.add_opt( "clip",             sp.typ.FLAG )
ps.add_opt( "mask-color",       sp.typ.STRING )
ps.add_opt( "log-scaling",      sp.typ.FLAG )
ps.add_opt( "min-value",        sp.typ.FLOAT )
ps.add_opt( "max-value",        sp.typ.FLOAT )
ps.add_opt( "mask-value",       sp.typ.FLOAT )

# Output
ps.add_opt( "file-prefix" )
ps.add_opt( "file-suffix" )

# Tree Output
# --write-x-tree arguments handled seperatly (check above)

# SVG Tree Output
ps.add_opt( "svg-tree-shape",           sp.typ.IN(['circular','rectangular']) )
ps.add_opt( "svg-tree-type",            sp.typ.IN(['cladogram','phylogram']) )
ps.add_opt( "svg-tree-stroke-width",    sp.typ.FLOAT )
ps.add_opt( "svg-tree-ladderize",       sp.typ.FLAG )

# Global Options
ps.add_opt( "allow-file-overwriting",   sp.typ.FLAG )
ps.add_opt( "verbose",                  sp.typ.FLAG )
ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

# gappa always uses the same output file names, so for independent instances of this rule,
# we have to use temp directories where we can do our work.
with TemporaryDirectory() as tempdir:

    ps.add( tempdir, "--out-dir {}", sp.typ.DIR )

    shell( ps.get_shell_string() )

    # Move the result files that we care about to our actual output
    out_filename = os.path.join( outdir, filename( snakemake.output[0] ) )
    for ext in exts:
        shell( f"mv {tempdir}/tree.{ext} {out_filename}.{ext}" )
