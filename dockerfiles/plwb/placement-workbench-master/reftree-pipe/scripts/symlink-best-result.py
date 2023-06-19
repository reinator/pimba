from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
from util import trim_path, dirname

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell.executable("bash")

treedir = "raxml-ng/tree"
postdir = "raxml-ng/post"

# So turns out symlinking like I wanted to doesnt really work. The output to the rule has to be the newick file
# such that the rules using that file can actually find it. But that means that snakemake "helpfully" pre creates
# the best_result/raxml... directory structure it sees in the output, which means that making the symlink fails.
# dirtiest hack I can think of to still make it work: delete that directory before the symlink is created

def best_of_the_best( wildcards, input, output ):
    """Compares a set of runs and returns the folder of the one with the highest loglh"""
    import re
    dirs = [trim_path( dirname(f), [treedir, postdir]) for f in input]

    # keep only unique dir paths, such that we can allow multiple redundant inputs
    # (for correct dependency resolution in snakemake)
    dirs = list(set(dirs))

    result_dict = dict()
    for rundir in dirs:
        with open( os.path.join( rundir, treedir, "search.raxml.log" ), 'r') as logfile:
            for line in logfile:
                if re.search( "Final LogLikelihood: ", line ):
                    loglh = line.split(' ')[-1]
                    result_dict[rundir] = float(loglh)
    best_dir = max( result_dict, key=result_dict.get )
    return best_dir

source = str(best_of_the_best(snakemake.wildcards, snakemake.input, snakemake.output))

# the output of the calling rule is expected to be the path to the actual best tree of the best run
# but we care only about the symlink to the best result dir, for which we need the symlink name
# and we get that from how the output file is specified
link            = str( trim_path( dirname( snakemake.output[0] ), treedir ) )
rel_source      = os.path.relpath( source, os.path.dirname(link) )

# snakemake will have created a directory here from the output, so we need to end that directories whole career
shell( f"rm -r {link}" )

shell( f"ln -s {rel_source} {link} {log}" )
