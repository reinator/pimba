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
outdir = util.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "raxml-ng", snakemake, ['params','raxml-ng','treesearch'] )

# select the run mode
ps.add( "--all" )

# Input and Output Options
ps.add( snakemake.input.msa,    "--msa {}",     sp.typ.FILE )
ps.add( snakemake.params.model, "--model {}" )

# optional tree constraint input file if specified
if 'tree_constraint' in snakemake.input.keys() and snakemake.input.tree_constraint != []:
    ps.add( snakemake.input.tree_constraint, "--tree-constraint {}", sp.typ.FILE )

# special handling of starting trees
pars_trees = snakemake.params.pars_trees
rand_trees = snakemake.params.rand_trees

trees = []
if pars_trees:
    trees.append( f"pars{{{{{pars_trees}}}}}" )
if rand_trees:
    trees.append( f"rand{{{{{rand_trees}}}}}" )
starting_trees = ",".join(trees)
if starting_trees:
    ps.add( starting_trees, "--tree {}" )

if "prefix" in snakemake.params.keys():
    prefix = snakemake.params.prefix
else:
    prefix = os.path.join( outdir, "search" )
ps.add( prefix,  "--prefix {}" )

ps.add_opt( "data-type" )
ps.add_opt( "log" )
ps.add_opt( "redo",         sp.typ.FLAG )
ps.add_opt( "precision",    sp.typ.UINT )
ps.add_opt( "site-weights", sp.typ.FILE )
ps.add_opt( "outgroup" )

# General options
def ONOFF():
    return sp.typ.IN( ['on', 'off'] )
ps.add_opt( "seed" )
ps.add_opt( "pat-comp",     ONOFF )
ps.add_opt( "tip-inner",    ONOFF )
ps.add_opt( "site-repeats", ONOFF )
ps.add_threads()
ps.add_opt( "workers",      sp.typ.UINT )
ps.add_opt( "simd",         sp.typ.IN(['none','sse3','avx','avx2']) )
ps.add_opt( "rate-scalers", ONOFF )
ps.add_opt( "force" )

# Model options
ps.add_opt( "brlen",        sp.typ.IN(['linked','scaled','unlinked']) )
ps.add_opt( "blmin",        sp.typ.FLOAT )
ps.add_opt( "blmax",        sp.typ.FLOAT )
ps.add_opt( "blopt",        sp.typ.IN(
    ['nr_fast','nr_safe','nr_oldfast','nr_oldsafe']) )
ps.add_opt( "opt-model",    ONOFF )
ps.add_opt( "opt-branches", ONOFF )
ps.add_opt( "prob-msa",     ONOFF )
ps.add_opt( "lh-epsilon",   sp.typ.FLOAT )

# Topology search options
ps.add_opt( "spr-radius",   sp.typ.UINT )
ps.add_opt( "spr-cutoff",   lambda a: sp.typ.FLOAT(a) if not a == "off" else True )

# Bootstrapping options
ps.add_opt( "bs-trees" )
ps.add_opt( "bs-cutoff" )
ps.add_opt( "bs-metric",    sp.typ.IN(['fbp','tbe']) )
ps.add_opt( "bs-write-msa", ONOFF )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

# =================================================================================================
#     Move over results
# =================================================================================================
def check_and_move( src, dest ):
    result_file = os.path.join( outdir, prefix + ".raxml.rfDistances" )
    util.expect_file_exists( src )
    # rename the result file
    shell( f"mv {src} {dest}" )

check_and_move( prefix + ".raxml.bestTree",     snakemake.output.best_tree )
check_and_move( prefix + ".raxml.bestModel",    snakemake.output.best_model )
check_and_move( prefix + ".raxml.support",      snakemake.output.support_tree )
check_and_move( prefix + ".raxml.mlTrees",      snakemake.output.ml_trees )
check_and_move( prefix + ".raxml.bootstraps",   snakemake.output.bs_trees )
