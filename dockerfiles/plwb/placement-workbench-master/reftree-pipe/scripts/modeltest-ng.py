# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import snakeparser as sp

# Taken from: https://github.com/BenoitMorel/ParGenes/blob/master/pargenes/pargenes_src/modeltest.py
def get_model_from_log(log_file, modeltest_criteria):
  with open(log_file) as reader:
    read_next_model = False
    for line in reader.readlines():
      if (line.startswith("Best model according to " + modeltest_criteria)):
          read_next_model = True
      if (read_next_model and line.startswith("Model")):
        model = line.split(" ")[-1][:-1]
        return model
  return None

def get_sites_from_log(log_file):
  with open(log_file) as reader:
    read_next_model = False
    for line in reader.readlines():
      if (line.startswith("  #sites:")):
        num_sites = line.split(" ")[-1][:-1]
        return num_sites
  return None

shell.executable("bash")

# Get the output directory
sample_outdir = util.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "modeltest-ng", snakemake )

# Required args
ps.add( snakemake.input[0], "--input {}", sp.typ.FILE )

# Optional args
ps.add_opt( "datatype",           sp.typ.IN(['nt','aa']) )
ps.add_opt( "partitions",         sp.typ.FILE )
ps.add_opt( "rngseed" )
ps.add_opt( "topology",
  sp.typ.IN(['ml', 'mp', 'fixed-ml-jc', 'fixed-ml-gtr', 'fixed-mp', 'random', 'user']) )
ps.add_opt( "utree",              sp.typ.FILE )
ps.add_opt( "force",              sp.typ.FLAG )
ps.add_opt( "disable-checkpoint", sp.typ.FLAG )
ps.add_opt( "asc-bias",           sp.typ.IN(['lewis','felsenstein','stamatakis']) )
ps.add_opt( "frequencies",        sp.typ.IN(['e','f']) )
ps.add_opt( "models-het",         sp.typ.IN(['u','i','g','f']) )
ps.add_opt( "models" )
ps.add_opt( "schemes",            sp.typ.IN(['3','5','7','11','203']) )
ps.add_opt( "template",           sp.typ.IN(['raxml','phyml','mrbayes','paup']) )
ps.add_opt( "eps",                sp.typ.FLOAT )
ps.add_opt( "tol",                sp.typ.FLOAT )
ps.add_opt( "smooth-frequencies", sp.typ.FLAG )
ps.add_opt( "gamma-rates",        sp.typ.IN(['a','g']) )
ps.add_opt( "no-compress",        sp.typ.FLAG )
ps.add_opt( "keep-params",        sp.typ.FLAG )
ps.add_opt( "verbose",            sp.typ.FLAG )

# Closing args
ps.add( snakemake.output[0], "--output {}" )
ps.add_threads( "--processes {}" )


# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

result_file = os.path.join( sample_outdir, "model.file.out" )
util.expect_file_exists( result_file )

model = get_model_from_log( result_file, "AIC" )
num_sites = get_sites_from_log( result_file )

if (not num_sites) or (not model):
    util.fail( "Modeltest output could not be parsed correctly." )

with open( snakemake.output[0], 'w+' ) as f:
    f.write( model + ", noname = 1-{}".format(num_sites) + '\n' )
