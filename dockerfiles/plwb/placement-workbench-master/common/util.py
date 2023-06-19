from shutil import rmtree, copyfile, copytree
import os
from os.path import join, isfile
import sys
import stat
import subprocess
from pathlib import Path
import re

# =================================================================================================
#     Error Handling
# =================================================================================================

def fail( msg ):
  raise Exception( "ERROR: " + msg )

def warn( msg ):
  print( "WARNING: " + msg )

# =================================================================================================
#     Input Validation
# =================================================================================================

def expect_dir_exists( dir_path ):
  if not os.path.isdir( dir_path ):
    fail( "Directory doesn't exist: " + dir_path )

def expect_file_exists( file_path ):
  if not os.path.isfile( file_path ):
    fail( "File doesn't exist: " + file_path )

def expect_executable_exists( executable ):
  import distutils.spawn
  if not distutils.spawn.find_executable( executable ):
    fail( "Executable not found: " + executable )

def parse_file_path( file_path ):
  file_path = os.path.normpath( file_path )
  expect_file_exists( file_path )
  return file_path

def parse_dir_path( dir_path ):
  dir_path = os.path.normpath( dir_path )
  expect_dir_exists( dir_path )
  return dir_path

def parse_executable_path( executable_path ):
  executable_path = os.path.normpath( executable_path )
  expect_executable_exists( executable_path )
  return executable_path

# =================================================================================================
#     File System Helpers
# =================================================================================================

def filename( path ):
  return os.path.splitext( os.path.basename( path ) )[0]

def dirname( path ):
  return os.path.dirname( path )

def extension( path ):
  parts = os.path.splitext( path )
  if( len(parts) == 1 ):
    fail( "file '{}' does not appear to have an extension, which is required.".format(path) )
  return parts[1]

def splitpath( path, maxdepth=20 ):
  path = os.path.normpath(path)
  ( head, tail ) = os.path.split(path)
  return splitpath(head, maxdepth - 1) + [ tail ] if maxdepth and head and head != path else [ head or tail ]

def num_dirs( path ):
  """Returns the number of directories in a path"""
  return len( Path( os.path.dirname( path ) ).parts )

def trim_path( dir_path, path_exts ):
  """
  Takes a directory path and removes a path from the end of it.
  Can take a list of possible matching path extensions, in whcih case the first matching
  extension is trimmed and the remaining path is returned
  """
  if type(path_exts) is not list: path_exts = [ path_exts ]
  expect_dir_exists( dir_path )
  path = os.path.normpath( dir_path )
  path_parts = Path( path ).parts
  for ext in path_exts:
    to_split = Path( ext ).parts
    n = len(to_split)
    if path_parts[-n:] == to_split:
      return Path(*(path_parts[:-n or None]))
  fail(f"dirpath '{dir_path}' did not match any given extensions '{path_exts}'")

def last_n_dirnames( path, n ):
  path = os.path.normpath(path)
  names = Path( os.path.dirname( path ) ).parts[ 1: ]
  # don't try to return more than there are
  n = min( n, len(names) )
  # return the last n parts
  return list( names[ -n: ] or None ) if n else []

def rstrip_gz( path ):
    """Strips .gz from a file ending, if it is there"""
    suffix = '.gz'
    if path.endswith( suffix ):
        return path[:-len(suffix)]
    else:
        return path

def ingest_paths( paths, extensions=None, allow_gz=False ):
  """Takes a list of paths, validates them to make sure they exist, and if a path is a directory
  globs all files in said directory that have the specified extension. If no specific extension
  is provided, returns all files in that directory (non-recursively)
  """
  file_list = []
  for path in paths:
    if os.path.isfile( path ):
      expect_file_exists( path )
      chk_path = rstrip_gz(path) if allow_gz else path
      if extensions and not extension(chk_path) in extensions:
        warn("file '{path}' does not have any of the expected file extensions: {extensions}")
      if not path in file_list:
        file_list.append( path )

    elif os.path.isdir( path ):
      expect_dir_exists( path )
      files = [join(path, f) for f in os.listdir( path ) if isfile( join(path, f) )]
      if extensions:
        file_list.extend( [f for f in files if extension( rstrip_gz(f) ) in extensions and not f in file_list] )
      else:
        file_list.extend( [f for f in files if not f in file_list] )
  # ensure duplicates are discarded
  file_list = list(set(file_list))
  return file_list

def get_unique_names( paths, allow_gz=False ):
  """
  We want to extract unique sample names from the file paths of the input files.
  To do so, we attempt to keep prepending directory names to the proposed sample
  names, until we have a unique set of names (the idea being that the user has encoded 
  the information about what the samples are called in their directory structure).
  Note that this can still fail, for example if the input files share name and directory,
  but not file extension ('x.fa x.fasta' for example).
  """

  # as we allow full duplicates, we need to track them and their indices
  # in 'paths', which we do here using a dict
  from collections import defaultdict
  paths_dict = defaultdict(list)
  for i, p in enumerate(paths):
    paths_dict[os.path.normpath(p)].append(i)
  
  # next we go through all paths, with an increasing number of path parts to go back
  names = []
  # extend by at most as long as the longest path
  for i in range(0, max( [num_dirs( f ) for f in paths_dict.keys()] )):

    failed=False
    for path in paths_dict.keys():
      
      # get at most i preceding dir names as prefixes
      prefixes = last_n_dirnames( path, i )
      # optional gzip ignore
      path = rstrip_gz( path ) if allow_gz else path
      prefixes.append( filename( path ) )
      new_name = '_'.join( prefixes )
      # ensure the name doesn't contain unescaped spaces
      new_name = re.sub( r'(?<=[^\\])\s', '_', new_name )
      if( new_name in names ):
        failed = True
        names = []
        break
      else:
        names.append( new_name )

    if not failed:
      break

  if failed:
    fail( f"Could not find assignment of unique names to list of input files. The list:\n{paths}" )

  assert( len(names) == len(paths_dict) )

  # finally, translate back from the paths dict to a name list that reflects the input paths
  names_of_paths = [''] * len(paths)
  for names_idx, paths_idxs in enumerate(paths_dict.values()):
    for i in paths_idxs:
      names_of_paths[ i ] = names[ names_idx ]

  return names_of_paths

def is_fastq( path ):
  return (extension( rstrip_gz(path) ) in [".fastq", ".fq"])

def is_fasta( path ):
  return (extension( rstrip_gz(path) ) in [".fasta", ".fa", ".afa"])

# =================================================================================================
#     File System Manipulation
# =================================================================================================

def copy( src, dest ):
  copyfile( src, dest )

def copy_dir( src, dest, ignore=None ):
  if ignore:
    ign_f = shutil.ignore_patterns(*ignore)
  else:
    ign_f = None
  copytree( src, dest, ignore=ign_f )

def clean_dir( path ):
  if os.path.exists( path ):
    rmtree( path, ignore_errors=True )

def clean_file ( path ):
  if os.path.exists( path ):
    os.remove( path )

def chmod_path( path ):
  if not os.path.isdir( path ):
    raise RuntimeError( "Directory doesn't exist: " + path )
  else:
    os.chmod( path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO )

def chmod_file( file_path ):
  if not os.path.isfile( file_path ):
    raise RuntimeError( "File doesn't exist: " + file_path )
  else:
    os.chmod( file_path, stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IROTH|stat.S_IWOTH )

def mkdirp( path ):
  if not os.path.exists( path ):
    os.mkdir( path )
    chmod_path( path )

def make_path( path ):
  if not os.path.exists( path ):
    os.makedirs( path )
    chmod_path( path )

def make_path_clean( path ):
  clean_dir( path )
  make_path( path )

# =================================================================================================
#     String Functions
# =================================================================================================

def find_string_between(input_str, marker1, marker2):
  """returns the first occurence of the string matching
    .*${marker1}${input_str}${marker2}.*  (bash notation)

    the function does not check that such a string exists
  """
  start = input_str.find(marker1) + len(marker1)
  end = input_str.find(marker2, start)
  return input_str[start:end]

def has_format_fields( s ):
  """Checks if a string has any formattable fields ('{}')"""
  from string import Formatter
  for _,a,b,c in Formatter().parse( s ):
    if a != None or b != None or c != None:
      return True
  return False

# =================================================================================================
#     System Functions
# =================================================================================================

# doesnt work on mac!
def num_physical_cores():
  out = subprocess.check_output(['lscpu', '--parse=Core,Socket'], encoding = 'utf-8')
  out = out.split('\n')
  num_cores = len(
     dict.fromkeys(
       [line for line in out if not line.startswith('#') and line != '']
  ))
  return num_cores

def is_tool(name):
  # via https://stackoverflow.com/a/34177358
  """Check whether `name` is on PATH and marked as executable."""

  # from whichcraft import which
  from shutil import which

  return which(name) is not None

# =================================================================================================
#     Snakemake-specific helper functions
# =================================================================================================
def config_to_file( config, out_dir ):
  import yaml
  with open( os.path.join( out_dir, "config.yaml" ), 'w' ) as outfile:
      yaml.dump( dict(config), outfile )

def listlen( obj ):
  """Returns the number of items in object, if obj is list, 1 otherwise"""
  if isinstance( obj, list ):
    return len(obj)
  else:
    return 1

def cluster_settings( clust_env, calling_dir ):
  # detect/set the job submission system to be used
  cluster_sys = None
  if clust_env == 'auto':
    if is_tool("sbatch"):
      cluster_sys = "slurm"
    elif is_tool("qsub"):
      cluster_sys = "sge"
    else:
      fail( "Could not autodetect job submission system." )
  else:
    cluster_sys = clust_env
  # set the correct cluster config and other needed settings
  # (that would otherwise be handled by --profile)
  profile = join( calling_dir, "profiles", cluster_sys )

  # jobscript = join( profile, f"{cluster_sys}-jobscript.py" )
  cluster_config = [join( profile, "config.yaml" ),
                    join( calling_dir, "profiles/cluster-config.yaml" )]
  cluster = join( profile, "slurm-submit.py" )

  return cluster, cluster_config
