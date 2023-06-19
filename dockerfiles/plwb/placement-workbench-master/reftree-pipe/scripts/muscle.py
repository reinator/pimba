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
ps = sp.Parser( "muscle", snakemake, arg_ident='-' )

# as per https://drive5.com/muscle5/manual/cmd_align.html
# OPTIONS
ps.add( snakemake.input[0],     "-align {}",    sp.typ.FILE )
ps.add( snakemake.output[0],    "-output {}" )

ps.add_opt( "perturb",      sp.typ.UINT )
ps.add_opt( "perm",         sp.typ.IN(['none','abc','acb','bca']) )
ps.add_opt( "stratified",   sp.typ.FLAG )
ps.add_opt( "diversified",  sp.typ.FLAG )
ps.add_opt( "replicates",   sp.typ.UINT )
ps.add_opt( "consiters",    sp.typ.UINT )
ps.add_opt( "refineiters",  sp.typ.UINT )

ps.add_opt( "nt",           sp.typ.FLAG )
ps.add_opt( "amino",        sp.typ.FLAG )

ps.add_threads( "-threads {}" )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
