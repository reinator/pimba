from os import access
import util
from collections.abc import Callable
from sys import maxsize
from operator import getitem
from functools import reduce

class typ:
    @staticmethod
    def NONE( arg ):
        return True
    @staticmethod
    def FILE( arg ):
        util.expect_file_exists( arg )
    @staticmethod
    def FILES( arg ):
        if type(arg) is not list: arg = list(arg)
        for f in arg:
            util.expect_file_exists( f )
    @staticmethod
    def EXEC( arg ):
        util.expect_executable_exists( arg )
    @staticmethod
    def DIR( arg ):
        util.expect_dir_exists( arg )
    @staticmethod
    def FLAG( arg: str ):
        if not isinstance( arg, bool ):
            util.fail( f"expected flag (True or False), but got '{arg}'")
    @staticmethod
    def IN( set: list ):
        def func( arg ):
            if not arg in set:
                util.fail( f"{arg} not in {set}" )
        return func
    @staticmethod
    def STRING( arg ):
        return True
    @staticmethod
    def FLOAT( lower: float = 0.0, upper: float = float("inf") ):
        """
        Returns validation function that checks if arg is within [lower,upper]
        """
        def func( arg ):
            f = float(arg)
            if not ( (f >= lower) and (f <= upper) ):
                util.fail( f"{arg} not within [{lower},{upper}]" )
        return func
    @staticmethod
    def UINT( lower: int = 0, upper: int = maxsize*2+1 ):
        """
        Returns validation function that checks if arg is within [lower,upper]
        """
        def func( arg ):
            f = int(arg)
            if not ( (f >= lower) and (f <= upper) ):
                util.fail( f"{arg} not within [{lower},{upper}]" )
        return func


def rule_key( s ):
    """Helper function that converts any '-' characters in a string to '-'.
    Needed to enable differing formatting between rules/params and config/params keys
    """
    return s.replace('-', '_')


class Parser:
    """Class for progressively building a validated shell command.
    Basically a fancy StringBuilder class with some input validation.
    """
    _shell_string: str
    _snakemake: object
    # list of accessors (lists of progressively deeper entries into config). In descending order
    # of preference (i.e., first try snakemake["config"]["params"][tool][key], then more if specified)
    # note that the params set in the calling rule always takes precedence
    _accessors: list[list] = []
    _arg_ident: str

    # custom constructor, to clarify inputs and to give direct control over accessor behaviour
    def __init__(   self,
                    tool: str,
                    snakemake: object,
                    accessors: list[list] = None,
                    arg_ident = '--' ):
        self._shell_string = tool
        self._snakemake = snakemake
        if not accessors and tool in snakemake.config['params'].keys():
            self._accessors = [['params',tool]] + self._accessors
        if accessors:
            if type(accessors) is not list[list]: accessors = [ accessors ]
            self._accessors = accessors + self._accessors
        self._arg_ident = arg_ident

    def add( self, arg,
            fstr: str = "{}",
            valid_func: Callable[[str], bool] = typ.NONE ):
        
        # Validate input
        valid_func( arg )

        # add to the shell string
        if valid_func is typ.FLAG:
            # differentiate between flags (they don't have arg values attached)
            assert( not util.has_format_fields( fstr ) )
            if arg:
                fstr = " " + fstr
            else:
                fstr = ""
        else:
            # ...and normal CLI arguments
            assert( util.has_format_fields( fstr ) )
            # surround by quotes if its a file, as the path may have spaces
            arg = f'\"{arg}\"' if (valid_func is typ.FILE) or (valid_func is typ.STRING)  else arg
            # same but for a list of files
            if valid_func is typ.FILES:
                if type(arg) is not list: arg = list(arg)
                arg = " ".join([f'\"{a}\"' for a in arg])
            fstr = " " + fstr.format( arg )

        # add to the complete shell string
        self._shell_string = self._shell_string + fstr

        # return the format string incase we want to check it
        return fstr
    
    def add_opt( self, key: str,
                valid_func: Callable[[str], bool] = typ.NONE,
                fstr: str = "",
                add_accessors=[] ):
        # we allow custom accessors to be added to the config accessor list
        if type(add_accessors) is not list[list]: add_accessors = [ add_accessors ]
        # if no format string was specified, build it from the key (--key)
        if not fstr:
            fstr = f"{self._arg_ident}{key}"
            if not valid_func is typ.FLAG:
                fstr = fstr + r" {}"

        # give absolute priority to the rule params, as they may have critical settings
        if rule_key( key ) in self._snakemake.params.keys():
            return self.add( self._snakemake.params[rule_key(key)], fstr, valid_func )

        # if the key isn't set in rule/params, look through the config accessors
        # and return the first hit
        for accessor in self._accessors + add_accessors:
            entry = reduce( getitem, accessor, self._snakemake.config )
            if key in entry.keys():
                return self.add( entry[key], fstr, valid_func )

    def add_threads( self, format_string: str = "--threads {}", valid_func: Callable[[str], bool] = typ.UINT ):
        return self.add( self._snakemake.threads, format_string, valid_func )

    def get_shell_string( self, do_log: bool = True, log_stdout: bool = True ):
        if do_log:
            return self._shell_string + f" {self._snakemake.log_fmt_shell(stdout=log_stdout, stderr=True)}"
        else:
            return self._shell_string 
