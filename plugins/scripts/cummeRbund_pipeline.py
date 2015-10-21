#*****************************************************************************
#  cummerbund_pipeline.py
#******************************************************************************

"""
####################
cummerbund.py
####################
Code defining functions to encapsulate cummeRbund functionalities.
"""

'''
Notes on what basic functions I want to inlude and what the default output should be:

readCufflinks


'''
import os
import sys
import argparse

try:
    from rpy2.robjects import r
    from rpy2.rinterface import RRuntimeError
except ImportError as ie:
    raise errors.BlacktieError('Unable to import required module: "rpy2".  Try installing it with "[sudo] pip install rpy2".')
except RuntimeError as rte:
    raise errors.BlacktieError('Importing required module "rpy2" failed because no R application could be found on your system. Try again after instaling R.')

def main():
    """
    The main loop.  Lets ROCK!
    """

    desc = """This script reads the files in a cuffdiff output directory into cummeRbund, generates some standard preliminary plots, and saves the output."""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--version', action='version', version='%(prog)s ' + blacktie.__version__,
                        help="""Print version number.""")    
    parser.add_argument('--cuffdiff-dir', type=str,
                        help="""Path to a cuffdiff output directory.""")
    #parser.add_argument('--cummerbund-db', type=str,
                        #help="""Path to a pre-built cummeRbund 'cuffData.db'. (this is rarely specified directly; usually --cuffdiff-dir works fine)""")
    parser.add_argument('--gtf-path', type=str, default='NULL',
                        help="""Path to gtf file used in cuffdiff analysis. This will provide transcript model information.""")
    parser.add_argument('--genome', type=str, default='NULL',
                        help="""String indicating which genome build the .gtf annotations are for (e.g. 'hg19' or 'mm9').""")
    parser.add_argument('--out', type=str, 
                        help="""A base directory to add to our saved plots into.""")
    parser.add_argument('--file-type', type=str, choices=['pdf','jpeg','png','ps'], default='pdf',
                        help="""The type of output file to use when saving our plots. (default: %(default)s)""")

