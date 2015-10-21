#!/usr/bin/python
"""
####################
Cufflinks_pipeline.py
####################
Code defining an object oriented python pipeline script to allow simplified
coordination of data through parts or all of the popular Cufflinks
RNA-seq analysis suite.
"""


import sys
import argparse
import re

import yaml
import json
try:
    import pprocess
except ImportError:
    pprocess = False
    pass

from bunch import bunchify
from externals import *
from errors import *
from call_utils import *
from misc import *

def start_prog(args_prog, test_vals):
    """
    Manage testing whether a program block should be entered based on the contents of args.prog.

    :param args_prog: `list` of program names.
    :param test_vals: `list` of strings that should trigger returning `True`
    :return: `bool`
    """
    return any((prog in test_vals) for prog in args_prog)

def main():
    """
    The main loop.  Lets ROCK!
    """

    desc = """This script reads options from a yaml formatted file and organizes the execution of tophat/cufflinks
    runs for multiple condition sets."""

    parser = argparse.ArgumentParser(description=desc)

    #parser.add_argument('--version', action='version', version='%(prog)s ' + blacktie.__version__,

    parser.add_argument('-input', type=str,required=True,nargs='+',help="Input GTF files ")
    parser.add_argument('-output', type=str,required=True,help = "Output directory name")
    parser.add_argument('-opts_dict', type=json.loads,help="""Options passed from the parameter file""")
    parser.add_argument('-gtf',type=str,help="gtf anotation for the reference in GTF format")
    parser.add_argument('-num-replicates',type=int,help="Number of replicates found in the data, 1 if no replicates are present")
   
    
    parser.add_argument('-labels',type=str,nargs='+',help="Lables for the different conditions in the data")
    parser.add_argument('-prog', type=str, nargs='+', choices=['tophat', 'cufflinks', 'cuffmerge', 'cuffdiff',
                                                                'cummerbund', 'all'], default='tophat',
                        help="""Which program do you want to run? (default: %(default)s)""")
    parser.add_argument('-hide-logs', action='store_true', default=False,
                        help="""Make your log directories hidden to keep a tidy 'looking' base directory. (default:
                        %(default)s)""")
    parser.add_argument('-base_dir',type=str, help= 'Option to set the base dir',required=True)
    parser.add_argument('-mode', type=str, choices=['analyze', 'dry_run'], default='analyze',
                        help="""1) 'analyze': run the analysis pipeline. 2) 'dry_run': walk through all steps that
                        would be run and print out the command lines; however, do not send the commands to the
                        system to be run. 3) 'qsub_script': generate bash scripts suitable to be sent to a compute
                        cluster's
                        SGE through the qsub command. (default: %(default)s)""")

    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    args = parser.parse_args()
    print args.output
    yargs = dict()
    yargs['run_options']=dict()
    yargs['run_options']['run_id'] = get_time()
    yargs['run_options']['base_dir'] = args.base_dir # set it to a default directory for the service /tmp
    yargs['cuffdiff_options'] = args.opts_dict
    yargs['cuffdiff_options']['o'] = args.output
    yargs['cuffdiff_options']['labels'] = args.labels
    yargs['cuffdiff_options']['positional_args']= dict()
    yargs['cuffdiff_options']['positional_args']['sample_bams'] = args.input
    yargs['cuffdiff_options']['positional_args']['transcripts_gtf'] = args.input

    yargs['condition'] = dict()
    yargs['condition']['name'] = 'name1'   # get the workspace object metadata name
    yargs['condition']['experiment_id'] =  'experiment1'# get the experiment id from the workspace object
    yargs['condition']['replicate_id'] =  1# get the replicate id from the workspace object
    if hasattr(args , 'gtf'):
        yargs['cuffdiff_options']['G'] = "from_conditions"
        yargs['condition']['gtf_annotation'] = args.gtf

    yargs= bunchify(yaml.load(yaml.dump(yargs)))
    # set up run_id, log files, and email info
    if yargs.run_options.run_id:
        run_id = yargs.run_options.run_id
    else:

        run_id = get_time()

    base_dir = yargs.run_options.base_dir.rstrip('/')
    if args.hide_logs:
        run_logs = '%s/.%s.logs' % (base_dir, run_id)
    else:
        run_logs = '%s/%s.logs' % (base_dir, run_id)

    if not args.mode == 'dry_run':
        mkdirp(run_logs)
    else:
        pass

    yaml_out = '%s/%s.yaml' % (run_logs, run_id)

    # copy yaml config file with run_id as name for records
    #if not args.mode == 'dry_run':
        #shutil.copyfile(args.config_file, yaml_out)
    # else:
    #    pass

    yargs.prgbar_regex = re.compile('>.+Processing.+\[.+\].+%\w*$')
    #print 'length is ' + str(len(yargs.condition_queue))
    #if len(yargs.condition_queue) > 1:
    # 	print 'inside the if conditions'
    # 	yargs.groups = map_condition_groups(yargs)
    #else:
    #	print 'inside the else contiions'
    #	yargs.groups = yargs.condition_queue[0]
    #yargs.call_records = {}

    #print ' yargs groups is ' + json.dumps(yargs.groups)
    # loop through the queued conditions and send reports for tophat

    if start_prog(args.prog, ['cuffdiff', 'all']):
        # attempt to run more than one cufflinks call in parallel since cufflinks
        # seems to use only one processor no matter the value of -p you give it and
        # doesn't seem to consume massive amounts of memory 
        print "[Note] Starting cuffdiff step.\n"
        
    	print yargs.condition
    	print yargs.cuffdiff_options
	cuffdiff_call = CuffdiffCall(yargs, run_id, run_logs, conditions=yargs.condition, mode=args.mode)
        cuffdiff_call.execute()

                # record the cufflinks_call object
        #yargs.call_records[cufflinks_call.call_id] = cufflinks_call
    else:
        print "[Note] Skipping cufflinks step.\n"

if __name__ == "__main__":
    main()
