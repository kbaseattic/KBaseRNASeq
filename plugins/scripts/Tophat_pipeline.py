#!/usr/bin/python

"""
####################
Tophat_pipeline.py

####################
Code defining an object oriented python pipeline script to allow simplified
coordination of data through parts or all of the popular Tophat/Cufflinks
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
from test_call_utils import *
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

    parser.add_argument('-input', type=str,required=True,nargs="+",help="Input fastq files ")
    parser.add_argument('-output', type=str,required=True,help = "Output directory name")
    parser.add_argument('-reference',type=str,required=True,help="Reference file for alignment")
    parser.add_argument('-opts_dict', type=json.loads,help="""Options passed from the parameter file""")
    parser.add_argument('-gtf',type=str,help="gtf anotation for the reference in GTF format")

    parser.add_argument('-prog', type=str, nargs='+', choices=['tophat', 'cufflinks', 'cuffmerge', 'cuffdiff',
                                                                'cummerbund', 'all'], default='tophat',
                        help="""Which program do you want to run? (default: %(default)s)""")
    parser.add_argument('-hide-logs', action='store_true', default=False,
                        help="""Make your log directories hidden to keep a tidy 'looking' base directory. (default:
                        %(default)s)""")
    parser.add_argument('-base_dir',type=str, help= 'Option to set the base dir',required=True)
    parser.add_argument('-library_type', type=str, choices=['SingleEnd','PairedEnd'],help = 'Library type for the input Fastq files',required=True)
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

############
############
#Things to do
############
############
# Check for Paired end data
# remove hard coded paths 
# identify a /tmp dir
# How to create bowtie indexes in a directory
# Fix bowtie index base name
# Fix required options

    #yargs = bunchify(yaml.load(open(args.config_file, 'rU')))
    #print 'yargs is ' + json.dumps(yargs)
    # set up run_id, log files, and email info
    #else:
    yargs = dict()
    print args.output
    yargs['run_options']=dict()	
    yargs['run_options']['run_id'] = get_time()
    yargs['run_options']['base_dir'] = args.base_dir # set it to a default directory for the service /tmp
    yargs['run_options']['bowtie_indexes_dir'] = "/home/srividya/KBaseRNASeq/t/data" # set it to a default directory for the service /tmp 
    yargs['tophat_options'] =args.opts_dict
    yargs['tophat_options']['o']=args.output
    yargs['tophat_options']['positional_args']= dict() 
    yargs['tophat_options']['positional_args']['bowtie2_index']= "from_conditions"
    yargs['condition'] = dict()
    yargs['condition']['name'] = 'name1'   # get the workspace object metadata name
    yargs['condition']['experiment_id'] =  'experiment1'# get the experiment id from the workspace object
    yargs['condition']['replicate_id'] =  1# get the replicate id from the workspace object
    yargs['condition']['genome_seq'] = args.reference    
    if hasattr(args , 'gtf'):
        yargs['tophat_options']['G'] = "from_conditions"
	yargs['condition']['gtf_annotation'] = args.gtf
    if args.library_type == "SingleEnd":
        yargs['tophat_options']['positional_args']['singleend_reads'] = "from_conditions"
    	yargs['condition']['singleend_reads'] = args.input[0]
    else:
	yargs['tophat_options']['positional_args']['left_reads'] = "from_conditions"
	yargs['tophat_options']['positional_args']['right_reads'] = "from_conditions"
	yargs['condition']['left_reads'] = args.input[0]
	yargs['condition']['right_reads'] = args.input[1]
    yargs['condition']['bowtie2_index'] = 'Ecoli_K12_MG1655_NC_000913.3' #### handle the bowtie index dir 
      
    yargs= bunchify(yaml.load(yaml.dump(yargs)))  
    print yargs  
    if yargs.run_options.run_id:
        run_id = yargs.run_options.run_id
		
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
    #else:
    #    pass

    yargs.prgbar_regex = re.compile('>.+Processing.+\[.+\].+%\w*$')
   # print 'length is ' + str(len(yargs.condition_queue))
   # if len(yargs.condition_queue) > 1:
   #	 print 'inside the if conditions'
   # 	yargs.groups = map_condition_groups(yargs)
   # else:
   # 	 print 'inside the else contiions'
   #	yargs.groups = yargs.condition_queue[0]
   # yargs.call_records = {}

    #print ' yargs groups is ' + json.dumps(yargs.groups)
    # loop through the queued conditions and send reports for tophat 
    if start_prog(args.prog, ['tophat', 'all']):
        print '[Note] Starting tophat step.\n'
        #for condition in yargs.condition_queue:
	    #print 'condition is ' + json.dumps(condition)
            # Prep Tophat Call
	try:
        	tophat_call = TophatCall(yargs, run_id, run_logs,conditions=yargs.condition, mode=args.mode)
      	    	#print 'tophat call is ' + json.dumps(tophat_call)
		tophat_call.execute()
	except:
		raise

	
            # record the tophat_call object
            #yargs.call_records[tophat_call.call_id] = tophat_call
    else:
        print "[Note] Skipping tophat step.\n"

if __name__ == "__main__":
    main()
