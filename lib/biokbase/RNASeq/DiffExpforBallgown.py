import simplejson, sys, shutil, os, ast , re
from mpipe import OrderedStage , Pipeline
import glob, json, uuid, logging  , time ,datetime 
import subprocess, threading,traceback
from collections import OrderedDict
from pprint import pprint, pformat
import parallel_tools as parallel
from mpipe import OrderedStage , Pipeline
import contig_id_mapping as c_mapping 
import script_util
import handler_utils as handler_util
from biokbase.workspace.client import Workspace
from biokbase.auth import Token
import multiprocessing as mp
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from biokbase.RNASeq.ExecutionBase import ExecutionBase
from biokbase.RNASeq import rnaseq_util
#import ExecutionBase.ExecutionBase as ExecutionBase

class DiffExpforBallgownException(Exception):
    pass

class DiffExpforBallgown(ExecutionBase): 

    def __init__(self, logger, directory, urls):
        logger.info( "in DiffExprforBallgown, type logger is " + pformat( type( logger ) ) )
        logger.info( " urls are " + pformat( urls ) )
        pprint(self.__class__)
        super(self.__class__, self).__init__(logger, directory, urls)

        # user defined shared variables across methods
        #self.num_threads = None
        self.num_threads = 1
        self.num_cores = 1
        self.tool_used = None
        self.tool_version = None

    def prepare(self):
        # for quick testing, we recover parameters here

        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        token = self.common_params['user_token']
        diffexp_dir = self.directory
        logger = self.logger
        logger.info( 'in DiffExpfoBallgown.prepare(), method params are')
        logger.info( pformat( self.method_params ) )

        #self.details = rnaseq_util.get_details_for_diff_exp(logger,ws_client,hs,params['ws_id'],self.urls,diffexp_dir,params['expressionset_id'],token)
        #logger.info( 'back from get_details_for_diff_exp(), details are')
        #logger.info( pformat( self.details ) )
        self.num_threads = mp.cpu_count()
        self.num_jobs = 1

        self.details = {}
        self.details["used_tool"] = "StringTie"    # Question: where does this really get set?
        #als = [] 
        #for l in self.details['labels']:
        #        rep_files=[ (os.path.join(diffexp_dir+'/'+l,sub+'/accepted_hits.bam'), os.path.join(diffexp_dir+'/'+l,sub+'/transcripts.gtf')) for sub in os.listdir(os.path.join(diffexp_dir,l)) if os.path.isdir(os.path.join(diffexp_dir,l+'/'+sub))]
        #        #rep_files=",".join([ os.path.join(diffexp_dir+'/'+l,sub+'/accepted_hits.bam') for sub in os.listdir(os.path.join(diffexp_dir,l)) if os.path.isdir(os.path.join(diffexp_dir,l+'/'+sub))])
        #        als += rep_files
        #### Call Cuffmerge function
        #used_tool = self.details['used_tool']
        #merge_dir = os.path.join(diffexp_dir,"merge")
        #if used_tool == 'StringTie':
        #   run_tool =  "StringTie"
        #   tool_version = "1.2.3"
        #   #merged_gtf = rnaseq_util.call_stringtiemerge(diffexp_dir,merge_dir,self.num_threads,self.details['gtf_file'],self.details['gtf_list_file'])
        #elif used_tool == 'Cufflinks':
        #   merged_gtf = rnaseq_util.call_cuffmerge(diffexp_dir,merge_dir,num_threads,gtf_file,self.details['gtf_list_file'])
        #   run_tool = "Tablemaker"
        #   tool_version = '2.0.9'
        #   merged_gtf = rnaseq_util.call_cuffmerge(diffexp_dir,merge_dir,self.num_threads,self.details['gtf_file'],self.details['gtf_list_file'])
#
        #self.bam_files = " ".join([i for i in als])
        #self.t_labels = ",".join(self.details['labels'])
        #ballgown_dir = os.path.join(diffexp_dir,"ballgown")
        #if not os.path.exists(ballgown_dir): os.mkdir(ballgown_dir)
        #### Make Input_dir from expression_file_name
        
        self.task_list = [self.__class__]
        logger.info( 'exiting ')


    def runEach(self,task_list):
         logger = self.logger
         ### Call Cuffmerge function
         used_tool = self.details['used_tool']
         logger.info(  'in DiffExpfoBallgown.runEach()' )
         if used_tool == 'StringTie':
           #merged_gtf = rnaseq_util.call_stringtiemerge(diffexp_dir,merge_dir,num_threads,self.details['gtf_file'],assembly_file)
           #run_tool =  "StringTie"
           #tool_version = "1.2.3"
           # For now, take no action for StringTie processing
           logger.info( 'Exiting immediately - StringTie case' )
           return
         elif used_tool == 'Cufflinks':
           merged_gtf = rnaseq_util.call_cuffmerge(diffexp_dir,merge_dir,num_threads,gtf_file,assembly_file)
           run_tool = "Tablemaker" 
           tool_version = '2.0.9'
         cuffmerge_dir = os.path.join(self.directory,"cuffmerge")
         merged_gtf = rnaseq_util.call_cuffmerge(self.directory,cuffmerge_dir,self.num_threads,self.details['gtf_file'],self.details['gtf_list_file'])
         ### Run DiffExpforBallgown
         output_dir = os.path.join(self.directory,self.method_params['output_obj_name'])
         diffexp_command = (' -p '+str(self.num_threads))

         ### Setting Advanced parameters for DiffExpforBallgown

         if('time_series' in self.method_params and self.method_params['time_series'] != 0) : diffexp_command += (' -T ')
         if('min_alignment_count' in self.method_params and self.method_params['min_alignment_count'] is not None ) : diffexp_command += (' -c '+str(self.method_params['min_alignment_count']))
         if('multi_read_correct' in self.method_params and self.method_params['multi_read_correct'] != 0 ): diffexp_command += (' --multi-read-correct ')
         if('library_type' in self.method_params and self.method_params['library_type'] is not None ) : diffexp_command += ( ' --library-type '+self.method_params['library_type'])
         if('library_norm_method' in self.method_params and self.method_params['library_norm_method'] is not None ) : diffexp_command += ( ' --library-norm-method '+self.method_params['library_norm_method'])
         try:
                diffexp_command += " -o {0} -L {1} -u {2} {3}".format(output_dir,self.t_labels,merged_gtf,self.bam_files)
                logger.info("Executing: diffexp {0}".format(diffexp_command))
                ret = script_util.runProgram(None,"diffexp",diffexp_command,None,self.directory)
                result = ret["result"]
                #error =  ret['stderr']
                #print result
                for line in result.splitlines(False):
                       logger.info(line)
                       stderr = ret["stderr"]
                       prev_value = ''
                       for line in stderr.splitlines(False):
                           if line.startswith('> Processing Locus'):
                                   words = line.split()
                                   cur_value = words[len(words) - 1]
                                   if prev_value != cur_value:
                                      prev_value = cur_value
                                      logger.info(line)
                                   else:
                                      prev_value = ''
                                      logger.info(line)
         except Exception,e:
                raise Exception(e)
                raise Exception("Error executing diffexp {0},{1}".format(diffexp_command,e))
         try:
                 logger.info("Zipping DiffExpforBallgown output")
                 out_file_path = os.path.join(self.directory,"{0}.zip".format(self.method_params['output_obj_name']))
                 script_util.zip_files(logger,output_dir,out_file_path)
         except Exception,e:
                 raise Exception("Error executing diffexp")
         try:
                 handle = self.common_params['hs_client'].upload(out_file_path)
         except Exception, e:
                 print " ".join(traceback.print_exc())
                 raise Exception("Failed to upload the DiffExpforBallgown output files: {0}".format(out_file_path))
         ## Save object to workspace
         try:
                 logger.info("Saving DiffExpforBallgown object to workspace")
                 self.cm_obj = { "tool_used" : self.tool_used,
                            "tool_version" : self.tool_version,
                            "condition" : self.details['labels'].split(","),
                            "genome_id" : self.details['genome_id'],
                            "expressionSet_id" : self.details['expressionset_id'],
                            "alignmentSet_id": self.details['alignmentset_id'],
                            "sampleset_id" : self.details['sampleset_id'],
                            "file" : handle
                           }
                 print self.cm_obj
         except Exception , e:
                raise Exception("Error Running DiffExpforBallgown {0} ".format(e))


    def collect(self):
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        #rscripts_dir = self.common_params['rscripts_dir']
        rscripts_dir = '/kb/module/rscripts'
        params = self.method_params
        token = self.common_params['user_token']
        diffexp_dir = self.directory
        logger = self.logger
        logger.info( 'in DiffExpforBallgown.collect')
        output_obj_name = self.method_params['output_obj_name']
        output_csv = "ballgown_diffexp.csv"
        stringtie_dir_prefix = "StringTie_outdir_"

        # 
        #  1) need a pattern RE to match all the StringTie subdirs, so prefix all
        #     unzipped dirs with "stringtie_out_"
        #  2) need a group identifier string i.e. "111000"
        #

        ballgown_sets = rnaseq_util.download_for_ballgown( logger, 
                                                           ws_client, 
                                                           hs, 
                                                           params['ws_id'],
                                                           self.urls,
                                                           diffexp_dir,
                                                           stringtie_dir_prefix,
                                                           params['expressionset_id'],
                                                           token
                                                          )
        logger.info( 'back from download_for_ballgown(), ballgown_sets are')
        logger.info( pformat( ballgown_sets ) )

        group_str = "1100"
        res = rnaseq_util.run_ballgown_diff_exp( logger, rscripts_dir, diffexp_dir, stringtie_dir_prefix, group_str, output_csv )

        raise Exception( "Debug quit")

        #load_diff_matrix( logger, ws_client, token, output_csv, output_object_name )    # ws id, etc

        ##################################        
        #try:
        #    e_sample = ws_client.get_objects( [ {'name' : params['expressionset_id'],
        #                                         'workspace' : params['ws_id']
        #                                         }
        #                                      ])[0]
        #except Exception,e:
        #    logger.exception("".join(traceback.format_exc()))
        #    raise Exception("Error Downloading objects from the workspace ")
        ### Get the Input object type and info #
        #e_sample_info = ws_client.get_object_info_new({"objects": [{'name': params['expressionset_id'], 'workspace': params['ws_id']}]})[0]
        #e_sample_type = e_sample_info[2].split('-')[0]
        #expressionset_id = str(e_sample_info[6]) + '/' + str(e_sample_info[0]) + '/' + str(e_sample_info[4])
        #alignmentset_id = e_sample['data']['alignmentSet_id'] 
        #sampleset_id = e_sample['data']['sampleset_id']
        #expression_ids = e_sample['data']['sample_expression_ids']
        #num_samples = len(expression_ids)
#
        #if num_samples < 2:
        #   raise ValueError("Please ensure you have atleast 2 expressions to run diffExpCallforBallgown in Set mode")
        ####################################

        res1 = self.common_params['ws_client'].save_objects(
                                     {"workspace":self.method_params['ws_id'],
                                      "objects": [{
                                      "type":"KBaseRNASeq.RNASeqDifferentialExpression",
                                      "data":self.cm_obj,
                                      "name":output_name}]})
        returnVal = { 'output'  : output_name ,'workspace' : self.method_params['ws_id']}
