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
        self.details["used_tool"] = "Ballgown (Bioconductor)"    # Question: where does this really get set?
        self.details["tool_version"] = "3.4"
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
         if used_tool == "Ballgown (Bioconductor)":
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
        params    = self.method_params
        ws_client = self.common_params['ws_client']
        hs_client = self.common_params['hs_client']

        ws_id = params['ws_id']
        #rscripts_dir = self.common_params['rscripts_dir']
        rscripts_dir = '/kb/module/rscripts'

        token = self.common_params['user_token']
        diffexp_dir = self.directory
        logger = self.logger
        logger.info( 'in DiffExpforBallgown.collect, method params (params) are')
        logger.info( pformat( params ) )
        output_object_name = params['output_obj_name']
        output_csv = "ballgown_diffexp.csv"
        stringtie_dir_prefix = "StringTie_outdir_"

        # 
        #  1) need a pattern RE to match all the StringTie subdirs, so prefix all
        #     unzipped dirs with "stringtie_out_"
        #  2) need a group identifier string i.e. "111000"
        #

        ballgown_set_info = rnaseq_util.get_info_and_download_for_ballgown( logger, 
                                                                            ws_client, 
                                                                            hs_client, 
                                                                            ws_id,
                                                                            self.urls,
                                                                            diffexp_dir,
                                                                            stringtie_dir_prefix,
                                                                            params['expressionset_id'],
                                                                            token
                                                                           )
        logger.info( 'back from download_for_ballgown(), ballgown_set_info are')
        logger.info( pformat( ballgown_set_info ) )

        # THIS IS TEMPORARY: REMOVE THIS WHEN group lists of expression object names are actually passed
        # incorporate params['group1_set'] and params['group2_set'] which are lists of sample names
        params['group1_name'] = "WT"
        params['group2_name'] = "exp"
        group1_set = []
        group2_set = []
        for subd in ballgown_set_info['subdirs']:
            if ( re.search( "(^WT|_WT|WT_|WT$)", subd, re.I ) ):
                group1_set.append( subd )
            else:
                group2_set.append( subd )
        params['group1_set'] = group1_set
        params['group2_set'] = group2_set
        # END OF TEMPORARY CODE 
        sample_dir_group_file = "sample_dir_group_table"  # output file
        group_list = rnaseq_util.create_sample_dir_group_file( ballgown_set_info['subdirs'], 
                                                               params['group1_name'],
                                                               params['group1_set'],
                                                               params['group2_name'],
                                                               params['group2_set'],
                                                               sample_dir_group_file )

        ballgown_output_dir = os.path.join( diffexp_dir, "ballgown_out" )
        logger.info( "ballgown output dir is {0}".format( ballgown_output_dir) )
        handler_util.setupWorkingDir( logger, ballgown_output_dir )

        logger.info( "about to run_ballgown_diff_exp" )
        rnaseq_util.run_ballgown_diff_exp( logger, rscripts_dir, diffexp_dir, sample_dir_group_file, ballgown_output_dir, output_csv )

        logger.info( "back from run_ballgown_diff_exp, about to load diff exp matrix file" )
        diff_expr_matrix = rnaseq_util.load_diff_expr_matrix( ballgown_output_dir, output_csv )    # read file before its zipped

        logger.info( "about to load ballgout output into workspace" )
        de_ws_save_obj_data = rnaseq_util.load_ballgown_output_into_ws( logger, 
                                                                        ws_id,
                                                                        ws_client, 
                                                                        hs_client,
                                                                        token, 
                                                                        diffexp_dir,
                                                                        ballgown_output_dir, 
                                                                        self.details["used_tool"],
                                                                        self.details["tool_version"],
                                                                        ballgown_set_info['sample_expression_ids'],  # for sample ids? Is this good?
                                                                        group_list,                                  # conditions
                                                                        ballgown_set_info['genome_id'],              # genome_id
                                                                        ballgown_set_info['expressionset_id'],       # expressionset_id
                                                                        ballgown_set_info['alignmentSet_id'],        # alignmentset_id
                                                                        ballgown_set_info['sampleset_id'],           # sampleset_id
                                                                        output_object_name 
                                                                      )
        logger.info( "back from loading ballgown output into workspace, object save data is " )
        logger.info( pformat( de_ws_save_obj_data ) )

        selected_gene_list = rnaseq_util.filter_genes_diff_expr_matrix( diff_expr_matrix, 
                                                                        params['fold_scale_type'], 
                                                                        params['alpha_cutoff'], 
                                                                        params['log2_fold_change_cutoff'],
                                                                        params['maximum_number_of_genes']
                                                                      )
        
        #  !!!!! IF selected_gene_list is empty print some kind of message, take no further action

        # get the unfiltered expression matrix
        em_name = params['expressionset_id'] + "_FPKM_ExpressionMatrix"
        logger.info( "about to fetch expression matrix  {0}".format( em_name ))
        try:
            emw = ws_client.get_objects( [ { "name": em_name, "workspace": ws_id } ] )[0]
        except:
            raise Exception( "unable to retrieve expression matrix object {0} from workspace {1}".format( em_name, ws_id  ))
        logger.info( pformat( emw ) )
        emo = emw["data"]
        # filter it
        filtered_emo = rnaseq_util.filter_expr_matrix_object( emo, selected_gene_list )
        # save it
        logger.info( "saving emo em_name {0}".format( em_name ))
        try:
            ret = ws_client.save_objects( { 'workspace' : ws_id,
                                            'objects' : [
                                                          { 'type'   : 'KBaseFeatureValues.ExpressionMatrix',
                                                            'data'   : filtered_emo,
                                                            'name'   : params["filtered_expr_matrix_name"]
                                                          }
                                                        ]
                                          }
                                        )
        except:
            raise Exception( "failed to save object " )
        logger.info( "ws save return:\n" + pformat(ret))

        # THIS NEEDS TO BE AN INPUT PARAMETER IN SPEC FILE
        #iltered_expr_matrix_name = expressionset_id + "_filtered_fpkm"
        #e_em_save_obj_data = created_and_save_filtered_expr_matrix( logger, 
        #                                                            ws_client, 
        #                                                            ws_id, 
        #                                                            token,
        #                                                            expression_set_name, 
        #                                                            fold_scale_type,      #"linear", "log2+1", "log10+1" 
        #                                                            alpha_cutoff,
        #                                                            q_value_cutoff,
        #                                                            log2_fold_change_cutoff,
        #                                                            maximum_num_genes,
        #                                                            filtered_expr_matrix_name
        #                                                           )

        returnVal = { 'output'  : output_object_name ,
                      #'filtered_expression_maxtrix': filtered_expr_matrix_name, 
                      'workspace' : ws_id }
