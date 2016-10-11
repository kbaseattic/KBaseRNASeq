import simplejson, sys, shutil, os, ast , re
from mpipe import OrderedStage , Pipeline
import glob, json, uuid, logging  , time ,datetime 
import subprocess, threading,traceback
from collections import OrderedDict
from pprint import pprint , pformat
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
#import ExecutionBase.ExecutionBase as ExecutionBase

class CufflinksException(Exception):
    pass

class Cufflinks(ExecutionBase): 

    def __init__(self, logger, directory, urls):
        pprint(self.__class__)
        super(Cufflinks, self).__init__(logger, directory, urls)

        # user defined shared variables across methods
        #self.sample = None
        #self.sampleset_info = None
        self.num_threads = None
	self.tool_used = "Cufflinks"
	self.tool_version = "1.2.3"

    def runEach(self,task_params):
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        
        s_alignment = task_params['job_id']
        gtf_file = task_params['gtf_file']
        directory = task_params['cufflinks_dir']
        genome_id = task_params['genome_id']
        annotation_id = task_params['annotation_id']
        sample_id = task_params['sample_id']
        alignmentset_id = task_params['alignmentset_id']
        ws_id = task_params['ws_id']
    
        print "Downloading Sample Alignment from workspace {0}".format(s_alignment)
        logger.info("Downloading Sample Alignment from workspace {0}".format(s_alignment))
        alignment_name = ws_client.get_object_info([{"ref" :s_alignment}],includeMetadata=None)[0][1]
        if not logger:
           logger = handler_util.create_logger(directory,"run_cufflinks_"+alignment_name)
        try:
           alignment = ws_client.get_objects(
                                        [{ 'ref' : s_alignment }])[0]
           input_direc = os.path.join(directory,alignment_name.split('_alignment')[0]+"_cufflinks_input")
           if not os.path.exists(input_direc) : os.mkdir(input_direc)
           output_name = alignment_name.split('_alignment')[0]+"_cufflinks_expression"
           output_dir = os.path.join(directory,output_name)
           #Download Alignment from shock
           a_file_id = alignment['data']['file']['id']
           a_filename = alignment['data']['file']['file_name']
           condition = alignment['data']['condition']
           try:
                script_util.download_file_from_shock(logger, shock_service_url=self.urls['shock_service_url'], shock_id=a_file_id,filename=a_filename,directory=input_direc,token=token)
           except Exception,e:
                raise Exception( "Unable to download shock file, {0},{1}".format(a_filename,"".join(traceback.format_exc())))
           try:
                input_dir = os.path.join(input_direc,alignment_name)
		if not os.path.exists(input_dir): os.mkdir(input_dir)
                script_util.unzip_files(logger,os.path.join(input_direc,a_filename), input_dir)
           except Exception, e:
		raise Exception(e)
                logger.error("".join(traceback.format_exc()))
                raise Exception("Unzip alignment files  error: Please contact help@kbase.us")

           input_file = os.path.join(input_dir,"accepted_hits.bam")
                ### Adding advanced options to tophat command
           tool_opts = { k:str(v) for k,v in params.iteritems() if not k in ('ws_id','alignmentset_id', 'num_threads') and v is not None  }
           cufflinks_command = (' -p '+str(self.num_threads))
	   if 'max_intron_length' in params and params['max_intron_length'] is not None:
               cufflinks_command += (' --max-intron-length '+str(params['max_intron_length']))
           if 'min_intron_length' in params and params['min_intron_length'] is not None:
               cufflinks_command += (' --min-intron-length '+str(params['min_intron_length']))
           if 'overhang_tolerance' in params  and params['overhang_tolerance'] is not None:
               cufflinks_command += (' --overhang-tolerance '+str(params['overhang_tolerance']))

           cufflinks_command += " -o {0} -G {1} {2}".format(output_dir,gtf_file,input_file)
           #cufflinks_command += " -o {0} -A {1} -G {2} {3}".format(t_file_name,g_output_file,gtf_file,input_file)
           logger.info("Executing: cufflinks {0}".format(cufflinks_command))
           print "Executing: cufflinks {0}".format(cufflinks_command)
           ret = script_util.runProgram(None,"cufflinks",cufflinks_command,None,directory)
           result = ret["result"]
	   for line in result.splitlines(False):
                       self.logger.info(line)
                       stderr = ret["stderr"]
                       prev_value = ''
                       for line in stderr.splitlines(False):
                           if line.startswith('> Processing Locus'):
                                   words = line.split()
                                   cur_value = words[len(words) - 1]
                                   if prev_value != cur_value:
                                      prev_value = cur_value
                                      self.logger.info(line)
                                   else:
                                      prev_value = ''
                                      self.logger.info(line)

	   ##Parse output files
           try:
	 	g_output_file = os.path.join(output_dir,"genes.fpkm_tracking")
                exp_dict = script_util.parse_FPKMtracking(g_output_file,'Cufflinks','FPKM')
                #tpm_exp_dict = script_util.parse_FPKMtracking(g_output_file,'Cufflinks','TPM')
           except Exception,e:
	        raise Exception(e)
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error parsing FPKMtracking")
        ##  compress and upload to shock
           try:
                logger.info("Zipping cufflinks output")
                print "Zipping cufflinks output"
                out_file_path = os.path.join(directory,"%s.zip" % output_name)
                script_util.zip_files(logger,output_dir,out_file_path)
           except Exception,e:
	        raise Exception(e)
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error executing cufflinks")
           try:
                handle = hs.upload(out_file_path)
           except Exception, e:
	        raise Exception(e)
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error while zipping the output objects: {0}".format(out_file_path))
                ## Save object to workspace
           try:
                logger.info("Saving cufflinks object to workspace")
                es_obj = { 'id' : output_name,
                           'type' : 'RNA-Seq',
                           'numerical_interpretation' : 'FPKM',
                           'expression_levels' : exp_dict,
                           #'tpm_expression_levels' : tpm_exp_dict,
                           'processing_comments' : "log2 Normalized",
                           'genome_id' : genome_id,
                           'annotation_id' : annotation_id,
                           'condition' : condition,
                           'mapped_rnaseq_alignment' : { sample_id : s_alignment },
                           'tool_used' : self.tool_used,
                           'tool_version' : self.tool_version,
                           'tool_opts' : tool_opts,
                           'file' : handle
                         }

                res= ws_client.save_objects(
                                   {"workspace":ws_id,
                                    "objects": [{
                                    "type":"KBaseRNASeq.RNASeqExpression",
                                    "data":es_obj,
                                    "name":output_name}
                                     ]})[0]
                expr_id = str(res[6]) + '/' + str(res[0]) + '/' + str(res[4])
           except Exception, e:
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Failed to upload the ExpressionSample: {0}".format(output_name))
   	except Exception,e:
       		logger.exception("".join(traceback.format_exc()))
       		raise Exception("Error executing cufflinks {0},{1}".format(cufflinks_command,directory))
   	finally:
                if os.path.exists(out_file_path): os.remove(out_file_path)
                if os.path.exists(output_dir): shutil.rmtree(output_dir)
                if os.path.exists(input_direc): shutil.rmtree(input_direc)
                ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqExpression",[output_name])
                if not ret is None:
                    return (alignment_name, output_name )
        return None

