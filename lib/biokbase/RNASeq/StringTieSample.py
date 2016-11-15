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
from biokbase.RNASeq import rnaseq_util
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from biokbase.RNASeq.StringTie import StringTie 

class StringTieSampleException(Exception):
    pass

class StringTieSample(StringTie): 

    def __init__(self, logger, directory, urls):
        super(self.__class__, self).__init__(logger, directory, urls)
        # user defined shared variables across methods
        self.alignment_info = None
        #self.sampleset_info = None
        self.num_threads = 1

    def prepare(self): 
        # for quick testing, we recover parameters here
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        stringtie_dir = self.directory
        try:
               a_sample = ws_client.get_objects(
                                        [{'name' : params['alignmentset_id'],'workspace' : params['ws_id']}])[0]
               a_alignment_info = ws_client.get_object_info_new({"objects" :
                                        [{'name' : params['alignmentset_id'],'workspace' : params['ws_id']}]})[0]
               self.alignment_info =  a_alignment_info
               s_alignment_id = str(a_alignment_info[6]) + '/' + str(a_alignment_info[0]) + '/' + str(a_alignment_info[4])
	except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise Exception("Error Downloading objects from the workspace ")
        read_sample_id = a_sample['data']['read_sample_id']
	### Check if the gtf file exists in the workspace. if exists download the file from that
        genome_id = a_sample['data']['genome_id']
        genome_name = ws_client.get_object_info([{"ref" :genome_id}],includeMetadata=None)[0][1]
        ws_gtf = genome_name+"_GTF_Annotation"
	gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,self.urls,params['ws_id'],ws_gtf,"KBaseRNASeq.GFFAnnotation",stringtie_dir,token)
        if gtf_file is None:
             rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],genome_id,genome_name,stringtie_dir,token)
	gtf_info = ws_client.get_object_info_new({"objects": [{'name': ws_gtf , 'workspace': params['ws_id']}]})[0]
        gtf_id = str(gtf_info[6]) + '/' + str(gtf_info[0]) + '/' + str(gtf_info[4])
	self.tool_opts = { k:str(v) for k,v in params.iteritems() if not k in ('ws_id','alignmentset_id', 'num_threads') and v is not None  }
	self.num_jobs =  1
	logger.info(" Number of threads used by each process {0}".format(self.num_threads))
        task_param = {'job_id' : s_alignment_id,
                      'gtf_file' : gtf_file,
                      'ws_id' : params['ws_id'],
                      'genome_id' : genome_id,
                      'stringtie_dir' : self.directory,
                      'annotation_id': gtf_id,
		      'sample_id' : read_sample_id , 
                      'alignmentset_id' : None
                      }
        self.task_list.append(task_param)

    def collect(self) :
        # do with 
        alignment_name = self.method_params['alignmentset_id']+"_stringtie_AlignmentSet"
        self.logger.info(" Creating Report for Alignment {0}".format(alignment_name))
	single_alignment , single_expression = self.results[0]
        # TODO: Split alignment set and report method
	sref = self.common_params['ws_client'].get_object_info_new({"objects": [{'name':single_expression, 'workspace': self.method_params['ws_id']}]})[0]
	self.returnVal = { 'output'  : single_expression ,'workspace' : self.method_params['ws_id']}
#        reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
#                                                 'description' : "RNA-seq Expression for read alignment : {0}".format(single_alignment)}],
#                                                 'text_message': "RNA-seq Alignment for reads alignment : {0}".format(single_alignment)}
#	reportName = 'Align_Reads_using_stringtie_'+str(hex(uuid.getnode()))
#        report_info = self.common_params['ws_client'].save_objects({
#                                                'id':self.alignment_info[6],
#                                                'objects':[
#                                                {
#                                                'type':'KBaseReport.Report',
#                                                'data':reportObj,
#                                                'name':reportName,
#                                                'meta':{},
#                                                'hidden':1, # important!  make sure the report is hidden
#                                                #'provenance':provenance
#                                                }
#                                                ]
#                                                })[0]
#
#        self.returnVal = { "report_name" : reportName,"report_ref" : str(report_info[6]) + '/' + str(report_info[0]) + '/' + str(report_info[4]) }
#
