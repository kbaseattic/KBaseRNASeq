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
from biokbase.RNASeq.HiSat2 import HiSat2

class HiSat2SampleException(Exception):
    pass

class HiSat2Sample(HiSat2): 

    def __init__(self, logger, directory, urls):
        super(HiSat2Sample, self).__init__(logger, directory, urls)
        # user defined shared variables across methods
        self.sample_info = None
        #self.sampleset_info = None
        self.num_threads = 1

    def prepare(self): 
        # for quick testing, we recover parameters here
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        hisat2_dir = self.directory

        try:
               sample,annotation_name = ws_client.get_objects(
                                        [{ 'name' : params['sampleset_id'], 'workspace' : params['ws_id']},
                                        { 'name' : params['genome_id'], 'workspace' : params['ws_id']}])
               self.sample =sample 
        except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise ValueError(" Error Downloading objects from the workspace ")
        ### Get object Info and  IDs
        sample_info,annotation_info = ws_client.get_object_info_new({"objects": [
                                           {'name': params['sampleset_id'], 'workspace': params['ws_id']},
                                           {'name': params['genome_id'], 'workspace': params['ws_id']}
                                           ]})
        self.sample_info = sample_info
        ### Get the workspace object ids for the objects ###
        sample_id = str(sample_info[6]) + '/' + str(sample_info[0]) + '/' + str(sample_info[4])
        annotation_id = str(annotation_info[6]) + '/' + str(annotation_info[0]) + '/' + str(annotation_info[4])
        sample_type = sample_info[2].split('-')[0]
	lib_types = ['KBaseAssembly.SingleEndLibrary', 'KBaseAssembly.PairedEndLibrary']
        ### Check if the Library objects exist in the same workspace
        if not sample_type in lib_types: #'KBaseAssembly.SingleEndLibrary' or sample_type != 'KBaseAssembly.PairedEndLibrary':
            raise HiSat2SampleException('Either of the Library typed objects SingleEndLibrary or PairedEndLibrary is required')
        r_label = 'Single'
        self.num_jobs = 1
	### Get the Genome Id for the genome selected and get fasta file
	ref_id , fasta_file =  rnaseq_util.get_fa_from_genome(logger,ws_client,self.urls,params['ws_id'],hisat2_dir,params['genome_id'])

	### Build Index for the fasta file
	hisat2base = os.path.basename(fasta_file) 
        #hisat2base =os.path.join(hisat2_dir,handler_util.get_file_with_suffix(hisat2_dir,".fa"))
        hisat2base_cmd = '{0} {1}'.format(fasta_file,hisat2base)
	try:
            logger.info("Building Index for Hisat2 {0}".format(hisat2base_cmd))
            cmdline_output = script_util.runProgram(logger,"hisat2-build",hisat2base_cmd,None,hisat2_dir)
        except Exception,e:
            raise Exception("Failed to run command {0}".format(hisat2base_cmd))
        ### Check if GTF object exists in the workspace pull the gtf
	ws_gtf = params['genome_id']+"_GTF"
	gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,self.urls,params['ws_id'],ws_gtf,"KBaseRNASeq.GFFAnnotation",hisat2_dir,token)
        if gtf_file is None:
             rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],ref_id,params['genome_id'],hisat2_dir,token)
	# Determine the num_threads provided by the user otherwise default the number of threads to 2
        logger.info(" Number of threads used by each process {0}".format(self.num_threads))
	task_param = {'job_id' : params['sampleset_id'],
                      'label' : r_label,
                      'ws_id' : params['ws_id'],
                      'reads_type' : sample_type,
                      'hisat2_dir' : self.directory,
                      'annotation_id': ref_id,
                      'sampleset_id' : None
                      }
	self.task_list.append(task_param)
	
		
    def collect(self) :
        # do with 
        alignment_name = self.method_params['sampleset_id']+"_hisat2_AlignmentSet"
        self.logger.info(" Creating Report for Alignment {0}".format(alignment_name))
	single_read , single_alignment = self.results[0]
        # TODO: Split alignment set and report method
	sref = self.common_params['ws_client'].get_object_info_new({"objects": [{'name':single_alignment, 'workspace': self.method_params['ws_id']}]})[0]
	self.returnVal = { 'output'  : single_alignment ,'workspace' : self.method_params['ws_id']}
#        reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
#                                                 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_read)}],
#                                                 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_read)}
#	reportName = 'Align_Reads_using_Hisat2_'+str(hex(uuid.getnode()))
#        report_info = self.common_params['ws_client'].save_objects({
#                                                'id':self.sample_info[6],
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
