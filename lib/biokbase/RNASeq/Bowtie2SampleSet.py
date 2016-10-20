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
from biokbase.RNASeq import rnaseq_util
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from biokbase.RNASeq.Bowtie2 import Bowtie2
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
#import ExecutionBase.ExecutionBase as ExecutionBase

class Bowtie2SampleSetException(Exception):
    pass

class Bowtie2SampleSet(Bowtie2): 

    def __init__(self, logger, directory, urls):
        super(self.__class__, self).__init__(logger, directory, urls)

        # user defined shared variables across methods
        self.sample = None
        self.sampleset_info = None
        #self.num_threads = None


    def prepare(self): 
        # for quick testing, we recover parameters here
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        bowtie2_dir = self.directory

        try:
	       sample,bowtie_index = ws_client.get_objects( 
                                        [{ 'name' : params['sampleset_id'], 'workspace' : params['ws_id']}, 
                                        { 'name' : params['bowtie_index'], 'workspace' : params['ws_id']}]) 
               self.sample =sample 
        except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise ValueError(" Error Downloading objects from the workspace ")
            ### Get obejct IDs
        sampleset_info,bowtie_index_info = ws_client.get_object_info_new({"objects": [
                                           {'name': params['sampleset_id'], 'workspace': params['ws_id']},
                                           {'name': params['bowtie_index'], 'workspace': params['ws_id']}
                                           ]})
        ### Get the workspace object ids for the objects ###
        sampleset_id = str(sampleset_info[6]) + '/' + str(sampleset_info[0]) + '/' + str(sampleset_info[4])
        bowtie_index_id = str(bowtie_index_info[6]) + '/' + str(bowtie_index_info[0]) + '/' + str(bowtie_index_info[4])
        self.sampleset_info = sampleset_info
        ### Get the workspace object ids for the objects ###
        sample_type = sampleset_info[2].split('-')[0]

        ### Check if the Library objects exist in the same workspace
        if sample_type != 'KBaseRNASeq.RNASeqSampleSet':
            raise Bowtie2SampleSetException('RNASeqSampleSet is required')
        logger.info("Check if the Library objects do exist in the current workspace")
        reads = sample['data']['sample_ids']
        r_label = sample['data']['condition']
        reads_type= sample['data']['Library_type']
        if reads_type == 'PairedEnd': r_type = 'KBaseAssembly.PairedEndLibrary'
        else: r_type = 'KBaseAssembly.SingleEndLibrary'
        e_ws_objs = script_util.if_ws_obj_exists(None,ws_client,params['ws_id'],r_type,reads)
        missing_objs = [i for i in reads if not i in e_ws_objs]
        if len(e_ws_objs) != len(reads):
            raise Bowtie2SampleSetException('Missing Library objects {0} in the {1}. please copy them and run this method'.format(",".join(missing_objs),params['ws_id']))
 
        self.num_jobs = len(reads)
	bw_index_files = script_util.check_and_download_existing_handle_obj(logger,ws_client,self.urls,params['ws_id'],params['bowtie_index'],"KBaseRNASeq.Bowtie2Indexes",bowtie2_dir,token)
        try:
                logger.info("Unzipping Bowtie2 Indices")
                script_util.unzip_files(logger,os.path.join(bowtie2_dir,bw_index_files),bowtie2_dir)
                mv_dir= handler_util.get_dir(bowtie2_dir)
                if mv_dir is not None:
                        script_util.move_files(logger,mv_dir,bowtie2_dir)
        except Exception, e:
                logger.error("".join(traceback.format_exc()))
                raise Exception("Unzip indexfile error: Please contact help@kbase.us")
        ### Build Index for the fasta file 
        fasta_file =os.path.join(bowtie2_dir,handler_util.get_file_with_suffix(bowtie2_dir,".fa")+".fa")
        bowtie2base =os.path.join(bowtie2_dir,handler_util.get_file_with_suffix(bowtie2_dir,".fa"))
        bowtie2base_cmd = '{0} {1}'.format(fasta_file,bowtie2base)
        try:
            logger.info("Building Index for Hisat2 {0}".format(bowtie2base_cmd))
            cmdline_output = script_util.runProgram(logger,"bowtie2-build",bowtie2base_cmd,None,bowtie2_dir)
        except Exception,e:
            raise Exception("Failed to run command {0}".format(bowtie2base_cmd))
        ### Check if GTF object exists in the workspace pull the gtf
        ref_id = bowtie_index['data']['genome_id']
        genome_name = ws_client.get_object_info_new({"objects": [{'ref' : ref_id }] })[0][1]
        ws_gtf = genome_name+"_GTF"
        gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,self.urls,params['ws_id'],ws_gtf,"KBaseRNASeq.GFFAnnotation",bowtie2_dir,token)
        if gtf_file is None:
	     rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],ref_id,genome_name,bowtie2_dir,token)
 
        count = 0
        logger.info(" Number of threads used by each process {0}".format(self.num_threads))
        for i in reads:
            try:
                    label = r_label[count]
                    task_param = {'job_id' : i,
                                  'label' : r_label[count],
                                  'ws_id' : params['ws_id'],
                                  'reads_type' : reads_type,
                                  'bowtie2_dir' : self.directory,
                                  'annotation_id': ref_id, # Changed annotation_id to ref_id for genome object 
                                  'sampleset_id' : sampleset_id
                                 }
                    self.task_list.append(task_param)
                    count = count + 1
            except Exception,e:
                    raise


    def collect(self) :
        # do with 
        alignmentSet_name = self.method_params['sampleset_id']+"_bowtie2_AlignmentSet"
        self.logger.info(" Creating AlignmentSet for the Alignments {0}".format(alignmentSet_name))
        # TODO: Split alignment set and report method
        reportObj=rnaseq_util.create_RNASeq_AlignmentSet_and_build_report(self.logger,self.common_params['ws_client'],self.method_params['ws_id'],self.sample['data']['sample_ids'],self.task_list[0]['sampleset_id'],self.task_list[0]['annotation_id'],None,self.results,alignmentSet_name)
	self.returnVal = { 'output'  : alignmentSet_name ,'workspace' : self.method_params['ws_id']}
#        reportName = 'Align_Reads_using_Hisat2_'+str(hex(uuid.getnode()))
#        report_info = self.common_params['ws_client'].save_objects({
#                                                'id':self.sampleset_info[6],
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
