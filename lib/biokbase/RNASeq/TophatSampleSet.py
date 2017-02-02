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
from biokbase.RNASeq.Tophat import Tophat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
#import ExecutionBase.ExecutionBase as ExecutionBase

class TophatSampleSetException(Exception):
    pass

class TophatSampleSet(Tophat): 

    def __init__(self, logger, directory, urls, max_cores):
        super(TophatSampleSet, self).__init__(logger, directory, urls, max_cores)

        # user defined shared variables across methods
        self.sample = None
	self.bowtie2index_id = None
        #self.num_threads = None


    def prepare(self): 
        # for quick testing, we recover parameters here
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        tophat_dir = self.directory

        try:
               #sample ,bowtie_index = ws_client.get_objects(
               #                         [{'name' : params['sampleset_id'],'workspace' : params['ws_id']},
               #                          { 'name' : params['bowtie_index'], 'workspace' : params['ws_id']}])
               sample = script_util.ws_get_obj(logger, ws_client, params['ws_id'],params['sampleset_id'])[0]
               bowtie_index = script_util.ws_get_obj(logger, ws_client, params['ws_id'],params['bowtie_index'])[0]
               self.sample = sample
        except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise ValueError(" Error Downloading objects from the workspace ")
        ### Get object Info and IDs
        sample_info = script_util.ws_get_obj_info(logger, ws_client, params['ws_id'], params['sampleset_id'])[0]
        sample_type = sample_info[2].split('-')[0]

        # SampleSet
        if not (sample_type == 'KBaseRNASeq.RNASeqSampleSet' or sample_type == 'KBaseSets.ReadsSet'):
            raise TophatSampleSetException('RNASeqSampleSet or ReadsSet is required')
        (reads, r_label) = rnaseq_util.get_reads_conditions(logger, sample, sample_type)
        #reads = sample['data']['sample_ids']
        #reads_type= sample['data']['Library_type']
        # Note: do not need the following as we support ws reference
        #e_ws_objs = script_util.if_ws_obj_exists_notype(None,ws_client,params['ws_id'],reads)
        #missing_objs = [i for i in reads if not i in e_ws_objs]
        #if len(e_ws_objs) != len(reads):
        #   raise ValueError('Missing Library objects {0} in the {1}. please copy them and run this method'.format(",".join(missing_objs),params['ws_id']))



	### Get obejct IDs
        #bowtie2_index_info,sampleset_info = ws_client.get_object_info_new({"objects": [{'name': params['bowtie_index'], 'workspace': params['ws_id']},{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
        #self.bowtie2index_id = str(bowtie2_index_info[6]) + '/' + str(bowtie2_index_info[0]) + '/' + str(bowtie2_index_info[4])  
        #sampleset_id = str(sampleset_info[6]) + '/' + str(sampleset_info[0]) + '/' + str(sampleset_info[4]) 
        self.bowtie2index_id = script_util.ws_get_ref(logger, ws_client, params['ws_id'], params['bowtie_index'])
        sampleset_id = script_util.ws_get_ref(logger, ws_client, params['ws_id'], params['sampleset_id'])
        bw_id = bowtie_index['data']['handle']['id'] 
        bw_name =  bowtie_index['data']['handle']['file_name']
        genome_id = bowtie_index['data']['genome_id']
        annotation_gtf = ws_client.get_object_info([{"ref" :genome_id}],includeMetadata=None)[0][1]
        shared_files={}
        shared_files[bw_name] = bw_id
        script_util.download_shock_files(logger,self.urls['shock_service_url'],tophat_dir,shared_files,token)
        try:
            logger.info("Unzipping Bowtie2 Indices")
            script_util.unzip_files(logger,os.path.join(tophat_dir,bw_name),tophat_dir)
            mv_dir= handler_util.get_dir(tophat_dir)
            if mv_dir is not None:
                    script_util.move_files(logger,mv_dir,tophat_dir)
        except Exception, e:
               logger.error("".join(traceback.format_exc()))
               raise Exception("Unzip indexfile error")
        fasta_file =os.path.join(tophat_dir,(handler_util.get_file_with_suffix(tophat_dir,".fa")+".fa"))
        bowtie2base =os.path.join(tophat_dir,handler_util.get_file_with_suffix(tophat_dir,".rev.1.bt2"))

	### Check if GTF annotation object exist or skip this step
	### Check if the gtf object exists in the workspace
        ### Only run create_gtf_annotation if object doesnt exist
	ws_gtf = annotation_gtf+"_GTF_Annotation"

        gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,self.urls,params['ws_id'],ws_gtf,"KBaseRNASeq.GFFAnnotation",tophat_dir,token)
        if gtf_file is None:
            gtf_file = rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],genome_id,genome_name,tophat_dir,token)
	#ret = script_util.if_obj_exists(None,ws_client,params['ws_id'],"KBaseRNASeq.GFFAnnotation",[ws_gtf]) # this line should be safe from reference
        #if not ret is None:
        #    logger.info("GFF Annotation Exist for Genome Annotation {0}.... Skipping step ".format(annotation_gtf))
	#    annot_name,annot_id = ret[0]
        #    gtf_obj=ws_client.get_objects([{'ref' : annot_id}])[0]
        #    gtf_id=gtf_obj['data']['handle']['id']
        #    gtf_name=gtf_obj['data']['handle']['file_name']
        #    try:
        #       script_util.download_file_from_shock(logger, shock_service_url=self.urls['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=tophat_dir,token=token)
        #       gtf_file = os.path.join(tophat_dir,gtf_name)
        #    except Exception,e:
	#       logger.exception(e)
        #       raise Exception( "Unable to download shock file, {0}".format(gtf_name))  
 	#else:		
	#    gtf_file =rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],genome_id,annotation_gtf,tophat_dir,token)		

	# Determine the num_threads provided by the user otherwise default the number of threads to 2
        #reads = sample['data']['sample_ids'] # duplicated lines
        #reads_type= sample['data']['Library_type'] # duplicated lines
        #r_label = sample['data']['condition']
        self.num_jobs =  len(reads)

        count = 0
        logger.info(" Number of threads used by each process {0}".format(self.num_threads))
        for i in reads:
            try:
                    label = r_label[count]
                    task_param = {'job_id' : i,
                                  'label' : r_label[count],
                                  'ws_id' : params['ws_id'],
                                  'tophat_dir' : self.directory,
				  'gtf_file' : gtf_file, # TODO: double check newly created file gtf_file is path or not
                                  'annotation_id': genome_id,
                                  'sampleset_id' : sampleset_id
                                 }
                    self.task_list.append(task_param)
                    count = count + 1
            except Exception,e:
                    raise


    def collect(self) :
        # do with 
        alignmentSet_name = script_util.ws_get_obj_name(self.logger, 
                                                        self.common_params['ws_client'], 
                                                        self.method_params['ws_id'], 
                                                        self.method_params['sampleset_id'])+"_tophat_AlignmentSet"
        self.logger.info(" Creating AlignmentSet for the Alignments {0}".format(alignmentSet_name))
        # TODO: Split alignment set and report method
        reportObj=rnaseq_util.create_RNASeq_AlignmentSet_and_build_report(self.logger,self.common_params['ws_client'],self.method_params['ws_id'],rnaseq_util.get_reads(self.logger, self.sample),self.task_list[0]['sampleset_id'],self.task_list[0]['annotation_id'],self.bowtie2index_id,self.results,alignmentSet_name)
        self.returnVal = { 'output'  : alignmentSet_name ,'workspace' : self.method_params['ws_id']}
#	reportName = 'Align_Reads_using_Tophat_'+str(hex(uuid.getnode()))
#        report_info = self.common_params['ws_client'].save_objects({
#                                                'workspace':self.method_params['ws_id'],
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
