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
            ### Get obejct IDs
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
	ref_info = ws_client.get_object_info_new({"objects": [{'name': params['genome_id'], 'workspace': params['ws_id']}]})[0]
        ref_id = str(ref_info[6]) + '/' + str(ref_info[0]) + '/' + str(ref_info[4])
        fasta_file = script_util.get_fasta_from_genome(logger,ws_client,self.urls,ref_id)
	## TODO Remove the below commented lines of code to complete get rid of genome annotation and data_api calls
	#fasta_file = script_util.generate_fasta(logger,self.urls,token,annotation_id,hisat2_dir,params['genome_id'])
        #logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
        #mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
        #c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
        #logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
        hisat2base =os.path.join(hisat2_dir,handler_util.get_file_with_suffix(hisat2_dir,".fasta"))
        hisat2base_cmd = '{0} {1}'.format(fasta_file,hisat2base)
	try:
            logger.info("Building Index for Hisat2 {0}".format(hisat2base_cmd))
            cmdline_output = script_util.runProgram(logger,"hisat2-build",hisat2base_cmd,None,hisat2_dir)
        except Exception,e:
            raise Exception("Failed to run command {0}".format(hisat2base_cmd))
        ws_gtf = params['genome_id']+"_GTF"
        ret = script_util.if_obj_exists(None,ws_client,params['ws_id'],"KBaseRNASeq.GFFAnnotation",[ws_gtf])
        if not ret is None:
            logger.info("GFF Annotation Exist for Genome Annotation {0}.... Skipping step ".format(params['genome_id']))
            annot_name,annot_id = ret[0]
            gtf_obj=ws_client.get_objects([{'ref' : annot_id}])[0]
            gtf_id=gtf_obj['data']['handle']['id']
            gtf_name=gtf_obj['data']['handle']['file_name']
 	    try:
                     script_util.download_file_from_shock(logger, shock_service_url=self.urls['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=hisat2_dir,token=token)
                     gtf_file = os.path.join(hisat2_dir,gtf_name)
            except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(gtf_name))
        else:
             script_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],ref_id,params['genome_id'],fasta_file,hisat2_dir,token)
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
	single_alignment , single_read = self.results[0]
        # TODO: Split alignment set and report method
	sref = self.common_params['ws_client'].get_object_info_new({"objects": [{'name':single_alignment, 'workspace': self.method_params['ws_id']}]})[0]
        reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
                                                 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_read)}],
                                                 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_read)}
	reportName = 'Align_Reads_using_Hisat2_'+str(hex(uuid.getnode()))
        report_info = self.common_params['ws_client'].save_objects({
                                                'id':self.sample_info[6],
                                                'objects':[
                                                {
                                                'type':'KBaseReport.Report',
                                                'data':reportObj,
                                                'name':reportName,
                                                'meta':{},
                                                'hidden':1, # important!  make sure the report is hidden
                                                #'provenance':provenance
                                                }
                                                ]
                                                })[0]

        returnVal = { "report_name" : reportName,"report_ref" : str(report_info[6]) + '/' + str(report_info[0]) + '/' + str(report_info[4]) }

