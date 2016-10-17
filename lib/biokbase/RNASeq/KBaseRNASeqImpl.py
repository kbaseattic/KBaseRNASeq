# -*- coding: utf-8 -*-
#BEGIN_HEADER
import simplejson
import sys
import shutil
import os
import ast
import glob
import json
import uuid
import numpy
import logging
import time
import subprocess
import threading, traceback
import multiprocessing
from collections import OrderedDict
from pprint import pprint,pformat
import parallel_tools as parallel
import script_util
from biokbase.RNASeq import rnaseq_util
import call_hisat2
import call_stringtie
import call_diffExpCallforBallgown
import contig_id_mapping as c_mapping
from biokbase.workspace.client import Workspace
import handler_utils as handler_util
from biokbase.auth import Token
from mpipe import OrderedStage , Pipeline
import multiprocessing as mp
import re
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
import datetime
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from biokbase.RNASeq.Bowtie2SampleSet import Bowtie2SampleSet
from biokbase.RNASeq.Bowtie2Sample import Bowtie2Sample
from biokbase.RNASeq.HiSat2SampleSet import HiSat2SampleSet
from biokbase.RNASeq.HiSat2Sample import HiSat2Sample
from biokbase.RNASeq.StringTieSampleSet import StringTieSampleSet
from biokbase.RNASeq.StringTieSample import StringTieSample
from biokbase.RNASeq.Cuffdiff import Cuffdiff
from biokbase.RNASeq.TophatSampleSet import TophatSampleSet
from biokbase.RNASeq.TophatSample import TophatSample
from biokbase.RNASeq.CufflinksSampleSet import CufflinksSampleSet
from biokbase.RNASeq.CufflinksSample import CufflinksSample

_KBaseRNASeq__DATA_VERSION = "0.2"

class KBaseRNASeqException(BaseException):
	def __init__(self, msg):
		self.msg = msg
	def __str__(self):
		return repr(self.msg)
#END_HEADER

class KBaseRNASeq:
    '''
    Module Name:
    KBaseRNASeq

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbase/KBaseRNASeq"
    GIT_COMMIT_HASH = "13e98933aea51ffac6b401451d8aebb874958614"
    
    #BEGIN_CLASS_HEADER
    __TEMP_DIR = 'temp'
    __PUBLIC_SHOCK_NODE = 'true'
    __ASSEMBLY_GTF_FN = 'assembly_GTF_list.txt'
    __GTF_SUFFIX = '_GTF_Annotation'
    __BOWTIE2_SUFFIX = '_bowtie2_Alignment'
    __TOPHAT_SUFFIX = '_tophat_Alignment'
    __STATS_DIR = 'stats'
    def generic_helper(self, ctx, params):
	  pass

    # we do normal help function call to parameter mapping
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
	 # This is where config variable for deploy.cfg are available
        #pprint(config)
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']
	if 'hs_url' in config:
	      self.__HS_URL = config['hs_url']
        if 'temp_dir' in config:
              self.__TEMP_DIR = config['temp_dir']
	if 'scratch' in config:
	      self.__SCRATCH= config['scratch']
	      #print self.__SCRATCH
        if 'svc_user' in config:
              self.__SVC_USER = config['svc_user']
        if 'svc_pass' in config:
              self.__SVC_PASS = config['svc_pass']
	if 'scripts_dir' in config:
	      self.__SCRIPTS_DIR = config['scripts_dir']
	if 'force_shock_node_2b_public' in config: # expect 'true' or 'false' string
	      self.__PUBLIC_SHOCK_NODE = config['force_shock_node_2b_public']
	self.__CALLBACK_URL = os.environ['SDK_CALLBACK_URL']	
	
	self.__SERVICES = { 'workspace_service_url' : self.__WS_URL,
			    'shock_service_url' : self.__SHOCK_URL,
			    'handle_service_url' : self.__HS_URL, 
			    'callback_url' : self.__CALLBACK_URL }
        # logging
        self.__LOGGER = logging.getLogger('KBaseRNASeq')
        if 'log_level' in config:
              self.__LOGGER.setLevel(config['log_level'])
        else:
              self.__LOGGER.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.__LOGGER.addHandler(streamHandler)
        self.__LOGGER.info("Logger was set")

        #END_CONSTRUCTOR
        pass
    

    def CreateRNASeqSampleSet(self, ctx, params):
        """
        :param params: instance of type "CreateRNASeqSampleSetParams"
           (FUNCTIONS used in the service) -> structure: parameter "ws_id" of
           String, parameter "sampleset_id" of String, parameter
           "sampleset_desc" of String, parameter "domain" of String,
           parameter "platform" of String, parameter "sample_ids" of list of
           String, parameter "condition" of list of String, parameter
           "source" of String, parameter "Library_type" of String, parameter
           "publication_id" of String, parameter "external_source_date" of
           String
        :returns: instance of type "RNASeqSampleSet" (Object to Describe the
           RNASeq SampleSet @optional platform num_replicates source
           publication_Id external_source_date sample_ids @metadata ws
           sampleset_id @metadata ws platform @metadata ws num_samples
           @metadata ws num_replicates @metadata ws length(condition)) ->
           structure: parameter "sampleset_id" of String, parameter
           "sampleset_desc" of String, parameter "domain" of String,
           parameter "platform" of String, parameter "num_samples" of Long,
           parameter "num_replicates" of Long, parameter "sample_ids" of list
           of String, parameter "condition" of list of String, parameter
           "source" of String, parameter "Library_type" of String, parameter
           "publication_Id" of String, parameter "external_source_date" of
           String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CreateRNASeqSampleSet
	
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    ### Create the working dir for the method; change it to a function call
	    out_obj = { k:v for k,v in params.iteritems() if not k in ('ws_id')}  	
	    sample_ids = params["sample_ids"]
	    out_obj['num_samples'] = len(sample_ids)
	    ## Validation to check if the Set contains more than one samples
	    if len(sample_ids) < 2:
		raise ValueError("This methods can only take 2 or more RNASeq Samples. If you have only one read sample, run either 'Align Reads using Tophat/Bowtie2' methods directly for getting alignment")

	    ## Validation to Check if the number of samples is equal to number of condition
	    if len(params["condition"]) != out_obj['num_samples']:
		raise ValueError("Please specify a treatment label for each sample in the RNA-seq SampleSet. Please enter the same label for the replicates in a sample type")
	    ## Validation to Check if the user is loading the same type as specified above
	    if params["Library_type"] == 'PairedEnd' : lib_type = 'KBaseAssembly.PairedEndLibrary'
	    else: lib_type = 'KBaseAssembly.SingleEndLibrary'
	    for i in sample_ids:
	    	s_info = ws_client.get_object_info_new({"objects": [{'name': i, 'workspace': params['ws_id']}]})
                obj_type = s_info[0][2].split('-')[0]
		if obj_type != lib_type:
			raise ValueError("Library_type mentioned : {0}. Please add only {1} typed objects in Reads fields".format(lib_type,lib_type)) 
	
   	    ## Code to Update the Provenance; make it a function later
            provenance = [{}]
            if 'provenance' in ctx:
                provenance = ctx['provenance']
            #add additional info to provenance here, in this case the input data object reference
            provenance[0]['input_ws_objects']=[ params['ws_id']+'/'+sample for sample in sample_ids]
	    
	    #Saving RNASeqSampleSet to Workspace
	    self.__LOGGER.info("Saving {0} object to workspace".format(params['sampleset_id']))
	    res= ws_client.save_objects(
                                {"workspace":params['ws_id'],
                                 "objects": [{
                                                "type":"KBaseRNASeq.RNASeqSampleSet",
                                                "data":out_obj,
                                                "name":out_obj['sampleset_id'],
						"provenance": provenance}]
                                })
            returnVal = out_obj
        except Exception,e:
                raise KBaseRNASeqException("Error Saving the object to workspace {0},{1}".format(out_obj['sampleset_id'],"".join(traceback.format_exc())))

        #END CreateRNASeqSampleSet

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method CreateRNASeqSampleSet return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def BuildBowtie2Index(self, ctx, params):
        """
        :param params: instance of type "Bowtie2IndexParams" -> structure:
           parameter "ws_id" of String, parameter "reference" of String,
           parameter "output_obj_name" of String
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BuildBowtie2Index
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
                bowtie_dir = os.path.join(self.__SCRATCH ,'tmp') 
	        handler_util.setupWorkingDir(self.__LOGGER,bowtie_dir)
		## Update the provenance
	     	provenance = [{}]
        	if 'provenance' in ctx:
            		provenance = ctx['provenance']
        	# add additional info to provenance here, in this case the input data object reference
        	provenance[0]['input_ws_objects']=[params['ws_id']+'/'+params['reference']]
		
		try:
			ref_id, outfile_ref_name = rnaseq_util.get_fa_from_genome(self.__LOGGER,ws_client,self.__SERVICES,params['ws_id'],bowtie_dir,params['reference'])
                except Exception, e:
			self.__LOGGER.exception("".join(traceback.format_exc()))
                        raise ValueError('Unable to get FASTA for object {}'.format("".join(traceback.format_exc())))
	        ## Run the bowtie_indexing on the  command line
		try:
	    		if outfile_ref_name:
				bowtie_index_cmd = "{0} {1}".format(outfile_ref_name,params['reference'])
			else:
				bowtie_index_cmd = "{0} {1}".format(params['reference'],params['reference']) 
	    	        self.__LOGGER.info("Executing: bowtie2-build {0}".format(bowtie_index_cmd))  	
			cmdline_output = script_util.runProgram(self.__LOGGER,"bowtie2-build",bowtie_index_cmd,None,bowtie_dir)
			if 'result' in cmdline_output:
				report = cmdline_output['result']
		except Exception,e:
			raise KBaseRNASeqException("Error while running BowtieIndex {0},{1}".format(params['reference'],e))
		
	    ## Zip the Index files
		try:
			script_util.zip_files(self.__LOGGER, bowtie_dir,os.path.join(self.__SCRATCH ,"%s.zip" % params['output_obj_name']))
			out_file_path = os.path.join(self.__SCRATCH,"%s.zip" % params['output_obj_name'])
        	except Exception, e:
			raise KBaseRNASeqException("Failed to compress the index: {0}".format(e))
	    ## Upload the file using handle service
		try:
			bowtie_handle = hs.upload(out_file_path)
		except Exception, e:
			raise KBaseRNASeqException("Failed to upload the Zipped Bowtie2Indexes file: {0}".format(e))
	    	bowtie2index = { "handle" : bowtie_handle ,"size" : os.path.getsize(out_file_path),'genome_id' : ref_id}   

	     ## Save object to workspace
	   	self.__LOGGER.info( "Saving bowtie indexes object to  workspace")
	   	res= ws_client.save_objects(
					{"workspace":params['ws_id'],
					 "objects": [{
					 "type":"KBaseRNASeq.Bowtie2Indexes",
					 "data":bowtie2index,
					 "name":params['output_obj_name'],
					 "provenance" : provenance}
					]})
		info = res[0]
	     ## Create report object:
                reportObj = {
                                'objects_created':[{
                                'ref':str(info[6]) + '/'+str(info[0])+'/'+str(info[4]),
                                'description':'Build Bowtie2 Index'
                                }],
                                'text_message':report
                            }

             # generate a unique name for the Method report
                reportName = 'Build_Bowtie2_Index_'+str(hex(uuid.getnode()))
                report_info = ws_client.save_objects({
                                                'id':info[6],
                                                'objects':[
                                                {
                                                'type':'KBaseReport.Report',
                                                'data':reportObj,
                                                'name':reportName,
                                                'meta':{},
                                                'hidden':1, # important!  make sure the report is hidden
                                                'provenance':provenance
                                                }
                                                ]
                                                })[0]

	    	returnVal = { "report_name" : reportName,"report_ref" : str(report_info[6]) + '/' + str(report_info[0]) + '/' + str(report_info[4]) }
	except Exception, e:
		raise KBaseRNASeqException("Build Bowtie2Index failed: {0}".format("".join(traceback.format_exc())))
	finally:
                handler_util.cleanup(self.__LOGGER,bowtie_dir)
        #END BuildBowtie2Index

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method BuildBowtie2Index return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def Bowtie2Call(self, ctx, params):
        """
        :param params: instance of type "Bowtie2Params" -> structure:
           parameter "ws_id" of String, parameter "sampleset_id" of String,
           parameter "genome_id" of String, parameter "bowtie_index" of
           String, parameter "phred33" of String, parameter "phred64" of
           String, parameter "local" of String, parameter "very-fast" of
           String, parameter "fast" of String, parameter "very-sensitive" of
           String, parameter "sensitive" of String, parameter
           "very-fast-local" of String, parameter "very-sensitive-local" of
           String, parameter "fast-local" of String, parameter
           "fast-sensitive" of String
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN Bowtie2Call
	
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        bowtie2_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,bowtie2_dir) 
        common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
        # Set the Number of threads if specified 

        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

        # Check to Call Bowtie2 in Set mode or Single mode
        wsc = common_params['ws_client']
        readsobj_info = wsc.get_object_info_new({"objects": [{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
        readsobj_type = readsobj_info[0][2].split('-')[0]
        if readsobj_type == 'KBaseRNASeq.RNASeqSampleSet':
                self.__LOGGER.info("Bowtie2 SampleSet Case")
                bw2ss = Bowtie2SampleSet(self.__LOGGER, bowtie2_dir, self.__SERVICES)
                returnVal = bw2ss.run(common_params, params)
        else:
                bw2ss = Bowtie2Sample(self.__LOGGER, bowtie2_dir, self.__SERVICES)
                returnVal = bw2ss.run(common_params,params)
	handler_util.cleanup(self.__LOGGER,bowtie2_dir)
        #END Bowtie2Call

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method Bowtie2Call return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def Hisat2Call(self, ctx, params):
        """
        :param params: instance of type "Hisat2Params" -> structure:
           parameter "ws_id" of String, parameter "sampleset_id" of String,
           parameter "genome_id" of String, parameter "num_threads" of Long,
           parameter "quality_score" of String, parameter "skip" of Long,
           parameter "trim3" of Long, parameter "trim5" of Long, parameter
           "np" of Long, parameter "minins" of Long, parameter "maxins" of
           Long, parameter "orientation" of String, parameter
           "min_intron_length" of Long, parameter "max_intron_length" of
           Long, parameter "no_spliced_alignment" of type "bool" (indicates
           true or false values, false <= 0, true >=1), parameter
           "transcriptome_mapping_only" of type "bool" (indicates true or
           false values, false <= 0, true >=1), parameter "tailor_alignments"
           of String
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN Hisat2Call
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        hisat2_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,hisat2_dir) 
	# Set the common Params
	common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
	# Set the Number of threads if specified 

        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

	# Check to Call HiSat2 in Set mode or Single mode
	wsc = common_params['ws_client']
	readsobj_info = wsc.get_object_info_new({"objects": [{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
        readsobj_type = readsobj_info[0][2].split('-')[0]
	if readsobj_type == 'KBaseRNASeq.RNASeqSampleSet':	
		self.__LOGGER.info("HiSat2 SampleSet Case")
        	hs2ss = HiSat2SampleSet(self.__LOGGER, hisat2_dir, self.__SERVICES)
        	returnVal = hs2ss.run(common_params, params)
	else:
		hs2ss = HiSat2Sample(self.__LOGGER, hisat2_dir, self.__SERVICES)
		returnVal = hs2ss.run(common_params,params)
	#finally:
        handler_util.cleanup(self.__LOGGER,hisat2_dir)
        #END Hisat2Call

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method Hisat2Call return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def TophatCall(self, ctx, params):
        """
        :param params: instance of type "TophatParams" -> structure:
           parameter "ws_id" of String, parameter "read_sample" of String,
           parameter "genome_id" of String, parameter "bowtie2_index" of
           String, parameter "read_mismatches" of Long, parameter
           "read_gap_length" of Long, parameter "read_edit_dist" of Long,
           parameter "min_intron_length" of Long, parameter
           "max_intron_length" of Long, parameter "num_threads" of Long,
           parameter "report_secondary_alignments" of String, parameter
           "no_coverage_search" of String, parameter "library_type" of
           String, parameter "annotation_gtf" of type
           "ws_referenceAnnotation_id" (Id for KBaseRNASeq.GFFAnnotation @id
           ws KBaseRNASeq.GFFAnnotation)
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN TophatCall
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        tophat_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,tophat_dir) 
	# Set the common Params
	common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
	# Set the Number of threads if specified 
        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

	# Check to Call Tophat in Set mode or Single mode
	wsc = common_params['ws_client']
	obj_info = wsc.get_object_info_new({"objects": [{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
        obj_type = obj_info[0][2].split('-')[0]
	if obj_type == 'KBaseRNASeq.RNASeqSampleSet':	
		self.__LOGGER.info("Tophat SampleSet Case")
        	tss = TophatSampleSet(self.__LOGGER, tophat_dir, self.__SERVICES)
        	returnVal = tss.run(common_params, params)
	else:
		self.__LOGGER.info("Tophat Sample Case")
		ts = TophatSample(self.__LOGGER, tophat_dir, self.__SERVICES)
		returnVal = ts.run(common_params,params)
        handler_util.cleanup(self.__LOGGER,tophat_dir)

        #END TophatCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method TophatCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def StringTieCall(self, ctx, params):
        """
        :param params: instance of type "StringTieParams" -> structure:
           parameter "ws_id" of String, parameter "sample_alignment" of
           String, parameter "num-threads" of Long, parameter "label" of
           String, parameter "min_isoform_abundance" of Double, parameter
           "a_juncs" of Long, parameter "min_length" of Long, parameter
           "j_min_reads" of Double, parameter "c_min_read_coverage" of
           Double, parameter "gap_sep_value" of Long, parameter
           "disable_trimming" of type "bool" (indicates true or false values,
           false <= 0, true >=1), parameter "ballgown_mode" of type "bool"
           (indicates true or false values, false <= 0, true >=1), parameter
           "skip_reads_with_no_ref" of type "bool" (indicates true or false
           values, false <= 0, true >=1), parameter "merge" of String
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN StringTieCall
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        stringtie_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,stringtie_dir) 
	# Set the common Params
	common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
	# Set the Number of threads if specified 
        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

	# Check to Call StringTie in Set mode or Single mode
	wsc = common_params['ws_client']
	obj_info = wsc.get_object_info_new({"objects": [{'name': params['alignmentset_id'], 'workspace': params['ws_id']}]})
        obj_type = obj_info[0][2].split('-')[0]
	if obj_type == 'KBaseRNASeq.RNASeqAlignmentSet':	
		self.__LOGGER.info("StringTie AlignmentSet Case")
        	sts = StringTieSampleSet(self.__LOGGER, stringtie_dir, self.__SERVICES)
        	returnVal = sts.run(common_params, params)
	else:
		sts = StringTieSample(self.__LOGGER, stringtie_dir, self.__SERVICES)
		returnVal = sts.run(common_params,params)
        handler_util.cleanup(self.__LOGGER,stringtie_dir)
        #END StringTieCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method StringTieCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def CufflinksCall(self, ctx, params):
        """
        :param params: instance of type "CufflinksParams" -> structure:
           parameter "ws_id" of String, parameter "sample_alignment" of
           String, parameter "num_threads" of Long, parameter
           "min-intron-length" of Long, parameter "max-intron-length" of
           Long, parameter "overhang-tolerance" of Long
        :returns: instance of type "ResultsToReport" (Object for Report type)
           -> structure: parameter "report_name" of String, parameter
           "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CufflinksCall
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        cufflinks_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,cufflinks_dir)
        # Set the common Params
        common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
        # Set the Number of threads if specified 
        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

        # Check to Call Cufflinks in Set mode or Single mode
        wsc = common_params['ws_client']
        obj_info = wsc.get_object_info_new({"objects": [{'name': params['alignmentset_id'], 'workspace': params['ws_id']}]})
        obj_type = obj_info[0][2].split('-')[0]
        if obj_type == 'KBaseRNASeq.RNASeqAlignmentSet':
                self.__LOGGER.info("Cufflinks AlignmentSet Case")
                sts = CufflinksSampleSet(self.__LOGGER, cufflinks_dir, self.__SERVICES)
                returnVal = sts.run(common_params, params)
        else:
		sts = CufflinksSample(self.__LOGGER, cufflinks_dir, self.__SERVICES)
                returnVal = sts.run(common_params,params)
        handler_util.cleanup(self.__LOGGER,cufflinks_dir)
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method CufflinksCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def CuffdiffCall(self, ctx, params):
        """
        :param params: instance of type "CuffdiffParams" -> structure:
           parameter "ws_id" of String, parameter "rnaseq_exp_details" of
           type "RNASeqSampleSet" (Object to Describe the RNASeq SampleSet
           @optional platform num_replicates source publication_Id
           external_source_date sample_ids @metadata ws sampleset_id
           @metadata ws platform @metadata ws num_samples @metadata ws
           num_replicates @metadata ws length(condition)) -> structure:
           parameter "sampleset_id" of String, parameter "sampleset_desc" of
           String, parameter "domain" of String, parameter "platform" of
           String, parameter "num_samples" of Long, parameter
           "num_replicates" of Long, parameter "sample_ids" of list of
           String, parameter "condition" of list of String, parameter
           "source" of String, parameter "Library_type" of String, parameter
           "publication_Id" of String, parameter "external_source_date" of
           String, parameter "output_obj_name" of String, parameter
           "time-series" of String, parameter "library-type" of String,
           parameter "library-norm-method" of String, parameter
           "multi-read-correct" of String, parameter "min-alignment-count" of
           Long, parameter "dispersion-method" of String, parameter
           "no-js-tests" of String, parameter "frag-len-mean" of Long,
           parameter "frag-len-std-dev" of Long, parameter
           "max-mle-iterations" of Long, parameter "compatible-hits-norm" of
           String, parameter "no-length-correction" of String
        :returns: instance of type "RNASeqDifferentialExpression" (Object
           RNASeqDifferentialExpression file structure @optional tool_opts
           tool_version sample_ids comments) -> structure: parameter
           "tool_used" of String, parameter "tool_version" of String,
           parameter "tool_opts" of list of mapping from String to String,
           parameter "file" of type "Handle" (@optional hid file_name type
           url remote_md5 remote_sha1) -> structure: parameter "hid" of type
           "HandleId" (Id for the handle object @id handle), parameter
           "file_name" of String, parameter "id" of String, parameter "type"
           of String, parameter "url" of String, parameter "remote_md5" of
           String, parameter "remote_sha1" of String, parameter "sample_ids"
           of list of String, parameter "condition" of list of String,
           parameter "genome_id" of String, parameter "expressionSet_id" of
           type "ws_expressionSet_id" (Id for expression sample set @id ws
           KBaseRNASeq.RNASeqExpressionSet), parameter "alignmentSet_id" of
           type "ws_alignmentSet_id" (The workspace id for a
           RNASeqAlignmentSet object @id ws KBaseRNASeq.RNASeqAlignmentSet),
           parameter "sampleset_id" of type "ws_Sampleset_id" (Id for
           KBaseRNASeq.RNASeqSampleSet @id ws KBaseRNASeq.RNASeqSampleSet),
           parameter "comments" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CuffdiffCall
		
	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        cuffdiff_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,cuffdiff_dir) 
	# Set the common Params
	common_params = {'ws_client' : Workspace(url=self.__WS_URL, token=ctx['token']),
                         'hs_client' : HandleService(url=self.__HS_URL, token=ctx['token']),
                         'user_token' : ctx['token']
                        }
	# Set the Number of threads if specified 
        if 'num_threads' in params and params['num_threads'] is not None:
            common_params['num_threads'] = params['num_threads']

	cuff = Cuffdiff(self.__LOGGER, cuffdiff_dir, self.__SERVICES)
        returnVal = cuff.run(common_params, params)

	#finally:
        handler_util.cleanup(self.__LOGGER,cuffdiff_dir)
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method CuffdiffCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def DiffExpCallforBallgown(self, ctx, params):
        """
        :param params: instance of type "DifferentialExpParams" -> structure:
           parameter "ws_id" of String, parameter "expressionset_id" of type
           "RNASeqExpressionSet" (Set object for RNASeqExpression objects
           @optional sample_ids condition tool_used tool_version tool_opts
           @metadata ws tool_used @metadata ws tool_version @metadata ws
           alignmentSet_id) -> structure: parameter "tool_used" of String,
           parameter "tool_version" of String, parameter "tool_opts" of
           mapping from String to String, parameter "alignmentSet_id" of type
           "ws_alignmentSet_id" (The workspace id for a RNASeqAlignmentSet
           object @id ws KBaseRNASeq.RNASeqAlignmentSet), parameter
           "sampleset_id" of type "ws_Sampleset_id" (Id for
           KBaseRNASeq.RNASeqSampleSet @id ws KBaseRNASeq.RNASeqSampleSet),
           parameter "genome_id" of String, parameter "sample_ids" of list of
           String, parameter "condition" of list of String, parameter
           "sample_expression_ids" of list of type "ws_expression_sample_id"
           (Id for expression sample @id ws KBaseRNASeq.RNASeqExpression),
           parameter "mapped_expression_objects" of list of mapping from
           String to String, parameter "mapped_expression_ids" of list of
           mapping from String to type "ws_expression_sample_id" (Id for
           expression sample @id ws KBaseRNASeq.RNASeqExpression), parameter
           "output_obj_name" of String, parameter "num_threads" of Long
        :returns: instance of type "RNASeqDifferentialExpression" (Object
           RNASeqDifferentialExpression file structure @optional tool_opts
           tool_version sample_ids comments) -> structure: parameter
           "tool_used" of String, parameter "tool_version" of String,
           parameter "tool_opts" of list of mapping from String to String,
           parameter "file" of type "Handle" (@optional hid file_name type
           url remote_md5 remote_sha1) -> structure: parameter "hid" of type
           "HandleId" (Id for the handle object @id handle), parameter
           "file_name" of String, parameter "id" of String, parameter "type"
           of String, parameter "url" of String, parameter "remote_md5" of
           String, parameter "remote_sha1" of String, parameter "sample_ids"
           of list of String, parameter "condition" of list of String,
           parameter "genome_id" of String, parameter "expressionSet_id" of
           type "ws_expressionSet_id" (Id for expression sample set @id ws
           KBaseRNASeq.RNASeqExpressionSet), parameter "alignmentSet_id" of
           type "ws_alignmentSet_id" (The workspace id for a
           RNASeqAlignmentSet object @id ws KBaseRNASeq.RNASeqAlignmentSet),
           parameter "sampleset_id" of type "ws_Sampleset_id" (Id for
           KBaseRNASeq.RNASeqSampleSet @id ws KBaseRNASeq.RNASeqSampleSet),
           parameter "comments" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN DiffExpCallforBallgown
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        #try:
        if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
        diffexp_dir = os.path.join(self.__SCRATCH,"tmp")
        handler_util.setupWorkingDir(self.__LOGGER,diffexp_dir)
        returnVal = call_diffExpCallforBallgown.runMethod(self.__LOGGER,user_token,ws_client,hs,self.__SERVICES,diffexp_dir,params)
	print returnVal
        #except Exception,e:
        #         self.__LOGGER.exception("".join(traceback.format_exc()))
         #        raise KBaseRNASeqException("Error Running StringTieCall")
        #finally:
        handler_util.cleanup(self.__LOGGER,stringtie_dir)
        #END DiffExpCallforBallgown

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method DiffExpCallforBallgown return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
