#BEGIN_HEADER

import simplejson
import sys
import os
import glob
import json
import logging
import time
import subprocess
from pprint import pprint
import script_util
from biokbase.workspace.client import Workspace
from biokbase.auth import Token
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

_KBaseRNASeq__DATA_VERSION = "0.2"

class KBaseRNASeqException(Exception):
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
    #BEGIN_CLASS_HEADER
    __TEMP_DIR = 'temp'
    __BOWTIE_DIR = 'bowtie'
    __TOPHAT_DIR = 'tophat'
    __CUFFLINKS_DIR = 'cufflinks'
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
	 # This is where config variable for deploy.cfg are available
        pprint(config)
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']
	if 'hs_url' in config:
	      self.__HS_URL = config['hs_url']
        if 'temp_dir' in config:
              self.__TEMP_DIR = config['temp_dir']
        if 'bowtie_dir' in config:
              self.__BOWTIE_DIR = config['bowtie_dir']
        if 'genome_input_fa' in config:
              self.__GENOME_FA = config['genome_input_fa']
        if 'svc_user' in config:
              self.__SVC_USER = config['svc_user']
        if 'svc_pass' in config:
              self.__SVC_PASS = config['svc_pass']
	if 'scripts_dir' in config:
	      self.__SCRIPTS_DIR = config['scripts_dir']
	
	self.__SCRIPT_TYPE = { 'ContigSet_to_fasta' : 'ContigSet_to_fasta.py',
			  	'RNASeqSample_to_fastq' : 'RNASeqSample_to_fastq',
			     } 

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

    def fastqcCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN fastqcCall
        #END fastqcCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method fastqcCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def associateReads(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN associateReads
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	out = dict()
	out['metadata'] = { k:v for k,v in params.iteritems() if not k in ('ws_id', "analysis_id", "genome_id","singleend_sample","pairedend_sample") and v is not None }
	self.__LOGGER.info( "Uploading RNASeqSample {0}".format(out['metadata']['sample_id']))
	if "genome_id" in params and params['genome_id'] is not None:
	    out["metadata"]["genome_id"] = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params["genome_id"]],params["ws_id"],user_token)[0]
	if "analysis_id" in params and params['analysis_id'] is not None:
            g_ref = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['analysis_id']],params['ws_id'],user_token)[0]
            out['analysis_id'] = g_ref
	if 'singleend_sample' in params and params['singleend_sample']  is not None:
	    try:
                s_res= ws_client.get_objects(
                                        [{'name' : params['singleend_sample'],
                                          'workspace' : params['ws_id']}])
               	out['singleend_sample'] = s_res[0]['data']
		print out['singleend_sample']
            except Exception,e:
                raise KBaseRNASeqException("Error Downloading SingleEndlibrary object from the workspace {0},{1}".format(params['singleend_sample'],e))

	if 'pairedend_sample' in params and params['pairedend_sample']  is not None:
 	    try:
		p_res= ws_client.get_objects(
                                         [{'name' : params['pairedend_sample'],
                                           'workspace' : params['ws_id']}])
                out['pairedend_sample'] = p_res[0]['data']

            except Exception,e:
                raise KBaseRNASeqException("Error Downloading PairedEndlibrary object from the workspace {0},{1}".format(params['pairedend_sample'],e))

	try:
        	res= ws_client.save_objects(
                                {"workspace":params['ws_id'],
                                 "objects": [{
                                                "type":"KBaseRNASeq.RNASeqSample",
                                                "data":out,
                                                "name":out['metadata']['sample_id']}]
                                })
	        returnVal = {"workspace": params['ws_id'],"output" : out['metadata']['sample_id'] }


	except Exception ,e:
		raise KBaseRNASeqException("Error Saving the object to workspace {0},{1} ".format(out['metadata']['sample_id'],e))
	

        #END associateReads

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method associateReads return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def SetupRNASeqAnalysis(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN SetupRNASeqAnalysis
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        out_obj = { k:v for k,v in params.iteritems() if not k in ('ws_id','genome_id','annotation_id') and v}
        pprint(out_obj)
        if "num_samples" in out_obj : out_obj["num_samples"] = int(out_obj["num_samples"])
        if "num_replicates" in out_obj : out_obj["num_replicates"] = int(out_obj["num_replicates"])
	if "genome_id" in params and params['genome_id'] is not None: out_obj["genome_id"] = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params["genome_id"]],params["ws_id"],user_token)[0]
	if "annotation_id" in params and params['annotation_id'] is not None: 
	    g_ref = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['annotation_id']],params['ws_id'],user_token)[0]
	    out_obj['annotation_id'] = g_ref
        self.__LOGGER.info( "Uploading RNASeq Analysis object to workspace {0}".format(out_obj['experiment_id']))
	try:
        	res= ws_client.save_objects(
                                {"workspace":params['ws_id'],
                                 "objects": [{
                                                "type":"KBaseRNASeq.RNASeqAnalysis",
                                                "data":out_obj,
                                                "name":out_obj['experiment_id']}]
                                })
		returnVal = {"workspace": params['ws_id'],"output" : out_obj['experiment_id'] }

	except Exception,e:
		raise KBaseRNASeqException("Error Saving the object to workspace {0},{1}".format(out_obj['experiment_id'],e))


        #END SetupRNASeqAnalysis

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method SetupRNASeqAnalysis return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def BuildBowtie2Index(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BuildBowtie2Index
        #END BuildBowtie2Index

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method BuildBowtie2Index return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def Bowtie2Call(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN Bowtie2Call
        #END Bowtie2Call

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method Bowtie2Call return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def TophatCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN TophatCall
	ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
	    ### Make a function to download the workspace object  and prepare dict of genome ,lib_type 

	    self.__LOGGER.info("Downloading RNASeq Sample file")
	    try:
                ret  = ws_client.get_objects(
                                        [{'name' : params['sample_id'],'workspace' : params['ws_id']},
					{ 'name' : params['reference'], 'workspace' : params['ws_id']},
					{ 'name' : params['bowtie_index'], 'workspace' : params['ws_id']},
					{ 'name' : params['annotation_gtf'] , 'workspace' : params['ws_id']}])
            except Exception,e:
		 raise KBaseRNASeqException("Error Downloading objects from the workspace ") 
	    
	    #Download reads from the JSON object
	    genome = params['reference']
	    if 'data' in reads: 
		#if 'metadata' in reads['data']:
		    	#genome = reads['data']['metadata']['ref_genome'] 
		if 'singleend_sample' in reads['data']:
			lib_type = "SingleEnd"
			#cmdstring =
	        elif 'pairedend_sample' in reads['data']:
			lib_type = "PairedEnd"
			#cmdstring =
	    ####Complete download reads
		#raise KBaseRNASeqException("Error Downloading FASTA object from the workspace {0}".format(params['reference']))
               
            #self.__LOGGER.info( "Downloading FASTA from ContigSet")
	
	   # cmdstring = "%s/%s --workspace_service_url %s --workspace_name %s --working_directory %s --output_file_name %s --object_name %s --shock_service_url %s" %( self.__SCRIPTS_DIR,self.__SCRIPT_TYPE['ContigSet_to_fasta'],self.__WS_URL,self.__TEMP_DIR,genome,genome,self.__SHOCK_URL)
            
	    #tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
	    #stdout, stderr = tool_process.communicate()
	    #if stdout is not None and len(stdout) > 0:
            #	self.__LOGGER.info(stdout)
	    #if stderr is not None and len(stderr) > 0:
	    #	self.__LOGGER.error("Unable to Download Fasta from ContigSet: " + stderr)
		

	    #try:
            #    assembly = ws_client.get_objects(
            #                            [{'name' : params['reference'],
            #                            'workspace' : params['ws_id']}])
            #except Exception,e:
            #    raise KBaseRNASeqException("Error Downloading FASTA object from the workspace {0}".format(params['reference']))
	# Move them to the temp tophat folder	 

        # Build Fasta object from the ContigSet reference
		
	# Call tophat script with options
	
        # run samtools bam file created 

        # create Json object for widget

	# Zip tophat folder
	
	# Upload to Shock
	
	# run samtools bam file created 
	
	# create Json object for widget

	# save to Tophat workspace
	except e:
		raise 
        job_id = "no_job_id"	
	     
        #END TophatCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method TophatCall return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def CufflinksCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CufflinksCall
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method CufflinksCall return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def CuffmergeCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CuffmergeCall
        #END CuffmergeCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method CuffmergeCall return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def CuffdiffCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CuffdiffCall
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method CuffdiffCall return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def getAlignmentStats(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN getAlignmentStats
	returnVal =  {	"dataset" : { 
			  "properly_paired" : 100 , 
			  "multiple_alignments": 50 , 
                          "singletons": 50 , 
	                  "alignment_rate": 80,
                          "mapped_reads" : 200,
                          "unmapped_reads"  :  50,
                          "total_reads" : 250,
  	                  "mapped_reads" : {
 	                  	"introns" : 50,
	                  	"exons" : 75,
	                  	"splice_junctions" : 100,
				"intergenic_regions" : 55
			  }
		      } }
        #END getAlignmentStats

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method getAlignmentStats return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def createExpressionHistogram(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN createExpressionHistogram
        #END createExpressionHistogram

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method createExpressionHistogram return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def cummeRbundCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN cummeRbundCall
        #END cummeRbundCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method cummeRbundCall return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def createExpressionSeries(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN createExpressionSeries
        #END createExpressionSeries

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method createExpressionSeries return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def createExpressionMatrix(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN createExpressionMatrix
        #END createExpressionMatrix

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method createExpressionMatrix return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]
