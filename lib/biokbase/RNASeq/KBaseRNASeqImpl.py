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
	
	__SCRIPT_TYPE = { 'ContigSet_to_fasta' : 'ContigSet_to_fasta.py',
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
	self.__LOGGER.info( "Uploading RNASeqSample {0}".format(out_obj['experiment_id']))
	out = dict()
	out['metadata'] = { k:v for k,v in params.iteritems() if not k in ('ws_id', "analysis_id") and v}

	if "analysis_id" in params:
		out['analysis_id'] = params['analysis_id']
	if 'singleend_sample' in params:
	    try:
               	out['singleend_sample'] = ws_client.get_objects(
                                        [{'name' : params['singleend_sample'],
                                          'workspace' : params['ws_id']}])
            except Exception,e:
                raise KBaseRNASeqException("Error Downloading SingleEndlibrary object from the workspace {0}".format(params['singleend_sample']))
	if 'pairedend_sample' in params:
 	    try:
                out['pairedend_sample'] = ws_client.get_objects(
                                         [{'name' : params['pairedend_sample'],
                                           'workspace' : params['ws_id']}])
            except Exception,e:
                raise KBaseRNASeqException("Error Downloading PairedEndlibrary object from the workspace {0}".format(params['pairedend_sample']))

	try:
        	res= ws_client.save_objects(
                                {"workspace":params['ws_id'],
                                 "objects": [{
                                                "type":"KBaseRNASeq.RNASeqSample",
                                                "data":out,
                                                "name":out_obj['output_obj_name']}]
                                })
	        returnVal = {"workspace": params['ws_id'],"output" : out_obj['output_obj_name'] }


	except Exception ,e:
		raise KBaseRNASeqException("Error Saving the object to workspace {0}".format(out_obj['output_obj_name']))
	

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
        out_obj = { k:v for k,v in params.iteritems() if not k in ('ws_id') and v}
        pprint(out_obj)
        if "num_samples" in out_obj : out_obj["num_samples"] = int(out_obj["num_samples"])
        if "num_replicates" in out_obj : out_obj["num_replicates"] = int(out_obj["num_replicates"])
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
		raise KBaseRNASeqException("Error Saving the object to workspace {0},{1} : {2}".format(out_obj['experiment_id'],e.errno,e.strerror))


        #END SetupRNASeqAnalysis

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method SetupRNASeqAnalysis return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def BuildBowtie2Index(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN BuildBowtie2Index
	user_token=ctx['token']
	print "start"
        #svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    self.__LOGGER.info( "Downloading KBaseGenome.ContigSet object from workspace")
            try:
	    	assembly = ws_client.get_objects(
                                	[{'name' : params['reference'],
                                  	'workspace' : params['ws_id']}])
	    except Exception,e:
		raise KBaseRNASeqException("Error Downloading FASTA object from the workspace {0}".format(params['reference']))	

	    ## Check if the bowtie_dir is present; remove files in bowtie_dir if exists ; create a new dir if doesnt exists	
	    bowtie_dir = self.__BOWTIE_DIR
	    if os.path.exists(bowtie_dir):
		files=glob.glob("%s/*" % bowtie_dir)
		for f in files: os.remove(f)
	    if not os.path.exists(bowtie_dir): os.makedirs(bowtie_dir)
	   
	    ## dump fasta object to a file in bowtie_dir
	       
	    ## Run the bowtie_indexing on the  command line

	    ## Zip the Index files 	
	    try:
		script_util.zip_files(self.__LOGGER, bowtie_dir, "%s.zip" % params['output_obj_name'])
	    except Exception, e:
		raise KBaseRNASeqException("Failed to compress the index: %s" %(e))
	    ## Upload the file using handle service
	    try:
		handle = hs.upload("%s.zip" % (params['output_obj_name']))
	    except Exception, e:
		raise KBaseRNASeqException("Failed to upload the index: %s" %(e))
	    ## Prepare handle dict object
	        	
	    ## Prepare bowtie indexes workspace object
	    self.__LOGGER.info( "Preparing the workspace object")

	    ## Save object to workspace
	    self.__LOGGER.info( "Saving bowtie indexes object to  workspace")
	    res= ws_client.save_objects(
					{"workspace":params['ws_id'],
					 "objects": [{
					 "type":"KBaseRNASeq.Bowtie2Indexes",
					 "data":bi,
					 "name":params['output_obj_name']}
					]})
	
	except Exception, e:
		raise 
    
	#if 'handle' in assembly and 'id' in assembly['handle']:
	#	shock_id = assembly['handle']['id']
	#script_name = "/kb/dev_container/modules/KBaseRNASeq/lib/biokbase/RNASeq/download_ContigSet.py"
	#res = subprocess.Popen(["python", script_name , "--workspace_service_url" ,self.__WS_URL,  "--workspace_name" , params['ws_id'] , --working_directory ,self.__TEMP_DIR , --output_file_name , params['output_obj_name'] , "--object_name" , params['reference'] , "--shock_service_url" , self.__SHOCK_URL ],stderr=subprocess.STDOUT,stdout = subprocess.PIPE )
	
	#out,err = res.communicate()

	pprint(params)
	job_id = "no_job_id"
	#pull data from shock

	#script_util.download_file_from_shock(self.__LOGGER,
        #                         shock_service_url = self.__SHOCK_URL, 
        #                         shock_id = shock_id,
        #                         filename = self.__INDEX_ZIP,
        #                         directory = self.__TEMP_DIR,
        #                         token = svc_token)


	#assembly_file = self.__TEMP_DIR + "/" + params["output_obj_name"] 
	#index = subprocess.Popen(["bowtie-build" , assembly_file , assembly_file], stderr=subprocess.STDOUT,stdout = subprocess.PIPE )
	#i_out, ierr = index.communicate()
	
	#upload file to shock
	#out_url= 
	
 	# create a dict for the handle object and wra Bowtie2Indexes

	#output =
	
	# save objec to workspace

	#ws_client.save_objects(
        #    {"workspace":params['ws_id'],
        #    "objects": [{
        #        "type":"KBaseRNASeq.Bowtie2Indexes",
        #        "data": output,
        #        "name":params['out_obj_name']}
        #    ]})
        #job_id = { "job_id" : "no_job" }	
        #END BuildBowtie2Index

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method BuildBowtie2Index return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

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
        # return variables are: job_id
        #BEGIN TophatCall
	ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
	    ### Make a function to download the workspace object  and prepare dict of genome ,lib_type 

	    self.__LOGGER.info("Downloading RNASeq Sample file")
	    try:
                reads = ws_client.get_objects(
                                        [{'name' : params['sample_id'],
                                        'workspace' : params['ws_id']}])
            except Exception,e:
		 raise KBaseRNASeqException("Error Downloading FASTQ object from the workspace {0}".format(params['sample_id'])) 
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
               
            self.__LOGGER.info( "Downloading FASTA from ContigSet")
	
	    cmdstring = "%s/%s --workspace_service_url %s --workspace_name %s --working_directory %s --output_file_name %s --object_name %s --shock_service_url %s" %( self.__SCRIPTS_DIR,self.__SCRIPT_TYPE['ContigSet_to_fasta'],self.__WS_URL,self.__TEMP_DIR,genome,genome,self.__SHOCK_URL)
            
	    tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
	    stdout, stderr = tool_process.communicate()
	    if stdout is not None and len(stdout) > 0:
		self.__LOGGER.info(stdout)
	    if stderr is not None and len(stderr) > 0:
		self.__LOGGER.error("Unable to Download Fasta from ContigSet: " + stderr)
		

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
        if not isinstance(job_id, basestring):
            raise ValueError('Method TophatCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CufflinksCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CufflinksCall
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CufflinksCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffmergeCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffmergeCall
        #END CuffmergeCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffmergeCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffdiffCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffdiffCall
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffdiffCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def getAlignmentStats(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN getAlignmentStats
        #END getAlignmentStats

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method getAlignmentStats return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionHistogram(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionHistogram
        #END createExpressionHistogram

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionHistogram return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def cummeRbundCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN cummeRbundCall
        #END cummeRbundCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method cummeRbundCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionSeries(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionSeries
        #END createExpressionSeries

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionSeries return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionMatrix(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionMatrix
        #END createExpressionMatrix

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionMatrix return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]
