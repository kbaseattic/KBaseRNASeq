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
import download_ContigSet as c_download 
from biokbase.auth import Token

_KBaseRNASeq__DATA_VERSION = "0.2"
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
        #pprint(config)
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

    def SetupRNASeqAnalysis(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN SetupRNASeqAnalysis
        #END SetupRNASeqAnalysis

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method SetupRNASeqAnalysis return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def BuildBowtie2Index(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN BuildBowtie2Index
	user_token=ctx['token']
	print "start"
        #svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	
	assembly = ws_client.get_objects(
				[{'name' : params['reference'],
				  'workspace' : params['ws_id']}])
	
	#if 'handle' in assembly and 'id' in assembly['handle']:
	#	shock_id = assembly['handle']['id']
	#script_name = "/kb/dev_container/modules/KBaseRNASeq/lib/biokbase/RNASeq/download_ContigSet.py"
	#res = subprocess.Popen(["python", script_name , "--workspace_service_url" ,self.__WS_URL,  "--workspace_name" , params['ws_id'] , --working_directory ,self.__TEMP_DIR , --output_file_name , params['output_obj_name'] , "--object_name" , params['reference'] , "--shock_service_url" , self.__SHOCK_URL ],stderr=subprocess.STDOUT,stdout = subprocess.PIPE )
	
	#out,err = res.communicate()

	pprint(params)

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
