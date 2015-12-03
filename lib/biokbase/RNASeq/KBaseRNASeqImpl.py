#BEGIN_HEADER

import simplejson
import sys
import os
import ast
import glob
import json
import logging
import time
import subprocess
import threading, traceback
from collections import OrderedDict
from pprint import pprint
import script_util
from biokbase.workspace.client import Workspace
import handler_utils as handler_util
from biokbase.auth import Token
from mpipe import OrderedStage , Pipeline
import  multiprocessing as mp

try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

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
    #BEGIN_CLASS_HEADER
    __TEMP_DIR = 'temp'
    __BOWTIE_DIR = 'bowtie'
    __BOWTIE2_DIR = 'bowtie2'
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
			  	'cufflinks' : 'cufflinks',
				'tophat_script' : 'Tophat_pipeline.py'
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

    #def parallel(function):                        
    #	 def apply(values,num_threads):
    #     	_pool = mp.Pool(num_threads)
    #    	_result = pool.map(function, values)
    #     	_pool.close()
    #     	_pool.join()         
    #		return _result    
    #	 return apply   

    #def multiprocess(processes, samples, x, widths):
    # 	 pool = mp.Pool(processes=processes)
    #	 results = [pool.apply_async(parzen_estimation, args=(samples, x, w)) for w in widths]
    #	 results = [p.get() for p in results]
    #	 results.sort() # to sort the results by input window width
    #return results	
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
        #pprint(out_obj)
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
	user_token=ctx['token']
    
        #svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	#hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	        self.__LOGGER.info( "Downloading KBaseGenome.ContigSet object from workspace")
	    ## Check if the bowtie_dir is present; remove files in bowtie_dir if exists ; create a new dir if doesnt exists	
	    	bowtie_dir = self.__BOWTIE_DIR
	    	if os.path.exists(bowtie_dir):
			#files=glob.glob("%s/*" % bowtie_dir)
			#for f in files: os.remove(f)
			handler_util.cleanup(self.__LOGGER,bowtie_dir)
	   	if not os.path.exists(bowtie_dir): os.makedirs(bowtie_dir)
	   
	   ## dump fasta object to a file in bowtie_dir
		try:
	   		dumpfasta= "--workspace_service_url {0} --workspace_name {1} --working_directory {2} --output_file_name {3} --object_name {4} --shock_service_url {5} --token \'{6}\'".format(self.__WS_URL , params['ws_id'],bowtie_dir,params['reference'],params['reference'],self.__SHOCK_URL,user_token)

            		script_util.runProgram(self.__LOGGER,self.__SCRIPT_TYPE['ContigSet_to_fasta'],dumpfasta,self.__SCRIPTS_DIR,os.getcwd())
		except Exception,e:
			raise KBaseRNASeqException("Error Creating  FASTA object from the workspace {0},{1},{2}".format(params['reference'],os.getcwd(),e))
		 
	   
	    ## Run the bowtie_indexing on the  command line
		try:
	    		bowtie_index_cmd = "{0} {1}".format(params['reference'],params['reference']) 
	    		script_util.runProgram(self.__LOGGER,"bowtie2-build",bowtie_index_cmd,None,bowtie_dir)
		except Exception,e:
			raise KBaseRNASeqException("Error while running BowtieIndex {0},{1}".format(params['reference'],e))
		
	    ## Zip the Index files
		try:
			script_util.zip_files(self.__LOGGER, bowtie_dir, "%s.zip" % params['output_obj_name'])
			out_file_path = os.path.join("%s.zip" % params['output_obj_name'])
        	except Exception, e:
			raise KBaseRNASeqException("Failed to compress the index: {0}".format(e))
	    ## Upload the file using handle service
		try:
			bowtie_handle = script_util.create_shock_handle(self.__LOGGER,"%s.zip" % params['output_obj_name'],self.__SHOCK_URL,self.__HS_URL,"Zip",user_token)
		except Exception, e:
			raise KBaseRNASeqException("Failed to upload the index: {0}".format(e))
	    	bowtie2index = { "handle" : bowtie_handle ,"size" : os.path.getsize(out_file_path)}	
	    ## Save object to workspace
	   	self.__LOGGER.info( "Saving bowtie indexes object to  workspace")
	   	res= ws_client.save_objects(
					{"workspace":params['ws_id'],
					 "objects": [{
					 "type":"KBaseRNASeq.Bowtie2Indexes",
					 "data":bowtie2index,
					 "name":params['output_obj_name']}
					]})
	    	returnVal = { "output" : params['output_obj_name'],"workspace" : params['ws_id'] }	
	except Exception, e:
		raise KBaseRNASeqException("Build Bowtie2Index failed: {0}".format(e))

        #END BuildBowtie2Index

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method BuildBowtie2Index return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def Bowtie2Call(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN Bowtie2Call
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
            bowtie2_dir = self.__BOWTIE2_DIR
            if os.path.exists(bowtie2_dir):
            #   files=glob.glob("%s/*" % bowtie2_dir)
            #    for f in files: os.remove(f)
                handler_util.cleanup(self.__LOGGER,bowtie2_dir)
            if not os.path.exists(bowtie2_dir): os.makedirs(bowtie2_dir)

            self.__LOGGER.info("Downloading RNASeq Sample file")
	    try:
                sample ,reference,bowtie_index = ws_client.get_objects(
                                        [{'name' : params['sample_id'],'workspace' : params['ws_id']},
                                        { 'name' : params['reference'], 'workspace' : params['ws_id']},
                                        { 'name' : params['bowtie2_index'], 'workspace' : params['ws_id']}])
            except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Downloading objects from the workspace ")
	    opts_dict = { k:v for k,v in params.iteritems() if not k in ('ws_id','sample_id','reference','bowtie_index','analysis_id','output_obj_name') and v is not None }

	    if 'data' in sample and sample['data'] is not None:
                self.__LOGGER.info("getting here")
                if 'metadata' in sample['data'] and sample['data']['metadata'] is not None:
                        genome = sample['data']['metadata']['genome_id']
                        self.__LOGGER.info(genome)
            if 'singleend_sample' in sample['data'] and sample['data']['singleend_sample'] is not None:
                lib_type = "SingleEnd"
                singleend_sample = sample['data']['singleend_sample']
                sample_shock_id = singleend_sample['handle']['id']
                sample_filename = singleend_sample['handle']['file_name']
                try:
                     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id,filename=singleend_sample['handle']['file_name'], directory=bowtie2_dir,token=user_token)
                except Exception,e:
                        raise Exception( "Unable to download shock file , {0}".format(e))
            if 'pairedend_sample' in sample['data'] and sample['data']['pairedend_sample'] is not None:
                lib_type = "PairedEnd"
                pairedend_sample = sample['data']['pairedend_sample']
                self.__LOGGER.info(lib_type)
                if "handle_1" in pairedend_sample and "id" in pairedend_sample['handle_1']:
                        sample_shock_id1  = pairedend_sample['handle_1']['id']
                if "handle_1" in pairedend_sample and "file_name" in pairedend_sample['handle_1']:
                        filename1 = pairedend_sample['handle']['file_name']
                if sample_shock_id1 is None:
                        raise Exception("Handle1 there was not shock id found.")
                if "handle_2" in pairedend_sample  and "id" in pairedend_sample['handle_2']:
                        sample_shock_id2  = pairedend_sample['handle_2']['id']
                if "handle_2" in pairedend_sample and "file_name" in pairedend_sample['handle_2']:
                        filename2 = pairedend_sample['handle']['file_name']

                if sample_shock_id2 is None:
                        raise Exception("Handle2 there was not shock id found.")
                try:
                        script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id1,filename=filename1,directory=bowtie2_dir, token=user_token)
                        script_util.download_file_from_shock(self.__LOGGER,shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id2,filename=filename2,directory=bowtie2_dir, token=user_token)
                except Exception,e:
                        raise Exception( "Unable to download shock file , {0}".format(e))

            if 'analysis_id' in sample['data'] and sample['data']['analysis_id'] is not None:
		# updata the analysis object with the alignment id
                analysis_id = sample['data']['analysis_id']
                self.__LOGGER.info("RNASeq Sample belongs to the {0}".format(analysis_id))
	    if 'handle' in bowtie_index['data'] and bowtie_index['data']['handle'] is not None:
                b_shock_id = bowtie_index['data']['handle']['id']
                b_filename = bowtie_index['data']['handle']['file_name']
                b_filesize = bowtie_index['data']['size']
            try:
                script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=b_shock_id,filename=b_filename,directory=bowtie2_dir,filesize=b_filesize,token=user_token)
            except Exception,e :
                self.__LOGGER.exception("".join(traceback.format_exc()))
                raise Exception( "Unable to download shock file , {0}".format(e))
	    try:
                script_util.unzip_files(self.__LOGGER,os.path.join(bowtie2_dir,b_filename),bowtie2_dir)
		script_util.move_files(self.__LOGGER,handler_util.get_dir(bowtie2_dir),bowtie2_dir)
            except Exception, e:
                   self.__LOGGER.error("Unzip indexfile error: Please contact help@kbase.us")
                   raise Exception("Unzip indexfile error: Please contact help@kbase.us")
            # Define the bowtie2 options
	    os.makedirs(os.path.join(bowtie2_dir,params['output_obj_name']))
	    output_dir = os.path.join(bowtie2_dir,params['output_obj_name'])
	    out_file = output_dir +"accepted_hits.sam"
	    bowtie2_base =os.path.join(bowtie2_dir,handler_util.get_file_with_suffix(bowtie2_dir,".1.bt2"))
	    if(lib_type == "SingleEnd"):
                sample_file = os.path.join(bowtie2_dir,sample_filename)
                bowtie2_cmd = "-U {0} -x {1} -S {2}".format(sample_file,bowtie2_base,out_file)
            elif(lib_type == "PairedEnd"):
                sample_file1 = os.path.join(bowtie2_dir,filename1)
                sample_file2 = os.path.join(bowtie2_dir,filename2)
                bowtie2_cmd = "-1 {0} -2 {1} -x {2} -S {3}".format(sample_file1,sample_file2,bowtie2_base,out_file)	
	    
            try:
		
                script_util.runProgram(self.__LOGGER,"bowtie2",bowtie2_cmd,None,os.getcwd())
		bam_file = os.path.join(output_dir,"accepted_hits.bam")
		sam_to_bam = "view -bS -o {0} {1}".format(bam_file,out_file)
		script_util.runProgram(self.__LOGGER,"samtools",sam_to_bam,None,os.getcwd())
                #script_util.runProgram(self.__LOGGER,self.__SCRIPT_TYPE['tophat_script'],tophat_cmd,self.__SCRIPTS_DIR,os.getcwd())
            except Exception,e:
                raise KBaseRNASeqException("Error Running the bowtie2 command {0},{1},{2}".format(bowtie2_cmd,bowtie2_dir,e))



        # Zip tophat folder
            try:
                script_util.zip_files(self.__LOGGER, output_dir, "%s.zip" % params['output_obj_name'])
                out_file_path = os.path.join("%s.zip" % params['output_obj_name'])
                #handler_util.cleanup(self.__LOGGER,tophat_dir)
            except Exception, e:
                raise KBaseRNASeqException("Failed to compress the index: {0}".format(e))
            ## Upload the file using handle service
            try:
                bowtie2_handle = script_util.create_shock_handle(self.__LOGGER,"%s.zip" % params['output_obj_name'],self.__SHOCK_URL,self.__HS_URL,"Zip",user_token)
            except Exception, e:
                raise KBaseRNASeqException("Failed to upload the index: {0}".format(e))
            bowtie2_out = { "file" : bowtie2_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "bowtie2" , "aligner_version" : "2.2.6","metadata" :  sample['data']['metadata']}
            #tophat_out = { "file" : tophat_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "tophat" , "aligner_version" : "3.1.0", "aligner_opts" : [ (k,v) for k,v in opts_dict.items()],"metadata" :  sample['data']['metadata']}
            returnVal = bowtie2_out

            ## Save object to workspace
            self.__LOGGER.info( "Saving bowtie2 object to  workspace")
            try:
                res= ws_client.save_objects(
                                        {"workspace":params['ws_id'],
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqSampleAlignment",
                                         "data":bowtie2_out,
                                         "name":params['output_obj_name']}
                                        ]})

            except Exception, e:
                raise KBaseRNASeqException("Failed to upload  the alignment: {0}".format(e))

	except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Running Bowtie2Call")

        #END Bowtie2Call

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method Bowtie2Call return value ' +
                             'returnVal is not type object as required.')
        # return the results
        return [returnVal]

    def TophatCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN TophatCall
	user_token=ctx['token']
	ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    ### Make a function to download the workspace object  and prepare dict of genome ,lib_type 
	    tophat_dir = self.__TOPHAT_DIR
            if os.path.exists(tophat_dir):
	    	handler_util.cleanup(self.__LOGGER,tophat_dir)
            if not os.path.exists(tophat_dir): os.makedirs(tophat_dir)

	    self.__LOGGER.info("Downloading RNASeq Sample file")
	    try:
		#values = [{'name' : params['sample_id'],'workspace' : params['ws_id']},
                #          {'name' : params['reference'], 'workspace' : params['ws_id']},
                #          {'name' : params['bowtie_index'], 'workspace' : params['ws_id']},
                #          {'name' : params['annotation_gtf'] , 'workspace' : params['ws_id']}]
	        	
		#def func_call():
        	#print current_thread()
		#print [[x] for x in values]
		
		#for x in values:
        	#	t1 = threading.Thread(target= ws_client.get_objects,args=([x]))
        	#	t1.start()
        	#	t1.join()
		#pool = mp.Pool(processes=4)
		#results = [pool.map(ws_client.get_objects,([x] for x in values))]
		#print results
		#output = [p.get() for p in results]
		#print(output)	
		#@parallel(ws_client.get_objects(values,2))

            	sample ,reference,bowtie_index,annotation = ws_client.get_objects(
                                        [{'name' : params['sample_id'],'workspace' : params['ws_id']},
					{ 'name' : params['reference'], 'workspace' : params['ws_id']},
					{ 'name' : params['bowtie_index'], 'workspace' : params['ws_id']},
					{ 'name' : params['annotation_gtf'] , 'workspace' : params['ws_id']}])
            except Exception,e:
		 self.__LOGGER.exception("".join(traceback.format_exc()))
		 raise KBaseRNASeqException("Error Downloading objects from the workspace ") 
                     
	    opts_dict = { k:v for k,v in params.iteritems() if not k in ('ws_id','sample_id','reference','bowtie_index','annotation_gtf','analysis_id','output_obj_name') and v is not None }
	    
 
            if 'data' in sample and sample['data'] is not None:
		self.__LOGGER.info("getting here")
		if 'metadata' in sample['data'] and sample['data']['metadata'] is not None:
			genome = sample['data']['metadata']['genome_id']
			self.__LOGGER.info(genome)
	    if 'singleend_sample' in sample['data'] and sample['data']['singleend_sample'] is not None:
		lib_type = "SingleEnd"
		singleend_sample = sample['data']['singleend_sample']
		sample_shock_id = singleend_sample['handle']['id']
		sample_filename = singleend_sample['handle']['file_name']
		try:
               	     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id,filename=singleend_sample['handle']['file_name'], directory=tophat_dir,token=user_token)
        	except Exception,e:
                	raise Exception( "Unable to download shock file , {0}".format(e))
	    if 'pairedend_sample' in sample['data'] and sample['data']['pairedend_sample'] is not None: 
		lib_type = "PairedEnd"
		pairedend_sample = sample['data']['pairedend_sample']
		self.__LOGGER.info(lib_type)
		if "handle_1" in pairedend_sample and "id" in pairedend_sample['handle_1']:
                	sample_shock_id1  = pairedend_sample['handle_1']['id']
        	if "handle_1" in pairedend_sample and "file_name" in pairedend_sample['handle_1']:
                	filename1 = pairedend_sample['handle']['file_name']
        	if sample_shock_id1 is None:
                	raise Exception("Handle1 there was not shock id found.")
        	if "handle_2" in pairedend_sample  and "id" in pairedend_sample['handle_2']:
                	sample_shock_id2  = pairedend_sample['handle_2']['id']
        	if "handle_2" in pairedend_sample and "file_name" in pairedend_sample['handle_2']:
                	filename2 = pairedend_sample['handle']['file_name']

        	if sample_shock_id2 is None:
                	raise Exception("Handle2 there was not shock id found.")
		try:
        		script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id1,filename=filename1,directory=tophat_dir, token=user_token)
        		script_util.download_file_from_shock(self.__LOGGER,shock_service_url=self.__SHOCK_URL, shock_id=sample_shock_id2,filename=filename2,directory=tophat_dir, token=user_token)
                except Exception,e:
                        raise Exception( "Unable to download shock file , {0}".format(e))

	    if 'analysis_id' in sample['data'] and sample['data']['analysis_id'] is not None:
		# updata the analysis object with the alignment id
		analysis_id = sample['data']['analysis_id']
	   	self.__LOGGER.info("RNASeq Sample belongs to the {0}".format(analysis_id)) 

            #self.__LOGGER.info("Tophat ran with the following options {0} ",format(str(opts_dict))) 
	    # Download bowtie_Indexes
	    if 'handle' in bowtie_index['data'] and bowtie_index['data']['handle'] is not None:
		b_shock_id = bowtie_index['data']['handle']['id']
		b_filename = bowtie_index['data']['handle']['file_name']
		b_filesize = bowtie_index['data']['size']
	    try:
		script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=b_shock_id,filename=b_filename,directory=tophat_dir,filesize=b_filesize,token=user_token)
	    except Exception,e :
		self.__LOGGER.exception("".join(traceback.format_exc()))
		raise Exception( "Unable to download shock file , {0}".format(e))

	    try:
		index_path = os.path.join(tophat_dir,b_filename)
                script_util.unzip_files(self.__LOGGER,index_path,tophat_dir)
		script_util.move_files(self.__LOGGER,handler_util.get_dir(tophat_dir),tophat_dir)
		#script_util.move_files(self.__LOGGER,os.path.join(tophat_dir,"kb_g.166828"),tophat_dir)
            except Exception, e:
                   self.__LOGGER.exception("".join(traceback.format_exc()))
                   raise Exception("Unzip indexfile error: Please contact help@kbase.us")
	    
            if 'handle' in annotation['data'] and annotation['data']['handle'] is not None:
                a_shock_id = annotation['data']['handle']['id']
                a_filename = annotation['data']['handle']['file_name']
		a_filesize = annotation['data']['size']
            try:
                script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=a_shock_id,filename=a_filename,directory=tophat_dir,filesize=a_filesize,token=user_token)
            except Exception,e :
		self.__LOGGER.exception("".join(traceback.format_exc()))
                raise Exception( "Unable to download shock file , {0}".format(e))
	    output_dir = os.path.join(tophat_dir,params['output_obj_name'])
	    gtf_file = os.path.join(tophat_dir,a_filename)
            #bowtie_base = os.path.join(tophat_dir,"kb_g.166828.fa")
	    bowtie_base =os.path.join(tophat_dir,handler_util.get_file_with_suffix(tophat_dir,".1.bt2"))
	    #topts_dict = { k:int(v) for (k,v) in opts_dict.items() if k in ('num_threads','read_mismatches','read_gap_length','read_edit_dist','min_intron_length','max_intron_length')}
	    #topts_dict['no_coverage_search'] = true
            #topts_dict['report_secondary_alignments'] = false
	    #print opts_dict
	    ### Build command line 
	    #tophat_cmd = "-input {0} -output {1} -reference {2} -opts_dict {3}  -gtf {4} -prog tophat -base_dir {5} -library_type {6} -mode dry_run".format(sample_file,output_dir,bowtie_base,ast.literal_eval(opts_dict),gtf_file,tophat_dir,lib_type)

	    if(lib_type == "SingleEnd"):
       		sample_file = os.path.join(tophat_dir,sample_filename)
            	tophat_cmd = "-o {0} -G {1} {2} {3}".format(output_dir,gtf_file,bowtie_base,sample_file)
	    elif(lib_type == "PairedEnd"):
		sample_file1 = os.path.join(tophat_dir,filename1)
		sample_file2 = os.path.join(tophat_dir,filename2) 
		tophat_cmd = "-o {0} -G {1} {2} {3} {4}".format(output_dir,gtf_file,bowtie_base,sample_file1,sample_file2)
	    #if('num_threads' in opts_dict ) : tophat_cmd = (' -p '+opts_dict['num_threads'])
	    #if('max_intron_length' in opts_dict ) : tophat_cmd += (' -I '+opts_dict['max_intron_length'])
	    #if('min_intron_length' in opts_dict ) : tophat_cmd += (' -i '+opts_dict['min_intron_length'])
	    #if('read_edit_dist' in opts_dict ) : tophat_cmd += (' --read-edit-dist '+opts_dict['read_edit_dist'])
	    #if('read_gap_length' in opts_dict ) : tophat_cmd += (' --read-gap-length '+opts_dict['read_gap_length'])
	    #if('read_mismatches' in opts_dict) : tophat_cmd += (' -N '+opts_dict['read_mismatches'])
	    #if('library_type' in opts_dict) : tophat_cmd += (' --library-type ' + opts_dict['library_type'])
	    #if('report_secondary_alignments' in opts_dict and int(opts_dict['report_secondary_alignments']) == 1) : tophat_cmd += ' --report-secondary-alignments'
	    if('no_coverage_search' in opts_dict and int(opts_dict['no_coverage_search']) == 1): tophat_cmd += ' --no-coverage-search'
#	    if(lib_type == "SingleEnd"):
#                sample_file = os.path.join(tophat_dir,sample_filename)
#                tophat_cmd += "-o {0} -G {1} {2} {3}".format(output_dir,gtf_file,bowtie_base,sample_file)
#            elif(lib_type == "PairedEnd"):
#                sample_file1 = os.path.join(tophat_dir,filename1)
#                sample_file2 = os.path.join(tophat_dir,filename2)
#                tophat_cmd += "-o {0} -G {1} {2} {3} {4}".format(output_dir,gtf_file,bowtie_base,sample_file1,sample_file2)

	    print tophat_cmd
	    try:  
            	script_util.runProgram(self.__LOGGER,"tophat",tophat_cmd,None,os.getcwd())
            	#script_util.runProgram(self.__LOGGER,self.__SCRIPT_TYPE['tophat_script'],tophat_cmd,self.__SCRIPTS_DIR,os.getcwd())
            except Exception,e:
                raise KBaseRNASeqException("Error Running the tophat command {0},{1},{2}".format(tophat_cmd,tophat_dir,e))



	# Zip tophat folder
            try:
                script_util.zip_files(self.__LOGGER, output_dir, "%s.zip" % params['output_obj_name'])
                out_file_path = os.path.join("%s.zip" % params['output_obj_name'])
		#handler_util.cleanup(self.__LOGGER,tophat_dir)
            except Exception, e:
                raise KBaseRNASeqException("Failed to compress the index: {0}".format(e))
            ## Upload the file using handle service
            try:
                tophat_handle = script_util.create_shock_handle(self.__LOGGER,"%s.zip" % params['output_obj_name'],self.__SHOCK_URL,self.__HS_URL,"Zip",user_token)
            except Exception, e:
                raise KBaseRNASeqException("Failed to upload the index: {0}".format(e))
            tophat_out = { "file" : tophat_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "tophat" , "aligner_version" : "3.1.0","metadata" :  sample['data']['metadata']}
            #tophat_out = { "file" : tophat_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "tophat" , "aligner_version" : "3.1.0", "aligner_opts" : [ (k,v) for k,v in opts_dict.items()],"metadata" :  sample['data']['metadata']}
            returnVal = tophat_out
	     
	    ## Save object to workspace
            self.__LOGGER.info( "Saving Tophat object to  workspace")
	    try:
            	res= ws_client.save_objects(
                                        {"workspace":params['ws_id'],
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqSampleAlignment",
                                         "data":tophat_out,
                                         "name":params['output_obj_name']}
                                        ]})
	       
	    except Exception, e:
                raise KBaseRNASeqException("Failed to upload  the alignment: {0}".format(e))

	except Exception,e:
            raise KBaseRNASeqException("Error Running Tophatcall {0}".format("".join(traceback.format_exc())))
	finally:
		handler_util.cleanup(self.__LOGGER,tophat_dir)
#	
	     
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
	user_token=ctx['token']
        self.__LOGGER.info("Started CufflinksCall")
        
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
            cufflinks_dir = self.__CUFFLINKS_DIR
            if os.path.exists(cufflinks_dir):
            #   files=glob.glob("%s/*" % tophat_dir)
            #    for f in files: os.remove(f)
                handler_util.cleanup(self.__LOGGER,cufflinks_dir)
            if not os.path.exists(cufflinks_dir): os.makedirs(cufflinks_dir)

            self.__LOGGER.info("Downloading Alignment Sample file")
	    try:
                sample,annotation_gtf = ws_client.get_objects(
                                        [{'name' : params['alignment_sample_id'],'workspace' : params['ws_id']},
                                         {'name' : params['annotation_gtf'], 'workspace' : params['ws_id']}])
            except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Downloading objects from the workspace ")
	    opts_dict = { k:v for k,v in params.iteritems() if not k in ('ws_id','alignment_sample_id','annotation_gtf','num_threads','min-intron-length','max-intron-length','overhang-tolerance','output_obj_name') and v is not None }

            ## Downloading data from shock
	    if 'data' in sample and sample['data'] is not None:
                self.__LOGGER.info("Downloading alignment sample")
                try:
                     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=sample['data']['file']['id'],filename=sample['data']['file']['file_name'], directory=cufflinks_dir,token=user_token)
                except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(e))
	        try:
                    script_util.unzip_files(self.__LOGGER,os.path.join(cufflinks_dir,sample['data']['file']['file_name']),cufflinks_dir)
                except Exception, e:
                       self.__LOGGER.error("Unzip indexfile error: Please contact help@kbase.us")
                       raise Exception("Unzip indexfile error: Please contact help@kbase.us")
            else:
                raise KBaseRNASeqException("No data was included in the referenced sample id");
#        typedef structure {
#                Handle handle;
#                int size;       
#                ws_genome_id genome_id;
#                string genome_scientific_name;
#        } ReferenceAnnotation;

	    if 'data' in annotation_gtf and annotation_gtf['data'] is not None:
                self.__LOGGER.info("Downloading ReferenceAnnotation")
                try:
                     agtf_fn = annotation_gtf['data']['handle']['file_name']
                     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=annotation_gtf['data']['handle']['id'],filename=agtf_fn, directory=".",token=user_token)
                except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(e))
            else:
                raise KBaseRNASeqException("No data was included in the referenced ReferenceAnnotation");

            ##  now ready to call
            try:
                command_list= ['cufflinks', '-G', cufflinks_dir, agtf_fn, "{0}/accepted_hits.bam".format(cufflinks_dir)]
                if 'num_threads' in params and params['num_threads'] is not None:
                     command_list.append('-p')
                     command_list.append(params['num_threads'])
                for arg in ['min-intron-length','max-intron-length','overhang-tolerance']:
                    if arg in params and params[arg] is not None:
                         command_list.append('--{0}'.format(arg))
                         command_list.append(params[arg])

                self.__LOGGER.info("Executing {0}".format(" ".join(command_list)))
            
                task = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                lines_iterator = iter(task.stdout.readline, b"")
                for line in lines_iterator:
                    self.callback(line)
  
                sub_stdout, sub_stderr = task.communicate()
  
                task_output = dict()
                task_output["stdout"] = sub_stdout
                task_output["stderr"] = sub_stderr
                
                if task.returncode != 0:
                    self.__LOGGER.info(sub_stdout)
                    self.__LOGGER.info(sub_stderr)
                    raise Exception(task_output["stdout"], task_output["stderr"])

            except Exception,e:
                raise KBaseRNASeqException("Error executing cufflinks {0},{1},{2}".format(" ".join(command_list),os.getcwd(),e))
            
            ##  compress and upload to shock
            try:
                self.__LOGGER.info("Ziping output")

                script_util.zip_files(self.__LOGGER,cufflinks_dir, "{0}.zip".format(params['output_obj_name']))
                handle = hs.upload("{0}.zip".format(params['output_obj_name']))
            except Exception,e:
                raise KBaseRNASeqException("Error executing cufflinks {0},{1},{2}".format(cl_args,os.getcwd(),e))
		

	    ## Save object to workspace
	    try:
                self.__LOGGER.info("Saving Cufflinks object to workspace")
                print annotation_gtf['data']['genome_id']
                es_obj = { 'id' : '1234',
                           'source_id' : 'source_id',
                           'type' : 'RNA-Seq',
                           'numerical_interpretation' : 'do_not_know',
                           'external_source_date' : 'external_source_date',
                           'expression_levels' : {},
                           #'genome_id' : 'kb.g.3472',
                           'genome_id' : annotation_gtf['data']['genome_id'],
                           'data_source' : 'data_source',
                           'shock_url' : "{0}/node/{1}".format(handle['url'],handle['id'])
                }

            	res= ws_client.save_objects(
                                        {"workspace":params['ws_id'],
                                         "objects": [{
                                         "type":"KBaseExpression.ExpressionSample",
                                         "data":es_obj,
                                         "name":params['output_obj_name']}
                                        ]})
	    except Exception, e:
                raise KBaseRNASeqException("Failed to upload the ExpressionSample: {0}".format(e))
            returnVal = es_obj
	except KBaseRNASeqException,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise
	except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Running Bowtie2Call")

        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, object):
            raise ValueError('Method CufflinksCall return value ' +
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

	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)        
	stats_data =  { 
    			"alignment_id": params['alignment_sample_id'], 
    			"alignment_rate": 80, 
    			"mapped_sections": {
        		"exons": 50, 
        		"intergenic_regions": 50, 
        		"introns": 50, 
        		"splice_junctions": 50
    			}, 
    			"multiple_alignments": 50, 
    			"properly_paired": 100, 
    			"singletons": 50, 
    			"total_reads": 250, 
    			"unmapped_reads": 50
			}
	
	## Save object to workspace
        self.__LOGGER.info( "Saving Alignment Statistics to the Workspace")
        try:
		res= ws_client.save_objects(
                                        {"workspace":params['ws_id'],
                                         "objects": [{
                                         "type":"KBaseRNASeq.AlignmentStatsResults",
                                         "data": stats_data,
                                         "name":params['output_obj_name']}
                                        ]})
                returnVal = stats_data
        except Exception, e:
                raise KBaseRNASeqException("get Alignment Statistics failed: {0}".format(e))


        #END getAlignmentStats

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method getAlignmentStats return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def createExpressionHistogram(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN createExpressionHistogram
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	
	try:
        	obj = ws_client.get_objects([{'name' : params['expression_sample'],'workspace' : params['ws_id'] }])[0]
        #return {"output" : str(status), "error": json_error}
    	except Exception as e:
        	raise FileNotFound("File Not Found: {}".format(e))
    	if 'expression_levels' in obj['data']:
        	hdict = obj['data']['expression_levels']
        	tot_genes =  len(hdict)
        	lmin = round(min([v for k,v in hdict.items()]))
        	lmax = round(max([v for k,v in hdict.items()]))
        	hist_dt = script_util.histogram(hdict.values(),lmin,lmax,int(params['number_of_bins']))
        	title = "Histogram  - " + params['expression_sample']
        	hist_json = {"title" :  title , "x_label" : "Gene Expression Level (FPKM)", "y_label" : "Number of Genes", "data" : hist_dt}
        	sorted_dt = OrderedDict({ "id" : "", "name" : "","row_ids" : [] ,"column_ids" : [] ,"row_labels" : [] ,"column_labels" : [] , "data" : [] })
        	sorted_dt["row_ids"] = [hist_json["x_label"]]
        	sorted_dt["column_ids"] = [hist_json["y_label"]]
        	sorted_dt['row_labels'] = [hist_json["x_label"]]
        	sorted_dt["column_labels"] =  [hist_json["y_label"]]
        	sorted_dt["data"] = [[float(i) for i in hist_json["data"]["x_axis"]],[float(j) for j in hist_json["data"]["y_axis"]]]
    		#sorted_dt["id"] = "kb|histogramdatatable."+str(idc.allocate_id_range("kb|histogramdatatable",1))
        	sorted_dt["id"] = params['output_obj_name']
        	sorted_dt["name"] = hist_json["title"]
        	res = ws_client.save_objects({"workspace": params['ws_id'],
                                  "objects": [{
                                                "type":"MAK.FloatDataTable",
                                                "data": sorted_dt,
                                                "name" : params['output_obj_name']}
                                            ]

                                 })
	returnVal = { "workspace" : params['ws_id']  , "output" :  params['output_obj_name'] }
        #END createExpressionHistogram

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method createExpressionHistogram return value ' +
                             'returnVal is not type dict as required.')
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
