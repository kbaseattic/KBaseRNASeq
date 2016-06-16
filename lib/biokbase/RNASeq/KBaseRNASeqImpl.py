#BEGIN_HEADER
import simplejson
import sys
import shutil
import os
import ast
import glob
import json
import uuid
import logging
import time
import subprocess
import threading, traceback
import multiprocessing
from collections import OrderedDict
from pprint import pprint,pformat
import parallel_tools as parallel
import script_util
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
import datetime
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
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
    

def parallelize(function,num_processors):
       def temp(_):
         def apply(args):
             final_result = []
             #num_processors=  mp.cpu_count()
             print "Number of processors running parallel jobs {0} ".format(num_processors)
             pool = mp.Pool(num_processors)
             result = pool.map_async(function, args)
             pool.close()
             pool.join()
             return result.get()
         return apply
       return temp

## Helper Function for Parallel call 
def CallBowtie2_helper(x):
    logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token = x
    return parallel._CallBowtie2(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token)

## Helper Function for Parallel call 
def CallTophat_helper(x):
	logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,gtf_file,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token = x
	return parallel._CallTophat(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,gtf_file,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token)

## Helper Function for Parallel call 
def CallCufflinks_helper(x):
	logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token = x
	return  parallel._CallCufflinks(logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token) 
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
    GIT_URL = "https://github.com/sjyoo/KBaseRNASeq"
    GIT_COMMIT_HASH = "9951f0d9b584861c8e828cdfd1c780a063184f9a"
    
    #BEGIN_CLASS_HEADER
    __TEMP_DIR = 'temp'
    __PUBLIC_SHOCK_NODE = 'true'
    __ASSEMBLY_GTF_FN = 'assembly_GTF_list.txt'
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
        if 'svc_user' in config:
              self.__SVC_USER = config['svc_user']
        if 'svc_pass' in config:
              self.__SVC_PASS = config['svc_pass']
	if 'scripts_dir' in config:
	      self.__SCRIPTS_DIR = config['scripts_dir']
	if 'force_shock_node_2b_public' in config: # expect 'true' or 'false' string
	      self.__PUBLIC_SHOCK_NODE = config['force_shock_node_2b_public']
	
	self.__SCRIPT_TYPE = { 'ContigSet_to_fasta' : 'ContigSet_to_fasta.py',
			     } 

	self.__SERVICES = { 'workspace_service_url' : self.__WS_URL,
			    'shock_service_url' : self.__SHOCK_URL,
			    'handle_service_url' : self.__HS_URL }
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
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CreateRNASeqSampleSet
	
	user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    ### Create the working dir for the method; change it to a function call
	    out_obj = { k:v for k,v in params.iteritems() if not k in ('ws_id', 'se_sample_ids', 'pe_sample_ids')}  	
	    sample_ids = params["sample_ids"]
	    out_obj['num_samples'] = len(sample_ids)
	    ## Validation to check if the Set contains more than one samples
	    if len(sample_ids) < 2:
		raise ValueError("This methods can only take 2 or more RNASeq Samples. If you have only one read sample, run either 'Align Reads using Tophat/Bowtie2' methods directly for getting alignment")

	    ## Validation to Check if the number of samples is equal to number of condition
	    if len(params["condition"]) != out_obj['num_samples']:
		raise ValueError("Please specify a treatment label for each sample in the RNA-seq SampleSet. Please enter the same label for the replicates in a sample type")
	    
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
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN BuildBowtie2Index
	user_token=ctx['token']
   	#pprint(params) 
        #svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
	hs = HandleService(url=self.__HS_URL, token=user_token)
	try:
	    	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
                bowtie_dir = os.path.join(self.__SCRATCH ,'/tmp') 
	        handler_util.setupWorkingDir(self.__LOGGER,bowtie_dir)
		## Update the provenance
	     	provenance = [{}]
        	if 'provenance' in ctx:
            		provenance = ctx['provenance']
        	# add additional info to provenance here, in this case the input data object reference
        	provenance[0]['input_ws_objects']=[params['ws_id']+'/'+params['reference']]
		
		ref_info = ws_client.get_object_info_new({"objects": [{'name': params['reference'], 'workspace': params['ws_id']}]})[0]
		genome_id = str(ref_info[6]) + '/' + str(ref_info[0]) + '/' + str(ref_info[4])
                self.__LOGGER.info( "Generating FASTA from Genome Annotation")
                outfile_ref_name = params['reference']+".fasta"
                try:
                    	output_file = script_util.generate_fasta(self.__LOGGER,self.__SERVICES,user_token,params['ws_id'],bowtie_dir,params['reference'])
			self.__LOGGER.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
        		mapping_filename = c_mapping.create_sanitized_contig_ids(output_file)
        		c_mapping.replace_fasta_contig_ids(output_file, mapping_filename, to_modified=True)
        		self.__LOGGER.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))	
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
	    	bowtie2index = { "handle" : bowtie_handle ,"size" : os.path.getsize(out_file_path),'genome_id' : genome_id}   

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

    def GetFeaturesToGTF(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GetFeaturesToGTF
        user_token=ctx['token']
        #pprint(params)
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
                self.__LOGGER.info( "Downloading Genome object from workspace")
            ## Check if the gtf_dir is present; remove files in gtf_dir if exists ; create a new dir if doesnt exists     
		#if os.path.exists(self.__SCRATCH):
                # 	handler_util.cleanup(self.__LOGGER,self.__SCRATCH)
            	if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
		gtf_dir = self.__SCRATCH+'/tmp'
                if os.path.exists(gtf_dir):
                        handler_util.cleanup(self.__LOGGER,gtf_dir)
                if not os.path.exists(gtf_dir): os.makedirs(gtf_dir)
                provenance = [{}]
                if 'provenance' in ctx:
                        provenance = ctx['provenance']
                # add additional info to provenance here, in this case the input data object reference
                provenance[0]['input_ws_objects']=[params['ws_id']+'/'+params['reference']]
		ref_info = ws_client.get_object_info_new({"objects": [{'name': params['reference'], 'workspace': params['ws_id']}]})
		#out_file_path = os.path.join(gtf_dir,params['output_obj_name']+'.gff')
		#output = open(out_file_path,'w')
		obj_type = ref_info[0][2].split('-')[0] 
		if obj_type == 'KBaseGenomeAnnotations.GenomeAnnotation':
			out_file_path = os.path.join(gtf_dir,params['output_obj_name']+'.gff')
			try:
				fasta_file= script_util.generate_fasta(self.__LOGGER,self.__SERVICES,user_token,params['ws_id'],gtf_dir,params['reference'])
                                self.__LOGGER.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
                                mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
                                c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
                                self.__LOGGER.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
				script_util.generate_gff(self.__LOGGER,self.__SERVICES,user_token,params['ws_id'],gtf_dir,params['reference'],out_file_path)
				c_mapping.replace_gff_contig_ids(out_file_path, mapping_filename, to_modified=True) 
			except Exception as e:
				self.__LOGGER.exception("".join(traceback.format_exc()))
				raise ValueError("Generating GFF file from Genome Annotation object Failed :  {}".format("".join(traceback.format_exc())))
		elif obj_type == 'KBaseGenomes.Genome':
		     try:
                	reference = ws_client.get_object_subset(
                                        [{ 'name' : params['reference'], 'workspace' : params['ws_id'],'included': ['features']}])
                	#reference = ws_client.get_objects(
                        #                [{ 'name' : params['reference'], 'workspace' : params['ws_id']}])
			out_file_path = os.path.join(gtf_dir,params['output_obj_name']+'.gtf')
                	output = open(out_file_path,'w')
			ref =reference[0]['data']
        		if "features" in ref:
                  		for f in ref['features']:
                     			if "type" in f and  f['type'] == 'CDS': f_type = f['type']
                     			if "id" in f: f_id =  f['id']
                     			if "location" in f:
                        			for contig_id,f_start,f_strand,f_len  in f['location']:
                                			f_end = script_util.get_end(int(f_start),int(f_len),f_strand)
			        			output.write(contig_id + "\tKBase\t" + f_type + "\t" + str(f_start) + "\t" + str(f_end) + "\t.\t" + f_strand + "\t"+ str(0) + "\ttranscript_id " + f_id + "; gene_id " + f_id + ";\n")
		     except Exception,e:
			raise KBaseRNASeqException("Failed to create Reference Annotation File: {0}".format(e))	
		     finally:
			output.close()
                try:
			#out_file_path = os.path.join(params['output_obj_name']+'.gtf')
                        gtf_handle = hs.upload(out_file_path)

                except Exception, e:
                        raise KBaseRNASeqException("Failed to create Reference Annotation: {0}".format(e))
                gtfhandle = { "handle" : gtf_handle ,"size" : os.path.getsize(out_file_path)}

             ## Save object to workspace
                self.__LOGGER.info( "Saving Reference Annotation object to  workspace")
                res= ws_client.save_objects(
                                        {"workspace":params['ws_id'],
                                         "objects": [{
                                         "type":"KBaseRNASeq.ReferenceAnnotation",
                                         "data":gtfhandle,
                                         "name":params['output_obj_name']}
                                        ]})
                info = res[0]
		report = "Extracting Features from {0}".format(params['reference'])
             ## Create report object:
                reportObj = {
                                'objects_created':[{
                                'ref':str(info[6]) + '/'+str(info[0])+'/'+str(info[4]),
                                'description':'Create Reference Annotation'
                                }],
                                'text_message':report
                            }
                reportName = 'Create_Reference_Annotation_'+str(hex(uuid.getnode()))
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

                print('saved Report: '+pformat(report_info))

		returnVal = { "report_name" : reportName,"report_ref" : str(report_info[6]) + '/' + str(report_info[0]) + '/' + str(report_info[4]) }
        except Exception, e:
                raise KBaseRNASeqException("Create Reference Annotation Failed: {0}".format(e))
        finally:
                handler_util.cleanup(self.__LOGGER,gtf_dir)
		#if os.path.exists(out_file_path): os.remove(out_file_path)
	
        #END GetFeaturesToGTF

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GetFeaturesToGTF return value ' +
                             'returnVal is not type dict as required.')
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
	    if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
            bowtie2_dir = os.path.join(self.__SCRATCH,'/tmp')
            handler_util.setupWorkingDir(self.__LOGGER,bowtie2_dir)
            self.__LOGGER.info("Downloading Bowtie2Indexes and GFF annotation objects")
	    try:
               sample,bowtie_index = ws_client.get_objects(
                                        [{ 'name' : params['sampleset_id'], 'workspace' : params['ws_id']},
                                        { 'name' : params['bowtie_index'], 'workspace' : params['ws_id']}])
            except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException(" Error Downloading objects from the workspace ")
	    ### Get obejct IDs
	    bowtie2_index_info,sampleset_info = ws_client.get_object_info_new({"objects": [{'name': params['bowtie_index'], 'workspace': params['ws_id']},{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
            bowtie2index_id = str(bowtie2_index_info[6]) + '/' + str(bowtie2_index_info[0]) + '/' + str(bowtie2_index_info[4]) 	
            sampleset_id = str(sampleset_info[6]) + '/' + str(sampleset_info[0]) + '/' + str(sampleset_info[4]) 	
#	    #annotation_name =  annotation['data']['handle']['file_name']
	    bw_id = bowtie_index['data']['handle']['id'] 
	    bw_name =  bowtie_index['data']['handle']['file_name']
            genome_id = bowtie_index['data']['genome_id']
	    annotation_gtf = ws_client.get_object_info([{"ref" :genome_id}],includeMetadata=None)[0][1]
	    shared_files={}
	    #shared_files[annotation_name] : annotation_id
	    shared_files[bw_name] = bw_id
	    script_util.download_shock_files(self.__LOGGER,self.__SHOCK_URL,bowtie2_dir,shared_files,user_token)
	    try:
                self.__LOGGER.info("Unzipping Bowtie2 Indices")
                script_util.unzip_files(self.__LOGGER,os.path.join(bowtie2_dir,bw_name),bowtie2_dir)
                mv_dir= handler_util.get_dir(bowtie2_dir)
                if mv_dir is not None:
                        script_util.move_files(self.__LOGGER,mv_dir,bowtie2_dir)
            except Exception, e:
                   self.__LOGGER.error("".join(traceback.format_exc()))
                   raise Exception("Unzip indexfile error: Please contact help@kbase.us")
	    fasta_file =os.path.join(bowtie2_dir,(handler_util.get_file_with_suffix(bowtie2_dir,".fasta")+".fasta"))
	    bowtie2base =os.path.join(bowtie2_dir,handler_util.get_file_with_suffix(bowtie2_dir,".rev.1.bt2"))
	    script_util.create_gtf_annotation(self.__LOGGER,ws_client,hs,self.__SERVICES,params['ws_id'],genome_id,annotation_gtf,fasta_file,bowtie2_dir,user_token)
	    #### Getting Samples info
	    sample_info = ws_client.get_object_info_new({"objects": [{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})[0]
            sample_type = sample_info[2].split('-')[0]
            shared_files = {}
	    # Determine the num_threads provided by the user otherwise default the number of threads to 2
	    if('num_threads' in params and params['num_threads'] is not None): 
			num_threads = int(params['num_threads']) 
	    else:
			num_threads = 2   
	    num_cores = mp.cpu_count()
	    self.__LOGGER.info("Number of available cores : {0}".format(num_cores))
	    b_tasks =[]
            if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
                reads = sample['data']['sample_ids']
		reads_type= sample['data']['Library_type']
                r_label = sample['data']['condition']
		num_samples =  len(reads)
		if num_cores != 1:
			pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
		else:
		   pool_size = 1
		   num_threads = 1
		count = 0 
		self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                for i in reads:
			try:
				label = r_label[count]
				b_tasks.append((None,self.__SERVICES,ws_client,hs,params['ws_id'],reads_type,num_threads,i,label,bowtie2_dir, bowtie2index_id,genome_id,sampleset_id,params,user_token))
				count = count + 1
				### Call multiprocessing of bowtie2 function
                        	#CallBowtie2(self.__LOGGER,self.__services,ws_client,params['ws_id'],params['Library_type'],i,bowtie2_dir,bowtie2_base,options,output_name,user_token)     
			except Exception,e:
				raise
            else:
   		try:
		    pool_size=1
		    num_threads = num_cores
		    self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                    b_tasks.append((self.__LOGGER,self.__SERVICES,ws_client,hs,params['ws_id'],sample_type,num_threads,params['sampleset_id'],'Single-Sample',bowtie2_dir, bowtie2index_id,genome_id,None,params,user_token))
		    #CallBowtie2(self.__LOGGER,self.__services,ws_client,params['ws_id'],params['Library_type'],i,bowtie2_dir,bowtie2_base,options,output_name,user_token)     
                except Exception,e:
                     raise

	    @parallelize(CallBowtie2_helper,pool_size)
	    def run_bowtie2_in_parallel(tasks):
        	pass
	    results=run_bowtie2_in_parallel(b_tasks)
	    reportName = 'Align_Reads_using_Bowtie2_'+str(hex(uuid.getnode()))
	    ### Create AlignmentSet object
	    if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
                alignmentSet_name = params['sampleset_id']+"_AlignmentSet"
                reportObj=script_util.create_RNASeq_AlignmentSet_and_build_report(self.__LOGGER,ws_client,params['ws_id'],reads,sampleset_id,genome_id,bowtie2index_id,results,alignmentSet_name)
 	    else:
		single_read, single_alignment = results
		single_align_obj = ws_client.get_objects(
                                        [{ 'name' : single_alignment, 'workspace' : params['ws_id']}])[0]['data'] 	
      		#returnVal = single_align_obj
                sref = ws_client.get_object_info_new({"objects": [{'name':single_alignment, 'workspace': params['ws_id']}]})[0]
                reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
                                                 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_read)}],
                                                 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_read)}         
            #### Save to report object #######
                #returnVal = single_align_obj   
            report_info = ws_client.save_objects({
                                                'id':sampleset_info[6],
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

	except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Running Bowtie2Call")
	finally:
                 handler_util.cleanup(self.__LOGGER,bowtie2_dir)
		 #if os.path.exists(out_file_path): os.remove(out_file_path)
        #END Bowtie2Call

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method Bowtie2Call return value ' +
                             'returnVal is not type dict as required.')
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
	    if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
	    tophat_dir = self.__SCRATCH +'/tmp'
	    handler_util.setupWorkingDir(self.__LOGGER,tophat_dir)
	    self.__LOGGER.info("Downloading Bowtie2Indexes and GFF annotation objects")
	    try:
                sample ,bowtie_index = ws_client.get_objects(
                                        [{'name' : params['sampleset_id'],'workspace' : params['ws_id']},
					{ 'name' : params['bowtie_index'], 'workspace' : params['ws_id']}])
            except Exception,e:
		 self.__LOGGER.exception("".join(traceback.format_exc()))
		 raise KBaseRNASeqException("Error Downloading objects from the workspace ") 
                     
	    ### Get obejct IDs
            bowtie2_index_info,sampleset_info = ws_client.get_object_info_new({"objects": [{'name': params['bowtie_index'], 'workspace': params['ws_id']},{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})
            bowtie2index_id = str(bowtie2_index_info[6]) + '/' + str(bowtie2_index_info[0]) + '/' + str(bowtie2_index_info[4])  
            sampleset_id = str(sampleset_info[6]) + '/' + str(sampleset_info[0]) + '/' + str(sampleset_info[4]) 
            bw_id = bowtie_index['data']['handle']['id'] 
            bw_name =  bowtie_index['data']['handle']['file_name']
            genome_id = bowtie_index['data']['genome_id']
            annotation_gtf = ws_client.get_object_info([{"ref" :genome_id}],includeMetadata=None)[0][1]
            shared_files={}
            shared_files[bw_name] = bw_id
            script_util.download_shock_files(self.__LOGGER,self.__SHOCK_URL,tophat_dir,shared_files,user_token)
            try:
                self.__LOGGER.info("Unzipping Bowtie2 Indices")
                script_util.unzip_files(self.__LOGGER,os.path.join(tophat_dir,bw_name),tophat_dir)
                mv_dir= handler_util.get_dir(tophat_dir)
                if mv_dir is not None:
                        script_util.move_files(self.__LOGGER,mv_dir,tophat_dir)
            except Exception, e:
                   self.__LOGGER.error("".join(traceback.format_exc()))
                   raise Exception("Unzip indexfile error: Please contact help@kbase.us")
            fasta_file =os.path.join(tophat_dir,(handler_util.get_file_with_suffix(tophat_dir,".fasta")+".fasta"))
            bowtie2base =os.path.join(tophat_dir,handler_util.get_file_with_suffix(tophat_dir,".rev.1.bt2"))

	    ### Check if GTF annotation object exist or skip this step
	    ### Check if the gtf object exists in the workspace
            ### Only run create_gtf_annotation if object doesnt exist
	    gtf_file = script_util.create_gtf_annotation(self.__LOGGER,ws_client,hs,self.__SERVICES,params['ws_id'],genome_id,annotation_gtf,fasta_file,tophat_dir,user_token)
            ### Need this code when if want to download exising gtf object
	    #gtf_obj = annotation_gtf+"_GTF_Annotation"
	    #gtf_obj=ws_client.get_objects([{'name' : params['annotation_gtf'],'workspace' : params['ws_id']])[0]
	    #gtf_id=gtf_obj['data']['handle']['id']
	    #gtf_name=gtf_obj['data']['handle']['name']
	    #try:
            #         script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=gtf_id,filename=gtf_name, directory=tophat_dir,token=user_token)
	    # 	     gtf_file = os.path.join(tophat_dir,gtf_name)
            #except Exception,e:
            #            raise Exception( "Unable to download shock file, {0}".format(gtf_name))  
	    
  	    #### Getting Samples info
            sample_info = ws_client.get_object_info_new({"objects": [{'name': params['sampleset_id'], 'workspace': params['ws_id']}]})[0]
            sample_type = sample_info[2].split('-')[0]
            shared_files = {}
	    # Determine the num_threads provided by the user otherwise default the number of threads to 2
            if('num_threads' in params and params['num_threads'] is not None): 
                        num_threads = int(params['num_threads']) 
            else:
                        num_threads = 2   
            num_cores = mp.cpu_count()
            self.__LOGGER.info("Number of available cores : {0}".format(num_cores))
            b_tasks =[]
            if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
                reads = sample['data']['sample_ids']
                reads_type= sample['data']['Library_type']
                r_label = sample['data']['condition']
                num_samples =  len(reads)
                if num_cores != 1:
                        pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
                else:
                   pool_size = 1
                   num_threads = 1
                count = 0 
                self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                for i in reads:
                        try:
                                label = r_label[count]
                                b_tasks.append((None,self.__SERVICES,ws_client,hs,params['ws_id'],reads_type,num_threads,i,gtf_file,label,tophat_dir, bowtie2index_id,genome_id,sampleset_id,params,user_token))
                                count = count + 1
                                ### Call multiprocessing of bowtie2 function
                                #CallBowtie2(self.__LOGGER,self.__services,ws_client,params['ws_id'],params['Library_type'],i,tophat_dir,bowtie2_base,options,output_name,user_token)     
                        except Exception,e:
                                raise
            else:
                try:
                    pool_size=1
  		    num_threads = num_cores
                    self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                    b_tasks.append((self.__LOGGER,self.__SERVICES,ws_client,hs,params['ws_id'],sample_type,num_threads,params['sampleset_id'],gtf_file,'Single-Sample',tophat_dir, bowtie2index_id,genome_id,None,params,user_token))
                except Exception,e:
                     raise

            @parallelize(CallTophat_helper,pool_size)
            def run_tophat_in_parallel(tasks):
                pass
            results=run_tophat_in_parallel(b_tasks)
	    ### Create report object 
	    #output_objs = {}
	    reportName = 'Align_Reads_using_Tophat_'+str(hex(uuid.getnode()))
            ### Create AlignmentSet object
            if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
		alignmentSet_name = params['sampleset_id']+"_AlignmentSet"
		reportObj=script_util.create_RNASeq_AlignmentSet_and_build_report(self.__LOGGER,ws_client,params['ws_id'],reads,sampleset_id,genome_id,bowtie2index_id,results,alignmentSet_name)
            else:
                single_read, single_alignment = results
                single_align_obj = ws_client.get_objects(
                                        [{ 'name' : single_alignment, 'workspace' : params['ws_id']}])[0]['data'] 
		sref = ws_client.get_object_info_new({"objects": [{'name':single_alignment, 'workspace': params['ws_id']}]})[0]
		reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
			    			 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_read)}],
						 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_read)}	
	    #### Save to report object #######
      		#returnVal = single_align_obj	
	    report_info = ws_client.save_objects({
                                                'id':sampleset_info[6],
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
        except Exception,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise KBaseRNASeqException("Error Running TophatCall")
        finally:
                 handler_util.cleanup(self.__LOGGER,tophat_dir)
                 #if os.path.exists(out_file_path): os.remove(out_file_path)
        
        #END TophatCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method TophatCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def CufflinksCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CufflinksCall
	user_token=ctx['token']
	#pprint(params)
        self.__LOGGER.info("Started CufflinksCall")
        
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
            if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH) 
	    cufflinks_dir = self.__SCRATCH +'/tmp'
	    handler_util.setupWorkingDir(self.__LOGGER,cufflinks_dir)
            self.__LOGGER.info("Downloading Alignment Sample file")
	    try:
                a_sample = ws_client.get_objects(
                                        [{'name' : params['alignmentset_id'],'workspace' : params['ws_id']}])[0]
            except Exception,e:
                self.__LOGGER.exception("".join(traceback.format_exc()))
                raise KBaseRNASeqException("Error Downloading objects from the workspace ")
	    ## Get the Input object type ##
	    a_sample_info = ws_client.get_object_info_new({"objects": [{'name': params['alignmentset_id'], 'workspace': params['ws_id']}]})[0]
            a_sample_type = a_sample_info[2].split('-')[0] 		
            alignmentset_id = str(a_sample_info[6]) + '/' + str(a_sample_info[0]) + '/' + str(a_sample_info[4])		
            #shared_files = {}

	    ### Check if there are existing cufflinks run on the alignments,
	    ### if exists : skip cufflinks run, 
		### get existing ids for set obj,
		### create new run list ,
	    ###  else : run on all alignments
	
	    ### Check if the gtf file exists in the workspace. if exists download the file from that
	    annotation_id = a_sample['data']['genome_id']
	    annotation_name = ws_client.get_object_info([{"ref" :annotation_id}],includeMetadata=None)[0][1]
	    gtf_obj_name = annotation_name+"_GTF_Annotation"
	    gtf_obj= ws_client.get_objects([{'name' : gtf_obj_name,'workspace' : params['ws_id']}])[0]
	    gtf_info = ws_client.get_object_info_new({"objects": [{'name': gtf_obj_name, 'workspace': params['ws_id']}]})[0]
            gtf_annotation_id = str(gtf_info[6]) + '/' + str(gtf_info[0]) + '/' + str(gtf_info[4])
            gtf_id=gtf_obj['data']['handle']['id']
            gtf_name=gtf_obj['data']['handle']['file_name']
            try:
                     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=gtf_id,filename=gtf_name, directory=cufflinks_dir,token=user_token)
                     gtf_file = os.path.join(cufflinks_dir,gtf_name)
            except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(gtf_name))  

	
	    # Determine the num_threads provided by the user otherwise default the number of threads to 2
            if('num_threads' in params and params['num_threads'] is not None): 
                        num_threads = int(params['num_threads']) 
            else:
                        num_threads = 2   
            num_cores = mp.cpu_count()
            self.__LOGGER.info("Number of available cores : {0}".format(num_cores))
            b_tasks =[]
            if a_sample_type == 'KBaseRNASeq.RNASeqAlignmentSet':
                alignment_ids = a_sample['data']['sample_alignments']
		m_alignment_names = a_sample['data']['mapped_rnaseq_alignments']
		sampleset_id =   a_sample['data']['sampleset_id']
		### Get List of Alignments Names
		align_names = []
		for d_align in m_alignment_names:
			for i , j  in d_align.items():
				align_names.append(j)

                m_alignment_ids = a_sample['data']['mapped_alignments_ids']	
                #reads_type= sample['data']['Library_type']
                #r_label = sample['data']['condition']
                num_samples =  len(alignment_ids)
                if num_cores != 1:
                        pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
                else:
                   pool_size = 1
                   num_threads = 1
                self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                for d_align in m_alignment_ids:
		  for s_name,a_id in d_align.items():
                        try:
                                b_tasks.append((None,self.__SERVICES,ws_client,hs,params['ws_id'],num_threads,a_id,gtf_file,cufflinks_dir,annotation_id,gtf_annotation_id,s_name,alignmentset_id,params,user_token))
                        except Exception,e:
                                raise
            else:
                try:
                    pool_size=1
  		    num_threads = num_cores
		    single_a_id = alignmentset_id
                    self.__LOGGER.info(" Number of threads used by each process {0}".format(num_threads)) 
                    b_tasks.append((None,self.__SERVICES,ws_client,hs,params['ws_id'],num_threads,single_a_id,gtf_file,cufflinks_dir,annotation_id,gtf_annotation_id,s_name,None,params,user_token))
                except Exception,e:
                     raise

            @parallelize(CallCufflinks_helper,pool_size)
            def run_cufflinks_in_parallel(tasks):
                pass
            results=run_cufflinks_in_parallel(b_tasks)
	    print results
	    ### Create report object 
	    #output_objs = {}
	    reportName = 'Assemble_Transcripts_Using_Cufflinks_'+str(hex(uuid.getnode()))
            ### Create ExpressionSet object
            if a_sample_type == 'KBaseRNASeq.RNASeqAlignmentSet':
		expressionSet_name = params['alignmentset_id']+"_ExpressionSet"
		reportObj=script_util.create_RNASeq_ExpressionSet_and_build_report(self.__LOGGER,ws_client,params['ws_id'],align_names,alignmentset_id,annotation_id,sampleset_id,results,expressionSet_name)
            else:
                single_alignment, single_expression = results
                single_expr_obj = ws_client.get_objects(
                                        [{ 'name' : single_expression, 'workspace' : params['ws_id']}])[0]['data'] 
		e_ref = ws_client.get_object_info_new({"objects": [{'name':single_expression, 'workspace': params['ws_id']}]})[0]
		reportObj = {'objects_created':[{'ref' :str(eref[6]) + '/' + str(eref[0]) + '/' + str(eref[4]),
			    			 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_alignment)}],
						 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_alignment)}	
	    #### Save to report object #######
      		#returnVal = single_align_obj	
	    report_info = ws_client.save_objects({
                                                'id':a_sample_info[6],
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
	except KBaseRNASeqException,e:
                self.__LOGGER.exception("".join(traceback.format_exc()))
                raise KBaseRNASeqException("Error Running Cufflinks ")
        finally:
                 handler_util.cleanup(self.__LOGGER,cufflinks_dir)
		 #if os.path.exists(out_file_path): os.remove(out_file_path)	
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method CufflinksCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def CuffdiffCall(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN CuffdiffCall
	user_token=ctx['token']
        self.__LOGGER.info("Started CuffdiffCall")
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        try:
            if not os.path.exists(self.__SCRATCH): os.makedirs(self.__SCRATCH)
	    cuffdiff_dir = self.__SCRATCH +'/tmp'
	    handler_util.setupWorkingDir(self.__LOGGER,cuffdiff_dir)
            self.__LOGGER.info("Downloading ExpressionSet  file")
            try:
                e_set = ws_client.get_objects(
                                        [{'name' : params['expressionset_id'],'workspace' : params['ws_id']}])[0]
            except Exception,e:
                self.__LOGGER.exception("".join(traceback.format_exc()))
                raise KBaseRNASeqException("Error Downloading ExpressionSet from the workspace ")
	    ### Getting all the set ids and genome_id
	    e_set_info = ws_client.get_object_info_new({"objects": [{'name' : params['expressionset_id'], 'workspace': params['ws_id']}]})[0] 
	    alignmentset_id = e_set['data']['alignmentSet_id'] 
	    expressionset_id =  str(e_set_info[6]) + '/' + str(e_set_info[0]) + '/' + str(e_set_info[4])
            sampleset_id =  e_set['data']['sampleset_id'] 
	    ### Check if there are existing cufflinks run on the alignments,
            ### if exists : skip cufflinks run, 
                ### get existing ids for set obj,
                ### create new run list ,
            ###  else : run on all alignments
            ### Check if the gtf file exists in the workspace. if exists download the file from that
            annotation_id = e_set['data']['genome_id']
            annotation_name = ws_client.get_object_info([{"ref" :annotation_id}],includeMetadata=None)[0][1]
            gtf_obj_name = annotation_name+"_GTF_Annotation"
            gtf_obj= ws_client.get_objects([{'name' : gtf_obj_name,'workspace' : params['ws_id']}])[0]
            gtf_info = ws_client.get_object_info_new({"objects": [{'name': gtf_obj_name, 'workspace': params['ws_id']}]})[0]
            gtf_annotation_id = str(gtf_info[6]) + '/' + str(gtf_info[0]) + '/' + str(gtf_info[4])
            gtf_id=gtf_obj['data']['handle']['id']
            gtf_name=gtf_obj['data']['handle']['file_name']
            try:
                     script_util.download_file_from_shock(self.__LOGGER, shock_service_url=self.__SHOCK_URL, shock_id=gtf_id,filename=gtf_name, directory=cuffdiff_dir,token=user_token)
                     gtf_file = os.path.join(cuffdiff_dir,gtf_name)
            except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(gtf_name))
            #### Getting the alignments and expression from the alignment set and expression set 
	    m_expr_ids = e_set['data']['mapped_expression_ids']
	    labels = []
	    alignments = []
	    counter = 0
	    assembly_file = os.path.join(cuffdiff_dir,self.__ASSEMBLY_GTF_FN)
            list_file = open(assembly_file,'w')
	    for i in m_expr_ids:
	    	for a_id ,e_id in i.items():
			print a_id  + ":" + e_id
			files = {}
			a_obj,e_obj = ws_client.get_objects(
                                        [{'ref' : a_id},{'ref': e_id}])
			### Get the condition name, replicate_id , shock_id and shock_filename
			condition = a_obj['data']['condition']
			if 'replicate_id' in a_obj['data'] : replicate_id = a_obj['data']['replicate_id']
			files[a_obj['data']['file']['file_name']] = a_obj['data']['file']['id']
			files[e_obj['data']['file']['file_name']] = e_obj['data']['file']['id']
                        if not condition in labels: labels.append(condition)
			else :  counter += 1 #### comment it when replicate_id is available from methods
			print condition
			s_path = os.path.join(cuffdiff_dir,condition+"/"+str(counter)) ### Comment this line when replicate_id is available from the methods
			#s_path = os.path.join(cuffdiff_dir,condition+"/"+replicate_id)
            	        if not os.path.exists(s_path): os.makedirs(s_path)
		        try:
				script_util.download_shock_files(self.__LOGGER,self.__SHOCK_URL,s_path,files,user_token)
                        except Exception,e:
                                raise Exception( "Unable to download shock file, {0}".format(e))
                        try:
                                script_util.unzip_files(self.__LOGGER,os.path.join(s_path,a_obj['data']['file']['file_name']),s_path)
                                script_util.unzip_files(self.__LOGGER,os.path.join(s_path,e_obj['data']['file']['file_name']),s_path)
				e_file_path =  os.path.join(s_path,"transcripts.gtf")
				a_file_path = os.path.join(s_path,"accepted_hits.bam")
				if os.path.exists(a_file_path) : print a_file_path
				if os.path.exists(e_file_path) : list_file.write("{0}\n".format(e_file_path))				
                        except Exception, e:
				self.__LOGGER.exception("".join(traceback.format_exc()))
                                raise Exception("Unzip file error: Please contact help@kbase.us")
	    list_file.close()
	    #output_dir = os.path.join(cuffdiff_dir, params['output_obj_name'])
            for l in labels:
                  rep_files=",".join([ os.path.join(cuffdiff_dir+'/'+l,sub+'/accepted_hits.bam') for sub in os.listdir(os.path.join(cuffdiff_dir,l)) if os.path.isdir(os.path.join(cuffdiff_dir,l+'/'+sub))])
                  alignments.append(rep_files)

            bam_files = " ".join([i for i in alignments])
            t_labels = ",".join(labels)
	   
	    num_threads = multiprocessing.cpu_count()

	    results = parallel.call_cuffmerge_and_cuffdiff(self.__LOGGER,ws_client,hs,params['ws_id'],num_threads,assembly_file,gtf_file,bam_files,t_labels,annotation_id,expressionset_id,alignmentset_id,sampleset_id,params,cuffdiff_dir,user_token)
	    expr_id, cuffdiff_obj = results
	    returnVal = { 'output'  : cuffdiff_obj ,'workspace' : params['ws_id']}
# 
#	   		merged_gtf = analysis['data']['transcriptome_id']
#	    		try:
#                		transcriptome = ws_client.get_objects([{ 'ref' : merged_gtf }])[0]
#            		except Exception,e:
#                 		self.__LOGGER.exception("".join(traceback.format_exc()))
#                 		raise KBaseRNASeqException("Error Downloading merged transcriptome ") 
#	    		t_url = transcriptome['data']['file']['url']
#	    		t_id = transcriptome['data']['file']['id']
#	    		t_name = transcriptome['data']['file']['file_name']
#	    		try:
#                 		script_util.download_file_from_shock(self.__LOGGER, shock_service_url=t_url, shock_id=t_id,filename=t_name, directory=cuffdiff_dir,token=user_token)
#
#            		except Exception,e:
#                 		raise Exception( "Unable to download transcriptome shock file, {0}".format(e))
#            		try:
#                 		script_util.unzip_files(self.__LOGGER,os.path.join(cuffdiff_dir,t_name),cuffdiff_dir)
#            		except Exception, e:
#                 		raise Exception("Unzip transcriptome zip file  error: Please contact help@kbase.us")
#            		gtf_file = os.path.join(cuffdiff_dir,"merged.gtf")
#	   
#            		### Adding advanced options
#	    		num_p = multiprocessing.cpu_count()
#            		#print 'processors count is ' +  str(num_p)
#	    		cuffdiff_command = (' -p '+str(num_p))
#            		#if('num-threads' in params and params['num-threads'] is not None) : cuffdiff_command += (' -p '+str(params['num-threads']))
#	    		if('time-series' in params and params['time-series'] != 0) : cuffdiff_command += (' -T ')
#	    		if('min-alignment-count' in params and params['min-alignment-count'] is not None ) : cuffdiff_command += (' -c '+str(params['min-alignment-count']))
#	    		if('multi-read-correct' in params and params['multi-read-correct'] != 0 ): cuffdiff_command += (' --multi-read-correct ')
#	    		if('library-type' in params and params['library-type'] is not None ) : cuffdiff_command += ( ' --library-type '+params['library-type'])
#	    		if('library-norm-method' in params and params['library-norm-method'] is not None ) : cuffdiff_command += ( ' --library-norm-method '+params['library-norm-method'])
# 
#	    		try:
#                		# TODO: add reference GTF later, seems googledoc command looks wrong
#                		cuffdiff_command += " -o {0} -L {1} -u {2} {3}".format(output_dir,labels,gtf_file,bam_files)
#				self.__LOGGER.info("Executing: cuffdiff {0}".format(cuffdiff_command))
#                		ret = script_util.runProgram(None,"cuffdiff",cuffdiff_command,None,cuffdiff_dir)
#				result = ret["result"]
#                		for line in result.splitlines(False):
#                    			self.__LOGGER.info(line)
#					stderr = ret["stderr"]
#                			prev_value = ''
#                			for line in stderr.splitlines(False):
#                    				if line.startswith('> Processing Locus'):
#                        				words = line.split()
#                        				cur_value = words[len(words) - 1]
#                        				if prev_value != cur_value:
#                            					prev_value = cur_value
#                            					self.__LOGGER.info(line)
#                    				else:
#                        				prev_value = ''
#                        				self.__LOGGER.info(line)
#        		except Exception,e:
#                		raise KBaseRNASeqException("Error executing cuffdiff {0},{1},{2}".format(cuffdiff_command,cuffdiff_dir,e))
#
#            		##  compress and upload to shock
#            		try:
#                		self.__LOGGER.info("Zipping Cuffdiff output")
#				out_file_path = os.path.join(self.__SCRATCH,"{0}.zip".format(params['output_obj_name']))
#                		script_util.zip_files(self.__LOGGER,output_dir,out_file_path)
#                		#handle = hs.upload("{0}.zip".format(params['output_obj_name']))
#            		except Exception,e:
#                		raise KBaseRNASeqException("Error executing cuffdiff {0},{1}".format(os.getcwd(),e))
#            		try:
#				#out_file_path = os.path.join("{0}.zip".format(params['output_obj_name']))
#                		handle = hs.upload(out_file_path)
#            		except Exception, e:
#                		raise KBaseRNASeqException("Failed to upload the Cuffdiff output files: {0}".format(e))
#			
#            		## Save object to workspace
#            		try:
#				self.__LOGGER.info("Saving Cuffdiff object to workspace")
#                		cm_obj = { 'file' : handle,
#                           		   'analysis' : analysis['data']
#                	  		 }
#                		res1= ws_client.save_objects(
#                                        		{"workspace":params['ws_id'],
#                                         		 "objects": [{
#                                         				"type":"KBaseRNASeq.RNASeqCuffdiffdifferentialExpression",
#                                         				"data":cm_obj,
#                                         				"name":params['output_obj_name']}
#                                        			]})	
#		
#                		analysis['data']['cuffdiff_diff_exp_id'] = "{0}/{1}".format(params['ws_id'],params['output_obj_name'])
#				res= ws_client.save_objects(
#                                        		{"workspace":params['ws_id'],
#                                         		 "objects": [{
#                                         		 "type":"KBaseRNASeq.RNASeqAnalysis",
#                                         		 "data":analysis['data'],
#                                         		 "name":params['rnaseq_exp_details']}
#                                        		]})
#            		except Exception, e:
#                		raise KBaseRNASeqException("Failed to upload the KBaseRNASeq.RNASeqCuffdiffdifferentialExpression and KBaseRNASeq.RNASeqAnalysis : {0}".format(e))
#		else:
#                        raise ValueError("Please run the methods 'Align Reads using Tophat/bowtie2' , 'Assemble Transcripts using Cufflinks' for all the RNASeqSamples. Also run the method 'Merge Trancripts using Cuffmerge' before running this step. The RNASeqAnalysis object has either of the missing tags 'alignments' , 'expression_values' , 'transcriptome_id' ");
#
#		returnVal = analysis['data']
        except KBaseRNASeqException,e:
                 self.__LOGGER.exception("".join(traceback.format_exc()))
                 raise
	#finally:
        #         handler_util.cleanup(self.__LOGGER,cuffdiff_dir)
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method CuffdiffCall return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
