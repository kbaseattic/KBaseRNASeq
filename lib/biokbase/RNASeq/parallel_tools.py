import simplejson , sys , shutil , os, ast , glob
import multiprocessing as mp
import json ,uuid , logging, time , re , datetime
import subprocess ,threading, traceback
from collections import OrderedDict
from pprint import pprint,pformat
from mpipe import OrderedStage , Pipeline
import script_util
import handler_utils as handler_util
import contig_id_mapping as c_mapping
from biokbase.workspace.client import Workspace
from biokbase.auth import Token
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

def create_logger(log_dir,name):
    """Create a logger

    args: name (str): name of logger

    returns: logger (obj): logging.Logger instance
    """
    logger = logging.getLogger(name)
    fmt = logging.Formatter('%(asctime)s - %(process)d - %(name)s - '
                            ' %(levelname)s -%(message)s')
    hdl = logging.FileHandler(os.path.join(log_dir,name+'.log'))
    hdl.setFormatter(fmt)

    logger.addHandler(hdl)

    return logger

def _CallBowtie2(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token):
	#logger.info("Downloading Read Sample{0}".format(read_sample))
	print "Downloading Read Sample{0}".format(read_sample)
	if not logger:
		logger = create_logger(directory,"run_Bowtie2_"+read_sample)
	
	logger.info("Downloading Read Sample{0}".format(read_sample))
	try:
		r_sample = ws_client.get_objects(
                                        [{ 'name' : read_sample, 'workspace' : ws_id}])[0]
		r_sample_info = ws_client.get_object_info_new({"objects": [{'name': read_sample, 'workspace': ws_id}]})[0]	
		sample_type = r_sample_info[2].split('-')[0]
		output_name = read_sample.split('.')[0]+"_bowtie2_alignment"
		output_dir = os.path.join(directory,output_name)
	        if not os.path.exists(output_dir): os.makedirs(output_dir)
            	out_file = output_dir +"/accepted_hits.sam"
            	bowtie2_base =os.path.join(directory,handler_util.get_file_with_suffix(directory,".rev.1.bt2"))
            	### Adding advanced options to Bowtie2Call
            	bowtie2_cmd = ''
		bowtie2_cmd += ( ' -p {0}'.format(num_threads))
            	if('quality_score' in params and params['quality_score'] is not None): bowtie2_cmd += ( ' --'+params['quality_score'])
            	if('alignment_type' in params and params['alignment_type'] is not None): bowtie2_cmd += ( ' --'+params['alignment_type'] )
            	if('preset_options' in params and params['preset_options'] is not None ) and ('alignment_type' in params and params['alignment_type'] is not None):
                	if (params['alignment_type'] == 'local'):
                        	 bowtie2_cmd += (' --'+params['preset_options']+'-local')
                	else: bowtie2_cmd += (' --'+params['preset_options'] )

		if sample_type  == 'KBaseAssembly.SingleEndLibrary':
			lib_type = 'SingleEnd'
			read_id = r_sample['data']['handle']['id']
			read_name =  r_sample['data']['handle']['file_name']
			try:
                     		script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read_id,filename=read_name, directory=directory,token=token)	
				bowtie2_cmd += " -U {0} -x {1} -S {2}".format(os.path.join(directory,read_name),bowtie2_base,out_file)
                	except Exception,e:
                        	raise Exception( "Unable to download shock file , {0}".format(read_name))
		if sample_type == 'KBaseAssembly.PairedEndLibrary':
			lib_type = 'PairedEnd'
			read1_id = r_sample['data']['handle_1']['id']
			read1_name = r_sample['data']['handle_1']['name']
			read2_id = r_sample['data']['handle_2']['id'] 
			read2_name = r_sample['data']['handle_2']['name']
			try:
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=directory,token=token)
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=directory,token=token)
				bowtie2_cmd += " -1 {0} -2 {1} -x {2} -S {3}".format(os.path.join(directory,read1_name),os.path.join(directory,read2_name),bowtie2_base,out_file)
			except Exception,e:
                        	raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
		try:
                	logger.info("Executing: bowtie2 {0}".format(bowtie2_cmd))
                	cmdline_output = script_util.runProgram(logger,"bowtie2",bowtie2_cmd,None,directory)
                 	#print cmdline_output
		except Exception,e:	
                	raise Exception("Failed to upload the index")
		try:
			#stats_obj_name = params['output_obj_name']+"_"+str(hex(uuid.getnode()))+"_AlignmentStats"
			stats_data = {}
                	stats_data = script_util.extractAlignmentStatsInfo(logger,"bowtie2",ws_client,ws_id,None,cmdline_output['stderr'],None)
                	bam_file = os.path.join(output_dir,"accepted_hits_unsorted.bam")
                	logger.info("Executing: sam_to_bam  {0}".format(bam_file))
			sam_to_bam = "view -bS -o {0} {1}".format(bam_file,out_file)
                	script_util.runProgram(logger,"samtools",sam_to_bam,None,directory)
                	final_bam_prefix = os.path.join(output_dir,"accepted_hits")
                	logger.info("Executing: Sorting bam file  {0}".format(bam_file))
                	sort_bam_cmd  = "sort {0} {1}".format(bam_file,final_bam_prefix)
                	script_util.runProgram(logger,"samtools",sort_bam_cmd,None,directory)
            	except Exception,e:
                	raise Exception("Error Running the bowtie2 command {0},{1} {2}".format(bowtie2_cmd,directory," ".join(traceback.print_exc())))
		
       	 	# Zip tophat folder
            	try:
                	out_file_path = os.path.join(directory,"%s.zip" % output_name)
                	logger.info("Zipping the output files".format(out_file_path))
                	script_util.zip_files(logger, output_dir,out_file_path)
            	except Exception, e:
                	raise Exception("Failed to compress the index: {0}".format(out_file_path))
            	## Upload the file using handle service
            	try:
                	bowtie2_handle = hs.upload(out_file_path)
            	except Exception, e:
                	raise Exception("Failed to upload zipped output file".format(out_file_path))
		#### Replace version with get_version command#####
            	bowtie2_out = { "file" : bowtie2_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "bowtie2" , "aligner_version" : "2.2.6" , 'library_type' : lib_type , 'condition' : condition ,'read_sample_id': read_sample, 'genome_id' : genome_id , 'bowtie2_index' : bowtie2index_id , "alignment_stats" : stats_data }
		
		#logger.info( "Saving bowtie2 object to  workspace")
            	if not sampleset_id is None: bowtie2_out['sampleset_id'] = sampleset_id
		try:
                	res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignment",
                                         "data":bowtie2_out,
                                         "name":output_name}
                                        ]})
                #map_key = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['sample_id']],params["ws_id"],user_token)[0]
                #map_value = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['output_obj_name']],params["ws_id"],user_token)[0]
            	except Exception, e:
                	raise Exception("Failed to save alignment to workspace")
	except Exception, e:
                        raise Exception("Failed to create bowtie2 Alignment {0}".format(" ".join(traceback.print_exc())))
    	return (read_sample,output_name )

def _CallTophat(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,gtf_file,condition,directory,bowtie2index_id,genome_id,sampleset_id,params,token):
	print "Downloading Read Sample{0}".format(read_sample)
	if not logger:
		logger = create_logger(directory,"run_Tophat_"+read_sample)	
	try:
		r_sample = ws_client.get_objects(
                                        [{ 'name' : read_sample, 'workspace' : ws_id}])[0]
		r_sample_info = ws_client.get_object_info_new({"objects": [{'name': read_sample, 'workspace': ws_id}]})[0]	
		sample_type = r_sample_info[2].split('-')[0]
		output_name = read_sample.split('.')[0]+"_tophat_alignment"
		output_dir = os.path.join(directory,output_name)
	        #if not os.path.exists(output_dir): os.makedirs(output_dir)
            	#out_file = output_dir +"/accepted_hits.sam"
            	bowtie2_base =os.path.join(directory,handler_util.get_file_with_suffix(directory,".rev.1.bt2"))

            	### Adding advanced options to tophat command
		tophat_cmd = (' -p '+str(num_threads))
            	if('max_intron_length' in params and params['max_intron_length'] is not None ) : tophat_cmd += (' -I '+str(params['max_intron_length']))
            	if('min_intron_length' in params and params['min_intron_length'] is not None ): tophat_cmd += (' -i '+str(params['min_intron_length']))
            	if('read_edit_dist' in params and params['read_edit_dist'] is not None ) : tophat_cmd += (' --read-edit-dist '+str(params['read_edit_dist']))
            	if('read_gap_length' in params and params['read_gap_length'] is not None) : tophat_cmd += (' --read-gap-length '+str(params['read_gap_length']))
            	if('read_mismatches' in params and params['read_mismatches'] is not None) : tophat_cmd += (' -N '+str(params['read_mismatches']))
            	if('library_type' in params and params['library_type']  is not None ) : tophat_cmd += (' --library-type ' + params['library_type'])
            	if('report_secondary_alignments' in params and int(params['report_secondary_alignments']) == 1) : tophat_cmd += ' --report-secondary-alignments'
            	if('no_coverage_search' in params and int(params['no_coverage_search']) == 1): tophat_cmd += ' --no-coverage-search'
		if sample_type  == 'KBaseAssembly.SingleEndLibrary':
			lib_type = 'SingleEnd'
			read_id = r_sample['data']['handle']['id']
			read_name =  r_sample['data']['handle']['file_name']
			try:
                     		script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read_id,filename=read_name, directory=directory,token=token)	
                		tophat_cmd += ' -o {0} -G {1} {2} {3}'.format(output_dir,gtf_file,bowtie2_base,os.path.join(directory,read_name))
                	except Exception,e:
                        	raise Exception( "Unable to download shock file , {0}".format(read_name))
		if sample_type == 'KBaseAssembly.PairedEndLibrary':
			lib_type = 'PairedEnd'
			read1_id = r_sample['data']['handle_1']['id']
			read1_name = r_sample['data']['handle_1']['name']
			read2_id = r_sample['data']['handle_2']['id'] 
			read2_name = r_sample['data']['handle_2']['name']
			try:
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=directory,token=token)
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=directory,token=token)
                		tophat_cmd += ' -o {0} -G {1} {2} {3} {4}'.format(output_dir,gtf_file,bowtie2_base,os.path.join(directory,read1_name),os.path.join(directory,read2_name))
			except Exception,e:
                        	raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
		try:
                	
                        script_util.runProgram(logger,"tophat",tophat_cmd,None,directory)	
            	except Exception,e:
                	raise Exception("Error Running the tophat command {0},{1} {2}".format(tophat_cmd,directory," ".join(traceback.print_exc())))
	 	try:
                	bam_file = output_dir+"/accepted_hits.bam"
                	align_stats_cmd="flagstat {0}".format(bam_file)
                	stats = script_util.runProgram(logger,"samtools",align_stats_cmd,None,directory)
			#print stats
			stats_data = {}
                	# Pass it to the stats['result']
                	#stats_obj_name = params['output_obj_name']+"_"+str(hex(uuid.getnode()))+"_AlignmentStats"
                	stats_data =script_util.extractAlignmentStatsInfo(logger,"samtools",ws_client,ws_id,None,stats['result'],None)
            	except Exception , e :
                	raise Exception("Failed to create RNASeqAlignmentStats: {0}".format(bam_file))
       	 	# Zip tophat folder
            	try:
                	out_file_path = os.path.join(directory,"%s.zip" % output_name)
                	script_util.zip_files(logger, output_dir,out_file_path)
            	except Exception, e:
                	raise Exception("Failed to compress the index: {0}".format(out_file_path))
            	## Upload the file using handle service
            	try:
                	tophat_handle = hs.upload(out_file_path)
            	except Exception, e:
                	raise Exception("Failed to upload the index")
		#### Replace version with get_version command#####
            	tophat_out = { "file" : tophat_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "tophat" , "aligner_version" : "2.1.1" , 'library_type' : lib_type , 'condition' : condition ,'read_sample_id': read_sample, 'genome_id' : genome_id , 'bowtie2_index' : bowtie2index_id , "alignment_stats" : stats_data}
		#logger.info( "Saving bowtie2 object to  workspace")
            	if not sampleset_id is None: tophat_out['sampleset_id'] = sampleset_id
		try:
                	res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignment",
                                         "data":tophat_out,
                                         "name":output_name}
                                        ]})
                #map_key = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['sample_id']],params["ws_id"],user_token)[0]
                #map_value = script_util.get_obj_info(self.__LOGGER,self.__WS_URL,[params['output_obj_name']],params["ws_id"],user_token)[0]
            	except Exception, e:
                	raise Exception("Failed to save alignment to workspace")
	except Exception, e:
                        raise Exception("Failed to create tophat Alignment {0}".format(" ".join(traceback.print_exc())))
    	return (read_sample,output_name )
