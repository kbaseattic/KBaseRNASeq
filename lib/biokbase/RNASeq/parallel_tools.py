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
	        if not os.path.exists(output_dir): os.mkdir(output_dir)
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
            	if('trim5' in params and params['trim5'] is not None): bowtie2_cmd += ( ' --trim5 '+str(params['trim5']))
            	if('trim3' in params and params['trim3'] is not None): bowtie2_cmd += ( ' --trim3 '+str(params['trim3']))
            	if('np' in params and params['np'] is not None): bowtie2_cmd += ( ' --np '+str(params['np']))
            	if('minins' in params and params['minins'] is not None): bowtie2_cmd += ( ' --minins '+str(params['minins']))
            	if('maxins' in params and params['maxins'] is not None): bowtie2_cmd += ( ' --maxins '+str(params['maxins']))
            	if('orientation' in params and params['orientation'] is not None): bowtie2_cmd += ( ' --'+params['orientation'])
		
		if sample_type  == 'KBaseAssembly.SingleEndLibrary':
			lib_type = 'SingleEnd'
			read_id = r_sample['data']['handle']['id']
			read_name =  r_sample['data']['handle']['file_name']
			try:
                     		script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read_id,filename=read_name, directory=directory,token=token)	
				bowtie2_cmd += " -U {0} -x {1} -S {2}".format(os.path.join(directory,read_name),bowtie2_base,out_file)
                	except Exception,e:
                        	#logger.exception( "Unable to download shock file , {0}".format(read_name))
                        	raise Exception( "Unable to download shock file , {0}".format(read_name))
		if sample_type == 'KBaseAssembly.PairedEndLibrary':
			lib_type = 'PairedEnd'
			read1_id = r_sample['data']['handle_1']['id']
			read1_name = r_sample['data']['handle_1']['file_name']
			read2_id = r_sample['data']['handle_2']['id'] 
			read2_name = r_sample['data']['handle_2']['file_name']
			try:
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=directory,token=token)
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=directory,token=token)
				bowtie2_cmd += " -1 {0} -2 {1} -x {2} -S {3}".format(os.path.join(directory,read1_name),os.path.join(directory,read2_name),bowtie2_base,out_file)
			except Exception,e:
                        	#logger.Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
                        	raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
		try:
                	logger.info("Executing: bowtie2 {0}".format(bowtie2_cmd))
                	cmdline_output = script_util.runProgram(logger,"bowtie2",bowtie2_cmd,None,directory)
                 	#print cmdline_output
		except Exception,e:
                	#logger.exception("Failed to upload the index")
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
                	#logger.exception("Error Running the bowtie2 command {0},{1} {2}".format(bowtie2_cmd,directory," ".join(traceback.print_exc())))
                	raise Exception("Error Running the bowtie2 command {0},{1} {2}".format(bowtie2_cmd,directory," ".join(traceback.print_exc())))
		
       	 	# Zip tophat folder
            	try:
                	out_file_path = os.path.join(directory,"%s.zip" % output_name)
                	logger.info("Zipping the output files".format(out_file_path))
                	script_util.zip_files(logger, output_dir,out_file_path)
            	except Exception, e:
                	#logger.exception("Failed to compress the index: {0}".format(out_file_path))
                	raise Exception("Failed to compress the index: {0}".format(out_file_path))
            	## Upload the file using handle service
            	try:
                	bowtie2_handle = hs.upload(out_file_path)
            	except Exception, e:
                	logger.exception("Failed to upload zipped output file".format(out_file_path))
                	#raise Exception("Failed to upload zipped output file".format(out_file_path))
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
            	except Exception, e:
                	#logger.exception("Failed to save alignment to workspace")
                	raise Exception("Failed to save alignment to workspace")
	except Exception, e:
                        #logger.exception("Failed to create bowtie2 Alignment {0}".format(" ".join(traceback.print_exc())))
                        raise Exception("Failed to create bowtie2 Alignment {0}".format(" ".join(traceback.print_exc())))
	finally:
		if os.path.exists(out_file_path): os.remove(out_file_path)
		if os.path.exists(output_dir): shutil.rmtree(output_dir)
		ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqAlignment",[output_name])	
		if not ret is None:
		    return (read_sample,output_name)
		else:
		    return None

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
            	if('min_anchor_length' in params and params['min_anchor_length'] is not None ): tophat_cmd += (' -a '+str(params['min_anchor_length']))
            	if('read_edit_dist' in params and params['read_edit_dist'] is not None ) : tophat_cmd += (' --read-edit-dist '+str(params['read_edit_dist']))
            	if('read_gap_length' in params and params['read_gap_length'] is not None) : tophat_cmd += (' --read-gap-length '+str(params['read_gap_length']))
            	if('read_mismatches' in params and params['read_mismatches'] is not None) : tophat_cmd += (' -N '+str(params['read_mismatches']))
            	if('library_type' in params and params['library_type']  is not None ) : tophat_cmd += (' --library-type ' + params['library_type'])
            	if('report_secondary_alignments' in params and int(params['report_secondary_alignments']) == 1) : tophat_cmd += ' --report-secondary-alignments'
            	if('no_coverage_search' in params and int(params['no_coverage_search']) == 1): tophat_cmd += ' --no-coverage-search'
            	if('preset_options' in params and params['preset_options'] is not None ): tophat_cmd += ' --'+params['preset_options']
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
			read1_name = r_sample['data']['handle_1']['file_name']
			read2_id = r_sample['data']['handle_2']['id'] 
			read2_name = r_sample['data']['handle_2']['file_name']
			try:
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=directory,token=token)
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=directory,token=token)
                		tophat_cmd += ' -o {0} -G {1} {2} {3} {4}'.format(output_dir,gtf_file,bowtie2_base,os.path.join(directory,read1_name),os.path.join(directory,read2_name))
			except Exception,e:
                        	raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
		try:
                	
                        eres , estrr = script_util.runProgram(logger,"tophat",tophat_cmd,None,directory)
			print eres
			print estrr
            	except Exception,e:
			raise Exception(e)
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
	finally:
		if os.path.exists(out_file_path): os.remove(out_file_path)
		if os.path.exists(output_dir): shutil.rmtree(output_dir)
		ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqAlignment",[output_name])	
		if not ret is None:
		    return (read_sample,output_name)
		else:
		    return None


def _CallCufflinks(logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token):
	print "Downloading Read Sample{0}".format(s_alignment)
	alignment_name = ws_client.get_object_info([{"ref" :s_alignment}],includeMetadata=None)[0][1]
	if not logger:
		logger = create_logger(directory,"run_Cufflinks_"+alignment_name)	
	try:
		alignment = ws_client.get_objects(
                                        [{ 'ref' : s_alignment }])[0]
		#alignment_info = ws_client.get_object_info_new({"objects": [{'name': read_sample, 'workspace': ws_id}]})[0]	
		#sample_type = r_sample_info[2].split('-')[0]
		output_name = alignment_name.split('_alignment')[0]+"_cufflinks_expression"
		output_dir = os.path.join(directory,output_name)
		#Download Alignment from shock
		a_file_id = alignment['data']['file']['id']
		a_filename = alignment['data']['file']['file_name']
		condition = alignment['data']['condition']
		#i_name = alignment_name+"_"+a_filename
		#if replicate_id in alignment['data'] : replicate_id = alignment['data']['replicate_id']
		try:
                     script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=a_file_id,filename=a_filename,directory=directory,token=token)
                except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(i_name))
                try:
		    input_dir = os.path.join(directory,alignment_name)
		    if not os.path.exists(input_dir): os.mkdir(input_dir)
                    script_util.unzip_files(logger,os.path.join(directory,a_filename), input_dir)
                except Exception, e:
                       logger.error("".join(traceback.format_exc()))
                       raise Exception("Unzip alignment files  error: Please contact help@kbase.us")

		input_file = os.path.join(input_dir,"accepted_hits.bam")
		### Adding advanced options to tophat command
		cufflinks_command = (' -p '+str(num_threads))
		if 'max_intron_length' in params and params['max_intron_length'] is not None:
                     cufflinks_command += (' --max-intron-length '+str(params['max_intron_length']))
                if 'min_intron_length' in params and params['min_intron_length'] is not None:
                     cufflinks_command += (' --min-intron-length '+str(params['min_intron_length']))
                if 'overhang_tolerance' in params  and params['overhang_tolerance'] is not None:
                     cufflinks_command += (' --overhang-tolerance '+str(params['overhang_tolerance']))

                cufflinks_command += " -o {0} -G {1} {2}".format(output_dir,gtf_file,input_file)
		logger.info("Executing: cufflinks {0}".format(cufflinks_command))
                ret = script_util.runProgram(None,"cufflinks",cufflinks_command,None,directory)
                result = ret["result"]
                for line in result.splitlines(False):
                    logger.info(line)
                stderr = ret["stderr"]
                prev_value = ''
                for line in stderr.splitlines(False):
                    if line.startswith('> Processing Locus'):
                        words = line.split()
                        cur_value = words[len(words) - 1]
                        if prev_value != cur_value:
				logger.info(line)
				#print line
                    else:
                        prev_value = ''
                        logger.info(line)
                        #print line

        	##Parse output files
        	exp_dict = script_util.parse_FPKMtracking(os.path.join(output_dir,"genes.fpkm_tracking"),'Cufflinks','FPKM')
        	##  compress and upload to shock
        	try:
                	logger.info("Zipping Cufflinks output")
                	out_file_path = os.path.join(directory,"%s.zip" % output_name)
                	script_util.zip_files(logger,output_dir,out_file_path)
        	except Exception,e:
                	logger.exception("".join(traceback.format_exc()))
                	raise Exception("Error executing cufflinks")
        	try:
                	handle = hs.upload(out_file_path)
        	except Exception, e:
                	logger.exception("".join(traceback.format_exc()))
                	raise Exception("Error while zipping the output objects: {0}".format(out_file_path))
        	## Save object to workspace
        	try:
                	logger.info("Saving Cufflinks object to workspace")
                	es_obj = { 'id' : output_name,
                           	'type' : 'RNA-Seq',
                           	'numerical_interpretation' : 'FPKM',
                           	'expression_levels' : exp_dict,
			   	'processing_comments' : "log2 Normalized",
                           	'genome_id' : genome_id,
			   	'annotation_id' : annotation_id,
			   	'condition' : condition, 
			   	'mapped_rnaseq_alignment' : { sample_id : s_alignment },
			   	'tool_used' : "Cufflinks",
			   	'tool_version' : "version",
			   	'file' : handle
                	   	}
		
        		res= ws_client.save_objects(
                                   {"workspace":ws_id,
                                    "objects": [{
                                    "type":"KBaseRNASeq.RNASeqExpression",
                                    "data":es_obj,
                                    "name":output_name}
                                     ]})[0]
			expr_id = str(res[6]) + '/' + str(res[0]) + '/' + str(res[4]) 
		except Exception, e:
                	logger.exception("".join(traceback.format_exc()))
                	raise Exception("Failed to upload the ExpressionSample: {0}".format(output_name))
        except Exception,e:
		logger.exception("".join(traceback.format_exc()))
                raise Exception("Error executing cufflinks {0},{1}".format(cufflinks_command,directory))
	finally:
		if os.path.exists(out_file_path): os.remove(out_file_path)
		if os.path.exists(output_dir): shutil.rmtree(output_dir)
		ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqExpression",[output_name])	
		if not ret is None:
    		    return (alignment_name, output_name )
		else:
		    return None

def call_cuffmerge_and_cuffdiff(logger,ws_client,hs,ws_id,num_threads,list_file,gtf_file,bam_files,t_labels,genome_id,expressionset_id,alignmentset_id,sampleset_id,params,directory,token):
	 ## Adding Advanced options for cuffmerge command
	 cuffmerge_dir = os.path.join(directory,"cuffmerge")	
	 cuffmerge_command = " -p {0} -o {1} -g {2} {3}".format(str(num_threads),cuffmerge_dir,gtf_file,list_file)
         try:
                logger.info("Executing: cuffmerge {0}".format(cuffmerge_command))
	        script_util.runProgram(logger,"cuffmerge",cuffmerge_command,None,directory)
		if os.path.exists(cuffmerge_dir+"/merged.gtf") : merged_gtf = os.path.join(cuffmerge_dir,"merged.gtf")
         except Exception,e:
                raise Exception("Error executing cuffmerge {0},{1}".format(cuffmerge_command,cuffmerge_dir))
             
	 ## Adding Advanced options for cuffdiff command
       	 output_dir = os.path.join(directory,params['output_obj_name']) 
	 cuffdiff_command = (' -p '+str(num_threads))
         if('time_series' in params and params['time_series'] != 0) : cuffdiff_command += (' -T ')
         if('min_alignment_count' in params and params['min_alignment_count'] is not None ) : cuffdiff_command += (' -c '+str(params['min_alignment_count']))
         if('multi_read_correct' in params and params['multi_read_correct'] != 0 ): cuffdiff_command += (' --multi-read-correct ')
         if('library_type' in params and params['library_type'] is not None ) : cuffdiff_command += ( ' --library-type '+params['library_type'])
         if('library_norm_method' in params and params['library_norm_method'] is not None ) : cuffdiff_command += ( ' --library-norm-method '+params['library_norm_method'])
	 try:
         	cuffdiff_command += " -o {0} -L {1} -u {2} {3}".format(output_dir,t_labels,gtf_file,bam_files)
                logger.info("Executing: cuffdiff {0}".format(cuffdiff_command))
                ret = script_util.runProgram(None,"cuffdiff",cuffdiff_command,None,directory)
                result = ret["result"]
                for line in result.splitlines(False):
                       logger.info(line)
                       stderr = ret["stderr"]
                       prev_value = ''
                       for line in stderr.splitlines(False):
                           if line.startswith('> Processing Locus'):
                                   words = line.split()
                                   cur_value = words[len(words) - 1]
                                   if prev_value != cur_value:
                                      prev_value = cur_value
                                      logger.info(line)
                                   else:
                                      prev_value = ''
                                      logger.info(line)
         except Exception,e:
                raise Exception("Error executing cuffdiff {0},{1}".format(cuffdiff_command,directory))

                        ##  compress and upload to shock
         try:
                 logger.info("Zipping Cuffdiff output")
                 out_file_path = os.path.join(directory,"{0}.zip".format(params['output_obj_name']))
                 script_util.zip_files(logger,output_dir,out_file_path)
         except Exception,e:
                 raise Exception("Error executing cuffdiff")
         try:
                 handle = hs.upload(out_file_path)
         except Exception, e:
		 print " ".join(traceback.print_exc())	
                 raise Exception("Failed to upload the Cuffdiff output files: {0}".format(out_file_path))
	 output_name = params['output_obj_name']
         ## Save object to workspace
         try:
                 logger.info("Saving Cuffdiff object to workspace")
                 cm_obj = { "tool_used" : "Cuffdiff",
			    "tool_version" : "2.2.1",
			    "condition" : t_labels.split(","),
			    "genome_id" : genome_id,
			    "expressionSet_id" : expressionset_id,
			    "alignmentSet_id":alignmentset_id,
			    "sampleset_id" : sampleset_id,
			    "file" : handle
                           }
		 print cm_obj
                 res1= ws_client.save_objects(
                                             {"workspace":params['ws_id'],
                                               "objects": [{
                                               "type":"KBaseRNASeq.RNASeqDifferentialExpression",
                                               "data":cm_obj,
                                               "name":output_name}]})
         except Exception, e:
                 raise Exception("Failed to upload the KBaseRNASeq.RNASeqDifferentialExpression : {0}".format(output_name))
	 
	 return ( expressionset_id , output_name)
