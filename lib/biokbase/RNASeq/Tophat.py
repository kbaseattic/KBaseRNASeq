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
from biokbase.RNASeq.ExecutionBase import ExecutionBase
#import ExecutionBase.ExecutionBase as ExecutionBase

class TophatException(Exception):
    pass

class Tophat(ExecutionBase): 

    def __init__(self, logger, directory, urls):
        pprint(self.__class__)
        super(Tophat, self).__init__(logger, directory, urls)

        # user defined shared variables across methods
        #self.sample = None
        #self.sampleset_info = None
        self.num_threads = None


    def runEach(self,task_params):
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        
        read_sample = task_params['job_id']
        condition = task_params['label']
        directory = task_params['tophat_dir']
        ws_id = task_params['ws_id']
        reads_type = task_params['reads_type']
        genome_id = task_params['annotation_id']
        sampleset_id = task_params['sampleset_id']
	gtf_file = task_params['gtf_file']

        print "Downloading Read Sample{0}".format(read_sample)
        logger.info("Downloading Read Sample{0}".format(read_sample))
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
                ### Adding advanced options to Bowtie2Call
		tophat_cmd = (' -p '+str(self.num_threads))
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
                #out_file = output_dir +"/accepted_hits.sam"
                if sample_type  == 'KBaseAssembly.SingleEndLibrary':
                        lib_type = 'SingleEnd'
                        read_id = r_sample['data']['handle']['id']
                        read_name =  r_sample['data']['handle']['file_name']
                        try:
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read_id,filename=read_name, directory=directory,token=token)
                		tophat_cmd += ' -o {0} -G {1} {2} {3}'.format(output_dir,gtf_file,bowtie2_base,os.path.join(directory,read_name))
                        except Exception,e:
                                self.logger.exception(e)
                                raise Exception( "Unable to download shock file , {0}".format(read_name))
                if sample_type == 'KBaseAssembly.PairedEndLibrary':
                        lib_type = 'PairedEnd'
                        if('orientation' in params and params['orientation'] is not None): tophat_cmd += ( ' --'+params['orientation'])
                        read1_id = r_sample['data']['handle_1']['id']
                        read1_name = r_sample['data']['handle_1']['file_name']
                        read2_id = r_sample['data']['handle_2']['id']
                        read2_name = r_sample['data']['handle_2']['file_name']
                        try:
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=directory,token=token)
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=directory,token=token)
                		tophat_cmd += ' -o {0} -G {1} {2} {3} {4}'.format(output_dir,gtf_file,bowtie2_base,os.path.join(directory,read1_name),os.path.join(directory,read2_name))
                        except Exception,e:
                                raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
                try:
                        self.logger.info("Executing: tophat {0}".format(tophat_cmd))
                        cmdline_output, cmd_err = script_util.runProgram(self.logger,"tophat",tophat_cmd,None,directory)
                except Exception,e:
                        raise Exception("Failed to run command {0}\n{1}\n{2}".format(tophat_cmd,cmdline_output,cmd_err))
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
                out_file_path = os.path.join(directory,"%s.zip" % output_name)
            	try:
                        logger.info("Zipping the output files".format(out_file_path))
                	script_util.zip_files(logger, output_dir,out_file_path)
            	except Exception, e:
                	raise Exception("Failed to compress the index: {0}".format(out_file_path))
                ## Upload the file using handle service
                try:
                        tophat_handle = hs.upload(out_file_path)
                except Exception, e:
                        raise Exception("Failed to upload zipped output file".format(out_file_path))
                #### Replace version with get_version command#####
		logger.info("Preparing output object")
                tophat_out = { "file" : tophat_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "tophat" , "aligner_version" : "2.2.1" , 'library_type' : lib_type , 'condition' : condition ,'read_sample_id': read_sample, 'genome_id' : genome_id , 'bowtie2_index': self.bowtie2index_id, "alignment_stats" : stats_data }
                if not sampleset_id is None: tophat_out['sampleset_id'] = sampleset_id
                pprint(tophat_out)
		try:
                        res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignment",
                                         "data":tophat_out,
                                         "name":output_name}
                                        ]})
                except Exception, e:
			raise Exception(e)
                        #logger.exception("Failed to save alignment to workspace")
                        raise Exception("Failed to save alignment to workspace")
        except Exception, e:
                        #logger.exception("Failed to create tophat Alignment {0}".format(" ".join(traceback.print_exc())))
                        raise Exception("Failed to create tophat Alignment {0}".format(" ".join(traceback.print_exc())))
        finally:
                if os.path.exists(out_file_path): os.remove(out_file_path)
                if os.path.exists(output_dir): shutil.rmtree(output_dir)
                ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqAlignment",[output_name])
                if not ret is None:
                    return (read_sample,output_name)
		else :
		    return None
