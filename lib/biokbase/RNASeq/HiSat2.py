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

class HiSat2Exception(Exception):
    pass

class HiSat2(ExecutionBase): 

    def __init__(self, logger, directory, urls):
        pprint(self.__class__)
        super(HiSat2, self).__init__(logger, directory, urls)

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
        directory = task_params['hisat2_dir']
        ws_id = task_params['ws_id']
        reads_type = task_params['reads_type']
        genome_id = task_params['annotation_id']
        sampleset_id = task_params['sampleset_id']

        print "Downloading Read Sample{0}".format(read_sample)
        logger.info("Downloading Read Sample{0}".format(read_sample))
        try:
                r_sample = ws_client.get_objects(
                                        [{ 'name' : read_sample, 'workspace' : ws_id}])[0]
                r_sample_info = ws_client.get_object_info_new({"objects": [{'name': read_sample, 'workspace': ws_id}]})[0]
                sample_type = r_sample_info[2].split('-')[0]
                input_direc = os.path.join(directory,read_sample.split('.')[0]+"_hisat2_input")
                if not os.path.exists(input_direc): os.mkdir(input_direc)
                output_name = read_sample.split('.')[0]+"_hisat2_alignment"
                output_dir = os.path.join(directory,output_name)
                if not os.path.exists(output_dir): os.mkdir(output_dir)
                print directory
                base = handler_util.get_file_with_suffix(directory,".1.ht2")
                print base
                hisat2_base =os.path.join(directory,base)
                ### Adding advanced options to Bowtie2Call
                hisat2_cmd = ''
                hisat2_cmd += ( ' -p {0}'.format(self.num_threads))
                if('quality_score' in params and params['quality_score'] is not None): hisat2_cmd += ( ' --'+params['quality_score'])
                if('alignment_type' in params and params['alignment_type'] is not None): hisat2_cmd += ( ' --'+params['alignment_type'] )
                if('trim5' in params and params['trim5'] is not None): hisat2_cmd += ( ' --trim5 '+str(params['trim5']))
                if('trim3' in params and params['trim3'] is not None): hisat2_cmd += ( ' --trim3 '+str(params['trim3']))
                if('np' in params and params['np'] is not None): hisat2_cmd += ( ' --np '+str(params['np']))
                if('minins' in params and params['minins'] is not None): hisat2_cmd += ( ' --minins '+str(params['minins']))
                if('maxins' in params and params['maxins'] is not None): hisat2_cmd += ( ' --maxins '+str(params['maxins']))
                #if('orientation' in params and params['orientation'] is not None): hisat2_cmd += ( ' --'+params['orientation'])
                if('min_intron_length' in params and params['min_intron_length'] is not None): hisat2_cmd += ( ' --min-intronlen '+str(params['min_intron_length']))
                if('max_intron_length' in params and params['max_intron_length'] is not None): hisat2_cmd += ( ' --max-intronlen '+str(params['max_intron_length']))
                if('no_spliced_alignment' in params and params['no_spliced_alignment'] != 0): hisat2_cmd += ( ' --no-spliced-alignment')
                if('transcriptome_mapping_only' in params and params['transcriptome_mapping_only'] != 0): hisat2_cmd += ( ' --transcriptome-mapping-only')
                if('tailor_alignments' in params and params['tailor_alignments'] is not None): 
                    hisat2_cmd += ( ' --'+params['tailor_alignments'])
                out_file = output_dir +"/accepted_hits.sam"
                if sample_type  == 'KBaseAssembly.SingleEndLibrary':
                        lib_type = 'SingleEnd'
                        read_id = r_sample['data']['handle']['id']
                        read_name =  r_sample['data']['handle']['file_name']
                        try:
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read_id,filename=read_name, directory=input_direc,token=token)
                                hisat2_cmd += " -U {0} -x {1} -S {2}".format(os.path.join(input_direc,read_name),hisat2_base,out_file)
                        except Exception,e:
                                self.logger.exception(e)
                                raise Exception( "Unable to download shock file , {0}".format(read_name))
                if sample_type == 'KBaseAssembly.PairedEndLibrary':
                        lib_type = 'PairedEnd'
                        if('orientation' in params and params['orientation'] is not None): hisat2_cmd += ( ' --'+params['orientation'])
                        read1_id = r_sample['data']['handle_1']['id']
                        read1_name = r_sample['data']['handle_1']['file_name']
                        read2_id = r_sample['data']['handle_2']['id']
                        read2_name = r_sample['data']['handle_2']['file_name']
                        try:
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=input_direc,token=token)
                                script_util.download_file_from_shock(self.logger, shock_service_url=self.urls['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=input_direc,token=token)
                                hisat2_cmd += " -1 {0} -2 {1} -x {2} -S {3}".format(os.path.join(input_direc,read1_name),os.path.join(input_direc,read2_name),hisat2_base,out_file)
                        except Exception,e:
                                self.logger.exception(e)
                                raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
                try:
                        self.logger.info("Executing: hisat2 {0}".format(hisat2_cmd))
                        cmdline_output = script_util.runProgram(self.logger,"hisat2",hisat2_cmd,None,directory)
                except Exception,e:
                        raise Exception("Failed to run command {0}".format(hisat2_cmd))
                try:
                        stats_data = {}
                        stats_data = script_util.extractAlignmentStatsInfo(self.logger,"bowtie2",ws_client,ws_id,None,cmdline_output['stderr'],None)
                        bam_file = os.path.join(output_dir,"accepted_hits_unsorted.bam")
                        logger.info("Executing: sam_to_bam  {0}".format(bam_file))
                        sam_to_bam = "view -bS -o {0} {1}".format(bam_file,out_file)
                        script_util.runProgram(self.logger,"samtools",sam_to_bam,None,directory)
                        final_bam_prefix = os.path.join(output_dir,"accepted_hits")
                        logger.info("Executing: Sorting bam file  {0}".format(bam_file))
                        sort_bam_cmd  = "sort {0} {1}".format(bam_file,final_bam_prefix)
                        script_util.runProgram(self.logger,"samtools",sort_bam_cmd,None,directory)
                except Exception,e:
                        raise Exception("Error Running the hisat2 command {0},{1} {2}".format(hisat2_cmd,directory," ".join(traceback.print_exc())))

                # Zip tophat folder
                out_file_path = os.path.join(directory,"%s.zip" % output_name)
                try:
                        logger.info("Zipping the output files".format(out_file_path))
                        script_util.zip_files(self.logger, output_dir,out_file_path)
                except Exception, e:
                        raise Exception("Failed to compress the index: {0}".format(out_file_path))
                ## Upload the file using handle service
                try:
                        hisat2_handle = hs.upload(out_file_path)
                except Exception, e:
                        raise Exception("Failed to upload zipped output file".format(out_file_path))
                #### Replace version with get_version command#####
		logger.info("Preparing output object")
                hisat2_out = { "file" : hisat2_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "hisat2" , "aligner_version" : "2.2.6" , 'library_type' : lib_type , 'condition' : condition ,'read_sample_id': read_sample, 'genome_id' : genome_id , "alignment_stats" : stats_data }
                if not sampleset_id is None: hisat2_out['sampleset_id'] = sampleset_id
                pprint(hisat2_out)
		try:
                        res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignment",
                                         "data":hisat2_out,
                                         "name":output_name}
                                        ]})
                except Exception, e:
			raise Exception(e)
                        #logger.exception("Failed to save alignment to workspace")
                        raise Exception("Failed to save alignment to workspace")
        except Exception, e:
                        #logger.exception("Failed to create hisat2 Alignment {0}".format(" ".join(traceback.print_exc())))
                        raise Exception("Failed to create hisat2 Alignment {0}".format(" ".join(traceback.print_exc())))
        finally:
                if os.path.exists(input_direc): shutil.rmtree(input_direc)
                if os.path.exists(out_file_path): os.remove(out_file_path)
                if os.path.exists(output_dir): shutil.rmtree(output_dir)
                ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqAlignment",[output_name])
                if not ret is None:
                    return (read_sample,output_name)
        return None
