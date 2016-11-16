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
from biokbase.RNASeq.KBParallelExecutionBase import KBParallelExecutionBase


class TophatException(Exception):
    pass

class Tophat(KBParallelExecutionBase): 

    def __init__(self, logger, directory, urls):
        pprint(self.__class__)
        super(Tophat, self).__init__(logger, directory, urls)

        # user defined shared variables across methods
        #self.sample = None
        #self.sampleset_info = None
        self.num_threads = None


    def runEach(self,task_params):

        print( "in Tophat.runEach(), task_params are")
        pprint( task_params )
        ws_client = self.common_params['ws_client']
        hs = self.common_params['hs_client']
        params = self.method_params
        logger = self.logger
        token = self.common_params['user_token']
        
        job_id = task_params['job_id']
        condition = task_params['label']
        directory = task_params['tophat_dir']
        ws_id = task_params['ws_id']
        reads_type = task_params['reads_type']
        genome_id = task_params['annotation_id']
        sampleset_id = task_params['sampleset_id']
        gtf_file = task_params['gtf_file']

        # fetch bowtie2 index file if necessary
        if ( handler_util.get_file_with_suffix(directory,".rev.1.bt2") == None ):
            try:
                   bowtie_index = ws_client.get_objects( [ { 'name' : task_params['bowtie_index'], 'workspace' : task_params['ws_id'] } ] )
            except Exception,e:
                   logger.exception("".join(traceback.format_exc()))
                   raise ValueError(" Error Downloading bowtie_index object from the workspace ")
    
            bowtie2_index_info = ws_client.get_object_info_new( {"objects": [{'name': task_params['bowtie_index'], 'workspace': task_params['ws_id']} ] } )
    
            self.bowtie2index_id = str(bowtie2_index_info[6]) + '/' + str(bowtie2_index_info[0]) + '/' + str(bowtie2_index_info[4])  
            bw_id = bowtie_index['data']['handle']['id'] 
            bw_name =  bowtie_index['data']['handle']['file_name']
            genome_id = bowtie_index['data']['genome_id']
            annotation_gtf = ws_client.get_object_info( [ {"ref" :genome_id}], includeMetadata=None )[0][1]
            shared_files={}
            shared_files[bw_name] = bw_id
            script_util.download_shock_files(logger,self.urls['shock_service_url'],directory,shared_files,token)
            try:
                logger.info("Unzipping Bowtie2 Indices")
                script_util.unzip_files(logger,os.path.join(directory,bw_name),directory)
                mv_dir= handler_util.get_dir(directory)
                if mv_dir is not None:
                        script_util.move_files(logger,mv_dir,directory)
            except Exception, e:
                   logger.error("".join(traceback.format_exc()))
                   raise Exception("Unzip indexfile error: Please contact help@kbase.us")
            #fasta_file =os.path.join(directory,(handler_util.get_file_with_suffix(directory,".fa")+".fa"))
            #bowtie2base =os.path.join(directory,handler_util.get_file_with_suffix(directory,".rev.1.bt2"))
            ##############

        # fetch GTF file if necessary  (right now this just repeats whats in prepare() - FIX)

        ### Check if GTF annotation object exist or skip this step
        ### Check if the gtf object exists in the workspace
        ### Only run create_gtf_annotation if object doesnt exist
        ws_gtf = annotation_gtf+"_GTF_Annotation"
        ret = script_util.if_obj_exists(None,ws_client,params['ws_id'],"KBaseRNASeq.GFFAnnotation",[ws_gtf])
        if not ret is None:
            logger.info("GFF Annotation Exist for Genome Annotation {0}.... Skipping step ".format(annotation_gtf))
            annot_name,annot_id = ret[0]
            gtf_obj=ws_client.get_objects([{'ref' : annot_id}])[0]
            gtf_id=gtf_obj['data']['handle']['id']
            gtf_name=gtf_obj['data']['handle']['file_name']
            try:
               script_util.download_file_from_shock(logger, shock_service_url=self.urls['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=tophat_dir,token=token)
               gtf_file = os.path.join(tophat_dir,gtf_name)
            except Exception,e:
               logger.exception(e)
               raise Exception( "Unable to download shock file, {0}".format(gtf_name))  
        else:
            gtf_file =rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],genome_id,annotation_gtf,tophat_dir,token)        


        print "Downloading Read Sample{0}".format(sampleset_id)
        logger.info("Downloading Read Sample{0}".format(sampleset_id))
        try:
                r_sample = ws_client.get_objects(
                                        [{ 'name' : sampleset_id, 'workspace' : ws_id}])[0]
                r_sample_info = ws_client.get_object_info_new({"objects": [{'name': sampleset_id, 'workspace': ws_id}]})[0]
                sample_type = r_sample_info[2].split('-')[0]
                output_name = sampleset_id.split('.')[0]+"_tophat_alignment"
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
                tophat_out = { "file" : tophat_handle ,
                               "size" : os.path.getsize(out_file_path),
                               "aligned_using" : "tophat" , 
                               "aligner_version" : "2.2.1" , 
                               'library_type' : lib_type , 
                               'condition' : condition ,
                               'sample_id': sample_id, 
                               'genome_id' : genome_id , 
                               'bowtie2_index': self.bowtie2index_id,
                               'alignment_stats' : stats_data }
                if not sampleset_id is None: tophat_out['sampleset_id'] = sampleset_id
                pprint(tophat_out)
                try:
                        res= ws_client.save_objects( { "workspace":ws_id,
                                                       "objects": [ { "type":"KBaseRNASeq.RNASeqAlignment",
                                                                      "data":tophat_out,
                                                                      "name":output_name } ]
                                                     } )
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
                    return { 'sampleset_id': sampleset_id, 'output_name': output_name }
                else :
                    return None
