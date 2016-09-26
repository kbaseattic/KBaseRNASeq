import simplejson, sys, shutil, os, ast , re
from mpipe import OrderedStage , Pipeline
import glob, json, uuid, logging  , time ,datetime
import subprocess, threading,traceback
from collections import OrderedDict
from pprint import pprint , pformat
import parallel_tools as parallel
from mpipe import OrderedStage , Pipeline
import contig_id_mapping as c_mapping
import numpy
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

TOOL_USED = 'StringTie'
TOOL_VERSION = '1.2.3'

try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

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
def CallStringtie_helper(x):
        logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token = x
        return  _CallStringtie(logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token)

def _CallStringtie(logger,services,ws_client,hs,ws_id,num_threads,s_alignment,gtf_file,directory,genome_id,annotation_id,sample_id,alignmentset_id,params,token):
        print "Downloading Read Sample{0}".format(s_alignment)
        alignment_name = ws_client.get_object_info([{"ref" :s_alignment}],includeMetadata=None)[0][1]
        if not logger:
                logger = handler_util.create_logger(directory,"run_Stringtie_"+alignment_name)
        try:
                alignment = ws_client.get_objects(
                                        [{ 'ref' : s_alignment }])[0]
                input_direc = os.path.join(directory,alignment_name.split('_alignment')[0]+"_stringtie_input")
                if not os.path.exists(input_direc) : os.mkdir(input_direc)
                output_name = alignment_name.split('_alignment')[0]+"_stringtie_expression"
                output_dir = os.path.join(directory,output_name)
                #Download Alignment from shock
                a_file_id = alignment['data']['file']['id']
                a_filename = alignment['data']['file']['file_name']
                condition = alignment['data']['condition']
                try:
                     script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=a_file_id,filename=a_filename,directory=input_direc,token=token)
                except Exception,e:
                        raise Exception( "Unable to download shock file, {0},{1}".format(a_filename,"".join.tracback.format_exc()))
                try:
                    input_dir = os.path.join(input_direc,alignment_name)
                    if not os.path.exists(input_dir): os.mkdir(input_dir)
                    script_util.unzip_files(logger,os.path.join(input_direc,a_filename), input_dir)
                except Exception, e:
                       logger.error("".join(traceback.format_exc()))
                       raise Exception("Unzip alignment files  error: Please contact help@kbase.us")

                input_file = os.path.join(input_dir,"accepted_hits.bam")
                ### Adding advanced options to tophat command
		tool_opts = { k:str(v) for k,v in params.iteritems() if not k in ('ws_id','alignmentset_id', 'num_threads') and v is not None  }
                stringtie_command = (' -p '+str(num_threads))
                if 'label' in params and params['label'] is not None:
                     stringtie_command += (' -l '+str(params['label']))
                if 'min_isoform_abundance' in params and params['min_isoform_abundance'] is not None:
                     stringtie_command += (' -f '+str(params['min_isoform_abundance']))
                if 'min_length' in params  and params['min_length'] is not None:
                     stringtie_command += (' -m '+str(params['min_length']))
                if 'a_juncs' in params  and params['a_juncs'] is not None:
                     stringtie_command += (' -a '+str(params['a_juncs']))
                if 'j_min_reads' in params  and params['j_min_reads'] is not None:
                     stringtie_command += (' -j '+str(params['j_min_reads']))
                if 'c_min_read_coverage' in params  and params['c_min_read_coverage'] is not None:
                     stringtie_command += (' -c '+str(params['c_min_read_coverage']))
                if 'gap_sep_value' in params  and params['gap_sep_value'] is not None:
                     stringtie_command += (' -g '+str(params['gap_sep_value']))
                if 'disable_trimming' in params  and params['disable_trimming'] != 0:
                     stringtie_command += (' -t ')
                if 'ballgown_mode' in params  and params['ballgown_mode'] != 0:
                     stringtie_command += (' -B ')
                if 'skip_reads_with_no_ref' in params  and params['skip_reads_with_no_ref'] != 0:
                     stringtie_command += (' -e ')
                t_file_name = os.path.join(output_dir,"transcripts.gtf")
		g_output_file = os.path.join(output_dir,"genes.fpkm_tracking")
                stringtie_command += " -o {0} -A {1} -G {2} {3}".format(t_file_name,g_output_file,gtf_file,input_file)
                logger.info("Executing: stringtie {0}".format(stringtie_command))
                ret = script_util.runProgram(None,"stringtie",stringtie_command,None,directory)
                ##Parse output files
		try:
                	exp_dict = script_util.parse_FPKMtracking(g_output_file,'StringTie','FPKM')
                	tpm_exp_dict = script_util.parse_FPKMtracking(g_output_file,'StringTie','TPM')
		except Exception,e:
                        logger.exception("".join(traceback.format_exc()))
                        raise Exception("Error parsing FPKMtracking")

                ##  compress and upload to shock
                try:
                        logger.info("Zipping Stringtie output")
                        out_file_path = os.path.join(directory,"%s.zip" % output_name)
                        script_util.zip_files(logger,output_dir,out_file_path)
                except Exception,e:
                        logger.exception("".join(traceback.format_exc()))
                        raise Exception("Error executing stringtie")
                try:
			handle = hs.upload(out_file_path)
                except Exception, e:
                        logger.exception("".join(traceback.format_exc()))
                        raise Exception("Error while zipping the output objects: {0}".format(out_file_path))
                ## Save object to workspace
                try:
                        logger.info("Saving Stringtie object to workspace")
                        es_obj = { 'id' : output_name,
                                'type' : 'RNA-Seq',
                                'numerical_interpretation' : 'FPKM',
                                'expression_levels' : exp_dict,
                                'tpm_expression_levels' : tpm_exp_dict,
                                'processing_comments' : "log2 Normalized",
                                'genome_id' : genome_id,
                                'annotation_id' : annotation_id,
                                'condition' : condition,
                                'mapped_rnaseq_alignment' : { sample_id : s_alignment },
                                'tool_used' : TOOL_USED,
                                'tool_version' : TOOL_VERSION,
                                'tool_opts' : tool_opts,
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
                raise Exception("Error executing stringtie {0},{1}".format(cufflinks_command,directory))
        finally:
                if os.path.exists(out_file_path): os.remove(out_file_path)
                if os.path.exists(output_dir): shutil.rmtree(output_dir)
                if os.path.exists(input_direc): shutil.rmtree(input_direc)
                ret = script_util.if_obj_exists(None,ws_client,ws_id,"KBaseRNASeq.RNASeqExpression",[output_name])
                if not ret is None:
                    return (alignment_name, output_name )
                else:
                    return None

def runMethod(logger,token,ws_client,hs,services,stringtie_dir,params):
	    try:
                a_sample = ws_client.get_objects(
                                        [{'name' : params['alignmentset_id'],'workspace' : params['ws_id']}])[0]
            except Exception,e:
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error Downloading objects from the workspace ")
            ## Get the Input object type ##
            a_sample_info = ws_client.get_object_info_new({"objects": [{'name': params['alignmentset_id'], 'workspace': params['ws_id']}]})[0]
            a_sample_type = a_sample_info[2].split('-')[0]
            alignmentset_id = str(a_sample_info[6]) + '/' + str(a_sample_info[0]) + '/' + str(a_sample_info[4])
	    ## Check if the Alignment objects exist in the same workspace
	    logger.info("Check if the Alignment objects do exist in the current workspace")
            if a_sample_type == 'KBaseRNASeq.RNASeqAlignmentSet':
                a_names = list(numpy.array([ i.values() for i in a_sample['data']['mapped_rnaseq_alignments']]).flatten())
                a_type = 'KBaseRNASeq.RNASeqAlignment'
                e_ws_objs = script_util.if_ws_obj_exists(None,ws_client,params['ws_id'],a_type,a_names)
                missing_objs = [i for i in a_names if not i in e_ws_objs]
                if len(e_ws_objs) != len(a_names):
                   raise ValueError('Missing Alignment objects {0} in the {1}. please copy them and run this method'.format(",".join(missing_objs),params['ws_id']))

            ### Check if the gtf file exists in the workspace. if exists download the file from that
            annotation_id = a_sample['data']['genome_id']
            annotation_name = ws_client.get_object_info([{"ref" :annotation_id}],includeMetadata=None)[0][1]
            gtf_obj_name = annotation_name+"_GTF_Annotation"
            ret = script_util.if_obj_exists(None,ws_client,params['ws_id'],"KBaseRNASeq.GFFAnnotation",[gtf_obj_name])
            if not ret is None:
                logger.info("GFF Annotation Exist for Genome Annotation {0}.... Skipping step ".format(annotation_name))
                gtf_obj= ws_client.get_objects([{'name' : gtf_obj_name,'workspace' : params['ws_id']}])[0]
                gtf_info = ws_client.get_object_info_new({"objects": [{'name': gtf_obj_name, 'workspace': params['ws_id']}]})[0]
                gtf_annotation_id = str(gtf_info[6]) + '/' + str(gtf_info[0]) + '/' + str(gtf_info[4])
                gtf_id=gtf_obj['data']['handle']['id']
                gtf_name=gtf_obj['data']['handle']['file_name']
                try:
                     script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=stringtie_dir,token=token)
                     gtf_file = os.path.join(stringtie_dir,gtf_name)
                except Exception,e:
                     raise Exception( "Unable to download shock file, {0}".format(gtf_name))
            else:
		rnaseq_util.create_gtf_annotation_from_genome(logger,ws_client,hs,self.urls,params['ws_id'],annotation_id,annotation_name,stringtie_dir,token)
                #fasta_file= script_util.generate_fasta(logger,services,token,annotation_id,stringtie_dir,annotation_name)
                #logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
                #mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
                #c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
                #logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))

                #gtf_file = script_util.create_gtf_annotation(logger,ws_client,hs,services['shock_service_url'],params['ws_id'],annotation_id,gtf_obj_name,fasta_file,stringtie_dir,token)
	    #Identify the tool options used
	    tool_opts = { k:str(v) for k,v in params.iteritems() if not k in ('ws_id','alignmentset_id', 'num_threads') and v is not None  }
            # Determine the num_threads provided by the user otherwise default the number of threads to 2
            if('num_threads' in params and params['num_threads'] is not None):
                        num_threads = int(params['num_threads'])
            else:
                        num_threads = 2
            num_cores = mp.cpu_count()
            logger.info("Number of available cores : {0}".format(num_cores))
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
                num_samples =  len(alignment_ids)
                if num_samples < 2:
                        raise ValueError("Please ensure you have atleast 2 alignments to run Stringtie in Set mode")
                if num_cores != 1:
                        pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
                else:
                   pool_size = 1
                   num_threads = 1
                logger.info(" Number of threads used by each process {0}".format(num_threads))
                for d_align in m_alignment_ids:
                  for s_name,a_id in d_align.items():
                        try:
                                b_tasks.append((None,services,ws_client,hs,params['ws_id'],num_threads,a_id,gtf_file,stringtie_dir,annotation_id,gtf_annotation_id,s_name,alignmentset_id,params,token))
                        except Exception,e:
                                raise

                @parallelize(CallStringtie_helper,pool_size)
                def run_stringtie_in_parallel(tasks):
                  pass
                results=run_stringtie_in_parallel(b_tasks)
                expressionSet_name = params['alignmentset_id']+"_StringTie_ExpressionSet"
                reportObj=script_util.create_RNASeq_ExpressionSet_and_build_report(logger,ws_client,TOOL_USED, TOOL_VERSION,tool_opts,params['ws_id'],align_names,alignmentset_id,annotation_id,sampleset_id,results,expressionSet_name)
            else:
                try:
                    pool_size=1
                    num_threads = num_cores
                    single_a_id = alignmentset_id
                    logger.info(" Number of threads used by each process {0}".format(num_threads))
                    results = _CallStringtie(None,services,ws_client,hs,params['ws_id'],num_threads,single_a_id,gtf_file,stringtie_dir,annotation_id,gtf_annotation_id,params['alignmentset_id'],None,params,token)
                except Exception,e:
                     raise
                single_alignment, single_expression = results
                single_expr_obj = ws_client.get_objects(
                                        [{ 'name' : single_expression, 'workspace' : params['ws_id']}])[0]['data']
                e_ref = ws_client.get_object_info_new({"objects": [{'name':single_expression, 'workspace': params['ws_id']}]})[0]
                reportObj = {'objects_created':[{'ref' :str(e_ref[6]) + '/' + str(e_ref[0]) + '/' + str(e_ref[4]),
                                                 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_alignment)}],
                                                 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_alignment)}
            #### Save to report object #######
                #returnVal = single_align_obj   
            reportName = 'Assemble_Transcripts_Using_Stringtie_'+str(hex(uuid.getnode()))
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
	    return returnVal
