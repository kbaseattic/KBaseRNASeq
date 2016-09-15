import simplejson, sys, shutil, os, ast , re
from mpipe import OrderedStage , Pipeline
import glob, json, uuid, logging  , time ,datetime
import subprocess, threading,traceback
from collections import OrderedDict
from pprint import pprint , pformat
import parallel_tools as parallel
import itertools
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

TOOL1_USED = 'StringTie'
TOOL2_USED = 'TableMaker'
TOOL1_VERSION = '1.2.3'
TOOL2_VERSION = '2.1.1'
ASSEMBLY_GTF_FN='assembly_GTF_list.txt'
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
def CalldiffExpCallforBallgown_helper(x):
        logger,services,ws_client,hs,ws_id,num_threads,alignment_file,transcripts_gtf,list_file,used_tool,directory,gtf_file = x
        return  _CalldiffExpCallforBallgown(logger,services,ws_client,hs,ws_id,num_threads,alignment_file,transcripts_gtf,list_file,used_tool,directory,gtf_file)

def call_cuffmerge(directory,num_threads,gtf_file,list_file):
	 #cuffmerge_dir = os.path.join(directory,"cuffmerge")
         cuffmerge_command = " -p {0} -o {1} -g {2} {3}".format(str(num_threads),directory,gtf_file,list_file)
         merged_gtf = None
	 try:
                logger.info("Executing: cuffmerge {0}".format(cuffmerge_command))
                script_util.runProgram(logger,"cuffmerge",cuffmerge_command,None,directory)
                if os.path.exists(cuffmerge_dir+"/merged.gtf") : merged_gtf = os.path.join(directory,"merged.gtf")
         except Exception,e:
                raise Exception("Error executing cuffmerge {0},{1}".format(cuffmerge_command,directory))
	 return merged_gtf

def call_stringtiemerge(directory,num_threads,gtf_file,list_file):
         #directory = os.path.join(directory,"cuffmerge")
         strmerge_command = " -p {0} -o {1} --merge -G {2} {3}".format(str(num_threads),directory,gtf_file,list_file)
         merged_gtf = None
         try:
                logger.info("Executing: stringtie {0}".format(strmerge_command))
                script_util.runProgram(logger,"stringtie",strmerge_command,None,directory)
                if os.path.exists(directory+"/merged.gtf") : merged_gtf = os.path.join(directory,"merged.gtf")
         except Exception,e:
                raise Exception("Error executing StringTie merge {0},{1}".format(strmerge_command,directory))
         return merged_gtf

def call_tablemaker(directory,num_threads,m_gtf_file,alignment_file):
         #cuffmerge_dir = os.path.join(directory,"cuffmerge")
         tm_command = " -p {0} -o {1} -q -W -G {2} {3}".format(str(num_threads),directory,m_gtf_file,alignment_file)
         try:
                logger.info("Executing: tablemaker {0}".format(tm_command))
                script_util.runProgram(logger,"tablemaker",tm_command,None,directory)
         except Exception,e:
                raise Exception("Error executing tablemaker {0},{1}".format(tm_command,directory))
         return directory


def call_stringtieBall(directory,num_threads,m_gtf_file,alignment_file):
         #directory = os.path.join(directory,"cuffmerge")
         strdiff_command = " -p {0} -o {1} -e -B -G {2} {3}".format(str(num_threads),directory,m_gtf_file,alignment_file)
         try:
                logger.info("Executing: stringtie {0}".format(strdiff_command))
                script_util.runProgram(logger,"stringtie",strdiff_command,None,directory)
         except Exception,e:
                raise Exception("Error executing StringTie differential expression {0},{1}".format(strdiff_command,directory))
         return directory

def _CalldiffExpCallforBallgown(logger,services,ws_client,hs,ws_id,num_threads,alignment_file,transcripts_gtf,list_file,used_tool,directory,gtf_file):
	### Create output directory name as ballgown/RNASeq_sample_name/ under diffexp_dir
	### Get i as  alignment_file
	### Get j as expression file
	### If tool is 'StringTie: Then call function call_stringtiemerge ; return ballgown/RNASeq_sample_name/merged.gtf ; Call function call_stringtieBall
        ### else if tool is 'TableMaker'; Then call function call_cuffmerge; return ballgown/RNASeq_sample_name/merged.gtf ; Call function call_tablemaker
	### return the  j and created paths. 
        print "Running Differential Expression steps for {0}".format(transcripts_gtf)
        if not logger:
                logger = handler_util.create_logger(directory,"run_diffExpCallforBallgown_"+str(hex(uuid.getnode())))
        try:
		merge_dir = os.path.join(directory,"merge") 
		if not os.path.exists(merge_dir): os.mkdir(merge_dir)
		print merge_dir
		ballgown_dir = os.path.join(directory,"ballgown")
		if not os.path.exists(ballgown_dir): os.mkdir(ballgown_dir)
		print ballgown_dir
		output_name = transcripts_gtf.split("_expression")[0]		
                output_dir = os.path.join(ballgown_dir,output_name)
		if not os.path.exists(output_dir): os.mkdir(output_dir)
		print output_dir
                #Download Alignment from shock
                #condition = expression['data']['condition']
		if used_tool == 'StringTie':
			print "Entering StringTie"
			merged_gtf = call_stringtiemerge(merge_dir,num_threads,gtf_file,list_file)
			call_stringtieBall(ballgown_dir,num_threads,merged_gtf,alignment_file)
                elif used_tool == 'Cufflinks':
			print "Entering Tablemaker"
			merged_gtf = call_cuffmerge(merge_dir,num_threads,gtf_file,list_file)	
			call_tablemaker(ballgown_dir,num_threads,m_gtf_file,alignment_file)
		if os.path.exists(ballgown_dir+"/t_data.ctab") :
			logger.info("Running Differential Expression for Sample {0} completed successfully".format(transcripts_gtf))
			print("Running Differential Expression for Sample {0} completed successfully".format(transcripts_gtf))
        except Exception,e:
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error executing stringtie {0},{1}".format(cufflinks_command,directory))
        finally:
                if os.path.exists(output_dir): shutil.rmtree(output_dir)
                if not ret is None:
                    return (transcripts_gtf, output_dir )
                else:
                    return None

def runMethod(logger,token,ws_client,hs,services,diffexp_dir,params):
	    try:
                e_sample = ws_client.get_objects(
                                        [{'name' : params['expressionset_id'],'workspace' : params['ws_id']}])[0]
            except Exception,e:
                logger.exception("".join(traceback.format_exc()))
                raise Exception("Error Downloading objects from the workspace ")
            ## Get the Input object type and info #
            e_sample_info = ws_client.get_object_info_new({"objects": [{'name': params['expressionset_id'], 'workspace': params['ws_id']}]})[0]
            e_sample_type = e_sample_info[2].split('-')[0]
            expressionset_id = str(e_sample_info[6]) + '/' + str(e_sample_info[0]) + '/' + str(e_sample_info[4])
	    alignmentset_id = e_sample['data']['alignmentSet_id'] 
	    sampleset_id = e_sample['data']['sampleset_id']
            expression_ids = e_sample['data']['sample_expression_ids']
            num_samples = len(expression_ids)
            if num_samples < 2:
               raise ValueError("Please ensure you have atleast 2 expressions to run diffExpCallforBallgown in Set mode")
            ### Check if the gtf file exists in the workspace. if exists download the file from that
            annotation_id = e_sample['data']['genome_id']
            logger.info("Check if the gtf file exists in the workspace".format(annotation_id))
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
                     script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=diffexp_dir,token=token)
                     gtf_file = os.path.join(diffexp_dir,gtf_name)
                except Exception,e:
                     raise Exception( "Unable to download shock file, {0}".format(gtf_name))
            else:
                fasta_file= script_util.generate_fasta(logger,services,token,annotation_id,diffexp_dir,annotation_name)
                logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
                mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
                c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
                logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
                gtf_file = script_util.create_gtf_annotation(logger,ws_client,hs,services,params['ws_id'],annotation_id,gtf_obj_name,fasta_file,diffexp_dir,token)
            m_expr_ids = e_sample['data']['mapped_expression_ids']
	    m_align_exp = []
            labels = []
            expressions = []
            counter = 0
            assembly_file = os.path.join(diffexp_dir,ASSEMBLY_GTF_FN)
            list_file = open(assembly_file,'w')
            for i in m_expr_ids:
                for a_id ,e_id in i.items():
                        #print a_id  + ":" + e_id
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
                        s_path = os.path.join(diffexp_dir,condition+"/"+str(counter)) ### Comment this line when replicate_id is available from the methods
                        if not os.path.exists(s_path): os.makedirs(s_path)
                        try:
                                script_util.download_shock_files(logger,services['shock_service_url'],s_path,files,token)
                        except Exception,e:
                                raise Exception( "Unable to download shock file, {0}".format(e))
                        try:
                                script_util.unzip_files(logger,os.path.join(s_path,a_obj['data']['file']['file_name']),s_path)
                                script_util.unzip_files(logger,os.path.join(s_path,e_obj['data']['file']['file_name']),s_path)
                                e_file_path =  os.path.join(s_path,"transcripts.gtf")
                                a_file_path = os.path.join(s_path,"accepted_hits.bam")
                                if os.path.exists(a_file_path) : print a_file_path
                                if os.path.exists(e_file_path) : 
					print e_file_path 
					list_file.write("{0}\n".format(e_file_path))
                        except Exception, e:
                                logger.exception("".join(traceback.format_exc()))
                                raise Exception("Unzip file error: Please contact help@kbase.us")
	    list_file.close()
	    print labels
            #output_dir = os.path.join(cuffdiff_dir, params['output_obj_name'])
            for l in labels:
                  #rep_files=",".join([ os.path.join(diffexp_dir+'/'+l,sub+'/accepted_hits.bam') for sub in os.listdir(os.path.join(diffexp_dir,l)) if os.path.isdir(os.path.join(diffexp_dir,l+'/'+sub))])
                  rep_files=[ (os.path.join(diffexp_dir+'/'+l,sub+'/accepted_hits.bam'), os.path.join(diffexp_dir+'/'+l,sub+'/transcripts.gtf')) for sub in os.listdir(os.path.join(diffexp_dir,l)) if os.path.isdir(os.path.join(diffexp_dir,l+'/'+sub))]
                  #itertools.chain(m_align_exp,rep_files)
		  m_align_exp += rep_files
	    print m_align_exp
	    #print "list of alignments and expression {0}".format(",".join(m_align_exp))
	    ### Get the tool_used from expressionset_obj
	    used_tool = e_sample['data']['tool_used']
	    if used_tool == 'StringTie':
		run_tool =  TOOL1_USED
		tool_version = TOOL1_VERSION
	    elif used_tool == 'Cufflinks':
		run_tool = TOOL2_USED
		tool_version = TOOL2_VERSION
            # Determine the num_threads provided by the user otherwise default the number of threads to 2
            if('num_threads' in params and params['num_threads'] is not None):
                        num_threads = int(params['num_threads'])
            else:
                        num_threads = 2
            num_cores = mp.cpu_count()
            logger.info("Number of available cores : {0}".format(num_cores))
            b_tasks =[]
            if num_cores != 1:
                pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
            else:
                pool_size = 1
                num_threads = 1
            count = 0
            logger.info(" Number of threads used by each process {0}".format(num_threads))
            for i,j in m_align_exp:
                        #try:
			print "Adding task {0} , {1} to task list".format(i,j)
                        b_tasks.append((None,services,ws_client,hs,params['ws_id'],num_threads,i,j,assembly_file,used_tool,diffexp_dir,gtf_file))
		

            @parallelize(CalldiffExpCallforBallgown_helper,pool_size)
            def run_diffexp_for_ballgown_in_parallel(tasks):
               pass
            results=run_diffexp_for_ballgown_in_parallel(b_tasks)
	    expr_file, single_ballgown_dir = results
	    #### Check if all the jobs passed
            ballgownobject_name = params['expressionset_id']+"_DifferentialExpression_Ballgown"
            ballgown_dir = os.path.join(diffexp_dir,"ballgown")
            #reportObj=script_util.create_RNASeq_ExpressionSet_and_build_report(logger,ws_client,TOOL_USED, TOOL_VERSION,tool_opts,params['ws_id'],align_names,expressionset_id,annotation_id,sampleset_id,results,expressionSet_name)
	    ### Save Ballgown differential Expression object to workspace
	    #except Exception,e:
            #    raise Exception("Error executing diffexp {0},{1}".format(cuffdiff_command,directory))

            ##  compress and upload to shock
            try:
                 logger.info("Zipping differential expression output for ballgown")
                 out_file_path = os.path.join(diffexp_dir,"{0}.zip".format(params['output_obj_name']))
                 script_util.zip_files(logger,ballgown_dir,out_file_path)
            except Exception,e:
                 raise Exception("Error zipping dir {0}".format(ballgown_dir)) 
            try:
                 handle = hs.upload(out_file_path)
            except Exception, e:
                 print " ".join(traceback.print_exc())
                 raise Exception("Failed to upload the diffexp output files: {0}".format(out_file_path))
            output_name = params['output_obj_name']
            ## Save object to workspace
            try:
                 logger.info("Saving diffexp object to workspace")
                 cm_obj = { "tool_used" : run_tool,
                            "tool_version" : tool_version,
                            "condition" : condition,
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

	    returnVal = { 'output'  : cuffdiff_obj ,'workspace' : params['ws_id']}
	    return returnVal
