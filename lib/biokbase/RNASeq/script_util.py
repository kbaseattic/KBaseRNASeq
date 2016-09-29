import sys
import os
import json
import re
import io
import urllib
import hashlib
import requests
import logging
import shutil
import time
import math
import traceback
from requests_toolbelt import MultipartEncoder
from multiprocessing import Pool
from collections import Counter
#from functools import partial
import subprocess
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join
import contig_id_mapping as c_mapping
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

from biokbase.workspace.client import Workspace
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI, GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI, AssemblyClientAPI
import datetime

def if_obj_exists(logger,ws_client,ws_id,o_type,obj_l):
    obj_list = ws_client.list_objects( {"workspaces" : [ws_id ] ,"type" : o_type,'showHidden' : 1})
    obj_names = [i[1] for i in obj_list]
    existing_names = [i for i in obj_l if i in obj_names]
    obj_ids = None
    if len(existing_names) != 0:
        e_queries = [{'name' : j , 'workspace' : ws_id } for j in existing_names]
        e_infos = ws_client.get_object_info_new({"objects": e_queries })      
	obj_ids =[ (str(k[1]) , (str(k[6]) + '/' + str(k[0]) + '/' + str(k[4])) ) for k in e_infos]
    return obj_ids

def check_and_download_existing_handle_obj(logger,ws_client,urls,ws_id,ws_object_name,ws_obj_type,directory,token):
	ret = if_obj_exists(logger,ws_client,ws_id,ws_obj_type,[ws_object_name])
        if not ret is None:
            logger.info("Object {0} exists in the workspace  {1}".format(ws_object_name,ws_id))
            obj_name,obj_id = ret[0]
            obj_info=ws_client.get_objects([{'ref' : obj_id}])[0]
            handle_id=obj_info['data']['handle']['id']
            handle_name=obj_info['data']['handle']['file_name']
            try:
                     download_file_from_shock(logger, shock_service_url=urls['shock_service_url'], shock_id=handle_id,filename=handle_name,directory=directory,token=token)
                     file_path = os.path.join(directory,handle_name)
            except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(handle_name))
	else:
	     return None
	return file_path

def if_ws_obj_exists(logger,ws_client,ws_id,o_type,obj_l):
    existing_names = None
    obj_list = ws_client.list_objects( {"workspaces" : [ws_id ] ,"type" : o_type,'showHidden' : 1})
    obj_names = [i[1] for i in obj_list]
    existing_names = [i for i in obj_l if i in obj_names]
    return existing_names

def find_read_objects(logger,ex_reads_alignments,suffix1,suffix2):
    #try:
    if len(ex_reads_alignments) != 0:
            objects = []
            for i in ex_reads_alignments:
                 objects.append(i.split(suffix1)[0]+suffix2)

            return objects
    else:
            return None
### TODO Remove this function from script_util , already moved to rnaseq_util
def get_fasta_from_genome(logger,ws_client,urls,genome_id):
    
    ref = ws_client.get_object_subset(
                                     [{ 'ref' : genome_id ,'included': ['contigset_ref']}])
    contig_id = ref[0]['data']['contigset_ref']
    logger.info( "Generating FASTA from Genome")
    try:
         ## get the FASTA
         assembly = AssemblyUtil(urls['callback_url'])
         ret = assembly.get_assembly_as_fasta({'ref':contig_id})
         output_file = ret['path']
         fasta_file = os.path.basename(output_file)
    	 return fasta_file
    except Exception, e:
	 raise Exception(e)
	 raise Exception("Unable to Create FASTA file from Genome : {0}".format(genome_id))
    return None
	
### TODO Remove this function from script_util , already moved to rnaseq_util
def generate_fasta(logger,internal_services,token,ref,output_dir,obj_name):
	try:
		ga = GenomeAnnotationAPI(internal_services,
                             token=token,
                             ref= ref)
	except Exception as e:
		raise Exception("Unable to Call GenomeAnnotationAPI : {0}".format("".join(traceback.format_exc())))
	logger.info("Generating FASTA file from Assembly for {}".format(obj_name))	
	fasta_start = datetime.datetime.utcnow()
	output_file = os.path.join(output_dir,'{}.fasta'.format(obj_name))
	fasta_file= io.open(output_file, 'wb')
    	try:
        	ga.get_assembly().get_fasta().to_file(fasta_file)
	except Exception as e:
		#raise Exception("Unable to Create FASTA file from Genome Annotation : {0}".format(obj_name))
		raise Exception("Unable to Create FASTA file from Genome Annotation : {0}".format("".join(traceback.format_exc())))
	finally:
		fasta_file.close()
    	fasta_end = datetime.datetime.utcnow()
	logger.info("Generating FASTA for {} took {}".format(obj_name, fasta_end - fasta_start))
	return output_file
	## Additional Step for sanitizing contig id
	#logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
	#mapping_filename = c_mapping.create_sanitized_contig_ids(output_file)
    	#c_mapping.replace_fasta_contig_ids(output_file, mapping_filename, to_modified=True)
	#logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))

### TODO Remove this function from script_util , already moved to rnaseq_util
def generate_gff(logger,internal_services,token,ref,output_dir,obj_name,output_file):
        try:
                ga = GenomeAnnotationAPI(internal_services,
                             token=token,
                             ref= ref)
        except:
                raise Exception("Unable to Call GenomeAnnotationAPI : {0}".format(("".join(traceback.format_exc()))))
        logger.info("Requesting GenomeAnnotation GFF for {}".format(obj_name))
    	gff_start = datetime.datetime.utcnow()
        gff_file= io.open(output_file, 'wb')
	#output_file = os.path.join(output_dir,'{}.gff'.format(obj_name))
	try:
        	ga.get_gff().to_file(gff_file)
	except Exception as e:
                #raise Exception("Unable to Create GFF  file from Genome Annotation : {0}: {1}".format(obj_name,e))
                raise Exception("Unable to Create GFF  file from Genome Annotation : {0}: {1}".format(obj_name,"".join(traceback.format_exc())))
        finally:
    		gff_file.close()
	gff_end = datetime.datetime.utcnow()
    	logger.info("Generating GFF for {} took {}".format(obj_name, gff_end - gff_start))
        ## Additional Step for sanitizing contig id
        #logger.info("Sanitizing the gff file to correct id names {}".format(datetime.datetime.utcnow()))

### TODO Remove this function from script_util , already moved to rnaseq_util
def create_gtf_annotation_from_genome(logger,ws_client,hs_client,urls,ws_id,genome_ref,genome_id,fasta_file,directory,token):
        try:
		#tmp_file = os.path.join(directory,genome_id + "_GFF.gff")
                ## get the GFF
		genome = GenomeFileUtil(urls['callback_url'])
		ret = genome.genome_to_gff({'genome_ref':genome_ref})
		file_path = ret['file_path']
		gtf_ext = ".gtf"
		if not file_path.endswith(gtf_ext): 
               		gtf_path = os.path.join(directory,genome_id+".gtf")
                	gtf_cmd = " -E {0} -T -o {1}".format(file_path,gtf_path)
                	try:
                   		logger.info("Executing: gffread {0}".format(gtf_cmd))
                   		cmdline_output = runProgram(None,"gffread",gtf_cmd,None,directory)
                	except Exception as e:
                   		raise Exception("Error Converting the GFF file to GTF using gffread {0},{1}".format(gtf_cmd,"".join(traceback.format_exc())))
		else:
			gtf_path = file_path
                if os.path.exists(gtf_path):
                               annotation_handle = hs_client.upload(gtf_path)
                               a_handle = { "handle" : annotation_handle ,"size" : os.path.getsize(gtf_path), 'genome_id' : genome_ref}
                ##Saving GFF/GTF annotation to the workspace
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.GFFAnnotation",
                                         "data":a_handle,
                                         "name":genome_id+"_GTF_Annotation",
                                        "hidden":1}
                                        ]})
        except Exception as e:
                raise ValueError("Generating GTF file from Genome Annotation object Failed :  {}".format("".join(traceback.format_exc())))
	return gtf_path

### TODO Remove this function from script_util , already moved to rnaseq_util
def create_gtf_annotation(logger,ws_client,hs_client,internal_services,ws_id,genome_ref,genome_id,fasta_file,directory,token):
        try:
		tmp_file = os.path.join(directory,genome_id + "_GFF.gff")
        	fasta_file= generate_fasta(logger,internal_services,token,genome_ref,directory,genome_id)
            	logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
                mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
                c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
                logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
                generate_gff(logger,internal_services,token,genome_ref,directory,genome_id,tmp_file)
                c_mapping.replace_gff_contig_ids(tmp_file, mapping_filename, to_modified=True)
                gtf_path = os.path.join(directory,genome_id+"_GTF.gtf")
                gtf_cmd = " -E {0} -T -o {1}".format(tmp_file,gtf_path)
                try:
                   logger.info("Executing: gffread {0}".format(gtf_cmd))
                   cmdline_output = runProgram(None,"gffread",gtf_cmd,None,directory)
                except Exception as e:
                   raise Exception("Error Converting the GFF file to GTF using gffread {0},{1}".format(gtf_cmd,"".join(traceback.format_exc())))
		#if os.path.exists(tmp_file): os.remove(tmp_file)
                if os.path.exists(gtf_path):
                               annotation_handle = hs_client.upload(gtf_path)
                               a_handle = { "handle" : annotation_handle ,"size" : os.path.getsize(gtf_path), 'genome_id' : genome_ref}
                ##Saving GFF/GTF annotation to the workspace
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.GFFAnnotation",
                                         "data":a_handle,
                                         "name":genome_id+"_GTF_Annotation",
                                        "hidden":1}
                                        ]})
        except Exception as e:
                raise ValueError("Generating GTF file from Genome Annotation object Failed :  {}".format("".join(traceback.format_exc())))
	return gtf_path
	

### TODO Remove this function from script_util , already moved to rnaseq_util
def create_RNASeq_AlignmentSet_and_build_report(logger,ws_client,ws_id,sample_list,sampleset_id,genome_id,bowtie2index_id,results,alignmentSet_name):
	 results =  [ ret for ret in results if not ret is None ]
	 if len(results) < 2:
	  	raise ValueError("Not enough alignments got created for a AlignmentSet obj")
	 set_obj = { 'sampleset_id' :sampleset_id ,'genome_id' : genome_id}
	 if not bowtie2index_id is None:
		set_obj['bowtie2_index'] = bowtie2index_id
         sids=[]
         m_alignments = []
         alignments = []
	 m_align_names = []
	 output_objs = []
	 num_samples = len(sample_list)
	 num_results = len(results)
	 num_failed = num_samples - num_results
	 run_list = [ k for (k,v) in results ]
	 print run_list
	 failed_list = [k for k in sample_list if k not in run_list ]
	 print  "\n".join(failed_list)
         for sid,s_alignments in results:
                    a_ref = ws_client.get_object_info_new({"objects": [{'name':s_alignments, 'workspace': ws_id}]})[0]
                    a_id = str(a_ref[6]) + '/' + str(a_ref[0]) + '/' + str(a_ref[4])
                    m_alignments.append({sid : a_id})
                    m_align_names.append({sid : s_alignments})
                    output_objs.append({'ref' : a_id , 'description': "RNA-seq Alignment for reads Sample :  {0}".format(sid)})
                    sids.append(sid)
                    alignments.append(a_id)
         set_obj['read_sample_ids']= sids
         set_obj['sample_alignments']= alignments
         set_obj['mapped_alignments_ids']=m_alignments
	 set_obj['mapped_rnaseq_alignments'] = m_align_names
         try:
        	logger.info( "Saving AlignmentSet object to  workspace")
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignmentSet",
                                         "data":set_obj,
                                         "name":alignmentSet_name}
                                        ]})[0]
                                                                
                output_objs.append({'ref': str(res[6]) + '/' + str(res[0]) + '/' + str(res[4]),'description' : "Set of Alignments for Sampleset : {0}".format(sampleset_id)})
	 except Exception as e:
                    logger.exception(e)
                    raise Exception("Failed Saving AlignmentSet to Workspace") 
	 ### Build Report obj ###
	 report = []
	 report.append("Total number of reads : {0}".format(str(num_samples)))
	 report.append("Number of reads ran successfully : {0}".format(str(num_results)))
	 report.append("Number of reads failed during this run : {0}".format(str(num_failed))) 
	 if len(failed_list) != 0:
		report.append("List of reads failed in this run : {0}".format("\n".join(failed_list)))
	 reportObj = {
                      'objects_created':output_objs,
                      'text_message':'\n'.join(report)
                     }
	 return reportObj

### TODO Remove this function from script_util , already moved to rnaseq_util
def create_RNASeq_ExpressionSet_and_build_report(logger,ws_client,tool_used, tool_version,tool_opts,ws_id,alignment_list,alignmentset_id,genome_id,sampleset_id,results,expressionSet_name):
	 results =  [ ret for ret in results if not ret is None ]
	 if len(results) < 2:
	  	raise ValueError("Not enough expression results to create a ExpressionSet object")
	 set_obj = { 'tool_used': tool_used, 'tool_version': tool_version,'alignmentSet_id' : alignmentset_id ,'genome_id' : genome_id,'sampleset_id' : sampleset_id }
	 if not tool_opts is None:
		set_obj['tool_opts'] = tool_opts
         sids=[]
         condition = []
	 expr_ids = []
         m_expr_names= []
	 m_expr_ids = []
	 output_objs = []
	 num_samples = len(alignment_list)
	 num_results = len(results)
	 num_failed = num_samples - num_results
	 run_list = [ k for (k,v) in results ]
	 failed_list = [k for k in alignment_list if k not in run_list ]
         for a_name, e_name in results:
                    a_ref,e_ref = ws_client.get_object_info_new({"objects": [{'name':a_name, 'workspace': ws_id},{'name':e_name, 'workspace': ws_id}]})
                    a_id = str(a_ref[6]) + '/' + str(a_ref[0]) + '/' + str(a_ref[4])
                    e_id = str(e_ref[6]) + '/' + str(e_ref[0]) + '/' + str(e_ref[4])
                    m_expr_ids.append({a_id : e_id})
                    m_expr_names.append({a_name : e_name})
                    output_objs.append({'ref' : e_id , 'description': "RNA-seq Alignment for reads Sample :  {0}".format(a_name)})
                    expr_ids.append(e_id)
         set_obj['sample_expression_ids']= expr_ids
         set_obj['mapped_expression_objects']= m_expr_names
         set_obj['mapped_expression_ids'] = m_expr_ids
         try:
        	logger.info( "Saving AlignmentSet object to  workspace")
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqExpressionSet",
                                         "data":set_obj,
                                         "name":expressionSet_name}
                                        ]})[0]
                                                                
                output_objs.append({'ref': str(res[6]) + '/' + str(res[0]) + '/' + str(res[4]),'description' : "Set of Expression objects for AlignmentSet : {0}".format(alignmentset_id)})
	 except Exception as e:
		    logger.exception("".join(traceback.format_exc()))
                    raise Exception("Failed Saving ExpressionSet to Workspace") 
	 ### Build Report obj ###
	 report = []
	 report.append("Total number of alignments given : {0}".format(str(num_samples)))
	 report.append("Number of assemblies ran successfully : {0}".format(str(num_results)))
	 report.append("Number of  assemblies failed during this run : {0}".format(str(num_failed))) 
	 if len(failed_list) != 0:
		report.append("List of reads failed in this run : {0}".format("\n".join(failed_list)))
	 reportObj = {
                      'objects_created':output_objs,
                      'text_message':'\n'.join(report)
                     }
	 return reportObj

		   	
### TODO Remove this function from script_util , already moved to rnaseq_util
def updateAnalysisTO(logger, ws_client, field, map_key, map_value, anal_ref, ws_id, objid):
    
        analysis = ws_client.get_objects([{'ref' : anal_ref}])[0]
        
        if field in analysis['data'] and analysis['data'][field] is not None:
                analysis['data'][field][map_key] = map_value
        else:
            analysis['data'][field] = {map_key : map_value}
	logger.info("Analysis object updated {0}".format(json.dumps(analysis['data'])))
        res1= ws_client.save_objects(
                            {"workspace":ws_id,
                             "objects": [{
                             "type":"KBaseRNASeq.RNASeqAnalysis",
                             "data":analysis['data'],
                             "objid":objid}
                            ]})


### TODO Remove this function from script_util , already moved to rnaseq_util
def extractStatsInfo(logger,ws_client,ws_id,sample_id,result,stats_obj_name):
	lines = result.splitlines()
        if  len(lines) != 11:
            raise Exception("Error not getting enough samtool flagstat information: {0}".format(result))
        # patterns
        two_nums  = re.compile(r'^(\d+) \+ (\d+)')
        two_pcts  = re.compile(r'\(([0-9.na\-]+)%:([0-9.na\-]+)%\)')
        # alignment rate
        m = two_nums.match(lines[0])
        total_qcpr = int(m.group(1))
        total_qcfr = int(m.group(2))
        total_read =  total_qcpr + total_qcfr
        m = two_nums.match(lines[2])
        mapped_r = int(m.group(1))
        unmapped_r = int(total_read - mapped_r)
	alignment_rate = float(mapped_r) / float(total_read)  * 100.0
        if alignment_rate > 100: alignment_rate = 100.0

        # singletons
        m = two_nums.match(lines[8])
        singletons = int(m.group(1))
	m = two_nums.match(lines[6])
        properly_paired = int(m.group(1))
        # Create Workspace object
        stats_data =  {
                       "alignment_id": sample_id,
                       "alignment_rate": alignment_rate,
                       #"multiple_alignments": 50, 
                       "properly_paired": properly_paired,
                       "singletons": singletons,
                       "total_reads": total_read,
                       "unmapped_reads": unmapped_r,
                       "mapped_reads": mapped_r
                       }
	logger.info(json.dumps(stats_data))
        ## Save object to workspace
        logger.info( "Saving Alignment Statistics to the Workspace")
        try:
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.AlignmentStatsResults",
                                         "data": stats_data,
				         "hidden" : 1,
                                         "name":stats_obj_name}
                                        ]})
                res = stats_data
        except Exception, e:
                raise Exception("get Alignment Statistics failed: {0}".format(e))

### TODO Remove this function from script_util , already moved to rnaseq_util
def extractAlignmentStatsInfo(logger,tool_used,ws_client,ws_id,sample_id,result,stats_obj_name):
        lines = result.splitlines()
	if tool_used == 'samtools':
        	if  len(lines) != 11:
            		raise Exception("Error not getting enough samtool flagstat information: {0}".format(result))
        	# patterns
        	two_nums  = re.compile(r'^(\d+) \+ (\d+)')
        	two_pcts  = re.compile(r'\(([0-9.na\-]+)%:([0-9.na\-]+)%\)')
        	# alignment rate
        	m = two_nums.match(lines[0])
        	total_qcpr = int(m.group(1))
        	total_qcfr = int(m.group(2))
        	total_read =  total_qcpr + total_qcfr
        	m = two_nums.match(lines[2])
        	mapped_r = int(m.group(1))
        	unmapped_r = int(total_read - mapped_r)
        	alignment_rate = float(mapped_r) / float(total_read)  * 100.0
        	if alignment_rate > 100: alignment_rate = 100.0

        	# singletons
       		m = two_nums.match(lines[8])
        	singletons = int(m.group(1))
        	m = two_nums.match(lines[6])
        	properly_paired = int(m.group(1))
		multiple_alignments = 0
	elif tool_used == 'bowtie2':
		if len(lines) not in [6,15]:
                        raise Exception("Error not getting enough bowtie2 alignment information: {0}".format(result))
		pattern1 = re.compile(r'^(\s*\d+)')
	        pattern2 = re.compile(r'^(\s*\d+.\d+)')	
		m =  pattern1.match(lines[0])
		total_read = int(m.group(1))
		m = pattern1.match(lines[2])
		unmapped_r =  int(m.group(1))
		mapped_r = total_read - unmapped_r
		m = pattern1.match(lines[4])
		multiple_alignments = int(m.group(1))
		if len(lines) == 6:
			m = pattern2.match(lines[5])
			alignment_rate = float(m.group(1))
			singletons = 0
			properly_paired = 0
		if len(lines) == 15:
			m =pattern1.match(lines[1])
			properly_paired = int(m.group(1))
			singletons = total_read - properly_paired	
			m = pattern2.match(lines[14])
			alignment_rate = float(m.group(1))
	elif tool_used == 'tophat':
		pass 
        # Create Workspace object
        stats_data =  {
                       #"alignment_id": sample_id,
                       "alignment_rate": alignment_rate,
                       "multiple_alignments": multiple_alignments, 
                       "properly_paired": properly_paired,
                       "singletons": singletons,
                       "total_reads": total_read,
                       "unmapped_reads": unmapped_r,
                       "mapped_reads": mapped_r
                       }
	#print stats_data
	return stats_data
        ## Save object to workspace
        #logger.info( "Saving Alignment Statistics to the Workspace")
#        try:
#                res= ws_client.save_objects(
#                                        {"workspace":ws_id,
#                                         "objects": [{
#                                         "type":"KBaseRNASeq.AlignmentStatsResults",
#                                         "data": stats_data,
#                                         "hidden" : 1,
#                                         "name":stats_obj_name}
#                                        ]})
#                res = stats_data
#        except Exception, e:
#                raise Exception("get Alignment Statistics failed: {0}".format(e))

def getExpressionHistogram(obj,obj_name,num_of_bins,ws_id,output_obj_name):
    if 'expression_levels' in obj['data']:
        hdict = obj['data']['expression_levels']
        tot_genes =  len(hdict)
        lmin = round(min([v for k,v in hdict.items()]))
        lmax = round(max([v for k,v in hdict.items()]))
        hist_dt = script_util.histogram(hdict.values(),lmin,lmax,int(num_of_bins))
        title = "Histogram  - " + obj_name
        hist_json = {"title" :  title , "x_label" : "Gene Expression Level (FPKM)", "y_label" : "Number of Genes", "data" : hist_dt}
        sorted_dt = OrderedDict({ "id" : "", "name" : "","row_ids" : [] ,"column_ids" : [] ,"row_labels" : [] ,"column_labels" : [] , "data" : [] })
        sorted_dt["row_ids"] = [hist_json["x_label"]]
        sorted_dt["column_ids"] = [hist_json["y_label"]]
        sorted_dt['row_labels'] = [hist_json["x_label"]]
        sorted_dt["column_labels"] =  [hist_json["y_label"]]
        sorted_dt["data"] = [[float(i) for i in hist_json["data"]["x_axis"]],[float(j) for j in hist_json["data"]["y_axis"]]]
        #sorted_dt["id"] = "kb|histogramdatatable."+str(idc.allocate_id_range("kb|histogramdatatable",1))
        sorted_dt["id"] = output_obj_name
        sorted_dt["name"] = hist_json["title"]
        res = ws_client.save_objects({"workspace": ws_id,
                                     "objects": [{
                                     "type":"MAK.FloatDataTable",
                                     "data": sorted_dt,
                                     "name" : output_obj_name}
                                     ]
                                     })
		

def stderrlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stderr handler attached and using a prefix
    format that will make logging consistent between scripts.
    """
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    
    return logger


def stdoutlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stdout handler attached and using a prefix
    format that will make logging consistent between scripts.
    """
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    
    return logger



def zip_files(logger, src_path, output_fn):
    """
    Compress all index files (not directory) into an output zip file on disk.
    """

    files = [ f for f in listdir(src_path) if isfile(join(src_path,f)) ]
    with ZipFile(output_fn, 'w', allowZip64=True) as izip:
        for f in files:
            izip.write(join(src_path,f),f)

def unzip_files(logger, src_fn, dst_path):
    """
    Extract all index files into an output zip file on disk.
    """

    with ZipFile(src_fn, 'r') as ozip:
        ozip.extractall(dst_path)

def move_files(logger, src, dest):
    """
    Move files from one folder to another.
    """
    
    src_files = os.listdir(src)
    for file_name in src_files:
       full_file_name = os.path.join(src, file_name)
       if (os.path.isfile(full_file_name)):
          shutil.copy(full_file_name, dest)

def download_file_from_shock(logger,
                             shock_service_url = None,
                             shock_id = None,
                             filename = None,
                             directory = None,
			     filesize= None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)
    #logger.info("Downloading shock node {0}/node/{1}".format(shock_service_url,shock_id))

    metadata_response = requests.get("{0}/node/{1}?verbosity=metadata".format(shock_service_url, shock_id), headers=header, stream=True, verify=True)
    shock_metadata = metadata_response.json()['data']
    if shock_metadata is not None:
    	shockFileName = shock_metadata['file']['name']
    	shockFileSize = shock_metadata['file']['size']
    metadata_response.close()
        
    download_url = "{0}/node/{1}?download_raw".format(shock_service_url, shock_id)
    try: 
    	data = requests.get(download_url, headers=header, stream=True, verify=True)
    except Exception,e:
	print(traceback.format_exc())
    if filename is not None:
        shockFileName = filename

    if directory is not None:
        filePath = os.path.join(directory, shockFileName)
    else:
        filePath = shockFileName

    if filesize is not None:
	shockFileSize = filesize

    chunkSize = shockFileSize/4
    
    maxChunkSize = 2**30
    
    if chunkSize > maxChunkSize:
        chunkSize = maxChunkSize
    
    f = io.open(filePath, 'wb')
    try:
        for chunk in data.iter_content(chunkSize):
            f.write(chunk)
    finally:
        data.close()
        f.close()      

def download_shock_files(logger,shock_url,directory,dict_of_files,token):
	for name, fid in dict_of_files.items():
		try:
			download_file_from_shock(logger, shock_service_url=shock_url, shock_id=fid,filename=name, directory=directory,token=token) 
	   	except Exception,e:
                        raise Exception( "Unable to download shock file , {0},{1}".format(name,fid))

def query_shock_node(logger,
                             shock_service_url = None,
                             condition = None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    query_str = urllib.urlencode(condition)
    
    logger.info("Querying shock node {0}/node/?query&{1}".format(shock_service_url,query_str))

    query_response = requests.get("{0}/node/?query&{1}".format(shock_service_url, query_str), headers=header, stream=True, verify=True)
    query_rst = query_response.json()['data']
    query_response.close()
    return query_rst


def upload_file_to_shock(logger,
                         shock_service_url = None,
                         filePath = None,
                         attributes = '{}',
                         ssl_verify = True,
                         token = None):
    """
    Use HTTP multi-part POST to save a file to a SHOCK instance.
    """

    if token is None:
        raise Exception("Authentication token required!")
    
    #build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(os.path.abspath(filePath), 'r')
    m = MultipartEncoder(fields={'attributes_str': json.dumps(attributes), 'upload': (os.path.split(filePath)[-1], dataFile)})
    header['Content-Type'] = m.content_type

    logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
    try:
        response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
        dataFile.close()
    except:
        dataFile.close()
        raise    

    if not response.ok:
        response.raise_for_status()

    result = response.json()

    if result['error']:            
        raise Exception(result['error'][0])
    else:
        return result["data"]    

def shock_node_2b_public(logger,
                         node_id = None,
                         shock_service_url = None,
                         ssl_verify = True,
                         token = None):
    """
    Ensure a node to be public
    """

    if token is None:
        raise Exception("Authentication token required!")

    if shock_service_url is None:
        raise Exception("Shock URL is required!")
    
    if node_id is None:
        raise Exception("Node ID is required!")

    #build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)


    logger.info("-X PUT {0}/node/{1}/acl/public_read".format(shock_service_url, node_id))
    try:
        response = requests.put("{0}/node/{1}/acl/public_read".format(shock_service_url, node_id), headers=header, allow_redirects=True, verify=ssl_verify)
    except Exception,e:
        raise Exception("Error making Shock ids Public{0}".format(e))   

    if not response.ok:
        response.raise_for_status()

    result = response.json()

    if result['error']:            
        raise Exception(result['error'][0])
    else:
        return result["data"]    


def getHandles(logger = None,
               shock_service_url = None,
               handle_service_url = None,
               shock_ids = None,
               handle_ids = None,
               token = None):
    """
    Retrieve KBase handles for a list of shock ids or a list of handle ids.
    """

    if token is None:
        raise Exception("Authentication token required!")

    hs = HandleService(url=handle_service_url, token=token)

    handles = list()
    if shock_ids is not None:
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        for sid in shock_ids:
            info = None

            try:
                logger.info("Found shock id {0}, retrieving information about the data.".format(sid))

                response = requests.get("{0}/node/{1}".format(shock_service_url, sid), headers=header, verify=True)
                info = response.json()["data"]
            except:
                logger.error("There was an error retrieving information about the shock node id {0} from url {1}".format(sid, shock_service_url))

            try:
                logger.info("Retrieving a handle id for the data.")
                handle = hs.persist_handle({"id" : sid,
                                           "type" : "shock",
                                           "url" : shock_service_url,
                                           "file_name": info["file"]["name"],
                                           "remote_md5": info["file"]["checksum"]["md5"]})
                handles.append(handle)
            except:
                try:
                    handle_id = hs.ids_to_handles([sid])[0]["hid"]
		    single_handle = hs.hids_to_handles([handle_id])

                    assert len(single_handle) != 0

                    if info is not None:
                        single_handle[0]["file_name"] = info["file"]["name"]
                        single_handle[0]["remote_md5"] = info["file"]["checksum"]["md5"]
                        logger.debug(single_handle)

                    handles.append(single_handle[0])
                except:
                    logger.error("The input shock node id {} is already registered or could not be registered".format(sid))

                    hs = HandleService(url=handle_service_url, token=token)
                    all_handles = hs.list_handles()

                    for x in all_handles:
                        if x[0] == sid:
                            logger.info("FOUND shock id as existing handle")
                            logger.info(x)
                            break
                    else:
                        logger.info("Unable to find a handle containing shock id")

                        logger.info("Trying again to get a handle id for the data.")
                        handle_id = hs.persist_handle({"id" : sid,
                                           "type" : "shock",
                                           "url" : shock_service_url,
                                           "file_name": info["file"]["name"],
                                           "remote_md5": info["file"]["checksum"]["md5"]})
                        handles.append(handle_id)

                    raise
    elif handle_ids is not None:
        for hid in handle_ids:
            try:
                single_handle = hs.hids_to_handles([hid])

                assert len(single_handle) != 0

                handles.append(single_handle[0])
            except:
                logger.error("Invalid handle id {0}".format(hid))
                raise
    return handles

def get_obj_info(logger,ws_url,objects,ws_id,token):
    """
    function to get the workspace object id from a object name
    """
    ret = []
    ws_client=Workspace(url=ws_url, token=token)
    for obj in  objects:
    	try:
            obj_infos = ws_client.get_object_info_new({"objects": [{'name': obj, 'workspace': ws_id}]})
            ret.append("{0}/{1}/{2}".format(obj_infos[0][6],obj_infos[0][0],obj_infos[0][4]))
        except Exception, e:
                     logger.error("Couldn't retrieve %s:%s from the workspace , %s " %(ws_id,obj,e))
    return ret

def whereis(program):
    """
    returns path of program if it exists in your ``$PATH`` variable or ``None`` otherwise
    """
    for path in os.environ.get('PATH', '').split(':'):
    	if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None

def runProgram(logger=None,
	       progName=None,
	       argStr=None,
	       script_dir=None,
	       working_dir=None):
        """
        Convenience func to handle calling and monitoring output of external programs.
    
        :param progName: name of system program command
        :param argStr: string containing command line options for ``progName``
    
        :returns: subprocess.communicate object
        """

        # Ensure program is callable.
        if script_dir is not None:
                progPath= os.path.join(script_dir,progName)
        else:
		progPath = progName
	#	progPath = whereis(progName)
        #       	if not progPath:
        #            raise RuntimeError(None,'{0} command not found in your PATH environmental variable. {1}'.format(progName,os.environ.get('PATH', '')))

        # Construct shell command
        cmdStr = "%s %s" % (progPath,argStr)
	print "Executing : "+cmdStr
	if logger is not None:
		logger.info("Executing : "+cmdStr)
        #if working_dir is None:
        #    logger.info("Executing: " + cmdStr + " on cwd")
        #else:
        #    logger.info("Executing: " + cmdStr + " on " + working_dir)

        # Set up process obj
        process = subprocess.Popen(cmdStr,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=working_dir)
        # Get results
        result,stderr  = process.communicate()
     	#print result
	#print stderr 
        # keep this until your code is stable for easier debugging
        if logger is not None and result is not None and len(result) > 0:
            logger.info(result)
	    print result
        if logger is not None and stderr is not None and len(stderr) > 0:
            logger.info(stderr)
	    print stderr

        # Check returncode for success/failure
        if process.returncode != 0:
		raise Exception("Command execution failed  {0}".format("".join(traceback.format_exc())))
                raise RuntimeError('Return Code : {0} , result {1} , progName {2}'.format(process.returncode,result,progName))

        # Return result
        return { "result" : result , "stderr" :stderr }

def hashfile(filepath):
       sha1 = hashlib.sha1()
       f = open(filepath, 'rb')
       try:
         sha1.update(f.read())
       finally:
         f.close()
       return sha1.hexdigest()


def create_shock_handle(logger=None,
			file_name=None,
			shock_url=None,
			handle_url=None,
			obj_type=None,
			token=None):
     	
        hs = HandleService(url=handle_url, token=token)
        f_shock = upload_file_to_shock(logger,shock_url,file_name,'{}',True,token)
        f_sha1 =  hashfile(file_name)
        hid = getHandles(logger,shock_url,handle_url,[f_shock['id']],None,token)[0]
        handle = { 'hid' : hid , 
	           "file_name" : f_shock['file']['name'] , 
                   "id" : f_shock['id'] , 
                   "type" : obj_type ,
	           "url" : shock_url,
	           "remote_md5" : f_shock['file']['checksum']['md5'],
	           "remote_sha1" : f_sha1 }	
   
        return handle

def parallel_function(f):
    def easy_parallize(f, sequence):
        pool = Pool(processes=8)
        # f is given sequence. guaranteed to be in order
        result = pool.map(f, sequence)
        cleaned = [x for x in result if not x is None]
        cleaned = asarray(cleaned)
        # not optimal but safe
        pool.close()
        pool.join()
        return cleaned
    from functools import partial
    # this assumes f has one argument, fairly easy with Python's global scope
    return partial(easy_parallize, f)


def histogram(iterable, low, high, bins):
    '''Count elements from the iterable into evenly spaced bins
        >>> scores = [82, 85, 90, 91, 70, 87, 45]
        >>> histogram(scores, 0, 100, 10)
        [0, 0, 0, 0, 1, 0, 0, 1, 3, 2]
    '''
    step = (high - low + 0.0) / bins
    ranges = range(int(round(low)),int(round(high)),int(round(step)))
    dist = Counter((float(x) - low) // step for x in iterable)
    return { "x_axis" : ranges , "y_axis" : [dist[b] for b in range(bins)] }



def parse_FPKMtracking(filename,tool,metric):
    result={}
    pos1= 0
    if tool == 'StringTie':
	if metric == 'FPKM': pos2 = 7
	if metric == 'TPM': pos2 = 8
    if tool == 'Cufflinks':
	pos2 = 9
    with open(filename) as f:
	next(f)
    	for line in f:
		larr = line.split("\t")
		if larr[pos1] != "":
			result[larr[pos1]] = math.log(float(larr[pos2])+1,2)
    return result

def get_end(start,leng,strand):
    stop = 0
    if strand == '+': 
	stop = start + ( leng - 1 )
    if strand == '-':
	stop = start - ( leng + 1)
    return stop
    
