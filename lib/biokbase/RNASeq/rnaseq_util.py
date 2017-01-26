import os
import re
import io
import csv
import math
import sys
import traceback
import datetime
import contig_id_mapping as c_mapping
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
import doekbase.data_api
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI, GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI, AssemblyClientAPI
from biokbase.RNASeq import handler_utils as handler_util
from biokbase.RNASeq import script_util
from pprint import pprint,pformat
from operator import itemgetter
import subprocess

def get_reads(logger, set_obj):
    if 'sample_ids' in set_obj['data']: # RNASeqSampleSet
        return set_obj['data']['sample_ids']
    if 'items' in set_obj['data']: # ReadsSet
        return [item['ref'] for item in set_obj['data']['items']]
    else:
        logger.error("Please use RNASeqSampleSet or ReadsSet :%s" % pformat(set_obj))
        raise Exception("Please use RNASeqSampleSet or ReadsSet")


def get_reads_conditions(logger, set_obj, set_type):
    if set_type == 'KBaseRNASeq.RNASeqSampleSet':
        reads = set_obj['data']['sample_ids']
        r_label = set_obj['data']['condition']
        return (reads, r_label)
    if set_type == 'KBaseSets.ReadsSet':
        reads = [item['ref'] for item in set_obj['data']['items']]
        r_label = [item['label'] for item in set_obj['data']['items']]
        return (reads, r_label)
    
    logger.error("Set object type is not one of  RNASeqSampleSet or ReadsSet : %s" % set_type) 
    raise Exception("Invalid set object type : %s" % set_type)
    

def get_fa_from_genome(logger,ws_client,urls,ws_id,directory,genome_name):
    #ref_info = ws_client.get_object_info_new({"objects": [{'name': genome_name, 'workspace': ws_id}]})[0]
    #ref_info = script_util.ws_get_object_info(logger, ws_client, ws_id, genome_name)[0]
    #genome_id = str(ref_info[6]) + '/' + str(ref_info[0]) + '/' + str(ref_info[4])
    genome_id = script_util.ws_get_ref(logger, ws_client, ws_id, genome_name)
    fasta_file =  os.path.join(directory,genome_name+ ".fa")
    ref = ws_client.get_object_subset(
                                     [{ 'ref' : genome_id ,'included': ['contigset_ref','assembly_ref']}])
    if 'contigset_ref' in ref[0]['data']:
            contig_id = ref[0]['data']['contigset_ref']
    elif 'assembly_ref' in ref[0]['data']:
        contig_id = ref[0]['data']['assembly_ref']
    if contig_id is None:
        raise ValueError("Genome {0} object does not have reference to the assembly object".format(genome_name))
    print contig_id
    logger.info( "Generating FASTA from Genome")
    try:
         ## get the FASTA
         assembly = AssemblyUtil(urls['callback_url'])
         ret = assembly.get_assembly_as_fasta({'ref':contig_id})
         output_file = ret['path']
         os.rename(output_file,fasta_file)
         logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
         mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
         c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
         logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
         #fasta_file = os.path.basename(fasta_file)
             #return (genome_id, fasta_file)
    except Exception, e:
         logger.exception(e)
         raise Exception("Unable to Create FASTA file from Genome : {0}".format(genome_name))
    finally:
         if os.path.exists(fasta_file): #os.remove(output_file)
             temp_fa = os.path.join(directory,handler_util.get_file_with_suffix(directory,"_temp.fa")+"._temp.fa")
             logger.debug (temp_fa)
             if os.path.exists(temp_fa): os.remove(temp_fa)
          
             return (genome_id, fasta_file)
    return None

def create_gtf_annotation_from_genome(logger,ws_client,hs_client,urls,ws_id,genome_ref,genome_name,directory,token):
    ref = ws_client.get_object_subset(
                                     [{ 'ref' : genome_ref ,'included': ['contigset_ref', 'assembly_ref']}])        
    if 'contigset_ref' in ref[0]['data']:
        contig_id = ref[0]['data']['contigset_ref']
    elif 'assembly_ref' in ref[0]['data']:
        contig_id = ref[0]['data']['assembly_ref']
    if contig_id is None:
        raise ValueError("Genome {0} object does not have reference to the assembly object".format(genome_name))
    print contig_id
    logger.info( "Generating GFF file from Genome")
    try:
                assembly = AssemblyUtil(urls['callback_url'])
                ret = assembly.get_assembly_as_fasta({'ref':contig_id})
                output_file = ret['path']
                mapping_filename = c_mapping.create_sanitized_contig_ids(output_file)
                os.remove(output_file)
                ## get the GFF
                genome = GenomeFileUtil(urls['callback_url'], service_ver="dev")
                ret = genome.genome_to_gff({'genome_ref':genome_ref})
                file_path = ret['file_path']
                c_mapping.replace_gff_contig_ids(file_path, mapping_filename, to_modified=True)
                gtf_ext = ".gtf"
                if not file_path.endswith(gtf_ext): 
                        gtf_path = os.path.join(directory,genome_name+".gtf")
                        gtf_cmd = " -E {0} -T -o {1}".format(file_path,gtf_path)
                        try:
                                   logger.info("Executing: gffread {0}".format(gtf_cmd))
                                   cmdline_output = script_util.runProgram(None,"gffread",gtf_cmd,None,directory)
                        except Exception as e:
                                   raise Exception("Error Converting the GFF file to GTF using gffread {0},{1}".format(gtf_cmd,"".join(traceback.format_exc())))
                else:
                        logger.info("GTF handled by GAU")
                        gtf_path = file_path
                logger.info("gtf file : " + gtf_path)
                if os.path.exists(gtf_path):
                               annotation_handle = hs_client.upload(gtf_path)
                               a_handle = { "handle" : annotation_handle ,"size" : os.path.getsize(gtf_path), 'genome_id' : genome_ref}
                ##Saving GFF/GTF annotation to the workspace
                res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.GFFAnnotation",
                                         "data":a_handle,
                                         "name":genome_name+"_GTF_Annotation",
                                         "hidden":1}
                                        ]})
    except Exception as e:
                raise ValueError("Generating GTF file from Genome Annotation object Failed :  {}".format("".join(traceback.format_exc())))
    return gtf_path

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
                    #a_ref = ws_client.get_object_info_new({"objects": [{'name':s_alignments, 'workspace': ws_id}]})[0]
                    #a_id = str(a_ref[6]) + '/' + str(a_ref[0]) + '/' + str(a_ref[4])
                    a_id = script_util.ws_get_ref(logger, ws_client, ws_id, s_alignments)
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


def store_table2exp_mat(logger, ws_client, tables, genome_ref, set_names, em_name, ws_id, hidden = 0):
    all_rows = {}    # build a dictionary of keys only which is a union of all row ids (gene_ids)
    #logger.info( "***** length of tables is {0}".format( len( tables )))
    for table in tables:
        for r in table.keys():
            all_rows[r] = []

    for gene_id in all_rows.keys():
        row = []
        for table in tables:
            if ( gene_id in table ):
                #logger.info( "append " + gene_id )
                #logger.info( pformat( table[gene_id]))
                           #all_rows[gene_id].append( table[gene_id] )
                row.append( table[gene_id] )
            else:
                #logger.info( "append  0" )
                row.append( 0 )
            all_rows[gene_id] = row
            #logger.info( all_rows[gene_id])

    emo= {
          "genome_ref" : genome_ref,
          "scale" : "log2",
          "type" : "level",
          "data" : {
                    "row_ids" : [],
                    "values" : [],
                    "col_ids" : set_names
                   },
          "feature_mapping" : {}
         }

    # we need to load row-by-row to preserve the order
    logger.info( "loading emo")

    for gene_id in all_rows.keys():
        emo["feature_mapping"][gene_id] = gene_id
        emo["data"]["row_ids"].append( gene_id )
        emo["data"]["values"].append( all_rows[gene_id] )
        emo["feature_mapping"][gene_id] = gene_id   # QUESTION: What to do here?

    try:
        logger.info( "saving emo em_name {0}".format( em_name ))
        ret = ws_client.save_objects( { 'workspace' : ws_id,
                                        'objects' : [
                                                      { 'type' : 'KBaseFeatureValues.ExpressionMatrix',
                                                        'data' : emo,
                                                        'name' : em_name,
                                                        'hidden' : hidden
                                                      }
                                                    ]
                                      }
                                    )
        logger.info( "ws save return:\n" + pformat(ret))
    except Exception as e:
        logger.exception(e)
        raise Exception( "Failed Saving Expression Matrix to Workspace" ) 

def create_and_load_expression_matrix( logger, ws_client, ws_id, expressionSet_name, genome_ref, hide_TPM ):

        try:
            logger.info( "*********getting expression set {0} from workspace*******".format( expressionSet_name ))
            expr_set = ws_client.get_objects( [ { 'name' : expressionSet_name, 'workspace': ws_id } ] )[0]
            logger.info( pformat( expr_set ))
        except Exception,e:
            logger.exception(e)
            raise Exception( "Unable to download expression set {0} from workspace {1}".format( expressionSet_name, ws_id ))

        eo_names = []
        for meo in expr_set["data"]["mapped_expression_objects"]: # better with ids
            eo_names.append( meo.values()[0] )

        fpkm_tables = []
        tpm_tables = []
        set_names = []
        tpm_table = None
        for eo_name in eo_names:
            try:
                logger.info( "*********getting expression set {0} from workspace*******".format( eo_name ))
                expr = ws_client.get_objects( [ { 'name' : eo_name, 'workspace': ws_id } ] )[0]
            except Exception,e:
                logger.exception(e)
                raise Exception( "Unable to download expression object {0} from workspace {1}".format( eo_name, ws_id ) )
            set_names.append( eo_name )

            num_interp = expr["data"]["numerical_interpretation"]
            if ( num_interp != "FPKM" ):
                raise Exception( "Did not get expected FPKM value from numerical interpretation key from Expression object {0}, instead got ".format(eo_name, num_interp) )

            pr_comments = expr["data"]["processing_comments"]     # log2 Normalized
            logger.info( "pr_comments are {0}".format( pr_comments ))

            fpkm_table = expr["data"]["expression_levels"]    # QUESTION: is this really FPKM levels?
            logger.info( "FPKM keycount: {0}".format( len(fpkm_table.keys()) ))
            fpkm_tables.append( fpkm_table )

            tpm_table = None                                  # Cufflinks doesn't generate TPM
            if "tpm_expression_levels" in expr["data"]:      # so we need to check for this key
                tpm_table = expr["data"]["tpm_expression_levels"]
                logger.info( "TPM keycount: {0}".format( len(tpm_table.keys()) ))
                tpm_tables.append( tpm_table )

        em_base_name = expressionSet_name 

        store_table2exp_mat(logger, ws_client, fpkm_tables, genome_ref, set_names, "{0}_FPKM_ExpressionMatrix".format(expressionSet_name), ws_id )
        if ( tpm_table != None ):
            store_table2exp_mat(logger, ws_client, tpm_tables, genome_ref, set_names, "{0}_TPM_ExpressionMatrix".format(expressionSet_name), ws_id, hide_TPM )


def create_RNASeq_ExpressionSet_and_build_report( logger,
                                                  ws_client,
                                                  tool_used,
                                                  tool_version,
                                                  tool_opts,
                                                  ws_id,
                                                  alignment_list,
                                                  alignmentset_id,
                                                  genome_id,
                                                  sampleset_id,
                                                  results,
                                                  expressionSet_name,
                                                  hide_TPM = 0 ):

        results =  [ ret for ret in results if not ret is None ]
        logger.info( "create_RNASeq_ExpressionSet, results:")
        logger.info( pformat(results ) )
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
               logger.info( "Saving ExpressionSet object to  workspace")
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

        try:
               logger.info("Creating ExpressionMatrix object")
               create_and_load_expression_matrix( logger, ws_client, ws_id, expressionSet_name, genome_id, hide_TPM )
        except Exception as e:
                   logger.exception(e)
                   raise e

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
        return stats_data


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

def get_details_for_diff_exp(logger,ws_client,hs,ws_id,urls,directory,expressionset_id,token):
        try:
           expression_set = ws_client.get_objects(
                                        [{'name' : expressionset_id,'workspace' : ws_id}])[0]
        except Exception,e:
           raise Exception("".join(traceback.format_exc()))
        ### Getting all the set ids and genome_id
        output_obj = {}
        expression_set_info = ws_client.get_object_info_new({"objects": [{'name' : expressionset_id, 'workspace': ws_id}]})[0] 
        output_obj['expressionset_id'] =  str(expression_set_info[6]) + '/' + str(expression_set_info[0]) + '/' + str(expression_set_info[4])
        output_obj['sampleset_id'] =  expression_set['data']['sampleset_id']
        output_obj['alignmentset_id'] = expression_set['data']['alignmentSet_id'] 
        output_obj['genome_id'] = expression_set['data']['genome_id']
        output_obj['genome_name'] = ws_client.get_object_info([{"ref" :output_obj['genome_id']}],includeMetadata=None)[0][1]
        ws_gtf = output_obj['genome_name']+"_GTF_Annotation"
    
        ### Check if GTF object exists in the workspace pull the gtf
        #ws_gtf = output_obj['genome_name']+"_GTF"
        gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,urls,ws_id,ws_gtf,"KBaseRNASeq.GFFAnnotation",directory,token)
        print 'GTF file is  ' +  gtf_file
        if gtf_file is None:
             create_gtf_annotation_from_genome(logger,ws_client,hs,urls,ws_id,output_obj['genome_id'],output_obj['genome_name'],directory,token)
        output_obj['gtf_file'] = gtf_file
        ### Getting the expression objects and alignment objects
        m_expr_ids = expression_set['data']['mapped_expression_ids']
        if len(m_expr_ids)  < 2:
           raise ValueError("Error the ExpressionSet object has less than 2 expression samples. Kindly check your reads files and repeat the previous step (Cufflinks)")
        output_obj['labels'] = []
        output_obj['alignments'] = []   
        counter = 0
        assembly_file = os.path.join(directory,"assembly_gtf.txt")
        list_file = open(assembly_file,'w')
        for i in m_expr_ids:
            for a_id ,e_id in i.items():
                   files = {}
                   a_obj,e_obj = ws_client.get_objects(
                                        [{'ref' : a_id},{'ref': e_id}])
                   ### Get the condition name, replicate_id , shock_id and shock_filename
                   condition = a_obj['data']['condition']
                   if 'replicate_id' in a_obj['data'] : replicate_id = a_obj['data']['replicate_id']
                   files[a_obj['data']['file']['file_name']] = a_obj['data']['file']['id']
                   files[e_obj['data']['file']['file_name']] = e_obj['data']['file']['id']
                   if not condition in output_obj['labels']: output_obj['labels'].append(condition)
                   else :  counter += 1 #### comment it when replicate_id is available from methods
                   s_path = os.path.join(directory,condition+"/"+str(counter)) ### Comment this line when replicate_id is available from the methods
                   if not os.path.exists(s_path): os.makedirs(s_path)
                   try:
                       script_util.download_shock_files(logger,urls['shock_service_url'],s_path,files,token)
                   except Exception,e:
                       raise Exception( "Unable to download shock file, {0}".format(e))
                   try:
                       script_util.unzip_files(logger,os.path.join(s_path,a_obj['data']['file']['file_name']),s_path)
                       script_util.unzip_files(logger,os.path.join(s_path,e_obj['data']['file']['file_name']),s_path)
                       e_file_path =  os.path.join(s_path,"transcripts.gtf")
                       a_file_path = os.path.join(s_path,"accepted_hits.bam")
                       if os.path.exists(a_file_path) : 
                                print a_file_path
                                output_obj['alignments'].append(a_file_path)
                       if os.path.exists(e_file_path) : list_file.write("{0}\n".format(e_file_path))
                   except Exception, e:
                       raise Exception("".join(traceback.format_exc()))
        list_file.close()
        output_obj['gtf_list_file'] = assembly_file
        print output_obj
        return output_obj        

# This downloads and unzips in a separate directory the stringtie output from each expression object
# in the given expression set, and verifies that it has the correct files for ballgown to run
# (*.ctab).  Each directory name is prefixed by the given dir_prefix so that the ballgown call can
# use it to specify the input directories to ballgown
#
# This also returns information extracted from the ExpressionSet object which is useful for further 
# workspace storage (genome_id, etc )
#
def get_info_and_download_for_ballgown( logger, ws_client, hs, ws_id, urls, directory, dir_prefix, expressionset_id, token ):
        try:
           expression_set = ws_client.get_objects( [ { 'name'      : expressionset_id,
                                                       'workspace' : ws_id }
                                                   ] )[0]

        except Exception,e:
           raise Exception( "".join(traceback.format_exc() ))
        ### Getting all the set ids and genome_id
        output_obj = {}
        expression_set_info = ws_client.get_object_info_new({"objects": [{'name' : expressionset_id, 'workspace': ws_id}]})[0] 
        output_obj['expressionset_id'] =  str(expression_set_info[6]) + '/' + str(expression_set_info[0]) + '/' + str(expression_set_info[4])
        #output_obj['sampleset_id'] =  expression_set['data']['sampleset_id']
        #output_obj['alignmentset_id'] = expression_set['data']['alignmentSet_id'] 
        #output_obj['genome_id'] = expression_set['data']['genome_id']
        #output_obj['genome_name'] = ws_client.get_object_info([{"ref" :output_obj['genome_id']}],includeMetadata=None)[0][1]
        #ws_gtf = output_obj['genome_name']+"_GTF_Annotation"

        # QUESTION:  whats the reason for getting this from object_info and not the object
        output_obj['genome_id']             = expression_set['data']['genome_id']
        output_obj['alignmentSet_id']       = expression_set['data']['alignmentSet_id']
        output_obj['sampleset_id']          = expression_set['data']['sampleset_id']
        output_obj['sample_expression_ids'] = expression_set['data']['sample_expression_ids']
        ### Check if GTF object exists in the workspace pull the gtf
        #ws_gtf = output_obj['genome_name']+"_GTF"
        #gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,urls,ws_id,ws_gtf,"KBaseRNASeq.GFFAnnotation",directory,token)
        #print 'GTF file is  ' +  gtf_file
        #if gtf_file is None:
        #     create_gtf_annotation_from_genome(logger,ws_client,hs,urls,ws_id,output_obj['genome_id'],output_obj['genome_name'],directory,token)
        #output_obj['gtf_file'] = gtf_file
        ### Getting the expression objects and alignment objects
        m_expr_ids = expression_set['data']['mapped_expression_ids']
        if len(m_expr_ids)  < 2:
           raise ValueError("Error the ExpressionSet object has less than 2 expression samples. Kindly check your reads files and repeat the previous step (Cufflinks)")
        output_obj['labels'] = []
        #output_obj['alignments'] = []   
        output_obj['subdirs'] = []
        counter = 0
        #assembly_file = os.path.join( directory, "assembly_gtf.txt" )
        #list_file = open( assembly_file, 'w' )
        for i in m_expr_ids:                  
            for a_id, e_id in i.items():
                #files = {}
                #a_obj,e_obj = ws_client.get_objects(
                #                     [{'ref' : a_id},{'ref': e_id}])
                e_obj= ws_client.get_objects( [ {'ref' : e_id} ] )[0]
                zipfile = e_obj['data']['file']['file_name']
                shock_id = e_obj['data']['file']['id']
                ### Get the condition name, replicate_id , shock_id and shock_filename
                condition = e_obj['data']['condition']                                            # Question: can we use this for group!?
                if 'replicate_id' in e_obj['data'] : replicate_id = e_obj['data']['replicate_id']
                #files[a_obj['data']['file']['file_name']] = a_obj['data']['file']['id']
                #files[e_obj['data']['file']['file_name']] = e_obj['data']['file']['id']
                if not condition in output_obj['labels']: 
                    output_obj['labels'].append( condition )
                else:
                    counter += 1 #### comment it when replicate_id is available from methods

                #subdir = os.path.join( directory, dir_prefix + "_" + condition + "_" + str(counter) ) ### Comment this line when replicate_id is available from the methods
                # Fix this - we're getting expression object name from the zip file
                if ( zipfile[-4:] != '.zip' ):
                    raise  Exception( "zip file {0} doesn't seem to have .zip extention, can't form a subdirectory name confidently" )

                subdir = os.path.join( directory, zipfile[0:-4] )
                logger.info( "subdir is {0}".format( subdir ) )
                output_obj['subdirs'].append( subdir )
                if not os.path.exists( subdir ): 
                    os.makedirs( subdir )
                try:
                    script_util.download_file_from_shock( logger = logger, 
                                                          shock_service_url = urls['shock_service_url'], 
                                                          shock_id = shock_id,
                                                          filename = zipfile,
                                                          directory = subdir, 
                                                          token = token )
                    #script_util.download_shock_files( logger, urls['shock_service_url'], s_path, files, token )
                except Exception,e:
                    raise Exception( "Unable to download shock file, {0}".format(e))
                try:
                    script_util.unzip_files( logger, 
                                             os.path.join( subdir, zipfile ),
                                             subdir )
                    logger.info( "listing of {0}".format( subdir ))
                    logger.info( os.listdir(subdir) )
                    #script_util.unzip_files(logger,os.path.join(s_path,e_obj['data']['file']['file_name']),s_path)
                    #e_file_path =  os.path.join(s_path,"transcripts.gtf")
                    #a_file_path = os.path.join( s_path, "accepted_hits.bam" )
                    #if os.path.exists(a_file_path) : 
                    #         print a_file_path
                    #         output_obj['alignments'].append(a_file_path)
                    #if os.path.exists(e_file_path) : list_file.write("{0}\n".format(e_file_path))
                except Exception, e:
                    raise Exception("".join(traceback.format_exc()))

                for f in [ 'e2t.ctab', 'e_data.ctab', 'i2t.ctab', 'i_data.ctab', 't_data.ctab' ]:
                    fullpath = os.path.join( subdir, f )
                    if ( not os.path.isfile( fullpath ) ):
                        raise Exception( "error: ballgown input file {0} not found. Can't proceed".format( fullpath ) )

        #list_file.close()
        #output_obj['gtf_list_file'] = assembly_file
        print output_obj
        return output_obj        

# this converts the input group set lists writes a file that can be loaded
# as a table by the ballgown_fpkmgenematrix.R  script to pass to
# ballgown to assign group ids to each setn

# this returns an ordered list of group names that corresponds to
# the input subdir_list

def create_sample_dir_group_file( subdir_list, 
                                  group1_name,
                                  group1_set,
                                  group2_name,
                                  group2_set,
                                  sample_dir_group_file ):
        try:
            f = open( sample_dir_group_file, "w")
        except Exception:
            raise Exception( "Can't open file {0} for writing {1}".format( sample_dir_group_file, traceback.format_exc() ) )
        group_name_list = []
        for subdir in subdir_list:
            exp = os.path.basename( subdir )
            if ( exp in group1_set ):
                if ( exp in group2_set ):
                    raise Exception( "group error - {0} is found in both group sets".format( exp ) )
                group = 0
                group_name_list.append( group1_name )
            elif ( exp in group2_set ):
                if ( exp in group1_set ):
                    raise Exception( "group error - {0} is found in both group sets".format( exp ) )
                group = 1
                group_name_list.append( group2_name )
            else:
                raise Exception( "group error - {0} is not found in either group set".format( exp ) )
            f.write( "{0}  {1}\n".format( subdir, group ))
        f.close()

        return( group_name_list )

def  run_ballgown_diff_exp( logger, rscripts_dir, directory, sample_dir_group_table_file, ballgown_output_dir, output_csv ):

        # sample_group_table is a listing of output Stringtie subdirectories,
        # (full path specification) paired with group label (0 or 1), ie
        #    /path/WT_rep1_stringtie    0
        #    /path/WT_rep2_stringtie    0 
        #    /path/EXP_rep1_stringtie   1
        #    /path/EXP_rep2_stringtie   1
        #  (order doesn't matter, but the directory-group correspondance does)

        #  1) Need to make a generic "Run rscript program"
        #  2) is this the best way to run the R script (subprocess Popen?)


        # Make call to execute the system.

        rcmd_list = [ 'Rscript', os.path.join( rscripts_dir, 'ballgown_fpkmgenematrix.R' ),
                      '--sample_dir_group_table', sample_dir_group_table_file,
                      '--output_dir',  ballgown_output_dir,
                      '--output_csvfile', output_csv
                    ] 
        rcmd_str = " ".join( str(x) for x in rcmd_list )
        logger.info( "rcmd_string is {0}".format( rcmd_str ) )
        openedprocess = subprocess.Popen( rcmd_str, shell=True )  # , stdout=subprocess.PIPE )
        openedprocess.wait()
        #Make sure the openedprocess.returncode is zero (0)
        if openedprocess.returncode != 0:
            logger.info( "R script did not return normally, return code - "
                + str(openedprocess.returncode))
            raise Exception( "Rscript failure" )

# reads csv diff expr matrix file from Ballgown and returns as a 
# dictionary of rows with the gene as key.  Each key gives a row of
# length three corresponding to fold_change, pval and qval in string form
# - can include 'NA'
#
def  load_diff_expr_matrix( ballgown_output_dir, output_csv ):

        diff_matrix_file = os.path.join( ballgown_output_dir, output_csv )

        if not os.path.isfile( diff_matrix_file ):
            raise Exception( "differential expression matrix csvfile {0} doesn't exist!".format(diff_matrix_file) )

        n = 0
        dm = {}
        with  open( diff_matrix_file, "r" ) as csv_file:
            csv_rows = csv.reader( csv_file, delimiter="\t", quotechar='"' )
            for row in csv_rows:
                n = n + 1
                if ( n == 1  ):
                    if ( row != ['id', 'fc', 'pval', 'qval'] ):
                        raise Exception( "did not get expected column heading from {0}".format( diff_matrix_file ) )
                else:
                    if ( len(row) != 4 ):
                        raise Exception( "did not get 4 elements in row {0} of csv file {1} ".format( n, diff_matrix_file ) )
                    key = row[0]
                    # put in checks for NA or numeric for row[1] through 4
                    if ( key in dm ):
                        raise Exception( "duplicate key {0} in row {1} of csv file {2} ".format( key, n, diff_matrix_file  ) )
                    dm[ key ] = row[1:5]

        return dm

# this takes the output_csv generated by run_ballgown_diff_exp() and loads it into
# an RNASeqDifferentialExpression named output_object_name

def  load_ballgown_output_into_ws( logger, 
                                   ws_id,
                                   ws_client, 
                                   hs_client, 
                                   token, 
                                   directory, 
                                   ballgown_output_dir, 
                                   tool_used,
                                   tool_version,
                                   sample_ids,
                                   conditions,
                                   genome_id,
                                   expressionset_id,
                                   alignmentset_id,
                                   sampleset_id,
                                   output_object_name
                                 ):

        logger.info( "Zipping ballgown output" )
        zip_file_path = os.path.join( directory, "{0}.zip".format( output_object_name ) )
        try:
            script_util.zip_files( logger, ballgown_output_dir, zip_file_path )
        except Exception,e:
            raise Exception( "Error creating zip file of ballgown output" )

        logger.info( "Creating handle from ballgown output zip file" )
        try:
            handle = hs_client.upload( zip_file_path )
        except Exception, e:
            raise Exception( "Error uploading ballgown output zip file to handle service: {0}".format( " ".join( traceback.print_exc() ) ) )

        # create object

        de_obj = { "tool_used"        : tool_used,
                   "tool_version"     : tool_version,
                   "sample_ids"       : sample_ids,
                   "condition"        : conditions,
                   "genome_id"        : genome_id,
                   "expressionSet_id" : expressionset_id,
                   "alignmentSet_id"  : alignmentset_id,
                   "sampleset_id"     : sampleset_id,
                   "file"             : handle
                  }

        objs_save_data = ws_client.save_objects(
                                      { "workspace":  ws_id,
                                        "objects":    [
                                                       {
                                                        "type":    "KBaseRNASeq.RNASeqDifferentialExpression",
                                                        "data":    de_obj,
                                                        "name":    output_object_name
                                                       } 
                                                      ]
                                      }
                                    )
        return objs_save_data[0]


#  returns a list of gene names from the keys of diff_expr_matrix 
#  assumed AND logic between all the tests

def  filter_genes_diff_expr_matrix( diff_expr_matrix, 
                                    scale_type,
                                    qval_cutoff,
                                    fold_change_cutoff,     # says this is log2 but we should make it match scale_type  
                                    max_genes ):

    # verify that qval_cutoff, fold change cutoff is None or numeric.
    # verify scale_type is an allowed value

    # !!!!! figure out what to do with scale conversion

    if  max_genes == None:
        max_genes = sys.maxint

    selected = []
    ngenes = 0
    for gene in diff_expr_matrix:

        fc, pval, qval = diff_expr_matrix[gene]

        # should NA, NaN fc, qval automatically cause exclusion?

        if ( qval_cutoff != None  ):                # user wants to apply qval_cutoff
            if qval == 'NA' or qval == "Nan":       # bad values automatically excluded
                continue 
            q = make_numeric( qval, "qval, gene {0}".format( gene ) )
            if ( q < qval_cutoff ):
                continue

        if ( fold_change_cutoff != None  ):         # user wants to apply fold_change_cutoff
            if fc == 'NA' or fc  == "Nan":          # bad values automatically excluded
                continue
            f = make_numeric( fc, "fc, gene {0}".format( gene ) )
            if ( f < fold_change_cutoff ):
                continue

        ngenes = ngenes + 1
        if  ngenes > max_genes:
            break

        # if we got here, it made the cut, so add the gene to the list
        selected.append( gene )

    return  selected


def  make_numeric( x, msg ):
    try:
        res = float( x )
    except ValueError:
        raise Exception( "illegal non-numeric, {0}".format( msg ) )
    return res

# this weeds out all the data (rows, mappings) in given expression matrix object
# for genes not in the given selected_list, and returns a new expression matrix object

def filter_expr_matrix_object( emo, selected_list ):

    fmo = {}
    fmo["type"]            = emo["type"]
    fmo["scale"]           = emo["scale"]
    fmo["genome_ref"]      = emo["genome_ref"]

    fmo["data"] = {}
    fmo["data"]["col_ids"] = emo["data"]["col_ids"]
    fmo["data"]["row_ids"] = []
    fmo["data"]["values"]  = []

    nrows = len( emo["data"]["row_ids"] )
    if nrows != len( emo["data"]["values"]):
        raise Exception( "filtering expression matrix:  row count mismatch in expression matrix" )

    for i in xrange( nrows ):
        if emo["data"]["row_ids"][i] in selected_list:
            fmo["data"]["row_ids"].append( emo["data"]["row_ids"][i] )
            fmo["data"]["values"].append( emo["data"]["values"][i] )

    fmo["feature_mapping"] = {}
    for g, v in emo["feature_mapping"].iteritems():
        if g in selected_list:
            fmo["feature_mapping"][g] = v

    return fmo


def call_cuffmerge(working_dir,directory,num_threads,gtf_file,list_file):
         #cuffmerge_dir = os.path.join(directory,"cuffmerge")
         print "Entering cuffmerge"
         print "Args passed {0},{1},{2},{3}".format(directory,num_threads,gtf_file,list_file)
         cuffmerge_command = " -p {0} -o {1} -g {2} {3}".format(str(num_threads),directory,gtf_file,list_file)
         merged_gtf = None
         try:
                #logger.info("Executing: cuffmerge {0}".format(cuffmerge_command))
                print  "Executing: cuffmerge {0}".format(cuffmerge_command)
                r,e = script_util.runProgram(None,"cuffmerge",cuffmerge_command,None,working_dir)
                print r + "\n" + e
                if os.path.exists(directory+"/merged.gtf") : merged_gtf = os.path.join(directory,"merged.gtf")
         except Exception,e:
                print "".join(traceback.format_exc())
                raise Exception("Error executing cuffmerge {0},{1}".format(cuffmerge_command,"".join(traceback.format_exc())))
         return merged_gtf

def call_stringtiemerge(working_dir,directory,num_threads,gtf_file,list_file):
         #directory = os.path.join(directory,"cuffmerge")
         strmerge_command = " -p {0} -o {1} --merge -G {2} {3}".format(str(num_threads),directory,gtf_file,list_file)
         merged_gtf = None
         try:
                print  "Executing: stringtie {0}".format(strmerge_command)
                script_util.runProgram(None,"stringtie",strmerge_command,None,working_dir)
                if os.path.exists(directory+"/merged.gtf") : merged_gtf = os.path.join(directory,"merged.gtf")
         except Exception,e:
                raise Exception("Error executing StringTie merge {0},{1}".format(strmerge_command,working_dir))
         return merged_gtf

def call_tablemaker(working_dir,directory,num_threads,m_gtf_file,alignment_file):
         print "Inside Tablemaker"
         print "Args passed : {0} , {1} , {2} , {3} , {4} ".format(working_dir,directory,num_threads,m_gtf_file,alignment_file)
         #cuffmerge_dir = os.path.join(directory,"cuffmerge")
         tm_command = " -p {0} -o {1} -q -W -G {2} {3}".format(str(num_threads),directory,m_gtf_file,alignment_file)
         try:
                print "Executing: tablemaker {0}".format(tm_command)
                script_util.runProgram(None,"tablemaker",tm_command,None,working_dir)
         except Exception,e:
                logger.exception(e)
                raise Exception("Error executing tablemaker {0},{1}".format(tm_command,working_dir))
         return directory

def call_stringtieBall(working_dir,directory,num_threads,m_gtf_file,alignment_file):
         #directory = os.path.join(directory,"cuffmerge")
         print "Inside stringtieBall"
         strdiff_command = " -p {0} -o {1} -e -B -G {2} {3}".format(str(num_threads),directory,m_gtf_file,alignment_file)
         try:
                print "Executing: stringtie {0}".format(strdiff_command)
                script_util.runProgram(None,"stringtie",strdiff_command,None,working_dir)
         except Exception,e:
                raise Exception("Error executing StringTie differential expression {0},{1}".format(strdiff_command,working_dir))
         return directory        

### TODO Function related to th Genome Annotation Changes and hence needs to be removed        
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
                raise Exception("Unable to Create FASTA file from Genome Annotation : {0}".format("".join(traceback.format_exc())))
        finally:
                fasta_file.close()
        fasta_end = datetime.datetime.utcnow()
        logger.info("Generating FASTA for {} took {}".format(obj_name, fasta_end - fasta_start))
        return output_file

### TODO Function related to the Genome Annotation Changes and hence needs to be removed        
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
        try:
                ga.get_gff().to_file(gff_file)
        except Exception as e:
                raise Exception("Unable to Create GFF  file from Genome Annotation : {0}: {1}".format(obj_name,"".join(traceback.format_exc())))
        finally:
                    gff_file.close()
        gff_end = datetime.datetime.utcnow()
        logger.info("Generating GFF for {} took {}".format(obj_name, gff_end - gff_start))


### TODO Function related to the Genome Annotation Changes and hence needs to be removed        
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
                   cmdline_output = script_util.runProgram(None,"gffread",gtf_cmd,None,directory)
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
        


### TODO Deprecated function  and hence needs to be removed        
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

### TODO Deprecated Function and hence needs to be removed        
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
                

# This was moved in from script_util.py (5 Jan 2017 mccorkle)
#
def parse_FPKMtracking( filename, tool, metric) :
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
                result[larr[pos1]] = math.log( float(larr[pos2]) + 1, 2 )
    return result

# 
# added to provide TPM expression for cufflinks output
# Cufflinks doesn't produce TPM, we infer from FPKM 
# (see discussion @ https://www.biostars.org/p/160989/)
#
def parse_FPKMtracking_calc_TPM( filename ) :

    fpkm_dict = {}
    tpm_dict = {}
    gene_col = 0
    fpkm_col = 9
    sum_fpkm = 0.0
    with open( filename ) as f:
        next( f )
        for line in f:
            larr = line.split("\t")
            gene_id = larr[gene_col]
            if gene_id != "":
                fpkm = float( larr[fpkm_col] )
                sum_fpkm = sum_fpkm + fpkm
                fpkm_dict[gene_id] = math.log( fpkm + 1, 2 )
                tpm_dict[gene_id] = fpkm

    for g in tpm_dict:
        tpm_dict[g] = math.log( ( tpm_dict[g] / sum_fpkm ) * 1e6 + 1, 2 )

    return  fpkm_dict, tpm_dict

