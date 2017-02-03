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

def parallelize1(function,num_processors):
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
def CallHisat2_helper(x):
    	logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,genome_id,sampleset_id,params,token = x
    	return _CallHisat2(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,genome_id,sampleset_id,params,token)

def _CallHisat2(logger,services,ws_client,hs,ws_id,sample_type,num_threads,read_sample,condition,directory,genome_id,sampleset_id,params,token):
        #logger.info("Downloading Read Sample{0}".format(read_sample))
        print "Downloading Read Sample{0}".format(read_sample)
        if not logger:
                logger = handler_util.create_logger(directory,"run_Hisat2_"+read_sample)
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
                hisat2_base =os.path.join(directory,handler_util.get_file_with_suffix(directory,".1.ht2"))
                ### Adding advanced options to Bowtie2Call
                hisat2_cmd = ''
                hisat2_cmd += ( ' -p {0}'.format(num_threads))
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
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read_id,filename=read_name, directory=input_direc,token=token)
                                hisat2_cmd += " -U {0} -x {1} -S {2}".format(os.path.join(input_direc,read_name),hisat2_base,out_file)
                        except Exception,e:
                                #logger.exception( "Unable to download shock file , {0}".format(read_name))
                                raise Exception( "Unable to download shock file , {0}".format(read_name))
                if sample_type == 'KBaseAssembly.PairedEndLibrary':
                        lib_type = 'PairedEnd'
                	if('orientation' in params and params['orientation'] is not None): hisat2_cmd += ( ' --'+params['orientation'])
                        read1_id = r_sample['data']['handle_1']['id']
                        read1_name = r_sample['data']['handle_1']['file_name']
                        read2_id = r_sample['data']['handle_2']['id']
                        read2_name = r_sample['data']['handle_2']['file_name']
                        try:
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read1_id,filename=read1_name, directory=input_direc,token=token)
                                script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=read2_id,filename=read2_name, directory=input_direc,token=token)
                                hisat2_cmd += " -1 {0} -2 {1} -x {2} -S {3}".format(os.path.join(input_direc,read1_name),os.path.join(output_dir,read2_name),hisat2_base,out_file)
                        except Exception,e:
                                #logger.Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
                                raise Exception( "Unable to download shock file , {0} or {1}".format(read1_name,read2_name))
                try:
                        logger.info("Executing: hisat2 {0}".format(hisat2_cmd))
                        cmdline_output = script_util.runProgram(logger,"hisat2",hisat2_cmd,None,directory)
  		except Exception,e:
                        raise Exception("Failed to run command {0}".format(hisat2_cmd))
                        #logger.exception("Failed to run command {0}".format(hisat2_cmd))
                try:
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
                        raise Exception("Error Running the hisat2 command {0},{1} {2}".format(hisat2_cmd,directory," ".join(traceback.print_exc())))
                        #logger.exception("Error Running the hisat2 command {0},{1} {2}".format(hisat2_cmd,directory," ".join(traceback.print_exc())))

                # Zip tophat folder
                try:
                        out_file_path = os.path.join(directory,"%s.zip" % output_name)
                        logger.info("Zipping the output files".format(out_file_path))
                        script_util.zip_files(logger, output_dir,out_file_path)
                except Exception, e:
                        raise Exception("Failed to compress the index: {0}".format(out_file_path))
                        #logger.exception("Failed to compress the index: {0}".format(out_file_path))
                ## Upload the file using handle service
                try:
                        hisat2_handle = hs.upload(out_file_path)
                except Exception, e:
                        raise Exception("Failed to upload zipped output file".format(out_file_path))
                        #logger.exception("Failed to upload zipped output file".format(out_file_path))
                #### Replace version with get_version command#####
                hisat2_out = { "file" : hisat2_handle ,"size" : os.path.getsize(out_file_path), "aligned_using" : "hisat2" , "aligner_version" : "2.2.6" , 'library_type' : lib_type , 'condition' : condition ,'read_sample_id': read_sample, 'genome_id' : genome_id , "alignment_stats" : stats_data }
                if not sampleset_id is None: hisat2_out['sampleset_id'] = sampleset_id
                try:
                        res= ws_client.save_objects(
                                        {"workspace":ws_id,
                                         "objects": [{
                                         "type":"KBaseRNASeq.RNASeqAlignment",
                                         "data":hisat2_out,
                                         "name":output_name}
                                        ]})
                except Exception, e:
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
                #else:
        return None
    	
def runMethod(logger,token,ws_client,hs,services,hisat2_dir,params):
	try:
               sample,annotation_name = ws_client.get_objects(
                                        [{ 'name' : params['sampleset_id'], 'workspace' : params['ws_id']},
                                        { 'name' : params['genome_id'], 'workspace' : params['ws_id']}])
        except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise ValueError(" Error Downloading objects from the workspace ")
            ### Get obejct IDs
        sampleset_info,annotation_info = ws_client.get_object_info_new({"objects": [
										   {'name': params['sampleset_id'], 'workspace': params['ws_id']},
										   {'name': params['genome_id'], 'workspace': params['ws_id']}
										   ]})
	### Get the workspace object ids for the objects ###
        sampleset_id = str(sampleset_info[6]) + '/' + str(sampleset_info[0]) + '/' + str(sampleset_info[4])
        annotation_id = str(annotation_info[6]) + '/' + str(annotation_info[0]) + '/' + str(annotation_info[4])
	sample_type = sampleset_info[2].split('-')[0]
	### Check if the Library objects exist in the same workspace
	logger.info("Check if the Library objects do exist in the current workspace")
        if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
        	reads = sample['data']['sample_ids']
        	reads_type= sample['data']['Library_type']
        	if reads_type == 'PairedEnd': r_type = 'KBaseAssembly.PairedEndLibrary'
        	else: r_type = 'KBaseAssembly.SingleEndLibrary'
        	e_ws_objs = script_util.if_ws_obj_exists(None,ws_client,params['ws_id'],r_type,reads)
        	missing_objs = [i for i in reads if not i in e_ws_objs]
        	if len(e_ws_objs) != len(reads):
            		raise Exception('Missing Library objects {0} in the {1}. please copy them and run this method'.format(",".join(missing_objs),params['ws_id']))

	### Build Hisat2 index
	fasta_file = script_util.generate_fasta(logger,services,token,annotation_id,hisat2_dir,params['genome_id'])
        logger.info("Sanitizing the fasta file to correct id names {}".format(datetime.datetime.utcnow()))
        mapping_filename = c_mapping.create_sanitized_contig_ids(fasta_file)
        c_mapping.replace_fasta_contig_ids(fasta_file, mapping_filename, to_modified=True)
        logger.info("Generating FASTA file completed successfully : {}".format(datetime.datetime.utcnow()))
        hisat2base =os.path.join(hisat2_dir,handler_util.get_file_with_suffix(hisat2_dir,".fasta"))
        hisat2base_cmd = '{0} {1}'.format(fasta_file,hisat2base)
	try:
            logger.info("Building Index for Hisat2 {0}".format(hisat2base_cmd))
            cmdline_output = script_util.runProgram(logger,"hisat2-build",hisat2base_cmd,None,hisat2_dir)
        except Exception,e:
            raise Exception("Failed to run command {0}".format(hisat2base_cmd))
        ws_gtf = params['genome_id']+"_GTF"
        ret = script_util.if_obj_exists(None,ws_client,params['ws_id'],"KBaseRNASeq.GFFAnnotation",[ws_gtf])
        print ret
        if not ret is None:
            logger.info("GFF Annotation Exist for Genome Annotation {0}.... Skipping step ".format(params['genome_id']))
            annot_name,annot_id = ret[0]
            gtf_obj=ws_client.get_objects([{'ref' : annot_id}])[0]
            gtf_id=gtf_obj['data']['handle']['id']
            gtf_name=gtf_obj['data']['handle']['file_name']
 	    try:
                     script_util.download_file_from_shock(logger, shock_service_url=services['shock_service_url'], shock_id=gtf_id,filename=gtf_name, directory=hisat2_dir,token=token)
                     gtf_file = os.path.join(hisat2_dir,gtf_name)
            except Exception,e:
                        raise Exception( "Unable to download shock file, {0}".format(gtf_name))
        else:
             script_util.create_gtf_annotation(logger,ws_client,hs,services,params['ws_id'],annotation_id,params['genome_id'],fasta_file,hisat2_dir,token)
	# Determine the num_threads provided by the user otherwise default the number of threads to 2
        if('num_threads' in params and params['num_threads'] is not None):
            num_threads = int(params['num_threads'])
        else:
            num_threads = 2
        num_cores = mp.cpu_count()
        logger.info("Number of available cores : {0}".format(num_cores))
        b_tasks =[]
        if sample_type == 'KBaseRNASeq.RNASeqSampleSet':
              	reads = sample['data']['sample_ids']
                reads_type= sample['data']['Library_type']
                r_label = sample['data']['condition']
                num_samples =  len(reads)
                if num_cores != 1:
                        pool_size,num_threads=handler_util.optimize_parallel_run(num_samples,num_threads,num_cores)
			pool_size = 2
	 	else:
                   	pool_size = 1
                   	num_threads = 1
                count = 0
                logger.info(" Number of threads used by each process {0}".format(num_threads))
                for i in reads:
                        try:
                                label = r_label[count]
                                b_tasks.append((None,services,ws_client,hs,params['ws_id'],reads_type,num_threads,i,label,hisat2_dir,annotation_id,sampleset_id,params,token))
                                count = count + 1
                        except Exception,e:
                                raise

                @parallelize1(CallHisat2_helper,pool_size)
                def run_hisat2_in_parallel(tasks):
                        pass
                results=run_hisat2_in_parallel(b_tasks)
		print results
		alignmentSet_name = params['sampleset_id']+"_hisat2_AlignmentSet"
                logger.info(" Creating AlignmentSet for the Alignments {0}".format(alignmentSet_name))
                reportObj=script_util.create_RNASeq_AlignmentSet_and_build_report(logger,ws_client,params['ws_id'],reads,sampleset_id,annotation_id,None,results,alignmentSet_name)
        else:
                try:
                    pool_size=1
                    num_threads = num_cores
                    logger.info(" Number of threads used by each process {0}".format(num_threads))
                    results = _CallHisat2(logger,services,ws_client,hs,params['ws_id'],sample_type,num_threads,params['sampleset_id'],'Single-Sample',hisat2_dir,annotation_id,None,params,token)
                except Exception,e:
                     raise
                single_read, single_alignment = results
                single_align_obj = ws_client.get_objects(
                                        [{ 'name' : single_alignment, 'workspace' : params['ws_id']}])[0]['data']
                sref = ws_client.get_object_info_new({"objects": [{'name':single_alignment, 'workspace': params['ws_id']}]})[0]
                reportObj = {'objects_created':[{'ref' :str(sref[6]) + '/' + str(sref[0]) + '/' + str(sref[4]),
                                                 'description' : "RNA-seq Alignment for reads Sample: {0}".format(single_read)}],
                                                 'text_message': "RNA-seq Alignment for reads Sample: {0}".format(single_read)}
            #### Save to report object #######
                #returnVal = single_align_obj   
        reportName = 'Align_Reads_using_Hisat2_'+str(hex(uuid.getnode()))
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
	return returnVal

