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
from biokbase.RNASeq import rnaseq_util
from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI , GenomeAnnotationClientAPI
from doekbase.data_api.sequence.assembly.api import AssemblyAPI , AssemblyClientAPI
import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()
try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from biokbase.RNASeq.Cufflinks import Cufflinks

class CufflinksSampleException(Exception):
    pass

class CufflinksSampleSet(Cufflinks): 

    def __init__(self, logger, directory, urls):
        pprint(self.__class__)
        super(self.__class__, self).__init__(logger, directory, urls)
        #super(StringtTieSampleSet, self).__init__(logger, directory, urls)
        # user defined shared variables across methods
        #self.sample_info = None
        self.alignmentset_info = None
        self.num_threads = 1

    def prepare( self, common_params, method_params ): 
        self.logger.info( "in CufflinksSampleSet.collect(), common_params are " )
        self.logger.info( pformat( common_params ) )
        self.logger.info( " and method_params are" )
        self.logger.info( pformat( method_params ) )

        # for quick testing, we recover parameters here
        ws_client = common_params['ws_client']
        hs = common_params['hs_client']
        params = method_params
        logger = self.logger
        token = common_params['user_token']
        cufflinks_dir = self.directory
        try:
               a_sampleset = ws_client.get_objects(
                                        [ {'name' : params['alignmentset_id'], 'workspace' : params['ws_id']} ] )[0]
               self.logger.info( "a_sampleset")
               self.logger.info( pformat( a_sampleset ) )
               a_sampleset_info = ws_client.get_object_info_new({"objects" : 
                                        [ {'name' : params['alignmentset_id'],'workspace' : params['ws_id']} ] } )[0]
               self.alignmentset_info =  a_sampleset_info    # QUESTION:  WHY IS THIS NEEDED?
               a_sampleset_id = str(a_sampleset_info[6]) + '/' + str(a_sampleset_info[0]) + '/' + str(a_sampleset_info[4])	
               alignmentset_id = a_sampleset_id
        except Exception,e:
               logger.exception("".join(traceback.format_exc()))
               raise Exception("Error Downloading objects from the workspace ")
        #read_sample_id']
        ### Check if the gtf file exists in the workspace. if exists download the file from that
        genome_id = a_sampleset['data']['genome_id']
        genome_name = ws_client.get_object_info([{"ref" :genome_id}],includeMetadata=None)[0][1]
        ws_gtf = genome_name+"_GTF_Annotation"
        gtf_file = script_util.check_and_download_existing_handle_obj( logger,
                                                                       ws_client,
                                                                       self.urls,
                                                                       params['ws_id'],
                                                                       ws_gtf,
                                                                       "KBaseRNASeq.GFFAnnotation",
                                                                       cufflinks_dir,
                                                                       token )
        if gtf_file is None:
             rnaseq_util.create_gtf_annotation_from_genome( logger,
                                                            ws_client,
                                                            hs,
                                                            self.urls,
                                                            params['ws_id'],
                                                            genome_id,
                                                            genome_name,
                                                            cufflinks_dir,
                                                            token )
        gtf_info = ws_client.get_object_info_new({"objects": [{'name': ws_gtf , 'workspace': params['ws_id']}]})[0]
        gtf_id = str(gtf_info[6]) + '/' + str(gtf_info[0]) + '/' + str(gtf_info[4])

        tool_opts = { k:str(v) for k,v in params.iteritems() if not k in ('ws_id','alignmentset_id', 'num_threads') and v is not None  }
        self.logger.info( "tool_opts:")
        self.logger.info( pformat( tool_opts ))

        alignment_ids = a_sampleset['data']['sample_alignments']
        self.logger.info( "alignment_ids:" )
        self.logger.info( pformat( alignment_ids ) )

        m_alignment_names = a_sampleset['data']['mapped_rnaseq_alignments']
        self.logger.info( "m_alignment_names:" )
        self.logger.info( pformat( m_alignment_names ) )

        sampleset_id = a_sampleset['data']['sampleset_id']
        logger.info( "sampleset_id:" )
        logger.info( sampleset_id )

        ### Get List of Alignments Names
        align_names = []
        for d_align in m_alignment_names:
            for i , j  in d_align.items():
                     align_names.append(j)
        self.logger.info( "align_names:" )
        self.logger.info( pformat( align_names ) )

        m_alignment_ids = a_sampleset['data']['mapped_alignments_ids']
        self.logger.info( "m_alignment_ids:" )
        self.logger.info( pformat( m_alignment_ids ) )

        self.num_jobs =  len(alignment_ids)
        if self.num_jobs < 2:
            raise ValueError("Please ensure you have atleast 2 alignments to run cufflinks in Set mode")

        logger.info(" Number of threads used by each process {0}".format(self.num_threads))
        count = 0
        for i in m_alignment_ids:
            for sample_name, alignment_id in i.items():
                    task_param = { "input_arguments":
                                    [
                                      { 'job_id'          : alignment_id,
                                        'gtf_file'        : gtf_file,
                                        'ws_gtf'          : ws_gtf,
                                        'ws_id'           : params['ws_id'],
                                        'genome_id'       : genome_id,
                                        'cufflinks_dir'   : self.directory,
                                        'annotation_id'   : gtf_id,
                                        'sample_id'       : sample_name,
                                        'alignmentset_id' : alignmentset_id 
                                      }
                                    ]
                                  }
                    self.task_list.append( task_param )
                    count = count + 1

        return self.task_list


    def collect( self, common_params, collect_params ):
        self.logger.info( "in CufflinksSample.collect(), common_params are " )
        self.logger.info( pformat( common_params ) )
        self.logger.info( " and collect_params are" )
        self.logger.info( pformat( collect_params ) )

        # do with 
        global_params = collect_params['global_params']
        input_result_pairs = collect_params['input_result_pairs']
        expressionSet_name = global_params['alignmentset_id']+"_cufflinks_ExpressionSet"
        self.logger.info(" Creating ExpressionSet for the Assemblies {0}".format(expressionSet_name))
        
        tool_used = 'Cufflinks'
        tool_version = "v2.2.1"
        tool_opts = None
        sampleset_id = '11600/6/1'  # not sure where to get this from
        results = []
        align_names = []
        for i in range( 0, len( input_result_pairs ) ):
             results.append(
                             ( input_result_pairs[i]['result']['alignmentset_id'],
                               input_result_pairs[i]['result']['output_name'] )
                           )
             align_names.append( input_result_pairs[i]['input']['input_arguments'][0]['sample_id'] )

        self.logger.info( 'results:')
        self.logger.info( pformat( results ) )
        self.logger.info( 'align_names:' )
        self.logger.info( pformat( align_names ) )
        # TODO: Split alignment set and report method

        reportObj = rnaseq_util.create_RNASeq_ExpressionSet_and_build_report( self.logger,
                                                                              common_params['ws_client'],
                                                                              tool_used,
                                                                              tool_version,
                                                                              tool_opts,
                                                                              global_params['ws_id'],
                                                                              align_names,
                                                                              input_result_pairs[0]['input']['input_arguments'][0]['alignmentset_id'],
                                                                              input_result_pairs[0]['input']['input_arguments'][0]['genome_id'],
                                                                              sampleset_id,
                                                                              results,
                                                                              expressionSet_name)
        returnVal = { 'output' : expressionSet_name , 'workspace' : global_params['ws_id']}

        return( returnVal )
#	reportName = 'Align_Reads_using_Hisat2_'+str(hex(uuid.getnode()))
#        report_info = self.common_params['ws_client'].save_objects({
#                                                'id':self.alignmentset_info[6],
#                                                'objects':[
#                                                {
#                                                'type':'KBaseReport.Report',
#                                                'data':reportObj,
#                                                'name':reportName,
#                                                'meta':{},
#                                                'hidden':1, # important!  make sure the report is hidden
#                                                #'provenance':provenance
#                                                }
#                                                ]
#                                               })[0]
#        self.returnVal = { "report_name" : reportName,"report_ref" : str(report_info[6]) + '/' + str(report_info[0]) + '/' + str(report_info[4]) }
