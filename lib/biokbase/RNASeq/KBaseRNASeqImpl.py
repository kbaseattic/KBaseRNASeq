#BEGIN_HEADER
import simplejson
import sys
import os
import glob
import json
import logging
import time
from pprint import pprint
from biokbase.genome_util import script_util
from biokbase.workspace.client import Workspace
from biokbase.auth import Token
#END_HEADER


class KBaseRNASeq:
    '''
    Module Name:
    KBaseRNASeq

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    __TEMP_DIR = 'temp_dir'
    __WS_URL = 'https://ci.kbase.us/services/ws'
    __SHOCK_URL = 'https://ci.kbase.us/services/shock-api/'
    __LOGGER = None
    __ERR_LOGGER = None
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        # This is where config variable for deploy.cfg are available
        pprint(config)
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']

        # logging
        self.__LOGGER = logging.getLogger('KBaseRNASeq')
        if 'log_level' in config:
              self.__LOGGER.setLevel(config['log_level'])
        else:
              self.__LOGGER.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.__LOGGER.addHandler(streamHandler)
        self.__LOGGER.info("Logger was set")

        #END_CONSTRUCTOR
        pass

    def CallFastqc(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CallFastqc
        #END CallFastqc

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CallFastqc return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def TophatCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN TophatCall
        #END TophatCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method TophatCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CufflinksCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CufflinksCall
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CufflinksCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffmergeCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffmergeCall

        # parse your params here

        #script call to your cuffmerge 

        # make sure the output file is properly stored back to WS or SHock

        #END CuffmergeCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffmergeCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffdiffCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffdiffCall
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffdiffCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def getAlignmentStats(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN getAlignmentStats
        #END getAlignmentStats

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method getAlignmentStats return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionHistogram(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionHistogram
        #END createExpressionHistogram

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionHistogram return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionSeries(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionSeries
        #END createExpressionSeries

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionSeries return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionMatrix(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionMatrix

        user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)


        # pull data from ws example
        #rnaseq_exp_details=ws_client.get_objects([{'name':params['in_id'], 
        #                                      'workspace': params['ws_id']}])

        # pull data from shock example
        #script_util.download_file_from_shock(self.__LOGGER,
        #                         shock_service_url = self.__SHOCK_URL, 
        #                         shock_id = query_rst[0]['id'],
        #                         filename = self.__INDEX_ZIP,
        #                         directory = self.__TEMP_DIR,
        #                         token = svc_token)
        #script_util.unzip_files(self.__LOGGER, zip_fn, blast_dir)

        pprint(params)


        # actual working code logic to be HERE

        # fill template empty expression matrix
        expr = {'type'  : 'level',
                'scale' : 'raw',
                'data'  : {
                            'row_ids' : [],
                            'col_ids' : [],
                            'values' : [[]]
                         }
               }
        
        # save back to workspace
        ws_client.save_objects(
            {"workspace":params['ws_id'],
            "objects": [{
                "type":"KBaseFeatureValues.ExpressionMatrix",
                "data":expr,
                "name":params['out_id']}
            ]})

        job_id ="no_job_id"

        #END createExpressionMatrix

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionMatrix return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]
