#import multiprocessing as mp
import pathos.multiprocessing as mp

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

class KBParallelExecutionBase(ExecutionBase):

    # Overriding run() method:

    # method is "TophatCall", etc
 
    def run(self, method, common_params, run_params):

        self._checkCommonParams(common_params)
        #self._setCommonParams(common_params)
        #self.method_params = method_params
        #self.common_params = common_params

        kbp = KBParallel( os.environ['SDK_CALLBACK_URL'], token=common_params['user_token'])
        returnVal = kbp.run( { 'method': { 'module_name': "KBaseRNASeq",
                                           'method_name': method,
                                           'service_ver': ""
                                          }
                               'is_local': 1,
                               'global_params': run_params,  # NOTE: this is called global_input in KBParallel.spec:  FIX!
                               'time_limit': 1000000} )
        return returnVal

