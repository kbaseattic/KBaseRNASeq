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

class KBParallelExecutionBase(ExcecutionBase):

    # Overriding run() method:

    def run(self, common_params, method_params):

        self._checkCommonParams(common_params)
        #self._setCommonParams(common_params)
        #self.method_params = method_params
        #self.common_params = common_params
        method_params["user_token"] = common_params["user_token"]   # don't know a better way pass token to kbp.run()

        kbp = KBParallel( os.environ['SDK_CALLBACK_URL'], token=common_params['user_token'])
        returnVal = kbp.run( { 'method': { 'module_name': "KBaseRNASeq",
                                           'method_name': "TophatCall"
                                           'service_ver': ""
                                          }
                               'is_local': 1,
                               'global_params': method_params,   # common_params not appropriate for prepare,runeach
                                                                 # NOTE: this is called global_input in KBParallel.spec:  FIX!
                               'time_limit': 1000000} )
        return returnVal

