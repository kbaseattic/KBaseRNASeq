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

class ExecutionBaseException(Exception):
    pass

class ExecutionBaseParamException(ExecutionBaseException):
    pass


class ExecutionBase(object):
#need to implement runEach(), optionally prepare(), collect() and writeReport()
    # defaut common parameters
    pl = ['ws_client', 'hs_client', 'user_token'];
    def __init__(self, logger, directory, urls):
        self.logger = logger # default logger
        self.directory = directory # working directory
        self.urls = urls # service URLs
        self.results = [] # parallel execution results
        self.returnVal = {} # final execution results
        self.task_list =[] # parallel task list

        #self.clearCommonParams()
        # user defined parameters
        self.common_params = None
        self.method_params = None
        self.num_threads = None
        self.num_jobs = None
	#self.callback_url = os.environ['SDK_CALLBACK_URL']
        # default parameters
        self.num_cores = mp.cpu_count()

    def _checkCommonParams(self, common_params):
        for p in self.pl:
            if p not in common_params:
                raise ExecutionBaseParamException("Common parameter `{0}' was not defined!".format(p))

    def _clearCommonParams(self):
        for p in self.pl:
            self.__dict__[p] = None
    
    def _setCommonParams(self, common_params):
        for p in pl:
            if p in common_params:
                self.__dict__[p] = common_params[p]
        
    # TODO: remove argument passing as we defined them
    def _optimizeParallel(self, num_jobs,num_threads,num_cores):
        """
        Optimizes the pool_size and the number of threads for any parallel operation
        """
        print "Optimizing parallel jobs for number ot threads and samples run at a time"
        print "Parameters passed : no. jobs {0} ,  no. threads {1} , no. cores {2}".format(num_jobs,num_threads,num_cores)
        if num_jobs * num_threads < num_cores:
            while num_jobs * num_threads < num_cores:
                num_threads = num_threads + 1
                continue 
            pool_size = num_jobs
            return (pool_size , num_threads)
        elif num_jobs * num_threads == num_cores:
            if num_threads < 2:
                num_threads = 2 
                num_jobs = int(round(num_jobs/num_threads))
            pool_size = num_jobs
            return (pool_size, num_threads)
        else: #elif num_jobs * num_threads > num_cores:           
            while num_jobs * num_threads > num_cores:
                if num_threads > 2:
                    num_threads = num_threads - 1
                if num_jobs > 1:
                    num_jobs = num_jobs - 1 
                continue
            pool_size = num_jobs
            return (pool_size, num_threads)

    

    # Prepapre common files and necessary input files for parallel execution
    # common_params contains common parameters regardless of method such as user token, ws_client, hs_client, etc
    # method_params is the parmaters passed by Narrative method call
    def prepare(self, method_params):
        raise NotImplementedMethod() 
        # set num_jobs, 

    # generic implementation of parallel execution, method independent
    def run(self, common_params, method_params):

        self._checkCommonParams(common_params)
        #self._setCommonParams(common_params)
        self.method_params = method_params
        self.common_params = common_params

        # parallel execution common additional preparation
        self.num_threads = 2
        if('num_threads' in common_params and common_params['num_threads'] is not None):
            self.num_threads = int(common_params['num_threads'])


        # method specific execution preparation
        self.prepare()

        if self.num_cores != 1:
            pool_size,self.num_threads=self._optimizeParallel(self.num_jobs,self.num_threads,self.num_cores)
        else:
            pool_size = 1
            self.num_threads = 1
        # sanity check
        if self.num_jobs is None:
            raise ExecutionBaseException("Need to define num_jobs during prepare()")
        if len(self.task_list) < 1:
            raise ExecutionBaseException("Need to define task_list during prepare()")

        pool = mp.Pool(pool_size)
        self.results = pool.map(self.runEach, self.task_list)
        #pool.close()
        #pool.join()

        self.collect()
        self.writeReport()

        # clean-ups
        self.num_jobs = None
        self.task_list = []
        self.common_params = None
        self.method_params = None
        return self.returnVal

    # individual file processing code
    def runEach(self, task_params):
	raise NotImplementedMethod() 
	

    # collect parallel execution results
    # (optionally) it may call writeReport
    def collect(self):
        pass

    # will be called within collect if it has to
    def writeReport(self):
        pass
