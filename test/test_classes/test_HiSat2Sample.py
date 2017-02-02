# -*- coding: utf-8 -*-
#BEGIN_HEADER
import unittest
import os
import json
import time
import requests
from ConfigParser import ConfigParser
from pprint import pprint
from biokbase.workspace.client import Workspace
requests.packages.urllib3.disable_warnings()
from biokbase.RNASeq.KBaseRNASeqImpl import KBaseRNASeq
from biokbase.RNASeq.HiSat2Sample import HiSat2Sample

class KBaseRNASeq_Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
	callback_url = 'http://172.17.0.1:54375'
        cls.ctx = {'token': token}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('KBaseRNASeq'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['ws_url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = KBaseRNASeq(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_KBaseRNASeq_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_your_method(self):
	## Test Create RNASeqSampleSet ###	
	## Test Align Reads Using Tophat - Set mode ###
	## Test Align Reads Using Tophat - Single mode ###
	
	## Test Align Reads Using Bowtie2 - Set mode ###
	## Test Align Reads Using Bowtie2 - Single mode ###

	## Test Align Reads Using HiSat2 - Set mode  ###
	## Test Align Reads Using HiSat2 - Single mode  ###
	hisat2params_singlemode = {
	 			"ws_id" : "srividya22:1474317595261",
	 			"sampleset_id" : "Ecoli_SE_8083_R1.fastq" ,
	 			"genome_id" : "Escherichia_coli_K12",
	 			"alignment_type" : null,
	 			"quality_score" : null,
	 			"np" : 1,
	 			"minins" : 0,
	 			"tailor_alignments" : "dta-cufflinks",
	 			"maxins" : 500}
	pprint(dir(self.__class__))
	try:
        	ret = self.getImpl().Hisat2Call(self.getContext(),hisat2params_singlemode)
	except Exception, e:
		raise Exception(e)
	## Test Assemble Transcripts Using Cufflinks - Set mode  ###
	## Test Assemble Transcripts  Using Cufflinks - Set mode  ###
	
	## Test Assemble Transcripts Using StringTie - Set mode  ###
	## Test Assemble Transcripts  Using StringTie - Set mode  ###
	
	## Test Identify Differential Expression Using Cuffdiff ###
	## Test Identify Differential Expression Using DiffExpforBallgown  ###

        pass 
