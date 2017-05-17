# Test script for genome_util package - it should be launched from
# the root of the genome_util module, ideally just with 'make test', as
# it looks for a hardcoded relative path to find the 'test.cfg' file
import unittest
import json
import ConfigParser

from pprint import pprint
from string import Template

from subprocess import call

import sys
from biokbase.auth import Token
from os import environ
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsUtils.ReadsUtilsClient import ReadsUtils
import shutil
import os


class TestRNASeqMethodsSetupUserIndependent(unittest.TestCase):
    '''
    These tests set up a private workspace for the user that is running the 
    tests. The test data is uploaded to the workspace and shock before the
    tests are executed. This allows any user (with a valid test_token) to run 
    the tests.
    '''

    @classmethod
    def setUpClass(cls):
      super(TestRNASeqMethodsSetupUserIndependent, cls).setUpClass()

      print('_______BEGIN_TEST_SETUP_________')
      token = environ.get('KB_AUTH_TOKEN', None)

      if token is None:
          sys.stderr.write(
              "Error: Unable to run tests without authentication token!\n")
          sys.exit(1)

      token_file = open('test/script_test/token.txt', 'w')
      token_file.write(token)

      from biokbase.RNASeq.authclient import KBaseAuth as _KBaseAuth
      from biokbase.workspace.client import Workspace as Workspace
      try:
          from ConfigParser import ConfigParser  # py2
      except:
          from configparser import ConfigParser  # py3

      config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
      cls.cfg = {}
      config = ConfigParser()
      config.read(config_file)
      for nameval in config.items('KBaseRNASeq'):
          cls.cfg[nameval[0]] = nameval[1]
      auth_service_url = cls.cfg.get('auth-service-url',
                                     "https://kbase.us/services/authorization/Sessions/Login")
      ws_url = cls.cfg['ws_url']
      auth_service_url_allow_insecure = cls.cfg[
          'auth-service-url-allow-insecure']
      auth_client = _KBaseAuth(auth_service_url)
      user_id = auth_client.get_user(token)

      # update workspace ids in input json
      print('updating workspace id template in input json files...\n')

      INPUT_META_DATA_DIR = "test/downsized_test_metadata/"
      input_meta_data_files = ['create_rnaseq_sample_set_input.json',
                               'build_bowtie2_index_input.json',
                               'hisat2_input.json',
                               'tophat2_input.json',
                               'stringtie_input.json',
                               'cufflinks_input.json',
                               'cuffdiff_input.json',
                               'ballgown_input.json']
      for input_meta_data_file in input_meta_data_files:
          with open(INPUT_META_DATA_DIR + input_meta_data_file,
                    'r') as infile:
              input_meta_data = json.load(infile)

          # update workspace name in input.json files and write to work dir
          ws_id_t = Template(input_meta_data['params'][0]['ws_id'])
          cls.ws_id = ws_id_t.substitute(user_id=user_id)
          input_meta_data['params'][0]['ws_id'] = cls.ws_id

          with open('work/' + input_meta_data_file, 'w') as outfile:
              json.dump(input_meta_data, outfile)

      print('workspace_name: ' + cls.ws_id)

      # create workspace that is local to the user if it does not exist
      cls.ws = Workspace(url=ws_url, token=token, auth_svc=auth_service_url,
                         trust_all_ssl_certificates=auth_service_url_allow_insecure)
      try:
          ws_info = cls.ws.get_workspace_info({'workspace': cls.ws_id})
          print("workspace already exists: " + str(ws_info))
      except:
          ws_info = cls.ws.create_workspace(
              {'workspace': cls.ws_id, 'description': 'Workspace for ' + str(
                  input_meta_data['method'])})
          print("Created new workspace: " + str(ws_info))

      # upload genbank file
      print('uploading input data to workspace...')
      INPUT_DATA_DIR = "/kb/module/test/downsized_test_data/"
      TMP_INPUT_DATA_DIR = "/kb/module/work/tmp/"
      input_file_name = 'at_chrom1_section.gbk'
      input_data_path = INPUT_DATA_DIR + input_file_name

      print('input data path: ' + input_data_path)

      with open(INPUT_META_DATA_DIR + 'build_bowtie2_index_input.json', 'r') as infile:
          input_meta_data = json.load(infile)

      # data has to be copied to tmp dir so it can be seen by
      # GenomeFileUtil subjob running in a separate docker container
      tmp_input_data_path = TMP_INPUT_DATA_DIR + input_file_name
      shutil.copy(input_data_path, tmp_input_data_path)

      genbankToGenomeParams = {"file": {"path": tmp_input_data_path},
                               "genome_name": str(
                                   input_meta_data['params'][0]['reference']),
                               "workspace_name": cls.ws_id,
                               "source": str(
                                   input_meta_data['params'][0]['source']),
                               "release": str(
                                   input_meta_data['params'][0]['version']),
                               "generate_ids_if_needed": True,
                               "type": "User upload"
                               }
      gfu = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'], token=token,
                           auth_svc=auth_service_url)
      save_result = gfu.genbank_to_genome(genbankToGenomeParams)
      print('genbank_to_genome save result: ' + str(save_result))

      # upload downsized single reads
      ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'], token=token)
      input_reads = ['extracted_hy5_rep1.fastq',
                     'extracted_hy5_rep2.fastq',
                     'extracted_WT_rep1.fastq',
                     'extracted_WT_rep2.fastq'
                     ]
      for input_file_name in input_reads:
          input_data_path = INPUT_DATA_DIR+input_file_name
          tmp_input_data_path = TMP_INPUT_DATA_DIR + input_file_name
          shutil.copy(input_data_path, tmp_input_data_path)
          print('input data path: ' + input_data_path)
          result = ru.upload_reads({"fwd_file":tmp_input_data_path,
                                   "sequencing_tech":"Illumina",
                                   "wsname":cls.ws_id,
                                   "name":input_file_name})
          print('reads upload save result: '+str(result))

      print('_______END_TEST_SETUP_________')


    def test_a_build_bowtie2_index(self):
          print("\n\n----------- test BuildBowtie2index ----------")

          # make call
          out =call(["./bin/run_KBaseRNASeq.sh",
          "work/build_bowtie2_index_input.json",
          "work/build_bowtie2_index_output.json",
          "test/script_test/token.txt"])

          print("Error code: "+str(out));
  
          with open('work/build_bowtie2_index_output.json') as o:
                  output =json.load(o)
          
          pprint("Output: "+str(output))
          
          self.assertTrue(os.path.isfile('work/at_chrom1_section_index.zip'), 
                                                'did not find output index file')
          self.assertEqual(9078678, os.path.getsize("work/at_chrom1_section_index.zip"),
                                                'output index file size did not match')
    

    def test_b_CreateRNASeqSampleSet(self):
        print("\n\n----------- test CreateRNASeqSampleSet ----------")

        out = call(["./bin/run_KBaseRNASeq.sh",
                    "work/create_rnaseq_sample_set_input.json",
                    "work/create_rnaseq_sample_set.output.json",
                    "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: '+str(out))

        with open('work/create_rnaseq_sample_set.output.json') as o:
            output = json.load(o)
        pprint('Output: '+str(output))

        reads = self.__class__.ws.get_objects2(
                    {'objects': [{
                          'workspace': self.__class__.ws_id,
                          'name': output['result'][0]['sampleset_id']}]})
        self.assertEqual('KBaseRNASeq.RNASeqSampleSet-4.0', reads['data'][0]['info'][2],
                                         "output sampleset object type did not match")


    def test_c_hisat2(self):
        print("\n\n----------- test HiSat2 ----------")

        # dependency: create a reads sampleset from the reads
        # that were uploaded in the setup
        #self.test_CreateRNASeqSampleSet()

        out =call(["./bin/run_KBaseRNASeq.sh",
        "work/hisat2_input.json",
        "work/hisat2_output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: '+str(out))

        with open('work/hisat2_output.json') as o:
                output =json.load(o)
        pprint('Output: '+str(output))

        alignment_set = self.__class__.ws.get_objects2({'objects': [
            {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqAlignmentSet-9.0', 
                                        alignment_set['data'][0]['info'][2],
                                        "output alignment set object type did not match")


    def test_d_tophat2(self):
        print("\n\n----------- test TophatCall ----------")

        # dependency: create a genome index from the genome
        # that was uploaded in the setup
        #self.test_build_bowtie2_index()

        out =call(["run_KBaseRNASeq.sh",
        "work/tophat2_input.json",
        "work/tophat2_output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: '+str(out))

        with open('work/tophat2_output.json') as o:
                output =json.load(o)
        pprint('Output: '+str(output))

        alignment_set = self.__class__.ws.get_objects2({'objects': [
            {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqAlignmentSet-9.0',
                         alignment_set['data'][0]['info'][2],
                         "output alignment set object type did not match")

    def test_e_stringtie(self):
        print("\n\n----------- test Stringtie ----------")

        # dependency: create a reads alignment set
        #elf.test_hisat2()

        out = call(["run_KBaseRNASeq.sh",
                    "work/stringtie_input.json",
                    "work/stringtie_output.json",
                    "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: ' + str(out))

        with open('work/stringtie_output.json') as o:
            output = json.load(o)
        pprint('Output: ' + str(output))

        expression_set = self.__class__.ws.get_objects2({'objects': [
           {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqExpressionSet-6.0',
                         expression_set['data'][0]['info'][2],
                         "output expression set object type did not match")

    def test_f_cufflinks(self):
        print("\n\n----------- test CufflinksCall ----------")

        # dependency: create a reads alignment set
        #self.test_hisat2()

        out =call(["run_KBaseRNASeq.sh",
        "work/cufflinks_input.json",
        "work/cufflinks_output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: ' + str(out))

        with open('work/stringtie_output.json') as o:
            output = json.load(o)
        pprint('Output: ' + str(output))

        expression_set = self.__class__.ws.get_objects2({'objects': [
            {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqExpressionSet-6.0',
                         expression_set['data'][0]['info'][2],
                         "output expression set object type did not match")


    def test_g_cuffdiff(self):
        print("\n\n----------- test CuffdiffCall ----------")

        out =call(["run_KBaseRNASeq.sh",
        "work/cuffdiff_input.json",
        "work/cuffdiff_output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: ' + str(out))

        with open('work/cuffdiff_output.json') as o:
            output = json.load(o)
        pprint('Output: ' + str(output))

        differential_expression = self.__class__.ws.get_objects2({'objects': [
            {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqDifferentialExpression-5.0',
                         differential_expression['data'][0]['info'][2],
                         "output differential expression object type did not match")
    '''

    def test_h_ballgown(self):
        print("\n\n----------- test DiffExpCallforBallgown ----------")

        out = call(["run_KBaseRNASeq.sh",
                    "work/ballgown_input.json",
                    "work/ballgown_output.json",
                    "test/script_test/token.txt"])

        # print error code of Implementation
        print('ErrorCode: ' + str(out))

        with open('work/cuffdiff_output.json') as o:
            output = json.load(o)
        pprint('Output: ' + str(output))

        differential_expression = self.__class__.ws.get_objects2({'objects': [
            {'workspace': self.__class__.ws_id,
             'name': output['result'][0]['output']}]})
        self.assertEqual('KBaseRNASeq.RNASeqDifferentialExpression-5.0',
                         differential_expression['data'][0]['info'][2],
                         "output differential expression object type did not match")



'''
# These tests use private workspaces that is accessible by 'kbasetest' user
# Before all the tests, read the config file and get a user token and
# save it to a file used by the main service script
class TestRNASeqMethodsSetup(unittest.TestCase):

  @classmethod
  def setUp(cls):
    token = environ.get('KB_AUTH_TOKEN', None)

    if token is None:
        sys.stderr.write("Error: Unable to run tests without authentication token!\n")
        sys.exit(1)

    token_file = open('test/script_test/token.txt', 'w')
    token_file.write(token)


# Define all our other test cases here
class TestRNASeqMethods(TestRNASeqMethodsSetup): 
# def test_0(self):
#        print("\n\n----------- test HiSat2 ----------")
#
#        out =call(["./bin/run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis2_hisat2_single_input.json",
#        "test/script_test/test_analysis2_hisat2_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis2_hisat2_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
# def test_a(self):
#        print("\n\n----------- test SetupRNASeqAnalysis ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_SetupRNASeq_analysis.json",
#        "test/script_test/test_SetupRNASeq_analysis_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_SetupRNASeq_analysis_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
# def test_b(self):
#       print("\n\n----------- test BuildBowtie2index ----------")
#
#       out =call(["run_KBaseRNASeq.sh",
#       "test/script_test/build_bowtie2_index_input.json",
#       "test/script_test/new_genome_build_bowtie_output.json",
#       "test/script_test/token.txt"])
#
#       # print error code of Implementation
#       print(out);
#
#       with open('test/script_test/new_genome_build_bowtie_output.json') as o:
#               output =json.load(o)
#       pprint(output)

# def test_Bowtie2Call(self):
#        print("\n\n----------- test Bowtie2Call ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis2_bowtie2_input.json",
#        "test/script_test/test_analysis2_bowtie2_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis2_bowtie2_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
## def test_Bowtie2Call(self):
##        print("\n\n----------- test Bowtie2Call ----------")
##
##        out =call(["run_KBaseRNASeq.sh",
##        "test/script_test/test_analysis1_bowtie2_input.json",
##        "test/script_test/test_analysis1_bowtie2_output.json",
##        "test/script_test/token.txt"])
##
##        # print error code of Implementation
##        print(out);
##
##        with open('test/script_test/test_analysis1_bowtie2_output.json') as o:
##                output =json.load(o)
##        pprint(output)
#
# def test_c(self):
#        print("\n\n----------- test TophatCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis_tophat_single.json",
#        "test/script_test/tophat_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/tophat_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
## def test_d(self):
#        print("\n\n----------- test TophatCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis_tophat1.json",
#        "test/script_test/tophat_output1.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/tophat_output1.json') as o:
#                output =json.load(o)
#        pprint(output)
#
# def test_e(self):
#        print("\n\n----------- test CufflinksCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis1_cufflinks_input.json",
#        "test/script_test/test_analysis1_cufflinks_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis1_cufflinks_output.json') as o:
#                output =json.load(o)
#        pprint(output)

# def test_ee(self):
#        print("\n\n----------- test StringtieCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis1_stringtie_input.json",
#        "test/script_test/test_analysis1_stringtie_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis1_stringtie_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
#
# def test_f(self):
#        print("\n\n----------- test CufflinksCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis2_cufflinks_input.json",
#        "test/script_test/test_analysis2_cufflinks_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis2_cufflinks_output.json') as o:
#                output =json.load(o)
#        pprint(output)
# 
# def test_g(self):
#        print("\n\n----------- test CuffmergeCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_cuffmerge_input.json",
#        "test/script_test/test_cuffmerge_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_cuffmerge_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#
#
# def test_h(self):
#        print("\n\n----------- test CuffdiffCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_cuffdiff_input.json",
#        "test/script_test/test_cuffdiff_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_cuffdiff_output.json') as o:
#                output =json.load(o)
#        pprint(output)
 def test_i(self):
        print("\n\n----------- test DiffExpCallforBallgown ----------")

        out =call(["run_KBaseRNASeq.sh",
        "test/script_test/test_ballgown_main.input.json",
        "test/script_test/test_ballgown_main.output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print(out);

        with open('test/script_test/test_ballgown_main.output.json') as o:
                output =json.load(o)
        pprint(output)
'''

#start the tests if run as a script
if __name__ == '__main__':
    unittest.main()
    


