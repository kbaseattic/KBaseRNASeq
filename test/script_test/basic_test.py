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

    '''
    #######################################
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
    auth_service_url_allow_insecure = cls.cfg['auth-service-url-allow-insecure']
    auth_client = _KBaseAuth(auth_service_url)
    print('>>>>>>>>>>>>>>>>>>>>AUTH_SERVICE_URL: ' + auth_service_url)
    print('>>>>>>>>>>>>>>>>>>>>TOKEN: ' + token)
    user_id = auth_client.get_user(token)
    print('>>>>>>>>>>>>>>>>>>>>USERID: '+user_id)

    ws = Workspace(url=ws_url, token=token, auth_svc=auth_service_url,
                   trust_all_ssl_certificates=auth_service_url_allow_insecure)

    filename = "test/script_test/input_meta_data/hisat2_input1.json"
    INPUT_DATA_DIR = "test/script_test/input_data"
    with open(filename, 'r') as infile:
        input_meta_data = json.load(infile)

    # create workspace that is local to the user if it does not exist
    workspace_name = str(input_meta_data['params'][0]['ws_id'])
    print('workspace_name: ' + workspace_name)

    try:
        ws_info = ws.get_workspace_info({'workspace': workspace_name})
        print("workspace already exists: " + str(ws_info))
    except:
        ws_info = ws.create_workspace(
            {'workspace': workspace_name, 'description': 'Workspace for ' + str(input_meta_data['method'])})
        print("Created new workspace: " + str(ws_info))

    input_data_filename = 'test/script_test/input_data/Athaliana_PhytozomeV11_TAIR10.json'
    print('input data filename: ' + input_data_filename)

    with open(input_data_filename, 'r') as infile:
        input_data = json.load(infile)

    # upload data (no effect if data already exists in workspace)
    print('uploading input data to workspace')
    ws.save_objects(
        {'workspace': workspace_name, 'objects': [{'type': "KBaseGenomes.Genome",
                                                   'data': input_data,
                                                   'name': 'Athaliana_PhytozomeV11_TAIR10'}]})

    input_data_filename = 'test/script_test/input_data/Ath_WT_Hy5_sampleset.json'
    print('input data filename: ' + input_data_filename)

    with open(input_data_filename, 'r') as infile:
        input_data = json.load(infile)

    # upload data (no effect if data already exists in workspace)
    print('uploading input data to workspace')
    ws.save_objects(
        {'workspace': workspace_name, 'objects': [{'type': 'KBaseRNASeq.RNASeqSampleSet',
                                                   'data': input_data,
                                                   'name': 'Ath_WT_Hy5_sampleset'}]})

    print('ws objects: ' + str(ws.list_objects({'workspaces': [workspace_name]})))
    #######################################
    '''


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
#       "test/script_test/new_genome_build_bowtie_input.json",
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

#start the tests if run as a script
if __name__ == '__main__':
    unittest.main()
