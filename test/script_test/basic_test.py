# Test script for genome_util package - it should be launched from
# the root of the genome_util module, ideally just with 'make test', as
# it looks for a hardcoded relative path to find the 'test.cfg' file
import unittest
import json
import ConfigParser

from pprint import pprint

from subprocess import call

from biokbase.auth import Token

# Before all the tests, read the config file and get a user token and
# save it to a file used by the main service script
class TestRNASeqMethodsSetup(unittest.TestCase):
  def setUp(self):
    config = ConfigParser.RawConfigParser()
    config.read('test/test.cfg')
    user_id = config.get('KBaseRNASeqTest','user')
    password = config.get('KBaseRNASeqTest','password')
    token = Token(user_id=user_id, password=password)
    token_file = open('test/script_test/token.txt','w')
    token_file.write(token.token)

# Define all our other test cases here
class TestRNASeqMethods(TestRNASeqMethodsSetup): 

# def test_SetupRNASeqAnalysis(self):
#        print("\n\n----------- test SetupRNASeqAnalysis ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/m_setupRNASeq.json",
#        "test/script_test/bowtie_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/bowtie_output.json') as o:
#                output =json.load(o)
#        pprint(output)
#  
# def test_associateReads(self):
#       print("\n\n----------- test associateReads ----------")
#
#       out =call(["run_KBaseRNASeq.sh",
#       "test/script_test/test_associateReads_input.json",
#       "test/script_test/test_associateReads_output.json",
#       "test/script_test/token.txt"])
#
#       # print error code of Implementation
#       print(out);
#
#       with open('test/script_test/test_associateReads_output.json') as o:
#               output =json.load(o)
#       pprint(output)
#
# def test_BuildBowtie2Index(self):
#       print("\n\n----------- test BuildBowtie2index ----------")
#
#       out =call(["run_KBaseRNASeq.sh",
#       "test/script_test/bowtie_input.json",
#       "test/script_test/bowtie_output.json",
#       "test/script_test/token.txt"])
#
#       # print error code of Implementation
#       print(out);
#
#       with open('test/script_test/bowtie_output.json') as o:
#               output =json.load(o)
#       pprint(output)
#
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

# def test_Bowtie2Call(self):
#        print("\n\n----------- test Bowtie2Call ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis1_bowtie2_input.json",
#        "test/script_test/test_analysis1_bowtie2_output.json",
#        "test/script_test/token.txt"])
#
#        # print error code of Implementation
#        print(out);
#
#        with open('test/script_test/test_analysis1_bowtie2_output.json') as o:
#                output =json.load(o)
#        pprint(output)

# def test_TophatCall(self):
#        print("\n\n----------- test TophatCall ----------")
#
#        out =call(["run_KBaseRNASeq.sh",
#        "test/script_test/test_analysis_tophat.json",
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
# def test_TophatCall1(self):
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
#
# def test_getAlignmentStats(self):
#       print("\n\n----------- test Get Alignment Statistics ----------")
#
#       out =call(["run_KBaseRNASeq.sh",
#       "test/script_test/test_getAlignmentStats_input.json",
#       "test/script_test/test_getAlignmentStats_output.json",
#       "test/script_test/token.txt"])
#
#       # print error code of Implementation
#       print(out);
#
#       with open('test/script_test/test_getAlignmentStats_output.json') as o:
#               output =json.load(o)
#       pprint(output)

 def test_CufflinksCall(self):
        print("\n\n----------- test CufflinksCall ----------")

        out =call(["run_KBaseRNASeq.sh",
        "test/script_test/test_analysis1_cufflinks_input.json",
        "test/script_test/test_analysis1_cufflinks_output.json",
        "test/script_test/token.txt"])

        # print error code of Implementation
        print(out);

        with open('test/script_test/test_analysis1_cufflinks_output.json') as o:
                output =json.load(o)
        pprint(output)

# def test_CufflinksCall1(self):
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
##  
# def test_CuffmergeCall(self):
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


# def test_CuffdiffCall(self):
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

#start the tests if run as a script
if __name__ == '__main__':
    unittest.main()
