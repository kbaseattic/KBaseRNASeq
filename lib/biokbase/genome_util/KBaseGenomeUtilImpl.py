#BEGIN_HEADER

import simplejson
import sys
import os
import glob
import json
import logging
import time
import subprocess
from pprint import pprint
import script_util
from biokbase.workspace.client import Workspace
from biokbase.auth import Token

try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

_KBaseGenomeUtil__DATA_VERSION = "0.5"


class KBaseGenomeUtilException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)



no_rst = """{
    "BlastOutput_db": "NoDB", 
    "BlastOutput_iterations": {
        "Iteration": [
            {
                "Iteration_hits": {
                    "Hit": []
                }, 
                "Iteration_iter-num": "1", 
                "Iteration_message": "ERR_MSG", 
                "Iteration_query-ID": "lcl|1_0", 
                "Iteration_query-len": "QRY_LNGTH", 
                "Iteration_stat": {
                    "Statistics": {
                        "Statistics_db-len": "1331648", 
                        "Statistics_db-num": "4280", 
                        "Statistics_eff-space": "2.6633e+06", 
                        "Statistics_entropy": "0.14", 
                        "Statistics_hsp-len": "0", 
                        "Statistics_kappa": "0.041", 
                        "Statistics_lambda": "0.267"
                    }
                }
            }, 
            {
                "Iteration_hits": {
                    "Hit": []
                }, 
                "Iteration_iter-num": "1", 
                "Iteration_stat": {
                    "Statistics": {
                        "Statistics_db-len": "1331648", 
                        "Statistics_db-num": "4280", 
                        "Statistics_eff-space": "2.6633e+06", 
                        "Statistics_entropy": "0.14", 
                        "Statistics_hsp-len": "0", 
                        "Statistics_kappa": "0.041", 
                        "Statistics_lambda": "0.267"
                    }
                }
            }
        ]
    }, 
    "BlastOutput_param": {
        "Parameters": {
            "Parameters_expect": "0.05", 
            "Parameters_filter": "F", 
            "Parameters_gap-extend": "1", 
            "Parameters_gap-open": "11", 
            "Parameters_matrix": "BLOSUM62"
        }
    }, 
    "BlastOutput_program": "error", 
    "BlastOutput_query-ID": "error", 
    "BlastOutput_query-def": "error", 
    "BlastOutput_query-len": "na", 
    "BlastOutput_reference": "error", 
    "BlastOutput_version": "error"
}"""


#END_HEADER


class KBaseGenomeUtil:
    '''
    Module Name:
    KBaseGenomeUtil

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    # Config variables that SHOULD get overwritten in the constructor
    __TEMP_DIR = 'index'
    __WS_URL = 'https://ci.kbase.us/services/ws'
    __HS_URL = 'https://ci.kbase.us/services/handle_service'
    __SHOCK_URL = 'https://ci.kbase.us/services/shock-api/'
    __BLAST_DIR = 'blast'
    __GENOME_FA = 'genome.fa'
    __ANNO_JSON = 'annotation.json'
    __QUERY_FA = 'query.fa'
    __INDEX_CMD = 'formatdb'
    __BLAST_CMD = 'blastall'
    __BLAST_OUT = 'result.txt'
    __INDEX_ZIP = 'index.zip'
    __SVC_USER = 'kbasetest'
    __SVC_PASS = ''
    __LOGGER = None
    __ERR_LOGGER = None


    __INDEX_TYPE = {'blastp' : 'protein_db',
                    'blastx' : 'protein_db',
                    'blastn' : 'transcript_db',
                    'tblastn' : 'transcript_db',
                    'tblastx' : 'transcript_db'}


    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        # This is where config variable for deploy.cfg are available
        #pprint(config)
        if 'ws_url' in config:
              self.__WS_URL = config['ws_url']
        if 'hs_url' in config:
              self.__HS_URL = config['hs_url']
        if 'shock_url' in config:
              self.__SHOCK_URL = config['shock_url']
        if 'temp_dir' in config:
              self.__TEMP_DIR = config['temp_dir']
        if 'blast_dir' in config:
              self.__BLAST_DIR = config['blast_dir']
        if 'genome_input_fa' in config:
              self.__GENOME_FA = config['genome_input_fa']
        if 'query_fa' in config:
              self.__QUERY_FA = config['query_fa']
        if 'svc_user' in config:
              self.__SVC_USER = config['svc_user']
        if 'svc_pass' in config:
              self.__SVC_PASS = config['svc_pass']

        # logging
        self.__LOGGER = logging.getLogger('GenomeUtil')
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

    def index_genomes(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN index_genomes
        user_token=ctx['token']
        svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        hs = HandleService(url=self.__HS_URL, token=user_token)
        gs = {'elements' : {}}
        try:
            self.__LOGGER.info( "Preparing Target FA")
         
            blast_dir =self.__BLAST_DIR
            if os.path.exists(blast_dir):
                files=glob.glob("%s/*" % blast_dir)
                for f in files: os.remove(f)
            if not os.path.exists(blast_dir): os.makedirs(blast_dir)
          
         
               
            target_nt_fn = "%s/%s_nt.fa" %( blast_dir, params['blastindex_name'])
            target_aa_fn = "%s/%s_aa.fa" %( blast_dir, params['blastindex_name'])
         
            try:
              target_nt=open(target_nt_fn,'w')
              target_aa=open(target_aa_fn,'w')
            except:
              self.__LOGGER.error("Couldn't open file")
              raise KBaseGenomeUtilException("Backend awe client error: Couldn't open files")
         
            have_nt_seq = False
            have_aa_seq = False
         
         
         
            # Iterate one at a time to cope with main memory limit for euk genomes
            for genome_id in params['genome_ids']: 
         
                try:
                    obj_infos = ws_client.get_object_info_new({"objects": [{'name':genome_id, # replace `0' with loop
                                                               'workspace': params['ws_id']}]})
                except:
                    self.__LOGGER.error("Couldn't retrieve %s:%s from the workspace" %(params['ws_id'],genome_id))
                    raise KBaseGenomeUtilException("Couldn't retrieve %s:%s from the workspace" %(params['ws_id'],genome_id))
                     
         
                if len(obj_infos) < 1:
                    self.__LOGGER.error("Couldn't find %s:%s from the workspace" %(params['ws_id'],genome_id))
                    continue
                    #err_msg += "Workspace error: Couldn't find %s:%s from the workspace\n" %(params['ws_id'],genome_id)                
                    # we can continue due to multiple genomes
                    #raise Exception("Couldn't find %s:%s from the workspace" %(params['ws_id'],genome_id)) 
         
                ref_id = "{0}/{1}/{2}".format(obj_infos[0][6],obj_infos[0][0],obj_infos[0][4])
                gs['elements'][genome_id] = [ref_id]
               
                self.__LOGGER.info( "Downloading genome object from workspace {0}".format(ref_id))
               
                # TODO: make the following procedures to be loop for each genome_ids 
                try:
                    genome_list=ws_client.get_object_subset([{'name':genome_id, # replace `0' with loop
                                                              'workspace': params['ws_id'], 
                                                              'included':['features']}])
                    #genome_list=ws_client.get_objects([{'name':genome_id, # replace `0' with loop
                    #                                          'workspace': params['ws_id']}])
                    genome = genome_list[0]
                except Exception, e:
                    raise KBaseGenomeUtilException("Failed to download genome object itself even though we got the object information")
     
               
               
                self.__LOGGER.info( "Dumping seq for %s" % genome_id)
                # Dump genome sequences
                check_seq=0
                #extract protein sequences from the genome object
                try:
                    for gene in genome['data']['features']:
                          #>kb.g.1234.CDS.1234#At1g3333 amalase...
                          function = "NA"
                          aliases = "NA"
                          if 'function' in gene: 
                              function = gene['function']
                          if 'aliases' in gene: aliases = ",".join(gene['aliases'])
                          if 'protein_translation' in gene:
                                target_aa.write(">%s#%s#%s#%s\n%s\n" % (gene['id'], ref_id, aliases, function, gene['protein_translation']))
                                have_aa_seq = True
                          if 'dna_sequence' in gene:
                                target_nt.write(">%s#%s#%s#%s\n%s\n" % (gene['id'], ref_id, aliases, function, gene['dna_sequence']))
                                have_nt_seq = True
                except Exception as e:
                    raise KBaseGenomeUtilException("Failed to dump target sequence for genome : %s" % genome_id)
            try:
                target_nt.close()
                target_aa.close()
            except Exception as e:
                raise KBaseGenomeUtilException("Failed to close sequence files")
                
                
               
            if not have_nt_seq :
                self.__LOGGER.info("The genome objects do not contain any dna sequences!")
            if not have_aa_seq :
                self.__LOGGER.info("The genome objects do not contain any amino acid sequences!")
         
            index_type = 'none'
               
            if have_nt_seq :
                try:
                    cmdstring="%s -i %s -p F" %(self.__INDEX_CMD, target_nt_fn)
                    # TODO: replace it to subprocess.Popen
                    tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
                    stdout, stderr = tool_process.communicate()
                    
                    if stdout is not None and len(stdout) > 0:
                        self.__LOGGER.info(stdout)
                    
                    if stderr is not None and len(stderr) > 0:
                        self.__LOGGER.error("Indexing error: " + stderr)
                        raise KBaseGenomeUtilException("Indexing error: " + stderr)
                except Exception, e:
                    raise KBaseGenomeUtilException("Failed to run indexing program (%s) : %s " %(self.__INDEX_CMD, e))
                   
                index_type = 'nucleotide'
               
            if have_aa_seq :
                try:
                    cmdstring="%s -i %s -p T" %(self.__INDEX_CMD, target_aa_fn)
                    # TODO: replace it to subprocess.Popen
                    tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
                    stdout, stderr = tool_process.communicate()
                    
                    if stdout is not None and len(stdout) > 0:
                        self.__LOGGER.info(stdout)
                    
                    if stderr is not None and len(stderr) > 0:
                        self.__LOGGER.error("Indexing error: " + stderr)
                        raise KBaseGenomeUtilException("Indexing error: " + stderr)
                except Exception, e:
                    raise KBaseGenomeUtilException("Failed to run indexing program (%s) : %s " %(self.__INDEX_CMD, e))
                if index_type == 'nucleotide': index_type = 'both'
                else: index_type = 'protein'
            
            #os.remove(target_nt_fn)
            #os.remove(target_aa_fn)
         
            # compress
            try: 
                script_util.zip_files(self.__LOGGER, blast_dir, "%s.zip" % params['blastindex_name'])
            except Exception, e:
                raise KBaseGenomeUtilException("Failed to compress the index: %s" %(e))
               
            try: 
                handle = hs.upload("%s.zip" % (params['blastindex_name']))
            except Exception, e:
                raise KBaseGenomeUtilException("Failed to upload the index: %s" %(e))

            bi = {'handle' : handle, 'genome_set' : gs, 'index_type' : index_type, 'index_program' : params['index_program']}
            if 'description' in params: bi['description'] = params['description']
         
            if index_type == 'none': 
                err_msg = 'No sequences were indexed'
                bi['description'] = err_msg
                res= ws_client.save_objects(
                    {"workspace":params['ws_id'],
                    "objects": [{
                        "type":"GenomeUtil.BlastIndex",
                        "data":bi,
                        "meta" : {'err_msg' : err_msg},
                        "name":params['blastindex_name']}
                    ]})
            else:
                res= ws_client.save_objects(
                    {"workspace":params['ws_id'],
                    "objects": [{
                        "type":"GenomeUtil.BlastIndex",
                        "data":bi,
                        "name":params['blastindex_name']}
                    ]})
            returnVal = { 'blastindex_ref' : "%s/%s" % (params['ws_id'], params['blastindex_name']) }
            if index_type == 'none':
                returnVal['err_msg'] = err_msg
        except MemoryError, e:
            handle = hs.new_handle()
            bi = {'handle' : handle, 'genome_set' : gs, 'index_type' : 'none', 'index_program' : params['index_program']}
            err_msg = 'Not enough main memory: please use smaller number of genomes only'
            bi['description'] = err_msg
            returnVal = {'err_msg' : err_msg }
            res= ws_client.save_objects(
                {"workspace":params['ws_id'],
                "objects": [{
                    "type":"GenomeUtil.BlastIndex",
                    "data":bi,
                    "meta" : {'err_msg' : err_msg},
                    "name":params['blastindex_name']}
                ]})
            
        except Exception, e:
            handle = hs.new_handle()
            bi = {'handle' : handle, 'genome_set' : gs, 'index_type' : 'none', 'index_program' : params['index_program']}
            err_msg = str(e)
            bi['description'] = err_msg
            returnVal = {'err_msg' : err_msg }
            res= ws_client.save_objects(
                {"workspace":params['ws_id'],
                "objects": [{
                    "type":"GenomeUtil.BlastIndex",
                    "data":bi,
                    "meta" : {'err_msg' : err_msg},
                    "name":params['blastindex_name']}
                ]})
        
        #END index_genomes

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method index_genomes return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def blast_against_genome(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN blast_against_genome

        # TODO: Rename blast_search

        try:
           self.__LOGGER.info( "Preparing FA")
           if len(params['query']) > 5:
               sequence=params['query']
           else:
               self.__LOGGER.error("The input sequence is too short!")
               raise KBaseGenomeUtilException("The input sequence is too short!")
        
           if not os.path.exists(self.__TEMP_DIR): os.makedirs(self.__TEMP_DIR)
         
           #print "generate input file for query sequence\n"
           query_fn = "%s/%s" %(self.__TEMP_DIR, self.__QUERY_FA)
           target=open(query_fn,'w')
           if sequence.startswith(">"):
             target.write(sequence)
           else:
             seqes = sequence.split("\n")
             for i in range(len(seqes)):
               target.write(">query_seq_%d\n" %(i))
               target.write(seqes[i])
           target.close()
         
           user_token=ctx['token']
           svc_token = Token(user_id=self.__SVC_USER, password=self.__SVC_PASS).token
           ws_client=Workspace(url=self.__WS_URL, token=user_token)
        
        
           err_msg = ""
        
           blast_dir =self.__BLAST_DIR
           if os.path.exists(blast_dir):
               files=glob.glob("%s/*" % blast_dir)
               for f in files: os.remove(f)
           if not os.path.exists(blast_dir): os.makedirs(blast_dir)
           target_fn = "%s/%s" %( blast_dir, self.__GENOME_FA)
           if 'target_seqs' in params:
               # let's build index directly and throw away
               sequence = params['target_seqs']
        
               target=open(target_fn,'w')
               if sequence.startswith(">"):
                 target.write(sequence)
               else:
                 seqes = sequence.split("\n")
                 for i in range(len(seqes)):
                   target.write(">target_seq_%d\n" %(i))
                   target.write(seqes[i])
               target.close()
            
               if(self.__INDEX_TYPE[params['blast_program']]  == 'protein_db'):
                   formatdb_type='T'
               elif(self.__INDEX_TYPE[params['blast_program']]  == 'transcript_db'):
                   formatdb_type='F'
               else:
                   self.__LOGGER.error("{0} is not yet supported".format(params['blast_program']))
                   raise KBaseGenomeUtilException("{0} is not yet supported".format(params['blast_program']))
               cmdstring="%s -i %s -p %s -o T" %(self.__INDEX_CMD, target_fn, formatdb_type)
               # TODO: replace it to subprocess.Popen
               tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
               stdout, stderr = tool_process.communicate()
   
               if stdout is not None and len(stdout) > 0:
                   self.__LOGGER.info(stdout)
   
               if stderr is not None and len(stderr) > 0:
                   self.__LOGGER.error("Index error: " + stderr)
                   raise KBaseGenomeUtilException("Index error: " + stderr)
        
           else:
               try:
                   blast_indexes=ws_client.get_object_subset([{'name':params['blastindex_name'],
                                                             'workspace': params['ws_id'], 
                                                             'included':['handle', 'index_type']}])
               except:
                   self.__LOGGER.error("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['blastindex_name']))
                   raise KBaseGenomeUtilException("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['genome_ids'][0]))
                   
               if len(blast_indexes) < 1:
                   self.__LOGGER.error("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['blastindex_name']))
                   raise KBaseGenomeUtilException("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['genome_ids'][0]))
        
               
               # TODO: Add err handling
               zip_fn = blast_indexes[0]['data']['handle']['file_name']
               target_fn = "%s/%s" %(blast_dir, zip_fn[:-4]) # remove '.zip'
        
               if(self.__INDEX_TYPE[params['blast_program']]  == 'protein_db'):
                   target_fn += '_aa.fa'
                   if blast_indexes[0]['data']['index_type'] == 'none' or blast_indexes[0]['data']['index_type'] == "nucleotide":
                       self.__LOGGER.error("The index object does not contain amino acid sequence indexes")
                       raise KBaseGenomeUtilException("The index object does not contain amino acid sequence indexes")                    
               elif(self.__INDEX_TYPE[params['blast_program']]  == 'transcript_db'):
                   target_fn += '_nt.fa'
                   if blast_indexes[0]['data']['index_type'] == 'none' or blast_indexes[0]['data']['index_type'] == "protein":
                       self.__LOGGER.error("The index object does not contain nucleotide sequence indexes")
                       raise KBaseGenomeUtilException("The index object does not contain nucleotide sequence indexes")                    
               else:
                   self.__LOGGER.error("{0} is not yet supported".format(params['blast_program']))
                   raise KBaseGenomeUtilException("{0} is not yet supported".format(params['blast_program']))
        
               # TODO: Add err handling
               zip_fn = blast_indexes[0]['data']['handle']['file_name']
               #pprint(blast_indexes[0])
              
               self.__LOGGER.info("Downloading the genome index")
               #hs = HandleService(url=self.__HS_URL, token=user_token)
               try:
                   script_util.download_file_from_shock(self.__LOGGER,
                                   shock_service_url= blast_indexes[0]['data']['handle']['url'],
                                   shock_id= blast_indexes[0]['data']['handle']['id'],
                                   filename= blast_indexes[0]['data']['handle']['file_name'],
                                   directory= '.',
                                   token = user_token)
               except Exception, e:
                   self.__LOGGER.error("Downloading error from shock: Please contact help@kbase.us")
                   raise KBaseGenomeUtilException("Downloading error from shock: Please contact help@kbase.us")
               try:
                   script_util.unzip_files(self.__LOGGER, zip_fn, blast_dir)
               except Exception, e:
                   self.__LOGGER.error("Unzip indexfile error: Please contact help@kbase.us")
                   raise KBaseGenomeUtilException("Unzip indexfile error: Please contact help@kbase.us")
        
        
           self.__LOGGER.info( "Searching...")
           #blast search
           cmdstring="%s -p %s -i %s -m 7 -o %s -d %s -e %s" % (self.__BLAST_CMD, params['blast_program'], query_fn, self.__BLAST_OUT, target_fn, params['e-value'])
        
           if 'gap_opening_penalty' in params:
             cmdstring += " -G %s" %(params['gap_opening_penalty'])
        
           if 'gap_extension_penalty' in params:
             cmdstring += " -E %s" %(params['gap_extension_penalty'])
        
           if 'nucleotide_match_reward' in params:
             cmdstring += " -r %s" %(params['nucleotide_match_reward'])
        
           if 'nucleotide_mismatch_penalty' in params:
             cmdstring += " -q %s" %(params['nucleotide_mismatch_penalty'])
        
           if 'word_size' in params:
             cmdstring += " -W %s" %(params['word_size'])
        
           if 'maximum_alignment_2show' in params:
             cmdstring += " -b %s" %(params['maximum_alignment_2show'])
        
           if 'substitution_matrix' in params and params['substitution_matrix'] != 'Default':
             cmdstring += " -M %s" %(params['substitution_matrix'])
        
           if 'mega_blast' in params:
             cmdstring += " -n %s" %(params['mega_blast'])
        
           if 'gapped_alignment' in params:
             cmdstring += " -g %s" %(params['gapped_alignment'])
        
           if 'filter_query_seq' in params:
             cmdstring += " -F %s" %(params['filter_query_seq'])
        
           if 'extending_hits' in params:
             cmdstring += " -f %s" %(params['extending_hits'])
        
           if 'maximum_seq_2show' in params:
             cmdstring += " -v %s" %(params['maximum_seq_2show'])
        
           # TODO: replace it to subprocess.Popen
           #print cmdstring
           try: 
               tool_process = subprocess.Popen(cmdstring, stderr=subprocess.PIPE, shell=True)
               stdout, stderr = tool_process.communicate()
   
               if stdout is not None and len(stdout) > 0:
                   self.__LOGGER.info(stdout)
   
               if stderr is not None and len(stderr) > 0:
                   self.__LOGGER.error("Search error: " + stderr)
                   raise KBaseGenomeUtilException("Search error: " + stderr)
               # TODO: Convert the following Perl script to python library code
               tool_process = subprocess.Popen("xml2kbaseblastjson result.txt > blastoutput_new.json", stderr=subprocess.PIPE, shell=True)
               stdout, stderr = tool_process.communicate()
   
               if stdout is not None and len(stdout) > 0:
                   self.__LOGGER.info(stdout)
   
               if stderr is not None and len(stderr) > 0:
                   self.__LOGGER.error("Output conversion error: " + stderr)
                   raise KBaseGenomeUtilException("Output conversion error: " + stderr)

               with open('blastoutput_new.json', 'r') as myfile:
                 	res1 = json.load(myfile)
           except Exception,e:
               self.__LOGGER.error("Search execution error: Please contact help@kbase.us")
               raise KBaseGenomeUtilException("Search execution error: Please contact help@kbase.us")
        
           #os.remove(query_fn)
         
           #extract the blast output
           # res=script_util.extract_blast_output(self.__BLAST_OUT, anno=g2f)
           
           #os.remove(self.__BLAST_OUT)
           #num_of_hits=len(res)
           
           #metadata=[{'input_genomes':params['genome_ids'][0],'input_sequence':sequence,'number_of_hits':float(num_of_hits)}]
 
           
           #res1={'hits' : res, 'info':metadata}
 
           self.__LOGGER.info( "Finished!!!")
           self.__LOGGER.debug( res1 )
          
 
           #store the BLAST output back into workspace
           res= ws_client.save_objects(
               {"workspace":params['ws_id'],
               "objects": [{
                   "type":"GenomeUtil.BlastOutput",
                   "data":res1,
                   "name":params['output_name']}
               ]})
           
           #print res1
        except KBaseGenomeUtilException, e:
           global no_rst
           res1 = json.loads(no_rst)
           res1['err_msg'] = str(e)
           res= ws_client.save_objects(
               {"workspace":params['ws_id'],
               "objects": [{
                   "type":"GenomeUtil.BlastOutput",
                   "data":res1,
                   "meta":{"err_msg": str(e)},
                   "name":params['output_name']}
               ]})
        except Exception, e:
           res1 = json.loads(no_rst)
           res1['err_msg'] = 'Contact help@kbase.us with the following messages: ' + str(e)
           res= ws_client.save_objects(
               {"workspace":params['ws_id'],
               "objects": [{
                   "type":"GenomeUtil.BlastOutput",
                   "data":res1,
                   "meta":{"err_msg": str(e)},
                   "name":params['output_name']}
               ]})
        finally:
           if not isinstance(res1, dict):
               res1 = json.loads(no_rst)
               res1['err_msg'] = 'Unable to store even the error message to workspace'
           returnVal = res1
        #END blast_against_genome

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method blast_against_genome return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def filter_BlastOutput(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN filter_BlastOutput
        user_token=ctx['token']
        ws_client=Workspace(url=self.__WS_URL, token=user_token)
        blast_outputs=ws_client.get_objects([{'name':params['in_id'], 
                                              'workspace': params['ws_id']}])

            

        fs ={'elements': {}}
        fs['description'] = "FeatureSet from BlastOutput by "
        printedEvalue = False
        printedEntries = False
        if 'evalue' in params and params['evalue'] != "":
            fs['description'] += " E-value:{0}".format(params['evalue'])
            printedEvalue = True
        if 'entries' in params and (params['entries'] != "" or params['entries'] > 0):
            if(printedEvalue): fs['description'] += ","
            fs['description'] += " # of entries :{0}".format(params['entries'])
            printedEntries = True
        if not printedEvalue and not printedEntries:
            fs['description'] += "no filtering"
        
        if len(blast_outputs) != 1:
            fs['description'] = "No such blast output object was found : {0}/{1}".format(param['workspace_name'], param['object_name'])
        else:
            fm = {}
            f2g = {}
            for boid in blast_outputs[0]['data']['BlastOutput_iterations']['Iteration']:
                for hitd in boid['Iteration_hits']['Hit']:
                    print hitd['Hit_def']
                    ali = hitd['Hit_def'].find('#')
                    if(ali < 0): next
                    fid = hitd['Hit_def'][0:ali]
                    gri = hitd['Hit_def'].find('#', ali+1)
                    if fid not in f2g: f2g[fid] = {}
                    if (gri >=  0 and not gri == (ali+1)): 
                        grid = hitd['Hit_def'][(ali+1):gri]
                        f2g[fid][grid] = 1
                    for hspd in hitd['Hit_hsps']['Hsp']:
                        if fid in fm:
                            if float(hspd['Hsp_evalue']) < fm[fid]:
                                fm[fid] = float(hspd['Hsp_evalue'])
                        else: fm[fid] = float(hspd['Hsp_evalue'])
           
            fms = sorted(fm.items(), key=lambda x: x[1], reverse=False)
            bol = len(fms)
            if params['entries'] != "" or int(params['entries']) > 0:
                if(int(params['entries']) < bol):
                    bol = int(params['entries'])
            for i in range(bol):
                if(fms[i][1] > float(params['evalue'])): break
                if fms[i][0] in f2g:
                    fs['elements'][fms[i][0]] = f2g[fms[i][0]].keys()
                else:
                    fs['elements'][fms[i][0]] = []

        ws_client.save_objects(
            {"workspace":params['ws_id'],
            "objects": [{
                "type":"KBaseCollections.FeatureSet",
                "data":fs,
                "name":params['out_id']}
            ]})

        #pprint(fs)
        returnVal = {'obj_name' : params['out_id'], 'ws_id' : params['ws_id']}

        #END filter_BlastOutput

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method filter_BlastOutput return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
