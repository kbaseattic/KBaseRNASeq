#BEGIN_HEADER

import simplejson
import sys
import os
import glob
import json
import logging
import time
from pprint import pprint
import script_util
from biokbase.workspace.client import Workspace
from biokbase.auth import Token

_KBaseGenomeUtil__DATA_VERSION = "0.5"
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
    __TEMP_DIR = 'temp_dir'
    __WS_URL = 'https://ci.kbase.us/services/ws'
    __SHOCK_URL = 'https://ci.kbase.us/services/shock-api/'
    __BLAST_DIR = 'blast_dir'
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

    def blast_against_genome(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN blast_against_genome

        self.__LOGGER.info( "Preparing FA")
        if len(params['query']) > 5:
            sequence=params['query']
        else:
            self.__LOGGER.error("The input sequence is too short!")
            raise Exception("The input sequence is too short!")

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

        if 'target_seqs' in params:
            # let's build index directly and throw away
            sequence = params['target_seqs']
            blast_dir = "%s/%s" %(self.__TEMP_DIR, self.__BLAST_DIR)
            if os.path.exists(blast_dir):
                files=glob.glob("%s/*" % blast_dir)
                for f in files: os.remove(f)
            if not os.path.exists(blast_dir): os.makedirs(blast_dir)
            target_fn = "%s/%s" %( blast_dir, self.__GENOME_FA)
            anno_fn = "%s/%s" %( blast_dir, self.__ANNO_JSON)

            target=open(target_fn,'w')
            if sequence.startswith(">"):
              target.write(sequence)
            else:
              seqes = sequence.split("\n")
              for i in range(len(seqes)):
                target.write(">target_seq_%d\n" %(i))
                target.write(seqes[i])
            target.close()
         
            with open(anno_fn,'w') as anno:
                json.dump({}, anno)

            if(self.__INDEX_TYPE[params['blast_program']]  == 'protein_db'):
                formatdb_type='T'
            elif(self.__INDEX_TYPE[params['blast_program']]  == 'transcript_db'):
                formatdb_type='F'
            else:
                self.__LOGGER.error("{0} is not yet supported".format(params['blast_program']))
                raise Exception("{0} is not yet supported".format(params['blast_program']))
            cmdstring="%s -i %s -p %s -o T" %(self.__INDEX_CMD, target_fn, formatdb_type)
            # TODO: replace it to subprocess.Popen
            os.system(cmdstring)

        else:
            obj_infos = ws_client.get_object_info_new({"objects": [{'name':params['genome_ids'][0],'workspace': params['ws_id']}]})
           
            if len(obj_infos) < 1:
                self.__LOGGER.error("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['genome_ids'][0]))
                raise Exception("Couldn't find %s:%s from the workspace" %(params['ws_id'],params['genome_ids'][0]))
            ref_id = "kb|ws.{0}.obj.{1}.ver.{2}".format(obj_infos[0][6],obj_infos[0][0],obj_infos[0][4])
            obj_chks = obj_infos[0][8]
           
            self.__LOGGER.info( "Querying genome object index from shock for workspace checksum {0}".format(obj_chks))
            query_rst = script_util.query_shock_node(self.__LOGGER, 
                shock_service_url = self.__SHOCK_URL,
                token = svc_token,
                condition = {'ws_obj_chks' : obj_chks, 
                             'service' : 'genome_util', 
                             'svc_data_version' : __DATA_VERSION, 
                             'index_type' : self.__INDEX_TYPE[params['blast_program']]})
           
            blast_dir = "%s/%s" %(self.__TEMP_DIR, self.__BLAST_DIR)
            if os.path.exists(blast_dir):
                files=glob.glob("%s/*" % blast_dir)
                for f in files: os.remove(f)
            if not os.path.exists(blast_dir): os.makedirs(blast_dir)
            zip_fn = "{0}/{1}".format(self.__TEMP_DIR,self.__INDEX_ZIP)
           
           
            target_fn = "%s/%s" %( blast_dir, self.__GENOME_FA)
            anno_fn = "%s/%s" %( blast_dir, self.__ANNO_JSON)
            g2f = {}
            if( query_rst is None or len(query_rst) == 0): # no index available
                self.__LOGGER.info( "Downloading genome object from workspace {0}".format(ref_id))
           
                # TODO: make the following procedures to be loop for each genome_ids 
                genome_list=ws_client.get_object_subset([{'name':params['genome_ids'][0], # replace `0' with loop
                                                          'workspace': params['ws_id'], 
                                                          'included':['features']}])
                genome = genome_list[0]
           
           
                self.__LOGGER.info( "Indexing genome\n")
                # Dump genome sequences
                check_seq=0
                if(self.__INDEX_TYPE[params['blast_program']]  == 'protein_db'):
                    formatdb_type='T'
                    #extract protein sequences from the genome object
                    target=open(target_fn,'w')
                    for gene in genome['data']['features']:
                          #>kb.g.1234.CDS.1234#At1g3333 amalase...
                          function = "NA"
                          aliases = "NA"
                          if 'function' in gene: 
                              g2f[gene['id']] = gene['function']
                              function = gene['function']
                          if 'aliases' in gene: aliases = ",".join(gene['aliases'])
                          if 'protein_translation' in gene.keys():
                                target.write(">%s#%s#%s\n%s\n" % (gene['id'], aliases, function, gene['protein_translation']))
                                check_seq=1
                    target.close()
              
              
                elif(self.__INDEX_TYPE[params['blast_program']]  == 'transcript_db'):
                    formatdb_type='F'
                    #extract dna sequence from the genome object
                    target=open(target_fn,'w')
                    for gene in genome['data']['features']:
                          function = ""
                          aliases = ""
                          if 'function' in gene: 
                              g2f[gene['id']] = gene['function']
                              function = gene['function']
                          if 'aliases' in gene: aliases = ",".join(gene['aliases'])
                          if 'dna_sequence' in gene.keys():
                                target.write(">%s#%s#%s\n%s\n" % (gene['id'], aliases, function, gene['dna_sequence']))
                                check_seq=1
                          if 'function' in gene: g2f[gene['id']] = gene['function']
                    target.close()
                else:
                    self.__LOGGER.error("{0} is not yet supported".format(params['blast_program']))
                    raise Exception("{0} is not yet supported".format(params['blast_program']))
           
                if check_seq == 0:
                    self.__LOGGER.error("The genome object does not contain any sequences!")
                    raise Exception("The genome object does not contain any sequences!")
           
           
                
           
           
                # dump function description
                with open(anno_fn,'w') as anno:
                    json.dump(g2f, anno)
                
              
                cmdstring="%s -i %s -p %s" %(self.__INDEX_CMD, target_fn, formatdb_type)
                # TODO: replace it to subprocess.Popen
                os.system(cmdstring)
           
                # compress
                script_util.zip_files(self.__LOGGER, blast_dir, zip_fn)
           
                # upload
                # TODO: Expand attribute information for future service expansion
                script_util.upload_file_to_shock(self.__LOGGER, 
                                     shock_service_url = self.__SHOCK_URL, 
                                     filePath = zip_fn,
                                     token = svc_token, attributes = {
                                       "ws_ref_id" : ref_id,
                                       "ws_obj_chks" : obj_chks,
                                       'service' : 'genome_util', 
                                       'svc_data_version' : __DATA_VERSION, 
                                       'index_type' : self.__INDEX_TYPE[params['blast_program']]
                                     })
           
            elif (len(query_rst) == 1): # index is available
                # download index
                self.__LOGGER.info("Downloading the genome index")
                script_util.download_file_from_shock(self.__LOGGER,
                                         shock_service_url = self.__SHOCK_URL, 
                                         shock_id = query_rst[0]['id'],
                                         filename = self.__INDEX_ZIP,
                                         directory = self.__TEMP_DIR,
                                         token = svc_token)
                script_util.unzip_files(self.__LOGGER, zip_fn, blast_dir)
                with open(anno_fn,'r') as anno:
                    g2f = json.load(anno)
           
            else:
                self.__LOGGER.error("There are multiple index for the same checksum!")
                raise

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
        os.system(cmdstring)
        # TODO: Convert the following Perl script to python library code
	os.system("xml2kbaseblastjson result.txt > blastoutput_new.json")
	with open('blastoutput_new.json', 'r') as myfile:
		res1 = json.load(myfile)

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
            for boid in blast_outputs[0]['data']['BlastOutput_iterations']['Iteration']:
                for hitd in boid['Iteration_hits']['Hit']:
                    ali = hitd['Hit_def'].find('#')
                    if(ali < 1): next
                    fid = hitd['Hit_def'][0:ali]
                    for hspd in hitd['Hit_hsps']['Hsp']:
                        if fid in fm:
                            if float(hspd['Hsp_evalue']) < fm[fid]:
                                fm[fid] = float(hspd['Hsp_evalue'])
                        else: fm[fid] = float(hspd['Hsp_evalue'])
           
            fms = sorted(fm.items(), key=lambda x: x[1], reverse=False)
            bol = len(fms)
            if 'entries' in params and (params['entries'] != "" or int(params['entries']) > 0):
                if(int(params['entries']) < bol):
                    bol = int(params['entries'])
            if 'evalue' not in params or params['evalue'] == "":
                evalue = 10
            else:
                evalue = float(params['evalue'])
            for i in range(bol):
                if(fms[i][1] > evalue): break
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
