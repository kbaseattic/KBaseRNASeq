import sys
import os
import json
import re
import io
import urllib
import hashlib
import requests
import logging
import shutil
import time
import traceback
from requests_toolbelt import MultipartEncoder
from multiprocessing import Pool
#from functools import partial
import subprocess
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join

try:
    from biokbase.HandleService.Client import HandleService
except:
    from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

from biokbase.workspace.client import Workspace

def stderrlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stderr handler attached and using a prefix
    format that will make logging consistent between scripts.
    """
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stderr)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    
    return logger


def stdoutlogger(name, level=logging.INFO):
    """
    Return a standard python logger with a stdout handler attached and using a prefix
    format that will make logging consistent between scripts.
    """
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # send messages to sys.stderr
    streamHandler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)

    logger.addHandler(streamHandler)
    
    return logger



def zip_files(logger, src_path, output_fn):
    """
    Compress all index files (not directory) into an output zip file on disk.
    """

    files = [ f for f in listdir(src_path) if isfile(join(src_path,f)) ]
    with ZipFile(output_fn, 'w') as izip:
        for f in files:
            izip.write(join(src_path,f),f)

def unzip_files(logger, src_fn, dst_path):
    """
    Extract all index files into an output zip file on disk.
    """

    with ZipFile(src_fn, 'r') as ozip:
        ozip.extractall(dst_path)

def move_files(logger, src, dest):
    """
    Move files from one folder to another.
    """
    
    src_files = os.listdir(src)
    for file_name in src_files:
       full_file_name = os.path.join(src, file_name)
       if (os.path.isfile(full_file_name)):
          shutil.copy(full_file_name, dest)

def download_file_from_shock(logger,
                             shock_service_url = None,
                             shock_id = None,
                             filename = None,
                             directory = None,
			     filesize= None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)
    print header

    logger.info("Downloading shock node {0}/node/{1}".format(shock_service_url,shock_id))

    metadata_response = requests.get("{0}/node/{1}?verbosity=metadata".format(shock_service_url, shock_id), headers=header, stream=True, verify=True)
    shock_metadata = metadata_response.json()['data']
    print "shock metadata is {0}".format(shock_metadata)
    if shock_metadata is not None:
    	shockFileName = shock_metadata['file']['name']
    	shockFileSize = shock_metadata['file']['size']
    metadata_response.close()
        
    download_url = "{0}/node/{1}?download_raw".format(shock_service_url, shock_id)
    print "download_url is {0}".format(download_url)
    try: 
    	data = requests.get(download_url, headers=header, stream=True, verify=True)
    except Exception,e:
	print(traceback.format_exc())
    if filename is not None:
        shockFileName = filename

    if directory is not None:
        filePath = os.path.join(directory, shockFileName)
    else:
        filePath = shockFileName

    if filesize is not None:
	shockFileSize = filesize

    chunkSize = shockFileSize/4
    
    maxChunkSize = 2**30
    
    if chunkSize > maxChunkSize:
        chunkSize = maxChunkSize
    
    f = io.open(filePath, 'wb')
    try:
        for chunk in data.iter_content(chunkSize):
            f.write(chunk)
    finally:
        data.close()
        f.close()      
    
def query_shock_node(logger,
                             shock_service_url = None,
                             condition = None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    query_str = urllib.urlencode(condition)
    print query_str
    
    logger.info("Querying shock node {0}/node/?query&{1}".format(shock_service_url,query_str))

    query_response = requests.get("{0}/node/?query&{1}".format(shock_service_url, query_str), headers=header, stream=True, verify=True)
    query_rst = query_response.json()['data']
    query_response.close()
    return query_rst


def upload_file_to_shock(logger,
                         shock_service_url = None,
                         filePath = None,
                         attributes = '{}',
                         ssl_verify = True,
                         token = None):
    """
    Use HTTP multi-part POST to save a file to a SHOCK instance.
    """

    if token is None:
        raise Exception("Authentication token required!")
    
    #build the header
    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(os.path.abspath(filePath), 'r')
    m = MultipartEncoder(fields={'attributes_str': json.dumps(attributes), 'upload': (os.path.split(filePath)[-1], dataFile)})
    header['Content-Type'] = m.content_type

    logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
    try:
        response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
        dataFile.close()
    except:
        dataFile.close()
        raise    

    if not response.ok:
        response.raise_for_status()

    result = response.json()

    if result['error']:            
        raise Exception(result['error'][0])
    else:
        return result["data"]    


def getHandles(logger = None,
               shock_service_url = None,
               handle_service_url = None,
               shock_ids = None,
               handle_ids = None,
               token = None):
    """
    Retrieve KBase handles for a list of shock ids or a list of handle ids.
    """

    if token is None:
        raise Exception("Authentication token required!")

    hs = HandleService(url=handle_service_url, token=token)

    handles = list()
    if shock_ids is not None:
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        for sid in shock_ids:
            info = None

            try:
                logger.info("Found shock id {0}, retrieving information about the data.".format(sid))

                response = requests.get("{0}/node/{1}".format(shock_service_url, sid), headers=header, verify=True)
                info = response.json()["data"]
            except:
                logger.error("There was an error retrieving information about the shock node id {0} from url {1}".format(sid, shock_service_url))

            try:
                logger.info("Retrieving a handle id for the data.")
                handle = hs.persist_handle({"id" : sid,
                                           "type" : "shock",
                                           "url" : shock_service_url,
                                           "file_name": info["file"]["name"],
                                           "remote_md5": info["file"]["checksum"]["md5"]})
                handles.append(handle)
            except:
                try:
                    handle_id = hs.ids_to_handles([sid])[0]["hid"]
		    single_handle = hs.hids_to_handles([handle_id])

                    assert len(single_handle) != 0

                    if info is not None:
                        single_handle[0]["file_name"] = info["file"]["name"]
                        single_handle[0]["remote_md5"] = info["file"]["checksum"]["md5"]
                        logger.debug(single_handle)

                    handles.append(single_handle[0])
                except:
                    logger.error("The input shock node id {} is already registered or could not be registered".format(sid))

                    hs = HandleService(url=handle_service_url, token=token)
                    all_handles = hs.list_handles()

                    for x in all_handles:
                        if x[0] == sid:
                            logger.info("FOUND shock id as existing handle")
                            logger.info(x)
                            break
                    else:
                        logger.info("Unable to find a handle containing shock id")

                        logger.info("Trying again to get a handle id for the data.")
                        handle_id = hs.persist_handle({"id" : sid,
                                           "type" : "shock",
                                           "url" : shock_service_url,
                                           "file_name": info["file"]["name"],
                                           "remote_md5": info["file"]["checksum"]["md5"]})
                        handles.append(handle_id)

                    raise
    elif handle_ids is not None:
        for hid in handle_ids:
            try:
                single_handle = hs.hids_to_handles([hid])

                assert len(single_handle) != 0

                handles.append(single_handle[0])
            except:
                logger.error("Invalid handle id {0}".format(hid))
                raise
    return handles

def get_obj_info(logger,ws_url,objects,ws_id,token):
    """
    function to get the workspace object id from a object name
    """
    ret = []
    ws_client=Workspace(url=ws_url, token=token)
    for obj in  objects:
    	try:
            obj_infos = ws_client.get_object_info_new({"objects": [{'name': obj, 'workspace': ws_id}]})
            print obj_infos
            ret.append("{0}/{1}/{2}".format(obj_infos[0][6],obj_infos[0][0],obj_infos[0][4]))
            print ret
        except Exception, e:
                     logger.error("Couldn't retrieve %s:%s from the workspace , %s " %(ws_id,obj,e))
    return ret

def whereis(program):
    """
    returns path of program if it exists in your ``$PATH`` variable or ``None`` otherwise
    """
    for path in os.environ.get('PATH', '').split(':'):
    	if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None

def runProgram(logger=None,
	       progName=None,
	       argStr=None,
	       script_dir=None,
	       working_dir=None):
        """
        Convenience func to handle calling and monitoring output of external programs.
    
        :param progName: name of system program command
        :param argStr: string containing command line options for ``progName``
    
        :returns: subprocess.communicate object
        """

        # Ensure program is callable.
        if script_dir is not None:
                progPath= os.path.join(script_dir,progName)
        else:
		progPath = progName
	#	progPath = whereis(progName)
        #       	if not progPath:
        #            raise RuntimeError(None,'{0} command not found in your PATH environmental variable. {1}'.format(progName,os.environ.get('PATH', '')))

        # Construct shell command
        cmdStr = "%s %s" % (progPath,argStr)

        # Set up process obj
        process = subprocess.Popen(cmdStr,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=working_dir)
        # Get results
        result,stderr  = process.communicate()
      
        # keep this until your code is stable for easier debugging
        if result is not None and len(result) > 0:
            logger.info(result)
        if stderr is not None and len(stderr) > 0:
            logger.info(stderr)

        # Check returncode for success/failure
        if process.returncode != 0:
                raise RuntimeError('Return Code : {0} , result {1} , progName {2}'.format(process.returncode,result[1],progName))

        # Return result
        return result

def hashfile(filepath):
       sha1 = hashlib.sha1()
       f = open(filepath, 'rb')
       try:
         sha1.update(f.read())
       finally:
         f.close()
       return sha1.hexdigest()


def create_shock_handle(logger=None,
			file_name=None,
			shock_url=None,
			handle_url=None,
			obj_type=None,
			token=None):
     	
        hs = HandleService(url=handle_url, token=token)
        f_shock = upload_file_to_shock(logger,shock_url,file_name,'{}',True,token)
        f_sha1 =  hashfile(file_name)
        hid = getHandles(logger,shock_url,handle_url,[f_shock['id']],None,token)[0]
        handle = { 'hid' : hid , 
	           "file_name" : f_shock['file']['name'] , 
                   "id" : f_shock['id'] , 
                   "type" : obj_type ,
	           "url" : shock_url,
	           "remote_md5" : f_shock['file']['checksum']['md5'],
	           "remote_sha1" : f_sha1 }	
   
        return handle

def parallel_function(f):
    def easy_parallize(f, sequence):
        pool = Pool(processes=8)
        # f is given sequence. guaranteed to be in order
        result = pool.map(f, sequence)
        cleaned = [x for x in result if not x is None]
        cleaned = asarray(cleaned)
        # not optimal but safe
        pool.close()
        pool.join()
        return cleaned
    from functools import partial
    # this assumes f has one argument, fairly easy with Python's global scope
    return partial(easy_parallize, f)
     
