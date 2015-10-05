import sys
import os
import json
import re
import io
import urllib
import requests
import logging
import time
from requests_toolbelt import MultipartEncoder
import subprocess
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join

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


def download_file_from_shock(logger,
                             shock_service_url = None,
                             shock_id = None,
                             filename = None,
                             directory = None,
                             token = None):
    """
    Given a SHOCK instance URL and a SHOCK node id, download the contents of that node
    to a file on disk.
    """

    header = dict()
    header["Authorization"] = "Oauth {0}".format(token)

    logger.info("Downloading shock node {0}/node/{1}".format(shock_service_url,shock_id))

    metadata_response = requests.get("{0}/node/{1}?verbosity=metadata".format(shock_service_url, shock_id), headers=header, stream=True, verify=True)
    shock_metadata = metadata_response.json()['data']
    shockFileName = shock_metadata['file']['name']
    shockFileSize = shock_metadata['file']['size']
    metadata_response.close()
        
    download_url = "{0}/node/{1}?download_raw".format(shock_service_url, shock_id)
        
    data = requests.get(download_url, headers=header, stream=True, verify=True)

    if filename is not None:
        shockFileName = filename

    if directory is not None:
        filePath = os.path.join(directory, shockFileName)
    else:
        filePath = shockFileName

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










