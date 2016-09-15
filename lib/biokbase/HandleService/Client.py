from biokbase.AbstractHandle.Client import AbstractHandle;
try:
    import json as _json
except ImportError:
    import sys
    sys.path.append('simplejson-2.3.3')
    import simplejson as _json

import requests as _requests
import urlparse as _urlparse
import random as _random
import base64 as _base64
from ConfigParser import ConfigParser as _ConfigParser
import os as _os
from requests_toolbelt import MultipartEncoder as _MultipartEncoder

##
# User defined import
import io



# I think these should be in the kb auth library function
def _get_token(user_id, password,
               auth_svc='https://nexus.api.globusonline.org/goauth/token?' +
                        'grant_type=client_credentials'):
    # This is bandaid helper function until we get a full
    # KBase python auth client released
    auth = _base64.encodestring(user_id + ':' + password)
    headers = {'Authorization': 'Basic ' + auth}
    ret = _requests.get(auth_svc, headers=headers, allow_redirects=True)
    status = ret.status_code
    if status >= 200 and status <= 299:
        tok = _json.loads(ret.text)
    elif status == 403:
        raise Exception('Authentication failed: Bad user_id/password ' +
                        'combination for user %s' % (user_id))
    else:
        raise Exception(ret.text)
    return tok['access_token']


def _read_rcfile(file=_os.environ['HOME'] + '/.authrc'):  # @ReservedAssignment
    # Another bandaid to read in the ~/.authrc file if one is present
    authdata = None
    if _os.path.exists(file):
        try:
            with open(file) as authrc:
                rawdata = _json.load(authrc)
                # strip down whatever we read to only what is legit
                authdata = {x: rawdata.get(x) for x in (
                    'user_id', 'token', 'client_secret', 'keyfile',
                    'keyfile_passphrase', 'password')}
        except Exception, e:
            print "Error while reading authrc file %s: %s" % (file, e)
    return authdata


def _read_inifile(file=_os.environ.get(  # @ReservedAssignment
                  'KB_DEPLOYMENT_CONFIG', _os.environ['HOME'] +
                  '/.kbase_config')):
    # Another bandaid to read in the ~/.kbase_config file if one is present
    authdata = None
    if _os.path.exists(file):
        try:
            config = _ConfigParser()
            config.read(file)
            # strip down whatever we read to only what is legit
            authdata = {x: config.get('authentication', x)
                        if config.has_option('authentication', x)
                        else None for x in ('user_id', 'token',
                                            'client_secret', 'keyfile',
                                            'keyfile_passphrase', 'password')}
        except Exception, e:
            print "Error while reading INI file %s: %s" % (file, e)
    return authdata

##
# User defined function 
def _upload_file_to_shock(node_url, filePath = None, ssl_verify = True, token = None): 

    if token is None:
        raise Exception("Authentication token required!")
    
    #build the header
    header = dict()
    header["Authorization"] = "OAuth %s" % token

    if filePath is None:
        raise Exception("No file given for upload to SHOCK!")

    dataFile = open(_os.path.abspath(filePath), 'r')
    m = _MultipartEncoder(fields={'upload': (_os.path.split(filePath)[-1], dataFile)})
    header['Content-Type'] = m.content_type

    try:
        response = _requests.put(node_url, headers=header, data=m, allow_redirects=True, verify=ssl_verify)
        dataFile.close()

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:            
            raise Exception(result['error'][0])
        else:
            return result["data"]    
    except:
        dataFile.close()
        raise    

class HandleService(object):
    def __init__(self, url=None, timeout=30 * 60, user_id=None,
                 password=None, token=None, ignore_authrc=False,
                 trust_all_ssl_certificates=False):

        self.url = url
        self.timeout = int(timeout)
        self._headers = dict()
        self.trust_all_ssl_certificates = trust_all_ssl_certificates
        self.token = token
        # token overrides user_id and password
        if token is not None: pass
        elif user_id is not None and password is not None:
            self.token = _get_token(user_id, password)
        elif 'KB_AUTH_TOKEN' in _os.environ:
            self.token = _os.environ.get('KB_AUTH_TOKEN')
        elif not ignore_authrc:
            authdata = _read_inifile()
            if authdata is None:
                authdata = _read_rcfile()
            if authdata is not None:
                if authdata.get('token') is not None:
                    self.token = authdata['token']
                elif(authdata.get('user_id') is not None
                     and authdata.get('password') is not None):
                    self.token = _get_token(
                        authdata['user_id'], authdata['password'])
        if self.timeout < 1:
            raise ValueError('Timeout value must be at least 1 second')

        self.dsi = AbstractHandle(url=url, token=self.token, trust_all_ssl_certificates=trust_all_ssl_certificates)



    def upload(self,infile) :

        handle = self.dsi.new_handle()
	url = "{}/node/{}".format(handle["url"], handle["id"]);
        ref_data = {}
        try:
            ref_data = _upload_file_to_shock(url,filePath=infile,token=self.token)
        except:
            raise

	remote_md5 = ref_data["file"]["checksum"]["md5"]
	#remote_sha1 = ref_data["file"]["checksum"]["sha1"] #SHOCK PUT would not create remote_sha1

        if remote_md5 is None: raise Exception("looks like upload failed and no md5 returned from remote server")
	
	handle["remote_md5"] = remote_md5
	#handle["remote_sha1"] = remote_sha1 	
        handle["file_name"] = _os.path.basename(infile)

	self.dsi.persist_handle(handle)
	return handle


    def download (self, handle, outfile):

        if( isinstance(handle, dict)): raise Exception("hande is not a dictionary")
        if "id" not in handle: raise Exception("no id in handle")
        if "url" not in handle: raise Exception("no url in handle")
        if outfile is None: raise Exception("outfile is not defined")
        raise Exception("Not implemented yet")

    def new_handle(self, *arg):
        return self.dsi.new_handle(*arg)

    def localize_handle(self, *arg):
        return self.dsi.localize_handle(*arg)

    def initialize_handle (self, *arg):
        return self.dsi.initialize_handle (*arg)

    def persist_handle(self, *arg):
        return self.dsi.persist_handle (*arg)

    def upload_metadata(self, handle, infile):
        raise Exception("Not implemented yet")

    def download_metadata (self, handle, outfile):
        raise Exception("Not implemented yet")

    def list_handles (self, *arg):
        return self.dsi.list_handles (*arg)

    def are_readable (self, *arg):
        return self.dsi.are_readable(*arg)

    def is_readable(self, *arg):
        return self.dsi.is_readable(*arg)

    def hids_to_handles(self, *arg):
        return self.dsi.hids_to_handles(*arg)

    def ids_to_handles(self, *arg):
        return self.dsi.ids_to_handles(*arg)
