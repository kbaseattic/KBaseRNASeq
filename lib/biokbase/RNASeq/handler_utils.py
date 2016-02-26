"""
Provides utility functions for handler/taskrunner level code to provide consistent
behavior across different taskrunners and reduce code repetition of similar tasks.
"""

import sys
import os
import shutil
import subprocess
import base64
import simplejson


def cleanup(logger=None, directory=None):
    """
    Clean up after the job.  At the moment this just means removing the working
    directory, but later could mean other things.
    """
    
    try:
        shutil.rmtree(directory, ignore_errors=True) # it would not delete if fold is not empty
        # need to iterate each entry
    except IOError, e:
        logger.error("Unable to remove working directory {0}".format(directory))
        raise


def gen_recursive_filelist(d):
    """
    Generate a list of all files present below a given directory.
    """
    
    for root, directories, files in os.walk(d):
        for file in files:
            yield os.path.join(root, file)

def get_dir(d):
    """
    Generate a list of all files present below a given directory.
    """
    
    dirs = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
    if len(dirs) > 1:
    	return dirs[0]
    else:
	return None
    #for directories, files in os.walk(d):
    #    return os.path(directories)

def get_file_with_suffix(d,suffix):
    """
    Generate a list of all files present below a given directory.
    """

    items = os.listdir(d)
    for file in items:
            if file.endswith(suffix):
		return file.split(suffix)[0]
    return None
