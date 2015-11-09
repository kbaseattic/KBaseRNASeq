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
        shutil.rmtree(directory)
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
