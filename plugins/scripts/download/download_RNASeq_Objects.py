# standard library imports
import os
import sys
import argparse
import logging
import string
import traceback
# 3rd party imports
# None

# KBase imports
import biokbase.workspace.client 
import biokbase.RNASeq.script_util as script_util

# conversion method that can be called if this module is imported
def transform(workspace_service_url=None, shock_service_url=None, handle_service_url=None, 
              workspace_name=None, object_name=None,version=None, working_directory=None, output_file_name=None, 
              logger=None):  
    """
    Converts KBaseRNASeq.RNASeqDifferentialExpression object to a Zipped file
    
    Args:
        workspace_service_url:  A url for the KBase Workspace service 
        shock_service_url: A url for the KBase SHOCK service.
        handle_service_url: A url for the KBase Handle Service.
        workspace_name: Name of the workspace
        object_name: Name of the object in the workspace 
        version : version of the object
        working_directory: The working directory where the output file should be stored.
        output_file_name: The name of the output file should be stored.
    
    Returns:
        Output folder or files from KBaseRNASeq.RNASeqDifferentialExpression.
    
    Authors:
        Srividya Ramakrishnan
    
    """ 
    if logger is None:
        logger = script_util.stderrlogger(__file__)
    
    logger.info("Starting conversion of KBaseRNASeq.RNASeqDifferentialExpression to output files")
    token = os.environ.get("KB_AUTH_TOKEN")
    
    if not os.path.isdir(args.working_directory): 
        raise Exception("The working directory does not exist {0} does not exist".format(working_directory)) 

    logger.info("Grabbing Data.")
 
    try:
        ws_client = biokbase.workspace.client.Workspace(workspace_service_url) 
        if object_name:
            cuffdiff_obj = ws_client.get_objects([{"workspace":workspace_name,"name":object_name}])[0]
    except Exception, e: 
        logger.exception("Unable to retrieve workspace object from {0}:{1}.".format(workspace_service_url,workspace_name))
        logger.exception(e)
        raise 

    shock_id = None
    file_name = None 
    if "file" in cuffdiff_obj["data"]: 
        shock_id = cuffdiff_obj["data"]["file"]["id"]
	if output_file_name:
		file_name = output_file_name
        file_name = cuffdiff_obj["data"]["file"]["file_name"] 
        logger.info("Trying to Retrieve data from Shock.")
        try:
            script_util.download_file_from_shock(logger = logger, 
                                                  shock_service_url = shock_service_url, 
                                                  shock_id = shock_id, 
						  filename = file_name,
                                                  directory = working_directory, 
                                                  token = token)
        except Exception, e:
	    raise Exception("Unable to download data from shock {0}".format("".join(traceback.format_exc())))
            logger.warning("Unable to retrive the file from shock.")

    logger.info("Conversion completed.")

if __name__ == "__main__":
    #script_details = script_utils.parse_docs(transform.__doc__)
    parser = argparse.ArgumentParser(description='Download script for KBaseRNASeq.RNASeqDifferentialExpression')
    #parser = argparse.ArgumentParser(prog=__file__, 
    #                                 description=script_details["Description"],
    #                                 epilog=script_details["Authors"])
    
    # The following 8 arguments should be fairly standard to all uploaders
    parser.add_argument("--workspace_service_url", 
                        #help=script_details["Args"]["workspace_service_url"], 
                        action="store", 
                        type=str, 
                        nargs="?", 
                        required=True)
    data_services = parser.add_mutually_exclusive_group(required=True)
    data_services.add_argument('--shock_service_url', 
				help='Shock url', 
				action='store', 
				type=str, 
				default='https://kbase.us/services/shock-api/', 
				nargs='?')
    data_services.add_argument('--handle_service_url', 
				help='Handle service url', 
				action='store', 
				type=str, 
				default='https://kbase.us/services/handle_service/', 
				nargs='?')
    parser.add_argument("--workspace_name", 
                        #help=script_details["Args"]["workspace_name"], 
                        action="store", 
                        type=str, 
                        nargs="?", 
                        required=True)

    parser.add_argument("--object_name", 
                        #help=script_details["Args"]["object_name"], 
                        action="store", 
                        type=str, 
                        nargs="?",
                        required=True) 

    parser.add_argument("--working_directory", 
                        #help=script_details["Args"]["working_directory"], 
                        action="store", 
                        type=str, 
                        nargs="?", 
                        required=True)

    parser.add_argument("--output_file_name", 
                        #help=script_details["Args"]["output_file_name"], 
                        action="store", 
                        type=str, 
                        nargs="?", 
                        required=False)

    parser.add_argument("--version", 
                        #help=script_details["Args"]["version"], 
                        action="store", 
                        type=int, 
                        nargs="?", 
                        required=False)
    
    args = parser.parse_args()

    logger = script_util.stderrlogger(__file__)
    logger.info("Starting download of KBaseRNASeq.RNASeqDifferentialExpression => zipped output")
    try:
        transform(workspace_service_url = args.workspace_service_url, 
                  shock_service_url=args.shock_service_url,
                  handle_service_url=args.handle_service_url,
		  workspace_name = args.workspace_name, 
                  object_name = args.object_name, 
                  version = args.version, 
                  output_file_name = args.output_file_name,
                  working_directory = args.working_directory, 
                  logger = logger)
    except Exception, e:
        logger.exception(e)
        sys.exit(1)
    
    sys.exit(0)
