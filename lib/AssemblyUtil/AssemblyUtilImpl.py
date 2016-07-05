#BEGIN_HEADER

import os
import sys
import shutil
import traceback
import uuid
from pprint import pprint, pformat

from biokbase.workspace.client import Workspace

from doekbase.data_api.sequence.assembly import api

import trns_transform_FASTA_DNA_Assembly_to_KBaseGenomeAnnotations_Assembly as uploader
from DataFileUtil.DataFileUtilClient import DataFileUtil

#END_HEADER


class AssemblyUtil:
    '''
    Module Name:
    AssemblyUtil

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""
    
    #BEGIN_CLASS_HEADER

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.sharedFolder = config['scratch']
        #END_CONSTRUCTOR
        pass
    

    def get_assembly_as_fasta(self, ctx, params):
        """
        Given a reference to an Assembly (or legacy ContigSet data object), along with a set of options,
        construct a local Fasta file with the sequence data.  If filename is set, attempt to save to the
        specified filename.  Otherwise, a random name will be generated.
        :param params: instance of type "GetAssemblyParams" (@optional
           filename) -> structure: parameter "ref" of String, parameter
           "filename" of String
        :returns: instance of type "FastaAssemblyFile" -> structure:
           parameter "path" of String
        """
        # ctx is the context object
        # return variables are: file
        #BEGIN get_assembly_as_fasta

        if 'ref' not in params:
            raise ValueError('Error calling get_assembly_as_fasta: "ref" parameter not defined.')

        obj = api.AssemblyAPI(
                token=ctx['token'],
                services={ 
                        'workspace_service_url': self.workspaceURL,
                        'shock_service_url': self.shockURL,
                        'handle_service_url': self.handleURL
                    },
                ref=params['ref'] )

        filename = str(uuid.uuid4()) + '.fasta'
        if 'filename' in params and params['filename'] is not None and len(params['filename'].strip())>0:
            filename = params['filename']

        filepath = os.path.join(self.sharedFolder, filename)
        f = open(filepath, 'w')
        obj.get_fasta().to_file(f)

        file = { 'path': filepath }

        #END get_assembly_as_fasta

        # At some point might do deeper type checking...
        if not isinstance(file, dict):
            raise ValueError('Method get_assembly_as_fasta return value ' +
                             'file is not type dict as required.')
        # return the results
        return [file]

    def save_assembly_from_fasta(self, ctx, params):
        """
        WARNING: has the side effect of moving the file to a temporary staging directory, because the upload
        script for assemblies currently requires a working directory, not a specific file.  It will attempt
        to upload everything in that directory.  This will move the file back to the original location, but
        if you are trying to keep an open file handle or are trying to do things concurrently to that file,
        this will break.  So this method is certainly NOT thread safe on the input file.
        :param params: instance of type "SaveAssemblyParams" (Options
           supported: workspace_name assembly_name Uploader options not yet
           supported taxon_reference: The ws reference the assembly points
           to.  (Optional) source: The source of the data (Ex: Refseq)
           date_string: Date (or date range) associated with data. (Optional)
           contig_information_dict: A mapping that has is_circular and
           description information (Optional)) -> structure: parameter "file"
           of type "FastaAssemblyFile" -> structure: parameter "path" of
           String, parameter "workspace_name" of String, parameter
           "assembly_name" of String
        :returns: instance of String
        """
        # ctx is the context object
        # return variables are: ref
        #BEGIN save_assembly_from_fasta

        if 'workspace_name' not in params:
            raise ValueError('workspace_name field was not defined')

        if 'assembly_name' not in params:
            raise ValueError('assembly_name field was not defined')

        # transform scripts assume the file gets dumped in its own unique directory, so create that special
        # directory and move the file there temporarily
        input_directory =  os.path.join(self.sharedFolder, 'assembly-upload-staging-'+str(uuid.uuid4()))
        os.makedirs(input_directory)

        orig_path = None
        temp_path = None
        if 'file' not in params:
            if 'shock_id' not in params:
                raise ValueError('Neither file field nor shock_id field was defined')
            data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=ctx['token'],
                                   service_ver='dev')
            file_name = data_file_cli.shock_to_file({'file_path': input_directory,
                'shock_id': params['shock_id']})['node_file_name']
            temp_path = os.path.join(input_directory, file_name)
            print("temp_path=" + temp_path)
        else:
            if 'path' not in params['file']:
                raise ValueError('file.path field was not defined')
            orig_path = params['file']['path']
            temp_path = os.path.join(input_directory, os.path.basename(orig_path))
            os.rename(orig_path, temp_path)

        
        # do the upload
        result = uploader.upload_assembly(
                logger=None,

                shock_service_url = self.shockURL,
                handle_service_url = self.handleURL,
                workspace_service_url = self.workspaceURL,

                input_directory=input_directory,
                workspace_name=params['workspace_name'],
                assembly_name=params['assembly_name']
            )

#                    taxon_reference = None, 
#                    source = None, 
#                    date_string = None,
#                    contig_information_dict = None,
#                    logger = None):


        if orig_path:
            # move the file back and trash that temporary directory
            os.rename(temp_path, orig_path)
        else:
            os.remove(temp_path)
        shutil.rmtree(input_directory)


        # get WS metadata to return the reference to the object
        ws = Workspace(url=self.workspaceURL)
        info = ws.get_object_info_new({'objects':[{'ref':params['workspace_name'] + '/' + result}],'includeMetadata':0, 'ignoreErrors':0})[0]

        ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])

        #END save_assembly_from_fasta

        # At some point might do deeper type checking...
        if not isinstance(ref, basestring):
            raise ValueError('Method save_assembly_from_fasta return value ' +
                             'ref is not type basestring as required.')
        # return the results
        return [ref]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
