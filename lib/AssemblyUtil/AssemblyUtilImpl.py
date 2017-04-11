# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
from pprint import pprint

from AssemblyUtil.FastaToAssembly import FastaToAssembly
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta

#END_HEADER


class AssemblyUtil:
    '''
    Module Name:
    AssemblyUtil

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.0"
    GIT_URL = "git@github.com:kbaseapps/AssemblyUtil"
    GIT_COMMIT_HASH = "b85b2ded513457f384868d0c324dcb7d42952d5b"

    #BEGIN_CLASS_HEADER

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.sharedFolder = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
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
           parameter "path" of String, parameter "assembly_name" of String
        """
        # ctx is the context object
        # return variables are: file
        #BEGIN get_assembly_as_fasta

        atf = AssemblyToFasta(self.callback_url, self.sharedFolder)
        file = atf.assembly_as_fasta(ctx, params)

        #END get_assembly_as_fasta

        # At some point might do deeper type checking...
        if not isinstance(file, dict):
            raise ValueError('Method get_assembly_as_fasta return value ' +
                             'file is not type dict as required.')
        # return the results
        return [file]

    def export_assembly_as_fasta(self, ctx, params):
        """
        A method designed especially for download, this calls 'get_assembly_as_fasta' to do
        the work, but then packages the output with WS provenance and object info into
        a zip file and saves to shock.
        :param params: instance of type "ExportParams" -> structure:
           parameter "input_ref" of String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_assembly_as_fasta

        atf = AssemblyToFasta(self.callback_url, self.sharedFolder)
        output = atf.export_as_fasta(ctx, params)

        #END export_assembly_as_fasta

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_assembly_as_fasta return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def save_assembly_from_fasta(self, ctx, params):
        """
        WARNING: has the side effect of moving the file to a temporary staging directory, because the upload
        script for assemblies currently requires a working directory, not a specific file.  It will attempt
        to upload everything in that directory.  This will move the file back to the original location, but
        if you are trying to keep an open file handle or are trying to do things concurrently to that file,
        this will break.  So this method is certainly NOT thread safe on the input file.
        :param params: instance of type "SaveAssemblyParams" (Options
           supported: file / shock_id / ftp_url - mutualy exclusive
           parameters pointing to file content workspace_name - target
           workspace assembly_name - target object name type - should be one
           of 'isolate', 'metagenome', (maybe 'transcriptome')
           min_contig_length - if set and value is greater than 1, this will
           only include sequences with length greater or equal to the
           min_contig_length specified, discarding all other sequences
           taxon_ref         - sets the taxon_ref if present contig_info     
           - map from contig_id to a small structure that can be used to set
           the is_circular and description fields for Assemblies (optional)
           Uploader options not yet supported taxon_reference: The ws
           reference the assembly points to.  (Optional) source: The source
           of the data (Ex: Refseq) date_string: Date (or date range)
           associated with data. (Optional)) -> structure: parameter "file"
           of type "FastaAssemblyFile" -> structure: parameter "path" of
           String, parameter "assembly_name" of String, parameter "shock_id"
           of type "ShockNodeId", parameter "ftp_url" of String, parameter
           "workspace_name" of String, parameter "assembly_name" of String,
           parameter "external_source" of String, parameter
           "external_source_id" of String, parameter "taxon_ref" of String,
           parameter "min_contig_length" of Long, parameter "contig_info" of
           mapping from String to type "ExtraContigInfo" (Structure for
           setting additional Contig information per contig is_circ - flag if
           contig is circular, 0 is false, 1 is true, missing indicates
           unknown) -> structure: parameter "is_circ" of Long, parameter
           "description" of String
        :returns: instance of String
        """
        # ctx is the context object
        # return variables are: ref
        #BEGIN save_assembly_from_fasta

        print('save_assembly_from_fasta -- paramaters = ')
        pprint(params)

        fta = FastaToAssembly(self.callback_url, self.sharedFolder)
        assembly_info = fta.import_fasta(ctx, params)
        ref = str(assembly_info[6]) + '/' + str(assembly_info[0]) + '/' + str(assembly_info[4])

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
