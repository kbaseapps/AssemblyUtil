# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os

from AssemblyUtil.FastaToAssembly import FastaToAssembly
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta
from AssemblyUtil.TypeToFasta import TypeToFasta
from installed_clients.WorkspaceClient import Workspace
from installed_clients.baseclient import ServerError

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
    VERSION = "1.2.3"
    GIT_URL = "git@github.com:kbaseapps/AssemblyUtil.git"
    GIT_COMMIT_HASH = "56ea7818e13cc38dadd74e52d4377070a488d74c"

    #BEGIN_CLASS_HEADER

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.sharedFolder = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.ws_url = config['workspace-url']
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
        file = atf.assembly_as_fasta(params)

        #END get_assembly_as_fasta

        # At some point might do deeper type checking...
        if not isinstance(file, dict):
            raise ValueError('Method get_assembly_as_fasta return value ' +
                             'file is not type dict as required.')
        # return the results
        return [file]

    def get_fastas(self, ctx, params):
        """
        Given a reference list of KBase objects constructs a local Fasta file with the sequence data for each ref.
        :param params: instance of type "KBaseOjbReferences" -> structure:
           parameter "ref_lst" of list of type "ref" (ref: workspace
           reference. KBaseOjbReferences: ref_lst: is an object wrapped array
           of KBase object references, which can be of the following types: -
           KBaseGenomes.Genome - KBaseSets.AssemblySet -
           KBaseMetagenome.BinnedContigs - KBaseGenomes.ContigSet -
           KBaseGenomeAnnotations.Assembly - KBaseSearch.GenomeSet -
           KBaseSets.GenomeSet ref_fastas paths - list of paths to fasta
           files associated with workspace object. type - workspace object
           type parent_refs - (optional) list of associated workspace object
           references if different from the output key)
        :returns: instance of mapping from type "ref" (ref: workspace
           reference. KBaseOjbReferences: ref_lst: is an object wrapped array
           of KBase object references, which can be of the following types: -
           KBaseGenomes.Genome - KBaseSets.AssemblySet -
           KBaseMetagenome.BinnedContigs - KBaseGenomes.ContigSet -
           KBaseGenomeAnnotations.Assembly - KBaseSearch.GenomeSet -
           KBaseSets.GenomeSet ref_fastas paths - list of paths to fasta
           files associated with workspace object. type - workspace object
           type parent_refs - (optional) list of associated workspace object
           references if different from the output key) to type "ref_fastas"
           -> structure: parameter "paths" of list of String, parameter
           "parent_refs" of list of type "ref" (ref: workspace reference.
           KBaseOjbReferences: ref_lst: is an object wrapped array of KBase
           object references, which can be of the following types: -
           KBaseGenomes.Genome - KBaseSets.AssemblySet -
           KBaseMetagenome.BinnedContigs - KBaseGenomes.ContigSet -
           KBaseGenomeAnnotations.Assembly - KBaseSearch.GenomeSet -
           KBaseSets.GenomeSet ref_fastas paths - list of paths to fasta
           files associated with workspace object. type - workspace object
           type parent_refs - (optional) list of associated workspace object
           references if different from the output key), parameter "type" of
           String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN get_fastas

        # Param check
        if not params:
            raise ValueError("Must provide params.")
        if not params.get("ref_lst"):
            raise ValueError("Must provice list of references 'ref_lst'")

        ref_lst = params.get("ref_lst")

        ws = Workspace(url=self.ws_url, token=ctx["token"])

        ttf = TypeToFasta(self.callback_url, self.sharedFolder, ws, ctx["token"])
        output = ttf.type_to_fasta(ref_lst)

        #END get_fastas

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method get_fastas return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

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
        output = atf.export_as_fasta(params)

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
           unknown description - if set, sets the description of the field in
           the assembly object which may override what was in the fasta file)
           -> structure: parameter "is_circ" of Long, parameter "description"
           of String
        :returns: instance of String
        """
        # ctx is the context object
        # return variables are: ref
        #BEGIN save_assembly_from_fasta

        print('save_assembly_from_fasta -- paramaters = ')
        #pprint(params)

        fta = FastaToAssembly(self.callback_url, self.sharedFolder, self.ws_url)
        try:
            assembly_info = fta.import_fasta(ctx, params)
        except ServerError as err:
            if 'Valid Content-Length header >= 0' in str(err):
                raise Exception('No output generated by the app.' 
                                  'The combination of data and input parameters '
                                  'resulted in no valid output.')
            else:
                raise err
        ref = f'{assembly_info[6]}/{assembly_info[0]}/{assembly_info[4]}'

        #END save_assembly_from_fasta

        # At some point might do deeper type checking...
        if not isinstance(ref, str):
            raise ValueError('Method save_assembly_from_fasta return value ' +
                             'ref is not type str as required.')
        # return the results
        return [ref]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
