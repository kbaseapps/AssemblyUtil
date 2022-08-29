# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
from pathlib import Path

from AssemblyUtil.FastaToAssembly import FastaToAssembly
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta
from AssemblyUtil.TypeToFasta import TypeToFasta
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

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
    VERSION = "2.0.0"
    GIT_URL = "https://github.com/kbaseapps/AssemblyUtil"
    GIT_COMMIT_HASH = "d3d3e8ac5f7d34cad858c8679ebd096680b95d58"

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
        :param params: instance of type "SaveAssemblyParams" (Required
           arguments: One, and only one, of: file - a pre-existing FASTA file
           to import. The 'assembly_name' field in the FastaAssemblyFile
           object is ignored. shock_id - an ID of a node in the Blobstore
           containing the FASTA file. workspace_name - target workspace
           assembly_name - target object name Optional arguments: type -
           should be one of 'isolate', 'metagenome', (maybe 'transcriptome')
           - defaults to 'Unknown min_contig_length - if set and value is
           greater than 1, this will only include sequences with length
           greater or equal to the min_contig_length specified, discarding
           all other sequences taxon_ref - sets the taxon_ref if present
           contig_info - map from contig_id to a small structure that can be
           used to set the is_circular and description fields for Assemblies
           (optional)) -> structure: parameter "file" of type
           "FastaAssemblyFile" -> structure: parameter "path" of String,
           parameter "assembly_name" of String, parameter "shock_id" of type
           "ShockNodeId", parameter "workspace_name" of String, parameter
           "assembly_name" of String, parameter "external_source" of String,
           parameter "external_source_id" of String, parameter "taxon_ref" of
           String, parameter "min_contig_length" of Long, parameter
           "contig_info" of mapping from String to type "ExtraContigInfo"
           (Structure for setting additional Contig information per contig
           is_circ - flag if contig is circular, 0 is false, 1 is true,
           missing indicates unknown description - if set, sets the
           description of the field in the assembly object which may override
           what was in the fasta file) -> structure: parameter "is_circ" of
           Long, parameter "description" of String
        :returns: instance of String
        """
        # ctx is the context object
        # return variables are: ref
        #BEGIN save_assembly_from_fasta

        print('save_assembly_from_fasta -- paramaters = ')
        #pprint(params)

        ref = FastaToAssembly(
            DataFileUtil(self.callback_url, token=ctx['token']),
            Workspace(self.ws_url, token=ctx['token']),
            Path(self.sharedFolder)
        ).import_fasta(params)

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
