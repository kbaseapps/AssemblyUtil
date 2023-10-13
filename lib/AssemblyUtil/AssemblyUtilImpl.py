# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os
from pathlib import Path

from AssemblyUtil.FastaToAssembly import FastaToAssembly, MAX_THREADS, THREADS_PER_CPU
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta
from AssemblyUtil.TypeToFasta import TypeToFasta
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


def _validate_max_threads_type(threads_count, var_name, default_val):
    if threads_count is None:
        print(f"Cannot retrieve {var_name} from the catalog, set {var_name}={default_val}")
        return default_val
    try:
        threads_count = int(threads_count)
    except ValueError as e:
        raise ValueError(f"{var_name} must be an integer") from e
    return threads_count


def _validate_threads_per_cpu_type(threads_count, var_name, default_val):
    if threads_count is None:
        print(f"Cannot retrieve {var_name} from the catalog, set {var_name}={default_val}")
        return default_val
    try:
        threads_count = float(threads_count)
    except ValueError as e:
        raise ValueError(f"{var_name} must be an integer or decimal") from e
    return threads_count
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
    VERSION = "3.0.0"
    GIT_URL = "https://github.com/kbaseapps/AssemblyUtil"
    GIT_COMMIT_HASH = "72183303c99ec8e88534173adbc47499da2e9b56"

    #BEGIN_CLASS_HEADER

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.sharedFolder = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.ws_url = config['workspace-url']

        max_threads = config.get("KBASE_SECURE_CONFIG_PARAM_MAX_THREADS")
        threads_per_cpu = config.get("KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU")
        self.max_threads = _validate_max_threads_type(max_threads, "MAX_THREADS", MAX_THREADS)
        self.threads_per_cpu = _validate_threads_per_cpu_type(threads_per_cpu, "THREADS_PER_CPU", THREADS_PER_CPU)
        #END_CONSTRUCTOR

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

    def save_assembly_from_fasta2(self, ctx, params):
        """
        Save a KBase Workspace assembly object from a FASTA file.
        :param params: instance of type "SaveAssemblyParams" (Required
           arguments: Exactly one of: file - a pre-existing FASTA file to
           import. The 'assembly_name' field in the FastaAssemblyFile object
           is ignored. shock_id - an ID of a node in the Blobstore containing
           the FASTA file. Exactly one of: workspace_id - the immutable,
           numeric ID of the target workspace. Always prefer providing the ID
           over the name. workspace_name - the name of the target workspace.
           assembly_name - target object name Optional arguments: type -
           should be one of 'isolate', 'metagenome', (maybe 'transcriptome').
           Defaults to 'Unknown' min_contig_length - if set and value is
           greater than 1, this will only include sequences with length
           greater or equal to the min_contig_length specified, discarding
           all other sequences contig_info - map from contig_id to a small
           structure that can be used to set the is_circular and description
           fields for Assemblies (optional)) -> structure: parameter "file"
           of type "FastaAssemblyFile" -> structure: parameter "path" of
           String, parameter "assembly_name" of String, parameter "shock_id"
           of type "ShockNodeId", parameter "workspace_id" of Long, parameter
           "workspace_name" of String, parameter "assembly_name" of String,
           parameter "type" of String, parameter "external_source" of String,
           parameter "external_source_id" of String, parameter
           "min_contig_length" of Long, parameter "contig_info" of mapping
           from String to type "ExtraContigInfo" (Structure for setting
           additional Contig information per contig is_circ - flag if contig
           is circular, 0 is false, 1 is true, missing indicates unknown
           description - if set, sets the description of the field in the
           assembly object which may override what was in the fasta file) ->
           structure: parameter "is_circ" of Long, parameter "description" of
           String
        :returns: instance of type "SaveAssemblyResult" (Results from saving
           an assembly. upa - the address of the resulting workspace object.
           filtered_input - the filtered input file if the minimum contig
           length parameter is present and > 0. null otherwise.) ->
           structure: parameter "upa" of type "upa" (A Unique Permanent
           Address for a workspace object, which is of the form W/O/V, where
           W is the numeric workspace ID, O is the numeric object ID, and V
           is the object version.), parameter "filtered_input" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN save_assembly_from_fasta2
        result = FastaToAssembly(
            DataFileUtil(self.callback_url, token=ctx['token']),
            Path(self.sharedFolder)
        ).import_fasta(params)
        #END save_assembly_from_fasta2

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method save_assembly_from_fasta2 return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def save_assembly_from_fasta(self, ctx, params):
        """
        @deprecated AssemblyUtil.save_assembly_from_fasta2
        :param params: instance of type "SaveAssemblyParams" (Required
           arguments: Exactly one of: file - a pre-existing FASTA file to
           import. The 'assembly_name' field in the FastaAssemblyFile object
           is ignored. shock_id - an ID of a node in the Blobstore containing
           the FASTA file. Exactly one of: workspace_id - the immutable,
           numeric ID of the target workspace. Always prefer providing the ID
           over the name. workspace_name - the name of the target workspace.
           assembly_name - target object name Optional arguments: type -
           should be one of 'isolate', 'metagenome', (maybe 'transcriptome').
           Defaults to 'Unknown' min_contig_length - if set and value is
           greater than 1, this will only include sequences with length
           greater or equal to the min_contig_length specified, discarding
           all other sequences contig_info - map from contig_id to a small
           structure that can be used to set the is_circular and description
           fields for Assemblies (optional)) -> structure: parameter "file"
           of type "FastaAssemblyFile" -> structure: parameter "path" of
           String, parameter "assembly_name" of String, parameter "shock_id"
           of type "ShockNodeId", parameter "workspace_id" of Long, parameter
           "workspace_name" of String, parameter "assembly_name" of String,
           parameter "type" of String, parameter "external_source" of String,
           parameter "external_source_id" of String, parameter
           "min_contig_length" of Long, parameter "contig_info" of mapping
           from String to type "ExtraContigInfo" (Structure for setting
           additional Contig information per contig is_circ - flag if contig
           is circular, 0 is false, 1 is true, missing indicates unknown
           description - if set, sets the description of the field in the
           assembly object which may override what was in the fasta file) ->
           structure: parameter "is_circ" of Long, parameter "description" of
           String
        :returns: instance of String
        """
        # ctx is the context object
        # return variables are: ref
        #BEGIN save_assembly_from_fasta

        print('save_assembly_from_fasta -- paramaters = ')
        #pprint(params)
        ref = self.save_assembly_from_fasta2(ctx, params)[0]['upa']

        #END save_assembly_from_fasta

        # At some point might do deeper type checking...
        if not isinstance(ref, str):
            raise ValueError('Method save_assembly_from_fasta return value ' +
                             'ref is not type str as required.')
        # return the results
        return [ref]

    def save_assemblies_from_fastas(self, ctx, params):
        """
        Save multiple assembly objects from FASTA files.
        WARNING: The code currently saves all assembly object data in memory before sending it
        to the workspace in a single batch. Since the object data doesn't include sequences,
        it is typically small and so in most cases this shouldn't cause issues. However, many
        assemblies and / or many contigs could conceivably cause memeory issues or could
        cause the workspace to reject the data package if the serialized data is > 1GB.
        TODO: If this becomes a common issue (not particularly likely?) update the code to
         Save assembly object data on disk if it becomes too large
         Batch uploads to the workspace based on data size
        :param params: instance of type "SaveAssembliesParams" (Input for the
           save_assemblies_from_fastas function. Required arguments:
           workspace_id - the numerical ID of the workspace in which to save
           the Assemblies. inputs - a list of FASTA files to import. All of
           the files must be from the same source - either all local files or
           all Blobstore nodes. Optional arguments: min_contig_length - an
           integer > 1. If present, sequences of lesser length will be
           removed from the input FASTA files.) -> structure: parameter
           "workspace_id" of Long, parameter "inputs" of list of type
           "FASTAInput" (An input FASTA file and metadata for import.
           Required arguments: Exactly one of: file - a path to an input
           FASTA file. Must be accessible inside the AssemblyUtil docker
           continer. node - a node ID for a Blobstore (formerly Shock) node
           containing an input FASTA file. assembly_name - the workspace name
           under which to save the Assembly object. Optional arguments: type
           - should be one of 'isolate', 'metagenome', (maybe
           'transcriptome'). Defaults to 'Unknown' external_source - the
           source of the input data. E.g. JGI, NCBI, etc. external_source_id
           - the ID of the input data at the source. contig_info - map from
           contig_id to a small structure that can be used to set the
           is_circular and description fields for Assemblies) -> structure:
           parameter "file" of String, parameter "node" of String, parameter
           "assembly_name" of String, parameter "type" of String, parameter
           "external_source" of String, parameter "external_source_id" of
           String, parameter "contig_info" of mapping from String to type
           "ExtraContigInfo" (Structure for setting additional Contig
           information per contig is_circ - flag if contig is circular, 0 is
           false, 1 is true, missing indicates unknown description - if set,
           sets the description of the field in the assembly object which may
           override what was in the fasta file) -> structure: parameter
           "is_circ" of Long, parameter "description" of String, parameter
           "min_contig_length" of Long
        :returns: instance of type "SaveAssembliesResults" (Results for the
           save_assemblies_from_fastas function. results - the results of the
           save operation in the same order as the input.) -> structure:
           parameter "results" of list of type "SaveAssemblyResult" (Results
           from saving an assembly. upa - the address of the resulting
           workspace object. filtered_input - the filtered input file if the
           minimum contig length parameter is present and > 0. null
           otherwise.) -> structure: parameter "upa" of type "upa" (A Unique
           Permanent Address for a workspace object, which is of the form
           W/O/V, where W is the numeric workspace ID, O is the numeric
           object ID, and V is the object version.), parameter
           "filtered_input" of String
        """
        # ctx is the context object
        # return variables are: results
        #BEGIN save_assemblies_from_fastas
        results = {
            'results': FastaToAssembly(
                DataFileUtil(self.callback_url, token=ctx['token']),
                Path(self.sharedFolder)
            ).import_fasta_mass(params, self.threads_per_cpu, self.max_threads)
        }
        #END save_assemblies_from_fastas

        # At some point might do deeper type checking...
        if not isinstance(results, dict):
            raise ValueError('Method save_assemblies_from_fastas return value ' +
                             'results is not type dict as required.')
        # return the results
        return [results]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
