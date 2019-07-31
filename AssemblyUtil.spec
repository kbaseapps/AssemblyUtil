/*

*/
module AssemblyUtil {


    typedef structure {
        string path;
        string assembly_name;
    } FastaAssemblyFile;

    /*
        @optional filename
    */
    typedef structure {
        string ref;
        string filename;
    } GetAssemblyParams;

    /*
        Given a reference to an Assembly (or legacy ContigSet data object), along with a set of options,
        construct a local Fasta file with the sequence data.  If filename is set, attempt to save to the
        specified filename.  Otherwise, a random name will be generated.
    */
    funcdef get_assembly_as_fasta(GetAssemblyParams params) 
                returns (FastaAssemblyFile file) authentication required;

    /*
        Structure for get_fasta function input and output:
            rKBaseOjbReferences - is an object wrapperd array of KBase object references, which can be of 
            the following types:
                - KBaseGenomes.Genome
                - KBaseSets.AssemblySet
                - KBaseMetagenome.BinnedContigs
                - KBaseGenomes.ContigSet
                - KBaseGenomeAnnotations.Assembly
                - KBaseSearch.GenomeSet
                - KBaseSets.GenomeSet

        ref - workspace reference.

        get_fastas_output_value
            paths - list of paths to fasta files associated with workspace object.
            type - workspace object type
            parent_refs - (optional) associated workspace object references if different from the output key
    */

    typedef string ref;

    typedef structure {
        list<ref> ref_lst;
    } KBaseOjbReferences;

    typedef structure {
        list<string> paths;
        list<ref> parent_refs;
        string type;
    } ref_fastas;
    /*
        Given a reference list of KBase objects constructs a local Fasta file with the sequence data for each ref.
    */
    funcdef get_fastas(KBaseOjbReferences params)
                returns (mapping<ref, ref_fastas> output) authentication required;

    typedef structure {
        string input_ref;
    } ExportParams;

    typedef structure {
        string shock_id;
    } ExportOutput;

    /*
        A method designed especially for download, this calls 'get_assembly_as_fasta' to do
        the work, but then packages the output with WS provenance and object info into
        a zip file and saves to shock.
    */
    funcdef export_assembly_as_fasta(ExportParams params)
                returns (ExportOutput output) authentication required;


    typedef string ShockNodeId;



    /*
        Structure for setting additional Contig information per contig
            is_circ - flag if contig is circular, 0 is false, 1 is true, missing
                      indicates unknown
            description - if set, sets the description of the field in the assembly object
                          which may override what was in the fasta file
    */
    typedef structure {
        int is_circ;
        string description;
    } ExtraContigInfo;


    /*
        Options supported:
            file / shock_id / ftp_url - mutualy exclusive parameters pointing to file content
            workspace_name - target workspace
            assembly_name - target object name

            type - should be one of 'isolate', 'metagenome', (maybe 'transcriptome')

            min_contig_length - if set and value is greater than 1, this will only include sequences
                                with length greater or equal to the min_contig_length specified, discarding
                                all other sequences

            taxon_ref         - sets the taxon_ref if present

            contig_info       - map from contig_id to a small structure that can be used to set the is_circular
                                and description fields for Assemblies (optional)

        Uploader options not yet supported
            taxon_reference: The ws reference the assembly points to.  (Optional)
            source: The source of the data (Ex: Refseq)
            date_string: Date (or date range) associated with data. (Optional)
    */
    typedef structure {
        FastaAssemblyFile file;
        ShockNodeId shock_id;
        string ftp_url;

        string workspace_name;
        string assembly_name;

        string external_source;
        string external_source_id;

        string taxon_ref;

        int min_contig_length;

        mapping<string,ExtraContigInfo> contig_info; 

    } SaveAssemblyParams;

    /*
        WARNING: has the side effect of moving the file to a temporary staging directory, because the upload
        script for assemblies currently requires a working directory, not a specific file.  It will attempt
        to upload everything in that directory.  This will move the file back to the original location, but
        if you are trying to keep an open file handle or are trying to do things concurrently to that file,
        this will break.  So this method is certainly NOT thread safe on the input file.
    */
    funcdef save_assembly_from_fasta(SaveAssemblyParams params) returns (string ref) authentication required;
};
