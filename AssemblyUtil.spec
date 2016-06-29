/*

*/
module AssemblyUtil {


    typedef structure {
        string path;
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
    funcdef get_assembly_as_fasta(GetAssemblyParams params) returns (FastaAssemblyFile file) authentication required;



    /*
        Options supported:
            workspace_name
            assembly_name

        Uploader options not yet supported
            taxon_reference: The ws reference the assembly points to.  (Optional)
            source: The source of the data (Ex: Refseq)
            date_string: Date (or date range) associated with data. (Optional)
            contig_information_dict: A mapping that has is_circular and description information (Optional)
    */
    typedef structure {
        FastaAssemblyFile file;

        string workspace_name;
        string assembly_name;

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
