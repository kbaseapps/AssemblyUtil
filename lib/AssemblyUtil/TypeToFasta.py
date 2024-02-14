

import os
import time as _time

from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta
from installed_clients.baseclient import ServerError as _MGUError


from shutil import copyfile

class TypeToFasta:

    def __init__(self, callback_url, scratch, wrkspc, token):
        self.ws = wrkspc
        self.scratch = scratch
        self.callback_url = callback_url
        self.mgu = MetagenomeUtils(callback_url, token=token)
        self.fasta_dict = {}

    def log(self, message, prefix_newline=False):
        print(('\n' if prefix_newline else '') + str(_time.time()) + ': ' + message)

    def add_to_dict(self, key, val):
        if key in self.fasta_dict:
            # if key is already dict, we want to add a field to the the 'parent_refs'
            if 'parent_refs' in self.fasta_dict[key]:
                self.fasta_dict[key]['parent_refs'] += val['parent_refs']
        else:
            self.fasta_dict[key] = val

    def genome_obj_to_fasta(self, ref, obj_type):

        # Initiate needed objects
        atf = AssemblyToFasta(self.callback_url, self.scratch)
        upas = []

        if 'KBaseSets.GenomeSet' in obj_type:
            obj_data = self.ws.get_objects2({'objects': [{"ref": ref}]})['data'][0]
            upas = [gsi['ref'] for gsi in obj_data['data']['items']]
        elif 'KBaseSearch.GenomeSet' in obj_type:
            obj_data = self.ws.get_objects2({'objects': [{"ref": ref}]})['data'][0]
            upas = [gse['ref'] for gse in obj_data['data']['elements'].values()]
        elif "KBaseGenomes.Genome" in obj_type:
            upas = [ref]

        if upas:
            for genome_upa in upas:
                # Get genome object assembly_ref or contigset_ref through subsetting object
                genome_data = self.ws.get_objects2({'objects': \
                            [{"ref": genome_upa, 'included' : ['/assembly_ref/','/contigset_ref/']}]}) \
                            ['data'][0]['data']

                # If genome object contains an assembly_ref or contigset_ref it will return a dictionary, genome_data.
                # If not an empty dictionary will be returned
                if genome_data:
                    # Get assembly_upa and fasta
                    assembly_upa = genome_upa + ';' + \
                                   str(genome_data.get('assembly_ref') or genome_data.get('contigset_ref'))

                    faf = atf.assembly_as_fasta({'ref': assembly_upa})
                    # Input data into object dict
                    self.add_to_dict(assembly_upa, {'paths' : [faf['path']], 'type': obj_type, 'parent_refs': [ref]})

                else:
                    raise TypeError("KBase object type %s does not contain an assembly reference or contig reference." % obj_type)


    def assembly_obj_to_fasta(self, ref, obj_type, input_ref=None, input_type=None):
        # Initiate needed objects
        atf = AssemblyToFasta(self.callback_url, self.scratch)
        obj = {"ref": ref}

        if "KBaseGenomes.ContigSet" in obj_type or "KBaseGenomeAnnotations.Assembly" in obj_type:
            # Get fasta
            faf = atf.assembly_as_fasta(obj)
            if input_ref and input_type:
                self.add_to_dict(input_ref, {'paths': [faf['path']], 'type': input_type, 'parent_refs': [input_ref, ref]})
            else:
                self.add_to_dict(ref, {'paths': [faf['path']], 'type': obj_type, 'parent_refs': [ref]})

        elif "KBaseSets.AssemblySet" in obj_type:
            # Get assembly set object
            obj_data = self.ws.get_objects2({'objects': [obj]})['data'][0]
            for item_upa in obj_data['data']['items']:
                # Get fasta
                faf = atf.assembly_as_fasta({"ref": item_upa['ref']})
                # Input data into object dict
                self.add_to_dict(item_upa['ref'], {'paths' : [faf['path']], 'type' : obj_type, 'parent_refs': [ref]})

    def metagenome_obj_to_fasta(self, ref, obj_type):

        if 'KBaseMetagenomes.BinnedContigs' in obj_type:
            fasta_paths = []
            try:
                # Binned_contigs_to_file saves fasta file to a directory in scratch.
                # Path: scratch/binned_contig_files_EXTENSION/Bin#.fasta
                bin_file_dir = self.mgu.binned_contigs_to_file({'input_ref': ref, 'save_to_shock': 0}) \
                    ['bin_file_directory']
                for (dirpath, dirnames, filenames) in os.walk(bin_file_dir):
                    for fasta_file in filenames:
                        # For fasta file in the binned contigs directory, copy fasta directly to scratch
                        # New path: scratch/Bin#.fasta
                        fasta_path = os.path.join(self.scratch, fasta_file)
                        copyfile(os.path.join(bin_file_dir, fasta_file), fasta_path)
                        fasta_paths.append(fasta_path)
                # Input data into object dict
                self.add_to_dict(ref, {'paths' : fasta_paths, 'type': obj_type})

            # Catch MetagenomeUtil Error
            except _MGUError as mgue:
                self.log('Logging exception loading binned contigs to file.')
                self.log(str(mgue))
                raise

        if 'KBaseMetagenomes.AnnotatedMetagenomeAssembly' in obj_type:
            ret = self.ws.get_objects2({'objects': [{'ref': ref, 'included': ['assembly_ref']}]})['data'][0]
            assembly_ref = ret['data']['assembly_ref']
            assembly_obj_type = self.ws.get_object_info3({'objects': [{'ref': assembly_ref}]})['infos'][0][2]
            self.assembly_obj_to_fasta(assembly_ref, assembly_obj_type, input_ref=ref, input_type=obj_type)

    def protein_seq_set_to_fasta(self, ref, obj_type):
        """
        Convert a KBaseSequences.ProteinSequenceSet object to a FASTA file.
        """

        # Check if the object type is ProteinSequenceSet
        if 'KBaseSequences.ProteinSequenceSet' in obj_type:
            # Fetch the ProteinSequenceSet object
            protein_seq_set = self.ws.get_objects2({'objects': [{'ref': ref}]})['data'][0]['data']
            
            # Create a FASTA file
            fasta_file_path = os.path.join(self.scratch, ref.replace('/', '_') + '.fasta')
            with open(fasta_file_path, 'w') as fasta_file:
                for protein_sequence in protein_seq_set['sequences']:
                    # Define a FASTA format record
                    fasta_record = f">{protein_sequence['id']}\n{protein_sequence['sequence']}\n"
                    fasta_file.write(fasta_record)

            # Add the path to the fasta_dict
            self.add_to_dict(ref, {'paths': [fasta_file_path], 'type': obj_type, 'parent_refs': [ref]})
    
    def type_to_fasta(self, ref_lst):
        """type_to_fasta takes in a list of KBase objects references. The object type of each reference
        is checked in functions: assembly_obj_to_fasta, metagenome_obj_to_fasta, and genome_obj_to_fasta. Depending
        on the type of KBase object input a fasta file is made through one of the functions mentioned above
        and a fasta object dictionary is created with structure: {ref: {'path' : fasta_paths, 'type': object type} }

        for objects of type AssemblySet and GenomeSet a parent ref key-value pair is added such that the structure is:
        {ref: {'path' : fasta_paths, 'type': object type, 'parent_refs': [ref]} }

        for objects of type KBaseMetagenomes.BinnedContigs a unique fasta path is made for each bin in binnedContigs
        Thus the output structure is: {ref: {'paths' : [fasta_contigbin1, fasta_contigbin2], 'type': object type} }

        where the key 'paths' points to an array of fasta paths for each contig bin in ascending order. """

        # Get type info for each ref in ref_lst
        for idx, ref in enumerate(ref_lst):

            # Get KBase object type with get_object_info3
            obj_info = self.ws.get_object_info3({"objects": [{"ref": ref}]})
            obj_type = obj_info["infos"][0][2]
            # Put object in object specific fasta dictionary by type
            self.genome_obj_to_fasta(ref, obj_type)
            self.assembly_obj_to_fasta(ref, obj_type)
            self.metagenome_obj_to_fasta(ref, obj_type)
            self.protein_seq_set_to_fasta(ref, obj_type)
            # Append all individual object dictionaries to complete fasta dictionary for reference list
        
        return self.fasta_dict
