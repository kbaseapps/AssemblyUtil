

import os
import time as _time

from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta
from installed_clients.baseclient import ServerError as _MGUError


from shutil import copyfile

class TypeToFasta:

    def __init__(self, callback_url, scratch, wrkspc):
        self.ws = wrkspc
        self.scratch = scratch
        self.callback_url = callback_url
        self.mgu = MetagenomeUtils(callback_url)

    def log(self, message, prefix_newline=False):
        print(('\n' if prefix_newline else '') + str(_time.time()) + ': ' + message)

    def genome_obj_to_fasta(self, ref, obj_type, fasta_array):

        upas = []
        atf = AssemblyToFasta(self.callback_url, self.scratch)

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
                genome_data = self.ws.get_objects2({'objects': [{"ref": genome_upa}]})['data'][0]['data']
                assembly_upa = genome_upa + ';' + str(genome_data.get('contigset_ref') \
                                                      or genome_data.get('assembly_ref'))

                faf = atf.assembly_as_fasta({'ref': assembly_upa})
                fasta_array.extend([faf['path'], assembly_upa])

        return fasta_array

    def assembly_obj_to_fasta(self, ref, obj_type, fasta_array):

        atf = AssemblyToFasta(self.callback_url, self.scratch)
        obj = {"ref": ref}

        if "KBaseGenomes.ContigSet" in obj_type or "KBaseGenomeAnnotations.Assembly" in obj_type:
            faf = atf.assembly_as_fasta(obj)
            fasta_array.extend([faf['path'], ref])

        elif "KBaseSets.AssemblySet" in obj_type:

            obj_data = self.ws.get_objects2({'objects': [{"ref": ref}]})['data'][0]

            for item_upa in obj_data['data']['items']:
                faf = atf.assembly_as_fasta({"ref": item_upa['ref']})
                fasta_array.extend([faf['path'], item_upa['ref']])

        return fasta_array

    def metagenome_obj_to_fasta(self, ref, obj_type, fasta_array):

        if 'KBaseMetagenomes.BinnedContigs' in obj_type:
            try:
                bin_file_dir = self.mgu.binned_contigs_to_file({'input_ref': ref, 'save_to_shock': 0}) \
                    ['bin_file_directory']
                for (dirpath, dirnames, filenames) in os.walk(bin_file_dir):
                    for fasta_file in filenames:
                        fasta_path = os.path.join(self.scratch, fasta_file)
                        copyfile(os.path.join(bin_file_dir, fasta_file), fasta_path)
                        fasta_array.extend([fasta_path, ref])
            except _MGUError as mgue:
                # not really any way to test this block
                self.log('Logging exception loading binned contigs to file.')
                self.log(str(mgue))
                raise

        return fasta_array

    def type_to_fasta(self, ref_lst):


        # Initiate objects
        fasta_dict, fasta_array = dict(), []

        # Get type info for each ref in ref_lst
        for idx, ref in enumerate(ref_lst):
            # Initiate objects and get KBase object type with get_object_info3

            obj_info = self.ws.get_object_info3({"objects": [{"ref": ref}]})
            obj_type = obj_info["infos"][0][2]

            fasta_array = self.genome_obj_to_fasta(ref, obj_type, fasta_array)
            fasta_array = self.assembly_obj_to_fasta(ref, obj_type, fasta_array)
            fasta_array = self.metagenome_obj_to_fasta(ref, obj_type, fasta_array)

        # return dictionary of FASTA
        fasta_dict["FASTA"] = fasta_array

        return fasta_dict
