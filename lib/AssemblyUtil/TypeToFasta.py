

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

    def genome_obj_to_fasta(self, ref, obj_type, fasta_dict):

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
                genome_data = self.ws.get_objects2({'objects': [{"ref": genome_upa, 'included' : ['/assembly_ref/']}]}) \
                                ['data'][0]['data']
                if genome_data:
                    assembly_upa = genome_upa + ';' + str(genome_data.get('assembly_ref'))
                    faf = atf.assembly_as_fasta({'ref': assembly_upa})
                    fasta_dict[assembly_upa] = faf['path']

                else:
                    raise TypeError("KBase object type %s does not contain an assembly reference." % obj_type)

        return fasta_dict

    def assembly_obj_to_fasta(self, ref, obj_type, fasta_dict):

        atf = AssemblyToFasta(self.callback_url, self.scratch)
        obj = {"ref": ref}

        if "KBaseGenomes.ContigSet" in obj_type or "KBaseGenomeAnnotations.Assembly" in obj_type:
            faf = atf.assembly_as_fasta(obj)
            fasta_dict[ref] = faf['path']

        elif "KBaseSets.AssemblySet" in obj_type:

            obj_data = self.ws.get_objects2({'objects': [{"ref": ref}]})['data'][0]

            for item_upa in obj_data['data']['items']:
                faf = atf.assembly_as_fasta({"ref": item_upa['ref']})
                fasta_dict[item_upa['ref']] = faf['path']

        return fasta_dict

    def metagenome_obj_to_fasta(self, ref, obj_type, fasta_dict):

        if 'KBaseMetagenomes.BinnedContigs' in obj_type:
            fasta_paths = []
            try:
                bin_file_dir = self.mgu.binned_contigs_to_file({'input_ref': ref, 'save_to_shock': 0}) \
                    ['bin_file_directory']
                for (dirpath, dirnames, filenames) in os.walk(bin_file_dir):

                    for fasta_file in filenames:
                        # Copy file to scratch
                        fasta_path = os.path.join(self.scratch, fasta_file)
                        fasta_path = os.path.splitext(fasta_path)[0] + ".fa"
                        copyfile(os.path.join(bin_file_dir, fasta_file), fasta_path)
                        fasta_paths.append(fasta_path)
                fasta_dict[ref] = fasta_paths

            except _MGUError as mgue:
                self.log('Logging exception loading binned contigs to file.')
                self.log(str(mgue))
                raise

        return fasta_dict

    def type_to_fasta(self, ref_lst):

        # Initiate objects
        fasta_dict = {}

        # Get type info for each ref in ref_lst
        for idx, ref in enumerate(ref_lst):
            # Initiate objects and get KBase object type with get_object_info3

            obj_info = self.ws.get_object_info3({"objects": [{"ref": ref}]})
            obj_type = obj_info["infos"][0][2]

            fasta_dict_genome_obj = self.genome_obj_to_fasta(ref, obj_type, fasta_dict)
            fasta_dict_assembly_obj = self.assembly_obj_to_fasta(ref, obj_type, fasta_dict_genome_obj)
            fasta_dict_metagenome_obj = self.metagenome_obj_to_fasta(ref, obj_type, fasta_dict_assembly_obj)

            fasta_dict = {**fasta_dict_genome_obj, **fasta_dict_assembly_obj, **fasta_dict_metagenome_obj}
        
        return fasta_dict
