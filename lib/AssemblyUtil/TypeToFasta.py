

import os

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from AssemblyUtil.AssemblyToFasta import AssemblyToFasta

from shutil import copyfile

class TypeToFasta:

    def __init__(self, callback_url, scratch, ws_url):
        self.ws_url = ws_url
        self.callback_url = callback_url
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)
        self.mgu = MetagenomeUtils(callback_url)

    def type_to_fasta(self, ctx, ref_lst):


        fasta_dict = dict()
        fasta_array = []
        atf = AssemblyToFasta(self.callback_url, self.scratch)

        # Get type info for each ref in ref_lst
        for idx, ref in enumerate(ref_lst):

            upas = []
            obj = {"ref": ref}
            obj_info = self.ws_url.get_object_info3({"objects": [obj]})
            obj_type = obj_info["infos"][0][2]


            # From type info get object
            if 'KBaseSets.GenomeSet' in obj_type:
                obj_data = self.dfu.get_objects({"object_refs": [ref]})['data'][0]
                upas = [gsi['ref'] for gsi in obj_data['data']['items']]
            elif 'KBaseSearch.GenomeSet' in obj_type:
                obj_data = self.dfu.get_objects({"object_refs": [ref]})['data'][0]
                upas = [gse['ref'] for gse in obj_data['data']['elements'].values()]
            elif "KBaseGenomes.Genome" in obj_type:
                upas = [ref]

            elif "KBaseGenomes.ContigSet" in obj_type or "KBaseGenomeAnnotations.Assembly" in obj_type:
                faf = [atf.assembly_as_fasta(ctx, obj)]
                fasta_array.extend([faf[0]['path'], ref])

            elif "KBaseSets.AssemblySet" in obj_type:
                fasta_paths = []

                obj_data = self.dfu.get_objects({"object_refs": [ref]})['data'][0]

                for item_upa in obj_data['data']['items']:
                    faf = [atf.assembly_as_fasta(ctx, {"ref": item_upa['ref']})]
                    fasta_paths.extend([faf[0]['path'], item_upa['ref']])
                    fasta_array = fasta_paths

            elif 'KBaseMetagenomes.BinnedContigs' in obj_type:
                fasta_paths = []

                bin_file_dir = self.mgu.binned_contigs_to_file({'input_ref': ref, 'save_to_shock': 0})['bin_file_directory']
                for (dirpath, dirnames, filenames) in os.walk(bin_file_dir):
                    for fasta_file in filenames:
                        fasta_path = os.path.join(self.scratch, fasta_file)
                        copyfile(os.path.join(bin_file_dir, fasta_file), fasta_path)
                        fasta_paths.extend([fasta_path, ref])
                    break
                fasta_array = fasta_paths

            if upas:
                for genome_upa in upas:
                    genome_data = self.ws_url.get_objects2({'objects': [{"ref": genome_upa}]})['data'][0]['data']
                    assembly_upa = genome_upa + ';' + str(genome_data.get('contigset_ref') or genome_data.get('assembly_ref'))
                    faf = [atf.assembly_as_fasta(ctx, {'ref': assembly_upa})]
                    fasta_array.extend([faf[0]['path'], assembly_upa])

        # return dictionary of FASTA
        fasta_dict["FASTA"] = fasta_array

        return fasta_dict
