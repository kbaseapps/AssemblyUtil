import json
import os
import os.path
import sys
import uuid
from collections import Counter
from hashlib import md5

from Bio import SeqIO

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


def sort_dict(in_struct):
    """Recursively sort a dictionary by dictionary keys. (saves WS the trouble)"""
    if isinstance(in_struct, dict):
        return {k: sort_dict(in_struct[k]) for k in sorted(in_struct)}
    elif isinstance(in_struct, list):
        return [sort_dict(k) for k in sorted(in_struct)]
    else:
        return in_struct


class FastaToAssembly:

    def __init__(self, callback_url, scratch, ws_url):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)
        self.ws = Workspace(ws_url)

        # Note added X due to kb|g.1886.fasta
        self.valid_chars = "-ACGTUWSMKRYBDHVNX"
        self.amino_acid_specific_characters = "PLIFQE"

    def import_fasta(self, ctx, params):
        print('validating parameters')
        self.validate_params(params)

        print('staging input files')
        fasta_file_path = self.stage_input(params)

        if 'min_contig_length' in params:
            min_contig_length = int(params['min_contig_length'])
            print(f'filtering FASTA file by contig length (min len={min_contig_length} bp)')
            fasta_file_path = self.filter_contigs_by_length(fasta_file_path, min_contig_length)

        print(f'parsing FASTA file: {fasta_file_path}')
        assembly_data = self.parse_fasta(fasta_file_path, params)
        print(f' - parsed {assembly_data["num_contigs"]} contigs,{assembly_data["dna_size"]} bp')
        print('saving assembly to KBase')

        # save file to shock and build handle
        fasta_file_handle_info = self.save_fasta_file_to_shock(fasta_file_path)
        # construct the output object
        assembly_object_to_save = self.build_assembly_object(assembly_data,
                                                             fasta_file_handle_info,
                                                             params)
        json.dump(assembly_object_to_save, open(self.scratch+"/example.json", 'w'))

        # save to WS and return
        if 'workspace_id' in params:
            workspace_id = int(params['workspace_id'])
        else:
            workspace_id = self.dfu.ws_name_to_id(params['workspace_name'])
        assembly_info = self.save_assembly_object(workspace_id,
                                                  params['assembly_name'],
                                                  assembly_object_to_save)

        return assembly_info

    def build_assembly_object(self, assembly_data, fasta_file_handle_info, params):
        """ construct the WS object data to save based on the parsed info and params """
        assembly_data['assembly_id'] = params['assembly_name']
        assembly_data['fasta_handle_ref'] = fasta_file_handle_info['handle']['hid']
        fasta_file_handle_info['handle'] = fasta_file_handle_info['handle']
        assembly_data['fasta_handle_info'] = fasta_file_handle_info

        assembly_data['type'] = 'Unknown'
        if 'type' in params:
            assembly_data['type'] = params['type']

        if 'taxon_ref' in params:
            info = self.ws.get_object_info3({'objects':[{'ref': params['taxon_ref']}]})['infos'][0]
            assembly_data['taxon_ref'] = f'{info[6]}/{info[0]}/{info[4]}'

        if 'external_source' in params:
            assembly_data['external_source'] = params['external_source']

        if 'external_source_id' in params:
            assembly_data['external_source_id'] = params['external_source_id']

        if 'external_source_origination_date' in params:
            assembly_data['external_source_origination_date'] = params['external_source_origination_date']

        return sort_dict(assembly_data)

    def parse_fasta(self, fasta_file_path, params):
        """ Do the actual work of inspecting each contig """

        # variables to store running counts of things
        total_length = 0
        base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        md5_list = []

        # map from contig_id to contig_info
        all_contig_data = {}
        extra_contig_info = {}
        if'contig_info' in params:
            extra_contig_info = params['contig_info']

        for record in SeqIO.parse(fasta_file_path, "fasta"):
            # SeqRecord(seq=Seq('TTAT...', SingleLetterAlphabet()),
            #           id='gi|113968346|ref|NC_008321.1|',
            #           name='gi|113968346|ref|NC_008321.1|',
            #           description='gi|113968346|ref|NC_008321.1| Shewanella sp. MR-4 chromosome, complete genome',
            #           dbxrefs=[])

            sequence = str(record.seq).upper()

            contig_info = {
                'contig_id': record.id,
                'name': record.id,
                'description': record.description[len(record.id):].strip(),
                'length': len(record.seq)
            }

            # 1) compute sequence character statistics running total
            total_length += contig_info['length']
            sequence_count_table = dict(Counter(sequence))
            for character in sequence_count_table:
                if character in base_counts:
                    base_counts[character] = base_counts[character] + sequence_count_table[character]
                else:
                    base_counts[character] = sequence_count_table[character]
                if character not in self.valid_chars:
                    if character in self.amino_acid_specific_characters:
                        raise ValueError('This FASTA file may have amino acids in it instead '
                                         'of the required nucleotides.')
                    raise ValueError(f"This FASTA file has non nucleic acid characters: "
                                     f"{character}")

            # 2) record number of 'N' characters (only set if there are some)
            Ncount = 0
            if 'N' in sequence_count_table:
                Ncount = sequence_count_table['N']
                contig_info['Ncount'] = Ncount

            # 2b) record if the contig is circular
            if record.id in extra_contig_info:
                if 'is_circ' in extra_contig_info[record.id]:
                    contig_info['is_circ'] = int(extra_contig_info[record.id]['is_circ'])
                if 'description' in extra_contig_info[record.id]:
                    contig_info['description'] = str(extra_contig_info[record.id]['description'])

            # 3) record md5 checksum
            contig_md5 = md5(sequence.encode()).hexdigest()
            contig_info['md5'] = contig_md5
            md5_list.append(contig_md5)

            # 4) record the all important GC to ~3 significant digits
            GC_count = 0
            for base in ['G', 'C']:
                if base in sequence_count_table:
                    GC_count += sequence_count_table[base]
            contig_info['gc_content'] = round(float(GC_count) / float(contig_info['length']), 5)

            # 5) add to contig list
            if contig_info['contig_id'] in all_contig_data:
                raise ValueError('The FASTA header key ' + contig_info['contig_id'] +
                                 'appears more than once in the file')

            all_contig_data[contig_info['contig_id']] = contig_info

        # Aggregate stats for the data
        total_gc_content = None
        if total_length > 0:
            total_gc_content = round(float(base_counts['G'] + base_counts['C']) / float(total_length), 5)
        assembly_data = {
            'md5': md5(",".join(sorted(md5_list)).encode()).hexdigest(),
            'base_counts': base_counts,
            'dna_size': total_length,
            'gc_content': total_gc_content,
            'contigs': all_contig_data,
            'num_contigs': len(all_contig_data)
        }
        return assembly_data

    @staticmethod
    def fasta_filter_contigs_generator(fasta_record_iter, min_contig_length):
        """ generates SeqRecords iterator for writing from a legacy contigset object """
        rows = 0
        rows_added = 0
        for record in fasta_record_iter:
            rows += 1
            if len(record.seq) >= min_contig_length:
                rows_added += 1
                yield record
        print(f' - filtered out {rows - rows_added} of {rows} contigs that were shorter '
              f'than {(min_contig_length)} bp.')

    def filter_contigs_by_length(self, fasta_file_path, min_contig_length):
        """ removes all contigs less than the min_contig_length provided """
        filtered_fasta_file_path = fasta_file_path + '.filtered.fa'

        fasta_record_iter = SeqIO.parse(fasta_file_path, 'fasta')
        SeqIO.write(self.fasta_filter_contigs_generator(fasta_record_iter, min_contig_length),
                    filtered_fasta_file_path, 'fasta')

        return filtered_fasta_file_path

    def save_assembly_object(self, workspace_id, assembly_name, obj_data):
        print('Saving Assembly to Workspace')
        sys.stdout.flush()
        if len(obj_data["contigs"]) == 0:
            raise ValueError('There are no contigs to save, thus there is no valid assembly.')
        obj_info = self.dfu.save_objects({'id': workspace_id,
                                          'objects': [{'type': 'KBaseGenomeAnnotations.Assembly',
                                                       'data': obj_data,
                                                       'name': assembly_name
                                                       }]
                                          })[0]
        return obj_info

    def save_fasta_file_to_shock(self, fasta_file_path):
        """ Given the path to the file, upload to shock and return Handle information
            returns:
                typedef structure {
                    string shock_id;
                    Handle handle;
                    string node_file_name;
                    string size;
                } FileToShockOutput;

        """
        print(f'Uploading FASTA file ({fasta_file_path}) to SHOCK')
        sys.stdout.flush()
        return self.dfu.file_to_shock({'file_path': fasta_file_path, 'make_handle': 1})

    def stage_input(self, params):
        """ Setup the input_directory by fetching the files and returning the path to the file"""
        file_path = None
        if 'file' in params:
            if not os.path.isfile(params['file']['path']):
                raise ValueError('KBase Assembly Utils tried to save an assembly, but the calling application specified a file ('+params['file']['path']+') that is missing. Please check the application logs for details.')
            file_path = os.path.abspath(params['file']['path'])
        elif 'shock_id' in params:
            print(f'Downloading file from SHOCK node: {params["shock_id"]}')
            sys.stdout.flush()
            input_directory = os.path.join(self.scratch, 'assembly-upload-staging-' + str(uuid.uuid4()))
            os.makedirs(input_directory)
            file_name = self.dfu.shock_to_file({'file_path': input_directory,
                                                'shock_id': params['shock_id']
                                                })['node_file_name']
            file_path = os.path.join(input_directory, file_name)
        elif 'ftp_url' in params:
            print(f'Downloading file from: {params["ftp_url"]}')
            sys.stdout.flush()
            file_path = self.dfu.download_web_file({'file_url': params['ftp_url'],
                                                    'download_type': 'FTP'
                                                    })['copy_file_path']

        # extract the file if it is compressed
        if file_path is not None:
            unpacked_file = self.dfu.unpack_file({'file_path': file_path})
            return unpacked_file['file_path']

        raise ValueError('No valid FASTA could be extracted based on the input parameters')


    @staticmethod
    def validate_params(params):
        for key in ('workspace_name', 'assembly_name'):
            if key not in params:
                raise ValueError('required "' + key + '" field was not defined')

        # one and only one of either 'file', 'shock_id', or ftp_url is required
        input_count = 0
        for key in ('file', 'shock_id', 'ftp_url'):
            if key in params and params[key] is not None:
                input_count = input_count + 1
                if key == 'file':
                    if not isinstance(params[key], dict) or 'path' not in params[key]:
                        raise ValueError('when specifying a FASTA file input, "path" field was not defined in "file"')

        if input_count == 0:
            raise ValueError('required FASTA file as input, set as either "file", "shock_id", or "ftp_url"')
        if input_count > 1:
            raise ValueError('required exactly one FASTA file as input source, you set more than one of ' +
                             'these fields: "file", "shock_id", or "ftp_url"')
