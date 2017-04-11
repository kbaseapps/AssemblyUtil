import os
import sys
import uuid

from hashlib import md5
from collections import Counter

from Bio import SeqIO

from DataFileUtil.DataFileUtilClient import DataFileUtil


class FastaToAssembly:

    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)

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
            print('filtering fasta file by contig length (min len=' + str(min_contig_length) + 'bp)')
            fasta_file_path = self.filter_contigs_by_length(fasta_file_path, min_contig_length)

        print('parsing FASTA file: ' + str(fasta_file_path))
        assembly_data = self.parse_fasta(fasta_file_path, params)
        print(' - parsed ' + str(assembly_data['num_contigs']) + ' contigs, ' +
              str(assembly_data['dna_size']) + 'bp')

        print('saving assembly to KBase')

        # save file to shock and build handle
        fasta_file_handle_info = self.save_fasta_file_to_shock(fasta_file_path)
        # construct the output object
        assembly_object_to_save = self.build_assembly_object(assembly_data,
                                                             fasta_file_handle_info,
                                                             params)

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
        ''' construct the WS object data to save based on the parsed info and params '''
        assembly_data['assembly_id'] = params['assembly_name']
        assembly_data['fasta_handle_ref'] = fasta_file_handle_info['handle']['hid']
        assembly_data['fasta_handle_info'] = fasta_file_handle_info

        assembly_data['type'] = 'Unknown'
        if 'type' in params:
            assembly_data['type'] = params['type']

        if 'taxon_ref' in params:
            assembly_data['taxon_ref'] = params['taxon_ref']

        if 'external_source' in params:
            assembly_data['external_source'] = params['external_source']

        if 'external_source_id' in params:
            assembly_data['external_source_id'] = params['external_source_id']

        if 'external_source_origination_date' in params:
            assembly_data['external_source_origination_date'] = params['external_source_origination_date']

        return assembly_data


    def parse_fasta(self, fasta_file_path, params):
        ''' Do the actual work of inspecting each contig '''

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
                        raise ValueError('This fasta file may have amino acids in it instead ' +
                                         'of the required nucleotides.')
                    raise ValueError("This FASTA file has non nucleic acid characters : {0}".format(character))

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
            contig_md5 = md5(sequence).hexdigest()
            contig_info['md5'] = contig_md5
            md5_list.append(contig_md5)

            # 4) record the all important GC to ~3 significant digits
            GC_count = 0
            for base in ['G', 'C']:
                if base in sequence_count_table:
                    GC_count += sequence_count_table[base]
            contig_info['gc_content'] = round(float(GC_count) / float(contig_info['length']), 3)

            # 5) add to contig list
            if contig_info['contig_id'] in all_contig_data:
                raise ValueError('The fasta header key ' + contig_info['contig_id'] +
                                 'appears more than once in the file')
            all_contig_data[contig_info['contig_id']] = contig_info

        # Aggregate stats for the data
        total_gc_content = None
        if total_length > 0:
            total_gc_content = round(float(base_counts['G'] + base_counts['C']) / float(total_length), 3)
        assembly_data = {
            'md5': md5(",".join(sorted(md5_list))).hexdigest(),
            'base_counts': base_counts,
            'dna_size': total_length,
            'gc_content': total_gc_content,
            'contigs': all_contig_data,
            'num_contigs': len(all_contig_data)
        }
        return assembly_data


    def fasta_filter_contigs_generator(self, fasta_record_iter, min_contig_length):
        ''' generates SeqRecords iterator for writing from a legacy contigset object '''
        rows = 0
        rows_added = 0
        for record in fasta_record_iter:
            rows += 1
            if len(record.seq) >= min_contig_length:
                rows_added += 1
                yield record
        print(' - filtered out ' + str(rows - rows_added) + ' of ' + str(rows) + ' contigs that were shorter than ' +
              str(min_contig_length) + 'bp.')


    def filter_contigs_by_length(self, fasta_file_path, min_contig_length):
        ''' removes all contigs less than the min_contig_length provided '''
        filtered_fasta_file_path = fasta_file_path + '.filtered.fa'

        fasta_record_iter = SeqIO.parse(fasta_file_path, 'fasta')
        SeqIO.write(self.fasta_filter_contigs_generator(fasta_record_iter, min_contig_length),
                    filtered_fasta_file_path, 'fasta')

        return filtered_fasta_file_path


    def save_assembly_object(self, workspace_id, assembly_name, obj_data):
        print('Saving Assembly to Workspace')
        sys.stdout.flush()
        obj_info = self.dfu.save_objects({'id': workspace_id,
                                          'objects': [{'type': 'KBaseGenomeAnnotations.Assembly',
                                                       'data': obj_data,
                                                       'name': assembly_name
                                                       }]
                                          })[0]
        return obj_info


    def save_fasta_file_to_shock(self, fasta_file_path):
        ''' Given the path to the file, upload to shock and return Handle information
            returns:
                typedef structure {
                    string shock_id;
                    Handle handle;
                    string node_file_name;
                    string size;
                } FileToShockOutput;

        '''
        print('Uploading fasta file (' + str(fasta_file_path) + ') to SHOCK')
        sys.stdout.flush()
        return self.dfu.file_to_shock({'file_path': fasta_file_path, 'make_handle': 1})


    def stage_input(self, params):
        ''' Setup the input_directory by fetching the files and returning the path to the file'''
        file_path = None
        if 'file' in params:
            file_path = os.path.abspath(params['file']['path'])
        elif 'shock_id' in params:
            print('Downloading file from SHOCK node: ' + str(params['shock_id']))
            sys.stdout.flush()
            input_directory = os.path.join(self.scratch, 'assembly-upload-staging-' + str(uuid.uuid4()))
            os.makedirs(input_directory)
            file_name = self.dfu.shock_to_file({'file_path': input_directory,
                                                'shock_id': params['shock_id']
                                                })['node_file_name']
            file_path = os.path.join(input_directory, file_name)
        elif 'ftp_url' in params:
            print('Downloading file from: ' + str(params['ftp_url']))
            sys.stdout.flush()
            file_path = self.dfu.download_web_file({'file_url': params['ftp_url'],
                                                    'download_type': 'FTP'
                                                    })['copy_file_path']

        # extract the file if it is compressed
        if file_path is not None:
            unpacked_file = self.dfu.unpack_file({'file_path': file_path})
            return unpacked_file['file_path']

        raise ValueError('No valid fasta could be extracted based on the input parameters')


    def validate_params(self, params):
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
                        raise ValueError('when specifying a fasta file input, "path" field was not defined in "file"')

        if input_count == 0:
            raise ValueError('required fasta file as input, set as either "file", "shock_id", or "ftp_url"')
        if input_count > 1:
            raise ValueError('required exactly one fasta file as input source, you set more than one of ' +
                             'these fields: "file", "shock_id", or "ftp_url"')
