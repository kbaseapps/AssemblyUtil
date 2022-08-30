import json
import os
import os.path
from pathlib import Path
import sys
import uuid
from typing import List
from collections import Counter
from hashlib import md5

from Bio import SeqIO

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


def _ref(object_info):
    return f'{object_info[6]}/{object_info[0]}/{object_info[4]}'

class FastaToAssembly:

    # Note added X due to kb|g.1886.fasta
    _VALID_CHARS = "-ACGTUWSMKRYBDHVNX"
    _AMINO_ACID_SPECIFIC_CHARACTERS = "PLIFQE"
    def __init__(self, dfu: DataFileUtil, ws: Workspace, scratch: Path):
        self._scratch = scratch
        self._dfu = dfu
        self._ws = ws

    def import_fasta(self, params):
        print('validating parameters')
        self._validate_single_params(params)
        inputs = dict(params)
        inputs.pop('workspace_name', None)
        inputs.pop('min_contig_length', None)
        inputs.pop('shock_id', None)
        mass_params = {
            'workspace_name': params['workspace_name'],
            'min_contig_length': params.get('min_contig_length'),
            'inputs': [inputs]
        }
        if 'file' in params:
            mass_params['inputs'][0]['file'] = params['file']['path']
        else:
            mass_params['inputs'][0]['node'] = params['shock_id']
        return self._import_fastas(mass_params)[0]

    def _import_fastas(self, params):
        # TODO TEST with muliple files
        # TODO check that inputs are all either files or shock nodes
        # TODO check inputs generally
        # TODO expose in API
        # For now this is completely serial, but theoretically we could start uploading
        # Blobstore nodes when some fraction of the initial checks are done, start uploading
        # Workspace obects when some fraction of the Blobstore nodes are done, parallelize
        # the file filtering / parsing, etc.
        # For now keep it simple
        # Also note all the assembly data is kept in memory once parsed, but it contains
        # no sequence info and so shouldn't be too huge. Could push to KBase in batches or
        # save to disk if that's an issue
        # We also probably want to add some retries for saving data, but that should go
        # in DataFileUtils if it's not there already
        # Finally, if more than 1G worth of assembly object data is sent to the workspace at once,
        # the call will fail. May need to add some checking / splitting code around this.
        if 'file' in params['inputs'][0]:
            input_files = self._stage_file_inputs(params['inputs'])
        else:
            input_files = self._stage_blobstore_inputs(params['inputs'])

        mcl = params.get('min_contig_length')
        mcl = int(mcl) if mcl else None  # TODO TEST no key and None
        assembly_data = []
        for i in range(len(input_files)):
            # Hmm, all through these printouts we should really put the blobstore node here as
            # well as the file if it exists... wait and see if that code path is still actually
            # used
            if mcl:
                print(f'filtering FASTA file {input_files[i]} by contig length '
                      + f'(min len={mcl} bp)')
                input_files[i] = self._filter_contigs_by_length(input_files[i], mcl)
            print(f'parsing FASTA file: {input_files[i]}')
            assdata = self._parse_fasta(
                input_files[i],
                params['inputs'][i].get('contig_info') or {})
            print(f' - parsed {assdata["num_contigs"]} contigs, {assdata["dna_size"]} bp')
            if not assdata["num_contigs"]:
                raise ValueError("Either the original FASTA file contained no sequences or they "
                                 + "were all filtered out based on the min_contig_length "
                                 + f"parameter for file {input_files[i]}")
            assembly_data.append(assdata)

        print('saving assemblies to KBase')
        file_handles = self._save_files_to_blobstore(input_files)
        assobjects = []
        for assdata, file_handle, inputs, sourcefile in zip(
                assembly_data, file_handles, params['inputs'], input_files):
            ao = self._build_assembly_object(assdata, file_handle, inputs)
            assobjects.append(ao)
            # this appears to be completely unused
            with open(sourcefile.parent / "example.json", "w") as f:
                json.dump(ao, f)

        # save to WS and return
        if 'workspace_id' in params:
            workspace_id = int(params['workspace_id'])
        else:
            workspace_id = self._dfu.ws_name_to_id(params['workspace_name'])
        assembly_infos = self._save_assembly_objects(
            workspace_id,
            [p['assembly_name'] for p in params['inputs']],
            assobjects
        )
        return [_ref(ai) for ai in assembly_infos]

    def _build_assembly_object(self, assembly_data, fasta_file_handle_info, params):
        """ construct the WS object data to save based on the parsed info and params """
        assembly_data['assembly_id'] = params['assembly_name']
        assembly_data['fasta_handle_ref'] = fasta_file_handle_info['handle']['hid']
        assembly_data['fasta_handle_info'] = fasta_file_handle_info

        assembly_data['type'] = 'Unknown'
        if 'type' in params:
            assembly_data['type'] = params['type']

        if 'taxon_ref' in params:
            info = self._ws.get_object_info3(
                {'objects':[{'ref': params['taxon_ref']}]})['infos'][0]
            assembly_data['taxon_ref'] = f'{info[6]}/{info[0]}/{info[4]}'

        if 'external_source' in params:
            assembly_data['external_source'] = params['external_source']

        if 'external_source_id' in params:
            assembly_data['external_source_id'] = params['external_source_id']

        if 'external_source_origination_date' in params:
            assembly_data['external_source_origination_date'] = params['external_source_origination_date']

        return assembly_data

    def _parse_fasta(self, fasta_file_path: Path, extra_contig_info):
        """ Do the actual work of inspecting each contig """

        # TODO TEST this needs more extensive unit testing
        # variables to store running counts of things
        total_length = 0
        base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        md5_list = []

        # map from contig_id to contig_info
        all_contig_data = {}

        for record in SeqIO.parse(str(fasta_file_path), "fasta"):
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
                if character not in self._VALID_CHARS:
                    if character in self._AMINO_ACID_SPECIFIC_CHARACTERS:
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
    def _fasta_filter_contigs_generator(fasta_record_iter, min_contig_length):
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

    def _filter_contigs_by_length(self, fasta_file_path: Path, min_contig_length) -> Path:
        """ removes all contigs less than the min_contig_length provided """
        filtered_fasta_file_path = Path(str(fasta_file_path) + '.filtered.fa')

        fasta_record_iter = SeqIO.parse(str(fasta_file_path), 'fasta')
        SeqIO.write(self._fasta_filter_contigs_generator(fasta_record_iter, min_contig_length),
                    str(filtered_fasta_file_path), 'fasta')

        return filtered_fasta_file_path

    def _save_assembly_objects(self, workspace_id, assembly_names, ass_data):
        print('Saving Assemblies to Workspace')
        sys.stdout.flush()
        ws_inputs = []
        for assname, assdata_singular in zip(assembly_names, ass_data):
            ws_inputs.append({
                'type': 'KBaseGenomeAnnotations.Assembly',  # This should really be versioned...
                'data': assdata_singular,
                'name': assname
            })
        return self._dfu.save_objects({'id': workspace_id, 'objects': ws_inputs})

    def _save_files_to_blobstore(self, files: List[Path]):
        print(f'Uploading FASTA files to the Blobstore')
        sys.stdout.flush()
        blob_input = [{'file_path': str(fp), 'make_handle': 1} for fp in files]
        return self._dfu.file_to_shock_mass(blob_input)

    def _stage_file_inputs(self, inputs) -> List[Path]:
        in_files = []
        for inp in inputs:
            if not os.path.isfile(inp['file']):
                raise ValueError(  # TODO TEST
                    "KBase Assembly Utils tried to save an assembly, but the calling "
                    + f"application specified a file ('{inp['file']}') that is missing. "
                    + "Please check the application logs for details.")
            # Ideally we'd have some sort of security check here but the DTN files could
            # be mounted anywhere...
            # TODO check with sysadmin about this - checked, waiting on clear list of safedirs
            fp = Path(inp['file']).resolve(strict=True)
            # make the downstream unpack call unpack into scratch rather than wherever the
            # source file might be
            file_path = self._create_temp_dir() / fp.name
            # symlink doesn't work, because in DFU filemagic doesn't follow symlinks, and so
            # DFU won't unpack symlinked files
            os.link(fp, file_path)
            in_files.append(file_path)
        # extract the file if it is compressed
        # could add a target dir argument to unpack_files, not sure how much work that might be
        fs = [{'file_path': str(fp), 'unpack': 'uncompress'} for fp in in_files]
        unpacked_files = self._dfu.unpack_files(fs)
        return [Path(uf['file_path']) for uf in unpacked_files]

    def _stage_blobstore_inputs(self, inputs) -> List[Path]:
        blob_params = []
        for inp in inputs:
            blob_params.append({
                'shock_id': inp['node'],
                'file_path': str(self._create_temp_dir()),
                'unpack': 'uncompress'  # Will throw an error for archives
            })
        dfu_res = self._dfu.shock_to_file_mass(blob_params)
        return [Path(dr['file_path']) for dr in dfu_res]

    def _create_temp_dir(self):
        tmpdir = self._scratch / ("import_fasta_" + str(uuid.uuid4()))
        os.makedirs(tmpdir)
        return tmpdir

    @staticmethod
    def _validate_single_params(params):
        for key in ('workspace_name', 'assembly_name'):
            if key not in params:
                raise ValueError('required "' + key + '" field was not defined')

        # one and only one of either 'file' or 'shock_id' is required
        input_count = 0
        for key in ('file', 'shock_id'):
            if key in params and params[key] is not None:
                input_count = input_count + 1
                if key == 'file':
                    if not isinstance(params[key], dict) or 'path' not in params[key]:
                        raise ValueError('when specifying a FASTA file input, "path" field was not defined in "file"')

        if input_count == 0:
            raise ValueError('required FASTA file as input, set as either "file" or "shock_id"')
        if input_count > 1:
            raise ValueError('required exactly one FASTA file as input source, you set more ' +
                             'than one of these fields: "file", "shock_id"')
