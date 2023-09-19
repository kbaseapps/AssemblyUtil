import itertools
import json
import math
import os
import sys
import uuid
from collections import Counter
from hashlib import md5
from pathlib import Path
from typing import Callable, List

from Bio import SeqIO
from installed_clients.baseclient import _JSONObjectEncoder
from installed_clients.DataFileUtilClient import DataFileUtil
from pathos.multiprocessing import ProcessingPool as Pool

_MAX_DATA_SIZE = 1024 * 1024 * 1024 # 1 GB
_SAFETY_FACTOR = 0.95
_SYSTEM_UTILIZATION = 1

_WSID = 'workspace_id'
_MCL = 'min_contig_length'
_INPUTS = 'inputs'
_FILE = 'file'
_NODE = 'node'
_ASSEMBLY_NAME = 'assembly_name'


def _upa(object_info):
    return f'{object_info[6]}/{object_info[0]}/{object_info[4]}'

def _get_serialized_object_size(assembly_object):
    arg_hash = {'params': assembly_object}
    serialized = json.dumps(arg_hash, cls=_JSONObjectEncoder)
    return sys.getsizeof(serialized)

class FastaToAssembly:

    # Note added X due to kb|g.1886.fasta
    _VALID_CHARS = "-ACGTUWSMKRYBDHVNX"
    _AMINO_ACID_SPECIFIC_CHARACTERS = "PLIFQE"
    def __init__(self,
             dfu: DataFileUtil,
             scratch: Path,
             uuid_gen: Callable[[], uuid.UUID] = lambda: uuid.uuid4()):
        self._scratch = scratch
        self._dfu = dfu
        self._uuid_gen = uuid_gen

    def import_fasta(self, params):
        print('validating parameters')
        mass_params = self._set_up_single_params(params)
        return self._import_fasta_mass(mass_params)[0]

    def import_fasta_mass(self, params, parallelize=True):
        print('validating parameters')
        self._validate_mass_params(params)
        if not parallelize or len(params[_INPUTS]) == 1:
            return self._import_fasta_mass(params)
        return self._run_parallel_import_fasta_mass(params)

    def _import_fasta_mass(self, params):
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
        if _FILE in params[_INPUTS][0]:
            input_files = self._stage_file_inputs(params[_INPUTS])
        else:
            input_files = self._stage_blobstore_inputs(params[_INPUTS])

        mcl = params.get(_MCL)
        assembly_data = []
        output = []
        for i in range(len(input_files)):
            # Hmm, all through these printouts we should really put the blobstore node here as
            # well as the file if it exists... wait and see if that code path is still actually
            # used
            if mcl:
                print(f'filtering FASTA file {input_files[i]} by contig length '
                      + f'(min len={mcl} bp)')
                input_files[i] = self._filter_contigs_by_length(input_files[i], mcl)
            output.append({'filtered_input': str(input_files[i]) if mcl else None})
            print(f'parsing FASTA file: {input_files[i]}')
            assdata = self._parse_fasta(
                input_files[i],
                params[_INPUTS][i].get('contig_info') or {})
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
                assembly_data, file_handles, params[_INPUTS], input_files):
            ao = self._build_assembly_object(assdata, file_handle, inputs)
            assobjects.append(ao)
            # this appears to be completely unused
            with open(sourcefile.parent / "example.json", "w") as f:
                json.dump(ao, f)

        # save to WS and return
        assembly_infos = []
        assembly_names = [p[_ASSEMBLY_NAME] for p in params[_INPUTS]]
        # generator to process files in bach. Avoid memory issues
        for assobject_batch, assembly_name_batch in self._assembly_objects_generator(
            assobjects, assembly_names
        ):
            assembly_infos.extend(
                self._save_assembly_objects(params[_WSID], assembly_name_batch, assobject_batch)
            )

        for out, ai in zip(output, assembly_infos):
            out['upa'] = _upa(ai)
        return output

    def _run_parallel_import_fasta_mass(self, params):
        threads = max(int(os.cpu_count() * min(_SYSTEM_UTILIZATION, 1)), 1)
        print(f' - running {threads} parallel threads')

        # distribute inputs evenly across threads
        param_inputs = params.pop(_INPUTS)
        chunk_size = math.ceil(len(param_inputs) / threads)
        batch_input = [
            (
                {
                    **params,
                    _INPUTS: param_inputs[i : i + chunk_size]
                }
            )
            for i in range(0, len(param_inputs), chunk_size)
        ]
        batch_result = Pool(threads).map(self._import_fasta_mass, batch_input)
        result = list(itertools.chain.from_iterable(batch_result))
        return result

    def _build_assembly_object(self, assembly_data, fasta_file_handle_info, params):
        """ construct the WS object data to save based on the parsed info and params """
        assembly_data['assembly_id'] = params[_ASSEMBLY_NAME]
        assembly_data['fasta_handle_ref'] = fasta_file_handle_info['handle']['hid']
        assembly_data['fasta_handle_info'] = fasta_file_handle_info

        assembly_data['type'] = 'Unknown'
        if 'type' in params:
            assembly_data['type'] = params['type']

        if 'external_source' in params:
            assembly_data['external_source'] = params['external_source']

        if 'external_source_id' in params:
            assembly_data['external_source_id'] = params['external_source_id']

        if 'external_source_origination_date' in params:
            # TODO this is an arbitrary string, which isn't useful. If this field is actually
            # used, make a new field with a standard timestamp format (epoch date?), validate that
            # format, and deprecate this field
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
            # TODO should throw an error if ECI has invalid record IDs
            if record.id in extra_contig_info:
                if 'is_circ' in extra_contig_info[record.id]:
                    # TODO supposed to be a boolean, should check for 1 or 0
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
    def _assembly_objects_generator(assembly_objects, assembly_names):
        """ generates assembly objects iterator for uploading to the target workspace"""
        start_idx = 0
        cumsize = 0
        # precalculate max_cumsize here to avoid memory issue
        max_cumsize = _MAX_DATA_SIZE * _SAFETY_FACTOR
        for idx, ao in enumerate(assembly_objects):
            aosize = _get_serialized_object_size(ao)
            if aosize + cumsize <= max_cumsize:
                cumsize += aosize
            else:
                yield assembly_objects[start_idx:idx], assembly_names[start_idx:idx]
                start_idx = idx
                cumsize = aosize
        # yield the last batch
        yield assembly_objects[start_idx:], assembly_names[start_idx:]

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
            if not os.path.isfile(inp[_FILE]):
                raise ValueError(
                    "KBase Assembly Utils tried to save an assembly, but the calling "
                    + f"application specified a file ('{inp[_FILE]}') that is missing. "
                    + "Please check the application logs for details.")
            # Ideally we'd have some sort of security check here but the DTN files could
            # be mounted anywhere...
            # TODO check with sysadmin about this - checked, waiting on clear list of safedirs
            fp = Path(inp[_FILE]).resolve(strict=True)
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
                'shock_id': inp[_NODE],
                'file_path': str(self._create_temp_dir()),
                'unpack': 'uncompress'  # Will throw an error for archives
            })
        dfu_res = self._dfu.shock_to_file_mass(blob_params)
        return [Path(dr['file_path']) for dr in dfu_res]

    def _create_temp_dir(self):
        tmpdir = self._scratch / ("import_fasta_" + str(self._uuid_gen()))
        os.makedirs(tmpdir, exist_ok=True)
        return tmpdir

    def _set_up_single_params(self, params):
        inputs = dict(params)
        ws_id = self._get_int(inputs.pop(_WSID, None), _WSID)
        ws_name = inputs.pop('workspace_name', None)
        if (bool(ws_id) == bool(ws_name)):  # xnor
            raise ValueError(f"Exactly one of a {_WSID} or a workspace_name must be provided")
        if not ws_id:
            print(f"Translating workspace name {ws_name} to a workspace ID. Prefer submitting "
                  + "a workspace ID over a mutable workspace name that may cause race conditions")
            ws_id = self._dfu.ws_name_to_id(params['workspace_name'])

        if not inputs.get(_ASSEMBLY_NAME):
            raise ValueError(f"Required parameter {_ASSEMBLY_NAME} was not defined")

        # one and only one of either 'file' or 'shock_id' is required
        file_ = inputs.pop(_FILE, None)
        shock_id = inputs.pop('shock_id', None)
        if (bool(file_) == bool(shock_id)):  # xnor
            raise ValueError(f"Exactly one of {_FILE} or shock_id is required")
        if file_:
            if not isinstance(file_, dict) or 'path' not in file_:
                raise ValueError('When specifying a FASTA file input, "path" field was '
                                 + f'not defined in "{_FILE}"')
        mass_params = {
            _WSID: ws_id,
            # Ideally set of minimum of 2 here, but left at 1 for backwards compatibility
            _MCL: self._get_int(inputs.pop(_MCL, None), f"If provided, {_MCL}"),
            _INPUTS: [inputs]
        }
        if file_:
            inputs[_FILE] = params[_FILE]['path']
        else:
            inputs[_NODE] = params['shock_id']
        return mass_params

    def _validate_mass_params(self, params):
        ws_id = self._get_int(params.get(_WSID), _WSID)
        if not ws_id:
            raise ValueError(f"{_WSID} is required")
        inputs = params.get(_INPUTS)
        if not inputs or type(inputs) != list:
            raise ValueError(f"{_INPUTS} field is required and must be a non-empty list")
        for i, inp in enumerate(inputs, start=1):
            if type(inp) != dict:
                raise ValueError(f"Entry #{i} in {_INPUTS} field is not a mapping as required")
        file_ = inputs[0].get(_FILE)
        if bool(file_) == bool(inputs[0].get(_NODE)):  # xnor
            raise ValueError(f"Entry #1 in {_INPUTS} field must have exactly one of "
                             + f"{_FILE} or {_NODE} specified")
        field = _FILE if file_ else _NODE
        for i, inp in enumerate(inputs, start=1):
            if not inp.get(field):
                raise ValueError(
                    f"Entry #{i} in {_INPUTS} must have a {field} field to match entry #1")
            if not inp.get(_ASSEMBLY_NAME):
                raise ValueError(f"Missing {_ASSEMBLY_NAME} field in {_INPUTS} entry #{i}")
        self._get_int(params.get(_MCL), f"If provided, {_MCL}", minimum=2)

    def _get_int(self, putative_int, name, minimum=1):
        if putative_int is not None:
            if type(putative_int) != int:
                raise ValueError(f"{name} must be an integer, got: {putative_int}")
            if putative_int < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        return putative_int
