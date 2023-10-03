'''
Unit tests for FastaToAssembly.py.
Integration tests are in the server test file.
'''

import os
import uuid
from pathlib import Path
from typing import Callable, Optional
from unittest.mock import create_autospec

from AssemblyUtil.FastaToAssembly import FastaToAssembly, _validate_threads_params
from conftest import assert_exception_correct
from installed_clients.DataFileUtilClient import DataFileUtil
from pytest import raises

# TODO Add more unit tests when changing things until entire file is covered by unit tests


def _set_up_mocks(
        path: str = 'fake_scratch',
        uuid_gen: Optional[Callable[[], uuid.UUID]] = None
        ):
    dfu = create_autospec(DataFileUtil, spec_set=True, instance=True)
    if uuid_gen:
        fta = FastaToAssembly(dfu, Path(path), uuid_gen=uuid_gen)
    else:
        fta = FastaToAssembly(dfu, Path(path))
    return fta, dfu


def _run_test_spec_fail(test_spec, mass=False, print_spec=False):
    fta, _ = _set_up_mocks()
    for err, params in test_spec:
        if print_spec:
            print(f"spec:\n{params}")
        if mass:
            _run_test_mass_fail(fta, params, err)
        else:
            _run_test_fail(fta, params, err)


def _run_test_fail(fta, params, expected):
        with raises(Exception) as got:
            fta.import_fasta(params)
        assert_exception_correct(got.value, ValueError(expected))

    
def _run_test_mass_fail(fta, params, expected):
        with raises(Exception) as got:
            fta.import_fasta_mass(params)
        assert_exception_correct(got.value, ValueError(expected))


def test_import_fasta_workspace_name_id_input_fail():
    '''
    Tests the cases where workspace identifiers are submitted incorrectly.
    '''
    err1 = 'Exactly one of a workspace_id or a workspace_name must be provided'
    err2 = 'workspace_id must be an integer >= 1'
    test_spec = [
        (err1, {}),
        (err1, {'werksce_nme': 'foo'}),
        (err1, {'workspace_name': 'bar', 'workspace_id': 1}),
        (err2, {'workspace_id': 0}),
        (err2, {'workspace_id': -1}),
        (err2, {'workspace_id': -100}),
        ('workspace_id must be an integer, got: foo', {'workspace_id': 'foo'}),
    ]
    _run_test_spec_fail(test_spec)


def test_import_fasta_no_assembly_name():
    fta, dfu = _set_up_mocks()
    dfu.ws_name_to_id.return_value = 34

    _run_test_fail(
        fta, {'workspace_name': 'whee'}, 'Required parameter assembly_name was not defined')

    dfu.ws_name_to_id.assert_called_once_with('whee')


def _update(d1, d2):
    d = dict(d1)
    d.update(d2)
    return d

def test_import_fasta_file_source_fail():
    '''
    Tests the case where file source (e.g. file path, shock ID) info is submitted incorrectly.
    '''
    err1 = 'Exactly one of file or shock_id is required'
    err2 = 'When specifying a FASTA file input, "path" field was not defined in "file"'
    b = {'workspace_id': 1, 'assembly_name': 'foo'}
    test_spec = [
        (err1, b),
        (err1, _update(b, {'filet': {'path': '/foo/bar'}})),
        (err1, _update(b, {'file':  {'path': '/foo/bar'}, 'shock_id': 'fakeid'})),
        (err2, _update(b, {'file': '/foo/bar'})),
        (err2, _update(b, {'file': {'perth': '/foo/bar'}})),
    ]
    _run_test_spec_fail(test_spec)


def test_import_fasta_min_contig_length_fail():
    '''
    Tests illegal min_contig_length values.
    '''
    err1 = 'If provided, min_contig_length must be an integer >= 1'
    b = {'workspace_id': 1, 'assembly_name': 'foo', 'shock_id': 'fake_id'}
    test_spec = [
        ('If provided, min_contig_length must be an integer, got: foo',
          _update(b, {'min_contig_length': 'foo'})),
        (err1, _update(b, {'min_contig_length': 0})),
        (err1, _update(b, {'min_contig_length': -1})),
        (err1, _update(b, {'min_contig_length': -100})),
        ('If provided, min_contig_length must be an integer, got: []',
          _update(b, {'min_contig_length': []})),
    ]
    _run_test_spec_fail(test_spec)


# https://docs.pytest.org/en/7.1.x/reference/reference.html#tmp-path
def test_import_fasta_mass_basic_file_no_mcl_key(tmp_path):
    '''
    Test of the mass importer with file inputs and no min_contig_length key.
    '''
    _test_import_fasta_mass_file(tmp_path, {})


def test_import_fasta_mass_basic_file_None_mcl_key(tmp_path):
    '''
    Test of the mass importer with file inputs and a min_contig_length key value of None.
    '''
    _test_import_fasta_mass_file(tmp_path, {'min_contig_length': None})


def _test_import_fasta_mass_file(tmp_path, params_root):
    ### Set up input files ###
    with open(tmp_path / 'f1.fasta', 'w') as f1:
        f1.writelines([
            '> contig1\n',
            'ANNTGGCC\n',
            '> contig2\n',
            'CCGNTTA\n',
        ])
    with open(tmp_path / 'f2.fasta', 'w') as f2:
        f2.writelines([
            '> contig3\n',
            'AATTGGC\n',
            '> contig4\n',
            'CCGGTAA\n',
        ])

    ### Set up FastaToAssembly ###
    scratch = tmp_path / 'scratch'
    uuid1 = uuid.uuid4()
    uuid2 = uuid.uuid4()
    dir1 = scratch / ('import_fasta_' + str(uuid1))
    dir2 = scratch / ('import_fasta_' + str(uuid2))
    fta, dfu = _set_up_mocks(
        path=str(scratch),
        uuid_gen=lambda i=iter([uuid1, uuid2]): next(i)
    )

    ### Set up mock return values ###
    dfu.unpack_files.return_value = [
        {'file_path': str(dir1 / 'f1.fasta')},
        {'file_path': str(dir2 / 'f2.fasta')}
    ]
    dfu.file_to_shock_mass.return_value = [
        {
            'shock_id': 'fake_id1',
            'handle': {
                'hid': 'KBH_24',
                'file_name': 'f1.fasta',
                'id': 'fake_id1',
                'url': 'https://kbase.us/services/shock-api',
                'type': 'shock',
                'remote_md5': 'fake_md5_1',
            },
            'node_file_name': 'f1.fasta',
            'size': 78
        },
        {
            'shock_id': 'fake_id2',
            'handle': {
                'hid': 'KBH_26',
                'file_name': 'f2.fasta',
                'id': 'fake_id2',
                'url': 'https://kbase.us/services/shock-api',
                'type': 'shock',
                'remote_md5': 'fake_md5_2',
            },
            'node_file_name': 'f2.fasta',
            'size': 90
        },
    ]
    dfu.save_objects.return_value = [
        [4, 'name', 'type', 'time', 75, 'user', 42, 'wsname', 'md5', 78, {}],
        [7, 'name', 'type', 'time', 1, 'user', 42, 'wsname', 'md5', 78, {}],
    ]

    ### Call the function ###
    params_root.update({
        'workspace_id': 42,
        'inputs': [
            {'file': str(tmp_path / 'f1.fasta'), 'assembly_name': 'foo1'},
            {
                'file': str(tmp_path / 'f2.fasta'),
                'assembly_name': 'foo2',
                'type': 'There should be a controlled vocabulary for this field',
                'external_source': 'ext source',
                'external_source_id': 'ext source id',
                'taxon_ref': 'MyWS/MyObj/3',  # Silently ignored, removed in 2.0.0
                'contig_info': {
                    'contig3': {'is_circ': 1, 'description': 'desc goes here'},
                    'contig4': {'is_circ': -2}
                }
            }
        ]
    })
    res = fta.import_fasta_mass(params_root, parallelize=False)
    assert res == [
        {'upa': '42/4/75', 'filtered_input': None},
        {'upa': '42/7/1', 'filtered_input': None}
    ]

    ### Check mock calls ###
    dfu.unpack_files.assert_called_once_with([
        {'file_path': str(dir1 / 'f1.fasta'), 'unpack': 'uncompress'},
        {'file_path': str(dir2 / 'f2.fasta'), 'unpack': 'uncompress'}
    ])    
    dfu.file_to_shock_mass.assert_called_once_with([
        {'file_path': str(dir1 / 'f1.fasta'), 'make_handle': 1},
        {'file_path': str(dir2 / 'f2.fasta'), 'make_handle': 1}
    ])
    dfu.save_objects.assert_called_once_with({
        'id': 42,
        'objects': [
            {
                'type': 'KBaseGenomeAnnotations.Assembly',
                'data': {
                    'md5': '1e3ac553a7d80d393859e867785cf164',
                    'base_counts': {'A': 2, 'G': 3, 'C': 4, 'T': 3, 'N': 3},
                    'dna_size': 15,
                    'gc_content': 0.46667,
                    'contigs': {
                        'contig1': {
                            'contig_id': 'contig1',
                            'name': 'contig1',
                            'description': '1',
                            'length': 8,
                            'Ncount': 2,
                            'md5': '95cc0920e348e08fdc8c0c36b27a8f25',
                            'gc_content': 0.5
                        },
                        'contig2': {
                            'contig_id': 'contig2',
                            'name': 'contig2',
                            'description': '2',
                            'length': 7,
                            'Ncount': 1,
                            'md5': '2b21f05114f53168aa9cef384ea24ce7',
                            'gc_content': 0.42857
                        }
                    },
                    'num_contigs': 2,
                    'assembly_id': 'foo1',
                    'fasta_handle_ref': 'KBH_24',
                    'fasta_handle_info': {
                        'shock_id': 'fake_id1',
                        'handle': {
                            'hid': 'KBH_24',
                            'file_name': 'f1.fasta',
                            'id': 'fake_id1',
                            'url': 'https://kbase.us/services/shock-api',
                            'type': 'shock',
                            'remote_md5': 'fake_md5_1'
                        },
                        'node_file_name': 'f1.fasta',
                        'size': 78},
                    'type': 'Unknown'
                },
                'name': 'foo1'
            },
            {
                'type': 'KBaseGenomeAnnotations.Assembly',
                'data': {
                    'md5': 'e988c591703e3fa68715af6ba56d11cf',
                    'base_counts': {'A': 4, 'G': 4, 'C': 3, 'T': 3},
                    'dna_size': 14,
                    'gc_content': 0.5,
                    'contigs': {
                        'contig3': {
                            'contig_id': 'contig3',
                            'name': 'contig3',
                            'description': 'desc goes here',
                            'length': 7,
                            'is_circ': 1,
                            'md5': 'f2111a004b33ba29e244df5fb2ed6ca8',
                            'gc_content': 0.42857
                        },
                        'contig4': {
                            'contig_id': 'contig4',
                            'name': 'contig4',
                            'description': '4',
                            'length': 7,
                            'is_circ': -2,
                            'md5': '6b7ab60e09e6af28b66a988f3ee3c807',
                            'gc_content': 0.57143
                        }
                    },
                    'num_contigs': 2,
                    'assembly_id': 'foo2',
                    'fasta_handle_ref': 'KBH_26',
                    'fasta_handle_info': {
                        'shock_id': 'fake_id2',
                        'handle': {
                            'hid': 'KBH_26',
                            'file_name': 'f2.fasta',
                            'id': 'fake_id2',
                            'url': 'https://kbase.us/services/shock-api',
                            'type': 'shock',
                            'remote_md5': 'fake_md5_2'
                        },
                        'node_file_name': 'f2.fasta',
                        'size': 90
                    },
                    'type': 'There should be a controlled vocabulary for this field',
                    'external_source': 'ext source',
                    'external_source_id': 'ext source id',
                },
                'name': 'foo2'
            }
        ]
    })


def test_import_fasta_mass_blobstore_min_contig_length(tmp_path):
    '''
    Test mass FASTA import with files sourced from the Blobstore and a > 0 min_contig_length
    parameter.
    '''
    ### Set up blobstore files ###
    scratch = tmp_path / 'scratch'
    uuid1 = uuid.uuid4()
    uuid2 = uuid.uuid4()
    dir1 = scratch / ('import_fasta_' + str(uuid1))
    dir2 = scratch / ('import_fasta_' + str(uuid2))
    os.makedirs(dir1)
    os.makedirs(dir2)
    file1 = dir1 / 'f1.blobstore.fasta'
    file2 = dir2 / 'f2.blobstore.fasta'
    with open(file1, 'w') as f1:
        f1.writelines([
            '> contig1\n',
            'AATTGGCC\n',
            '> contig2\n',
            'CCGNTTA\n',
            '> expect removal\n',
            'AAAAAA\n',
            '> and here\n',
            'A\n'
        ])
    with open(file2, 'w') as f2:
        f2.writelines([
            '> contig3\n',
            'AATTGGT\n',
            '> contig4\n',
            'CCGGTAA\n',
        ])

    ### Set up FastaToAssembly ###
    fta, dfu = _set_up_mocks(
        path=str(scratch),
        uuid_gen=lambda i=iter([uuid1, uuid2]): next(i)
    )

    ### Set up mock return values ###
    dfu.shock_to_file_mass.return_value = [{'file_path' : str(file1)}, {'file_path': str(file2)}]
    dfu.file_to_shock_mass.return_value = [
        {
            'shock_id': 'fake_id_1',
            'handle': {
                'hid': 'KBH_24',
                'file_name': 'f1.blobstore.fasta.filtered.fa',
                'id': 'fake_id_1',
                'url': 'https://kbase.us/services/shock-api',
                'type': 'shock',
                'remote_md5': 'fake_md5_1',
            },
            'node_file_name': 'f1.blobstore.fasta.filtered.fa',
            'size': 78
        },
        {
            'shock_id': 'fake_id_2',
            'handle': {
                'hid': 'KBH_26',
                'file_name': 'f2.blobstore.fasta.filtered.fa',
                'id': 'fake_id_2',
                'url': 'https://kbase.us/services/shock-api',
                'type': 'shock',
                'remote_md5': 'fake_md5_2',
            },
            'node_file_name': 'f2.blobstore.fasta.filtered.fa',
            'size': 90
        },
    ]
    dfu.save_objects.return_value = [
        [3, 'name', 'type', 'time', 75, 'user', 42, 'wsname', 'md5', 78, {}],
        [6, 'name', 'type', 'time', 1, 'user', 42, 'wsname', 'md5', 78, {}],
    ]

    ### Call the function ###
    res = fta.import_fasta_mass({
        'workspace_id': 42,
        'min_contig_length': 7,
        'inputs': [
            {'node': 'fake_id_1', 'assembly_name': 'foo1'},
            {'node': 'fake_id_2', 'assembly_name': 'foo2'}
        ]
    }, parallelize=False)
    assert res == [
        {'upa': '42/3/75', 'filtered_input': str(dir1 / 'f1.blobstore.fasta.filtered.fa')},
        {'upa': '42/6/1', 'filtered_input': str(dir2 / 'f2.blobstore.fasta.filtered.fa')}
    ]

    ### Check mock calls ###
    dfu.shock_to_file_mass.assert_called_once_with([
        {
            'shock_id': 'fake_id_1',
            'file_path': str(dir1),
            'unpack': 'uncompress'
        },
        {
            'shock_id': 'fake_id_2',
            'file_path': str(dir2),
            'unpack': 'uncompress'
        },
    ])
    dfu.file_to_shock_mass.assert_called_once_with([
        {'file_path': str(dir1 / 'f1.blobstore.fasta.filtered.fa'), 'make_handle': 1},
        {'file_path': str(dir2 / 'f2.blobstore.fasta.filtered.fa'), 'make_handle': 1}
    ])
    dfu.save_objects.assert_called_once_with({
        'id': 42,
        'objects': [
            {
                'type': 'KBaseGenomeAnnotations.Assembly',
                'data': {
                    'md5': 'a3445818c6c795c8d17d38ce1851b882',
                    'base_counts': {'A': 3, 'G': 3, 'C': 4, 'T': 4, 'N': 1},
                    'dna_size': 15,
                    'gc_content': 0.46667,
                    'contigs': {
                        'contig1': {
                            'contig_id': 'contig1',
                            'name': 'contig1',
                            'description': '1',
                            'length': 8,
                            'md5': 'bfa6af6c781dccc1453e2a52074de5c7',
                            'gc_content': 0.5
                        },
                        'contig2': {
                            'contig_id': 'contig2',
                            'name': 'contig2',
                            'description': '2',
                            'length': 7,
                            'Ncount': 1,
                            'md5': '2b21f05114f53168aa9cef384ea24ce7',
                            'gc_content': 0.42857
                        }
                    },
                    'num_contigs': 2,
                    'assembly_id': 'foo1',
                    'fasta_handle_ref': 'KBH_24',
                    'fasta_handle_info': {
                        'shock_id': 'fake_id_1',
                        'handle': {
                            'hid': 'KBH_24',
                            'file_name': 'f1.blobstore.fasta.filtered.fa',
                            'id': 'fake_id_1',
                            'url': 'https://kbase.us/services/shock-api',
                            'type': 'shock',
                            'remote_md5': 'fake_md5_1'
                        },
                        'node_file_name': 'f1.blobstore.fasta.filtered.fa',
                        'size': 78},
                    'type': 'Unknown'
                },
                'name': 'foo1'
            },
            {
                'type': 'KBaseGenomeAnnotations.Assembly',
                'data': {
                    'md5': 'e42898241f2337ecbd4fc71abe0863f3',
                    'base_counts': {'A': 4, 'G': 4, 'C': 2, 'T': 4},
                    'dna_size': 14,
                    'gc_content': 0.42857,
                    'contigs': {
                        'contig3': {
                            'contig_id': 'contig3',
                            'name': 'contig3',
                            'description': '3',
                            'length': 7,
                            'md5': '40693dc60253d4f100e2a8686620b907',
                            'gc_content': 0.28571
                        },
                        'contig4': {
                            'contig_id': 'contig4',
                            'name': 'contig4',
                            'description': '4',
                            'length': 7,
                            'md5': '6b7ab60e09e6af28b66a988f3ee3c807',
                            'gc_content': 0.57143
                        }
                    },
                    'num_contigs': 2,
                    'assembly_id': 'foo2',
                    'fasta_handle_ref': 'KBH_26',
                    'fasta_handle_info': {
                        'shock_id': 'fake_id_2',
                        'handle': {
                            'hid': 'KBH_26',
                            'file_name': 'f2.blobstore.fasta.filtered.fa',
                            'id': 'fake_id_2',
                            'url': 'https://kbase.us/services/shock-api',
                            'type': 'shock',
                            'remote_md5': 'fake_md5_2'
                        },
                        'node_file_name': 'f2.blobstore.fasta.filtered.fa',
                        'size': 90
                    },
                    'type': 'Unknown',
                },
                'name': 'foo2'
            }
        ]
    })


def test_import_fasta_mass_fail_workspace_id_input():
    '''
    Tests the cases where workspace identifiers are submitted incorrectly.
    '''
    err1 = 'workspace_id is required'
    err2 = 'workspace_id must be an integer >= 1'
    test_spec = [
        (err1, {}),
        (err1, {'werksce_nme': 'foo'}),
        (err1, {'workspace_name': 'bar', 'worcspace_id': 1}),
        (err2, {'workspace_id': 0}),
        (err2, {'workspace_id': -1}),
        (err2, {'workspace_id': -100}),
        ('workspace_id must be an integer, got: foo', {'workspace_id': 'foo'}),
    ]
    _run_test_spec_fail(test_spec, mass=True)


def test_import_fasta_mass_fail_bad_inputs_field():
    err1 = 'inputs field is required and must be a non-empty list'
    err2 = 'Entry #3 in inputs field is not a mapping as required'
    test_spec = [
        (err1, {'workspace_id': 3}),
        (err1, {'workspace_id': 3, 'inpt': [{}]}),
        (err1, {'workspace_id': 3, 'inputs': []}),
        (err2, {'workspace_id': 3, 'inputs': [{}, {}, [], {}]}),
    ]
    _run_test_spec_fail(test_spec, mass=True)


def test_import_fasta_mass_fail_mixed_input_types():
    err1 = 'Entry #1 in inputs field must have exactly one of file or node specified'
    test_spec = [
        (err1, {'workspace_id': 3, 'inputs': [{}, {'node': 'a'}]}),
        (err1, {'workspace_id': 3, 'inputs': [{'file': 'a', 'node': 'b'}, {'node': 'a'}]}),
        ('Entry #3 in inputs must have a node field to match entry #1',
         {'workspace_id': 3, 'inputs': [
            {'node': 'b', 'assembly_name': 'x'},
            {'node': 'a', 'assembly_name': 'x'},
            {'file': 'b', 'assembly_name': 'x'}
        ]}),
        ('Entry #2 in inputs must have a file field to match entry #1',
         {'workspace_id': 3, 'inputs': [
            {'file': 'b', 'assembly_name': 'x'},
            {'node': 'a', 'assembly_name': 'x'},
            {'file': 'b', 'assembly_name': 'x'}
        ]}),
    ]
    _run_test_spec_fail(test_spec, mass=True)


def test_import_fasta_mass_fail_missing_assembly_name():
    test_spec = [
        ('Missing assembly_name field in inputs entry #1',
         {'workspace_id': 3, 'inputs': [
            {'node': 'b',},
            {'node': 'a', 'assembly_name': 'x'},
            {'node': 'b', 'assembly_name': 'x'}
        ]}),
        ('Missing assembly_name field in inputs entry #2',
         {'workspace_id': 3, 'inputs': [
            {'file': 'b', 'assembly_name': 'x'},
            {'file': 'a'},
            {'file': 'b', 'assembly_name': 'x'}
        ]}),
    ]
    _run_test_spec_fail(test_spec, mass=True)


def test_import_fasta_mass_fail_min_contig_length():
    err1 = 'If provided, min_contig_length must be an integer >= 2'
    b = {'workspace_id': 1, 'inputs': [{'assembly_name': 'foo', 'node': 'fake_id'}]}
    test_spec = [
        ('If provided, min_contig_length must be an integer, got: foo',
          _update(b, {'min_contig_length': 'foo'})),
        (err1, _update(b, {'min_contig_length': 1})),
        (err1, _update(b, {'min_contig_length': 0})),
        (err1, _update(b, {'min_contig_length': -1})),
        (err1, _update(b, {'min_contig_length': -100})),
        ('If provided, min_contig_length must be an integer, got: {}',
          _update(b, {'min_contig_length': {}})),
    ]
    _run_test_spec_fail(test_spec, mass=True)


def test_assembly_objects_generator_with_overflow():
    fta, _ = _set_up_mocks()
    # each serialized assembly object size is roughtly 180 below
    assembly_objects = [
        {'md5': '1e007bad0811a6d6e09a882d3bf802ab', 'base_counts': {'A': 1200415, 'G': 852652, 'C': 846697, 'T': 1191689, 'N': 847}, 'dna_size': 4092300, 'gc_content': 0.41526},
        {'md5': 'a26f200923f8c860f86a8d728055fd02', 'base_counts': {'A': 1184671, 'G': 847148, 'C': 850586, 'T': 1196308, 'N': 218}, 'dna_size': 4078931, 'gc_content': 0.41622},
        {'md5': '8aa6b1244e18c4f93bb3307902bd3a4d', 'base_counts': {'A': 1199516, 'G': 850011, 'C': 845769, 'T': 1183802, 'N': 106}, 'dna_size': 4079204, 'gc_content': 0.41571}
    ]
    assembly_names = ["name_1", "name_2", "name_3"]
    res = list(fta._assembly_objects_generator(assembly_objects, assembly_names, 400))
    assert len(res) == 2
    assert res[0][1] == ["name_1", "name_2"]
    assert res[1][1] == ["name_3"]


def test_assembly_objects_generator_without_overflow():
    fta, _ = _set_up_mocks()
    # each serialized assembly object size is roughtly 180 below
    assembly_objects = [
        {'md5': '1e007bad0811a6d6e09a882d3bf802ab', 'base_counts': {'A': 1200415, 'G': 852652, 'C': 846697, 'T': 1191689, 'N': 847}, 'dna_size': 4092300, 'gc_content': 0.41526},
        {'md5': 'a26f200923f8c860f86a8d728055fd02', 'base_counts': {'A': 1184671, 'G': 847148, 'C': 850586, 'T': 1196308, 'N': 218}, 'dna_size': 4078931, 'gc_content': 0.41622},
        {'md5': '8aa6b1244e18c4f93bb3307902bd3a4d', 'base_counts': {'A': 1199516, 'G': 850011, 'C': 845769, 'T': 1183802, 'N': 106}, 'dna_size': 4079204, 'gc_content': 0.41571}
    ]
    assembly_names = ["name_1", "name_2", "name_3"]
    res = list(fta._assembly_objects_generator(assembly_objects, assembly_names, 600))
    assert len(res) == 1
    assert res[0][1] == assembly_names


def test_validate_threads_params():
    var_name = "MAX_THREADS"
    with raises(ValueError, match=f"{var_name} is required"):
        _validate_threads_params(None, var_name)

    with raises(ValueError, match=f"{var_name} must be an integer"):
        _validate_threads_params("abc", var_name)

    with raises(ValueError, match=f"{var_name} must be >= 1"):
        _validate_threads_params("0", var_name)
