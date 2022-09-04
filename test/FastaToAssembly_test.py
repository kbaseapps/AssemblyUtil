'''
Unit tests for FastaToAssembly.py.
Integration tests are in the server test file.
'''

from pathlib import Path
from pytest import raises
from unittest.mock import create_autospec

from AssemblyUtil.FastaToAssembly import FastaToAssembly
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from testing_helpers import assert_exception_correct


# TODO Add more unit tests when changing things until entire file is covered by unit tests


def _set_up_mocks(path: str = 'fake_scratch'):
    dfu = create_autospec(DataFileUtil, spec_set=True, instance=True)
    ws = create_autospec(Workspace, spec_set=True, instance=True)
    fta = FastaToAssembly(dfu, ws, Path(path))
    return fta, dfu, ws


def _run_test_spec(test_spec, print_spec=False):
    fta, _, _ = _set_up_mocks()
    for err, params in test_spec:
        if print_spec:
            print(f"spec:\n{params}")
        with raises(Exception) as got:
            fta.import_fasta(params)
        assert_exception_correct(got.value, ValueError(err))


def test_import_fasta_workspace_name_id_input_fail():
    '''
    Tests the cases where workspace identifiers are submitted incorrectly.
    '''
    err1 = 'Exactly one of a workspace_id or a workspace_name must be provided'
    err2 = 'workspace_id must be an integer > 0'
    test_spec = [
        (err1, {}),
        (err1, {'werksce_nme': 'foo'}),
        (err1, {'workspace_name': 'bar', 'workspace_id': 1}),
        (err2, {'workspace_id': 0}),
        (err2, {'workspace_id': -1}),
        (err2, {'workspace_id': -100}),
        ('workspace_id must be an integer, got: foo', {'workspace_id': 'foo'}),
    ]
    _run_test_spec(test_spec)


def test_import_fasta_no_assembly_name():
    fta, dfu, _ = _set_up_mocks()
    dfu.ws_name_to_id.return_value = 34

    with raises(Exception) as got:
        fta.import_fasta({'workspace_name': 'whee'})
    assert_exception_correct(got.value, ValueError(
        'Required parameter assembly_name was not defined'))

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
    _run_test_spec(test_spec)


def test_import_fasta_min_contig_length_fail():
    '''
    Tests illegal min_contig_length values.
    '''
    err1 = 'If provided, min_contig_length must be an integer > 0'
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
    _run_test_spec(test_spec)
