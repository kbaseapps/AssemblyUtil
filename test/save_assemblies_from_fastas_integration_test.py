'''
Minimal integration tests for the save_assemblies_from_fastas function. Detailed tests should
go in the unit tests. For instance, here we test only a couple of error types and do not test
error cases exhaustively.

These tests are designed to be run via the kb-sdk test command, and therefore have many
expectations about the environment in which they run.
'''

import dateutil.parser
import os
import re
import requests
import uuid
import shutil

from configparser import ConfigParser
from pathlib import Path
from pytest import fixture, raises

from AssemblyUtil.authclient import KBaseAuth
from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
from AssemblyUtil.AssemblyUtilServer import MethodContext
from AssemblyUtil.FastaToAssembly import FastaToAssembly
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

from conftest import assert_close_to_time, assert_exception_correct


def _ref(object_info):
    return f'{object_info[6]}/{object_info[0]}/{object_info[4]}'


@fixture(scope='module')
def config():
    cfg = {'token': os.environ['KB_AUTH_TOKEN']}
    parser = ConfigParser()
    parser.read(os.environ['KB_DEPLOYMENT_CONFIG'])
    cfg['file_config'] = {k: v for k, v in parser.items('AssemblyUtil')}

    # add catalog secure params
    cfg['file_config']["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10"
    cfg['file_config']["KBASE_SECURE_CONFIG_PARAM_THREADS_PER_CPU"] = "2"

    auth_client = KBaseAuth(cfg['file_config']['auth-service-url'])
    cfg['user'] = auth_client.get_user(cfg['token'])

    cfg['callback_url'] = os.environ['SDK_CALLBACK_URL']

    ws = _get_workspace(cfg)
    ws_info = ws.create_workspace({'workspace': f'{Path(__file__).stem}_{uuid.uuid4()}'})
    cfg['ws_id'] = ws_info[0]
    cfg['ws_name'] = ws_info[1]

    yield cfg

    _clean_up(cfg)


def _get_workspace(cfg):
    return Workspace(cfg['file_config']['workspace-url'], token=cfg['token'])


def _clean_up(cfg):
    ws = _get_workspace(cfg)
    # assume here that the tests don't create 10K+ objects
    objects = ws.list_objects({'ids': [cfg['ws_id']], 'showAllVersions': 1})
    if not objects:
        ws.delete_workspace({'id': cfg['ws_id']})
        return
    objdata = ws.get_objects2({
        'objects': [{'ref': _ref(oi)} for oi in objects],
        'no_data': 1,
        'skip_external_system_updates': 1,  # we own any created nodes so no need for ACL updates
        })
    hs = HandleService(url=cfg['file_config']['handle-service-url'], token=cfg['token'])
    hids = []
    for od in objdata['data']:
        if 'handle' in od['extracted_ids']:
            hids.extend(od['extracted_ids']['handle'])
    handles = hs.hids_to_handles(hids)
    header = {'Authorization': f'Oauth {cfg["token"]}'}
    for h in handles:
        requests.delete(  # a mass delete endpoint would be nice
            f'{cfg["file_config"]["shock-url"]}/node/{h["id"]}',
            headers=header,
            allow_redirects=True)
        print(('Deleted shock node ' + h["id"]))
    hs.delete_handles(handles)
    ws.delete_workspace({'id': cfg['ws_id']})
    

@fixture(scope='module')
def context(config) -> MethodContext:
    # WARNING: don't call any logging methods on the context object,
    # it'll result in a NoneType error
    ctx = MethodContext(None)
    ctx.update({'token': config['token'],
                'user_id': config['user'],
                'authenticated': 1
                })
    yield ctx


@fixture(scope='module')
def scratch(config) -> Path:
    yield Path(config['file_config']['scratch'])


@fixture(scope='module')  # no need to recreate I think?
def impl(config) -> AssemblyUtil:
    yield AssemblyUtil(config['file_config'])


def test_import_local_files(config, impl, context, scratch):
    tmp_dir = scratch / ("test_import_local_files" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    shutil.copy(data / 'legacy_test.fna', tmp_dir)
    shutil.copy(data / 'test2.fna.gz', tmp_dir)

    res = impl.save_assemblies_from_fastas(
        context,
        {
            'workspace_id': config['ws_id'],
            'min_contig_length': 9,
            'inputs': [
                {
                    'file': tmp_dir / 'legacy_test.fna',
                    'assembly_name': 'legacy',
                    'type': 'isolate',
                    'external_source': "yer mum's doody",
                    'external_source_id': '#458150',
                    'contig_info': {'s3': {'is_circ': 1, 'description': 'whee'}}
                },
                {
                    'file': tmp_dir / 'test2.fna.gz',
                    'assembly_name': 'test2',
                }
            ]
        }
    )[0]['results']
    # middle part of path is a UUID, so no way to test in an integration test
    assert res[0]['filtered_input'].endswith('/legacy_test.fna.filtered.fa')
    assert res[1]['filtered_input'].endswith('/test2.fna.filtered.fa')

    handleid, blobid = _check_object(
        config,
        res[0]['upa'],
        {
            'name': 'legacy',
            'meta': {
                'GC content': '0.44444',
                'MD5': 'eba4d1771060e19671a56832d159526e',
                'N Contigs': '1',
                'Size': '18'
            },
            'data': {
                'assembly_id': 'legacy',
                'base_counts': {'A': 0, 'C': 1, 'G': 7, 'T': 10},
                'contigs': {'s3': {'contig_id': 's3',
                                   'description': 'whee',
                                   'gc_content': 0.44444,
                                   'is_circ': 1,
                                   'length': 18,
                                   'md5': '4f339bd56e5f43ecb52e8682a790a111',
                                   'name': 's3'
                }},
                'dna_size': 18,
                'external_source': "yer mum's doody",
                'external_source_id': '#458150',
                'fasta_handle_info': {
                    'handle': {
                        'file_name': 'legacy_test.fna.filtered.fa',
                        'remote_md5': '74fd2a681507d651e772316c2f5abf36',
                        'type': 'shock',
                    },
                    'node_file_name': 'legacy_test.fna.filtered.fa',
                    'size': 23
                },
                'gc_content': 0.44444,
                'md5': 'eba4d1771060e19671a56832d159526e',
                'num_contigs': 1,
                'type': 'isolate'
            }
        }
    )
    _check_handle(
        config,
        {
            'hid': handleid,
            'id': blobid,
            'created_by': config['user'],
            'file_name': 'legacy_test.fna.filtered.fa',
            'remote_md5': '74fd2a681507d651e772316c2f5abf36',
            'remote_sha1': None,
            'type': 'shock',
        }
    )
    _check_blobstore_data(
        config,
        {
            'id': blobid,
            'attributes': None,
            'file': {
                'checksum': {'md5': '74fd2a681507d651e772316c2f5abf36'},
                'name': 'legacy_test.fna.filtered.fa',
                'size': 23,
            },
            'format': ''
        },
        "\n".join([
            ">s3",
            "ctgtgtttgtgtgtgtgt\n"
        ])
    )
    handleid, blobid = _check_object(
        config,
        res[1]['upa'],
        {
            'name': 'test2',
            'meta': {
                'GC content': '0.45595',
                'MD5': '666fba6a931ba3e490a3554e41032636',
                'N Contigs': '1',
                'Size': '2520'
            },
            'data': {
                'assembly_id': 'test2',
                'base_counts': {'A': 700, 'C': 559, 'G': 590, 'T': 671},
                'contigs': {'gi|113968346|ref|NC_008321.1|': {
                    'contig_id': 'gi|113968346|ref|NC_008321.1|',
                    'description': 'Shewanella sp. MR-4 chromosome, complete genome',
                    'gc_content': 0.45595,
                    'length': 2520,
                    'md5': '5262007827bef551836469d81057c1ed',
                    'name': 'gi|113968346|ref|NC_008321.1|'
                }},
                'dna_size': 2520,
                'fasta_handle_info': {
                    'handle': {
                        'file_name': 'test2.fna.filtered.fa',
                        'remote_md5': '7d1abc8233a710de84107c49deef9143',
                        'type': 'shock',
                    },
                    'node_file_name': 'test2.fna.filtered.fa',
                    'size': 2641
                },
                'gc_content': 0.45595,
                'md5': '666fba6a931ba3e490a3554e41032636',
                'num_contigs': 1,
                'type': 'Unknown'
            }
        }
    )
    _check_handle(
        config,
        {
            'hid': handleid,
            'id': blobid,
            'created_by': config['user'],
            'file_name': 'test2.fna.filtered.fa',
            'remote_md5': '7d1abc8233a710de84107c49deef9143',
            'remote_sha1': None,
            'type': 'shock',
        }
    )
    with open(data / 'test2.fna') as f:
        file_lines = f.readlines()
    file_txt = _reformat_1_line_fasta(file_lines, line_length=60)
    _check_blobstore_data(
        config,
        {
            'id': blobid,
            'attributes': None,
            'file': {
                'checksum': {'md5': '7d1abc8233a710de84107c49deef9143'},
                'name': 'test2.fna.filtered.fa',
                'size': 2641,
            },
            'format': ''
        },
        file_txt
    )


def test_import_blobstore_files(config, impl, context, scratch):
    # since DFU can't read the data dir
    source_tmp = scratch / ("test_import_blobstore_files_source" + str(uuid.uuid4()))
    os.makedirs(source_tmp)
    tmp_dir = scratch / ("test_import_blobstore_files" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    shutil.copy(data / 'legacy_test.fna', source_tmp)
    shutil.copy(data / 'test2.fna.gz', source_tmp)

    dfu = DataFileUtil(config['callback_url'], token=config['token'])
    nodes = dfu.file_to_shock_mass([
        {'file_path': str(source_tmp / 'legacy_test.fna')},
        {'file_path': str(source_tmp / 'test2.fna.gz')}
    ])

    res = impl.save_assemblies_from_fastas(
        context,
        {
            'workspace_id': config['ws_id'],
            'inputs': [
                {'node': nodes[0]['shock_id'], 'assembly_name': 'legacy2'},
                {'node': nodes[1]['shock_id'], 'assembly_name': 'test22'}
            ]
        }
    )[0]['results']
    assert res[0]['filtered_input'] is None
    assert res[1]['filtered_input'] is None

    handleid, blobid = _check_object(
        config,
        res[0]['upa'],
        {
            'name': 'legacy2',
            'meta': {
                'GC content': '0.5',
                'MD5': '961ca06ad332212d3797edd3a0515d84',
                'N Contigs': '3',
                'Size': '34'
            },
            'data': {
                'assembly_id': 'legacy2',
                'base_counts': {'A': 3, 'C': 2, 'G': 15, 'T': 14},
                'contigs': {
                    's1': {
                        'contig_id': 's1',
                        'description': 'something',
                        'gc_content': 0.75,
                        'length': 8,
                        'md5': 'db6906e8c40c19bcd822b92e82e5abe9',
                        'name': 's1'
                    },
                    's2': {
                        'contig_id': 's2',
                        'description': 'something2',
                        'gc_content': 0.375,
                        'length': 8,
                        'md5': '4c0e88f87fb2651722474f812b62507c',
                        'name': 's2'
                    },
                    's3': {
                        'contig_id': 's3',
                        'description': '',
                        'gc_content': 0.44444,
                        'length': 18,
                        'md5': '4f339bd56e5f43ecb52e8682a790a111',
                        'name': 's3'
                    }
                },
                'dna_size': 34,
                'fasta_handle_info': {
                    'handle': {
                        'file_name': 'legacy_test.fna',
                        'remote_md5': '0e9cb19516096741da1d55c6ecfb0cec',
                        'type': 'shock',
                    },
                    'node_file_name': 'legacy_test.fna',
                    'size': 70
                },
                'gc_content': 0.5,
                'md5': '961ca06ad332212d3797edd3a0515d84',
                'num_contigs': 3,
                'type': 'Unknown'
            }
        }
    )
    _check_handle(
        config,
        {
            'hid': handleid,
            'id': blobid,
            'created_by': config['user'],
            'file_name': 'legacy_test.fna',
            'remote_md5': '0e9cb19516096741da1d55c6ecfb0cec',
            'remote_sha1': None,
            'type': 'shock',
        }
    )
    with open(data / "legacy_test.fna") as f:
        file_ = f.read()
    _check_blobstore_data(
        config,
        {
            'id': blobid,
            'attributes': None,
            'file': {
                'checksum': {'md5': '0e9cb19516096741da1d55c6ecfb0cec'},
                'name': 'legacy_test.fna',
                'size': 70,
            },
            'format': ''
        },
        file_
    )
    handleid, blobid = _check_object(
        config,
        res[1]['upa'],
        {
            'name': 'test22',
            'meta': {
                'GC content': '0.45595',
                'MD5': '666fba6a931ba3e490a3554e41032636',
                'N Contigs': '1',
                'Size': '2520'
            },
            'data': {
                'assembly_id': 'test22',
                'base_counts': {'A': 700, 'C': 559, 'G': 590, 'T': 671},
                'contigs': {'gi|113968346|ref|NC_008321.1|': {
                    'contig_id': 'gi|113968346|ref|NC_008321.1|',
                    'description': 'Shewanella sp. MR-4 chromosome, complete genome',
                    'gc_content': 0.45595,
                    'length': 2520,
                    'md5': '5262007827bef551836469d81057c1ed',
                    'name': 'gi|113968346|ref|NC_008321.1|'
                }},
                'dna_size': 2520,
                'fasta_handle_info': {
                    'handle': {
                        'file_name': 'test2.fna',
                        'remote_md5': 'c3fb50316e293dcb937f9c73d09b443f',
                        'type': 'shock',
                    },
                    'node_file_name': 'test2.fna'   ,
                    'size': 2634
                },
                'gc_content': 0.45595,
                'md5': '666fba6a931ba3e490a3554e41032636',
                'num_contigs': 1,
                'type': 'Unknown'
            }
        }
    )
    _check_handle(
        config,
        {
            'hid': handleid,
            'id': blobid,
            'created_by': config['user'],
            'file_name': 'test2.fna',
            'remote_md5': 'c3fb50316e293dcb937f9c73d09b443f',
            'remote_sha1': None,
            'type': 'shock',
        }
    )
    with open(data / 'test2.fna') as f:
        file_ = f.read()
    _check_blobstore_data(
        config,
        {
            'id': blobid,
            'attributes': None,
            'file': {
                'checksum': {'md5': 'c3fb50316e293dcb937f9c73d09b443f'},
                'name': 'test2.fna',
                'size': 2634,
            },
            'format': ''
        },
        file_
    )


def _check_object(config, upa, expected):
    ws = _get_workspace(config)
    obj = ws.get_objects2({'objects': [{'ref': upa}]})['data'][0]

    # Check relevant object info fields
    info = obj['info']
    assert info[1] == expected['name']
    assert info[2].split('-')[0] == 'KBaseGenomeAnnotations.Assembly'
    assert info[6] == config['ws_id']
    # md5 and size will change as handle & shock data changes
    assert info[10] == expected['meta']

    # Check handle info and remove mutable data
    data = obj['data']
    handle_id = data.pop('fasta_handle_ref')
    assert handle_id.split('_')[0] == 'KBH'
    handle_info = data['fasta_handle_info']
    blobid = handle_info.pop('shock_id')
    handle = handle_info['handle']
    assert handle.pop('hid') == handle_id
    assert handle.pop('id') == blobid
    url = handle.pop('url')
    assert url.startswith('https://')
    assert url.endswith('kbase.us/services/shock-api')
    assert obj['extracted_ids']['handle'] == [handle_id]

    # Check cleaned object data
    assert obj['data'] == expected['data']

    # Doesn't seem useful to test the rest of the contents of obj, we'd just be testing
    # DFU and WS code

    return handle_id, blobid


def _check_handle(cfg, expected):
    hs = HandleService(cfg['file_config']['handle-service-url'], token=cfg['token'])
    handle = hs.hids_to_handles([expected['hid']])[0]
    url = handle.pop('url')
    assert url.startswith('https://')
    assert url.endswith('kbase.us/services/shock-api')
    # TODO why does this test take so long? The data sizes are tiny
    assert_close_to_time(handle.pop('creation_date'), seconds_slop=10)
    assert handle == expected


def _check_blobstore_data(cfg, expected_node, expected_file: str):
    node_id = expected_node['id']
    headers = {'Authorization': f'OAuth {cfg["token"]}'}
    url = f'{cfg["file_config"]["shock-url"]}/node/{node_id}'
    got_node = requests.get(url, headers=headers).json()

    created_str = got_node['data'].pop('created_on')
    created = dateutil.parser.isoparse(created_str).timestamp()
    assert_close_to_time(created, seconds_slop=10)
    assert created_str == got_node['data'].pop('last_modified')

    assert got_node == {
        'data': expected_node,
        'error': None,
        'status': 200
    }
    
    got_file = requests.get(f'{url}?download', headers=headers).text
    assert got_file == expected_file

def _reformat_1_line_fasta(file_lines, line_length=70):
    # Expand to multiple lines when needed
    lines = [file_lines[0].strip()]
    seq = "".join([l.strip() for l in file_lines[1:]])
    lines.extend([seq[i:i+line_length] for i in range(0, len(seq), line_length)])
    return '\n'.join(lines) + '\n'


def test_fail_no_such_file(config, impl, context, scratch):
    tmp_dir = scratch / ("test_fail_no_such_file" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    shutil.copy(data / 'legacy_test.fna', tmp_dir)

    with raises(Exception) as got:
        impl.save_assemblies_from_fastas(
            context,
            {
                'workspace_id': config['ws_id'],
                'inputs': [
                    {'file': tmp_dir / 'legacy_test.fna', 'assembly_name': 'legacy'},
                    {'file': tmp_dir / 'test2.fna.gz', 'assembly_name': 'test2'}
                ]
            }
        )
    assert_exception_correct(got.value, ValueError(
        "KBase Assembly Utils tried to save an assembly, but the calling application "
        + f"specified a file ('{tmp_dir / 'test2.fna.gz'}') that is missing. Please check the " +
        "application logs for details."))


def test_fail_mismatched_inputs(config, impl, context, scratch):
    tmp_dir = scratch / ("test_fail_no_such_file" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    shutil.copy(data / 'legacy_test.fna', tmp_dir)

    with raises(Exception) as got:
        impl.save_assemblies_from_fastas(
            context,
            {
                'workspace_id': config['ws_id'],
                'inputs': [
                    {'file': tmp_dir / 'legacy_test.fna', 'assembly_name': 'legacy'},
                    {'node': 'fake_node', 'assembly_name': 'test2'}
                ]
            }
        )
    assert_exception_correct(got.value, ValueError(
        'Entry #2 in inputs must have a file field to match entry #1'))


def test_parallelize_import_fasta_mass(config, impl, context, scratch):
    tmp_dir = scratch / ("test_parallelize_import_fasta_mass" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    file_names = [
        'GCA_002506415.1_ASM250641v1_genomic.fna.gz',
        'GCF_000007065.1_ASM706v1_genomic.fna.gz',
        'GCF_000970165.1_ASM97016v1_genomic.fna.gz',
        'GCF_000970185.1_ASM97018v1_genomic.fna.gz',
        'GCF_000970205.1_ASM97020v1_genomic.fna.gz',
        'GCF_000970245.1_ASM97024v1_genomic.fna.gz'
    ]

    object_metas = [
        {'GC content': '0.41675', 'Size': '3839682', 'N Contigs': '130', 'MD5': '8118290bf6d3369f78ebb70a59d85dd3'},
        {'GC content': '0.41485', 'Size': '4096345', 'N Contigs': '1', 'MD5': '418bf8e730b1b948ffb9cded9acfaf26'},
        {'GC content': '0.41457', 'Size': '4096482', 'N Contigs': '1', 'MD5': '949a0fe665048cb917c8cf74f75c74b7'},
        {'GC content': '0.41487', 'Size': '4066551', 'N Contigs': '1', 'MD5': 'd33802829ba0686714a5d74280527615'},
        {'GC content': '0.41421', 'Size': '4142816', 'N Contigs': '1', 'MD5': 'cf47d74f66a16dffcbaa7a05eb9eec70'},
        {'GC content': '0.41488', 'Size': '4166241', 'N Contigs': '1', 'MD5': '90178f629aa7bfbeea19bac8e616c467'}
    ]

    # copy 6 assembly files into the data dir
    for file_name in file_names:
        shutil.copy(data / file_name, tmp_dir)

    params = {
        'workspace_id': config['ws_id'],
        'inputs': [
            {'file': tmp_dir / 'GCA_002506415.1_ASM250641v1_genomic.fna.gz', 'assembly_name': 'GCA_002506415.1_ASM250641v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000007065.1_ASM706v1_genomic.fna.gz', 'assembly_name': 'GCF_000007065.1_ASM706v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970165.1_ASM97016v1_genomic.fna.gz', 'assembly_name': 'GCF_000970165.1_ASM97016v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970185.1_ASM97018v1_genomic.fna.gz', 'assembly_name': 'GCF_000970185.1_ASM97018v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970205.1_ASM97020v1_genomic.fna.gz', 'assembly_name': 'GCF_000970205.1_ASM97020v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970245.1_ASM97024v1_genomic.fna.gz', 'assembly_name': 'GCF_000970245.1_ASM97024v1_genomic.fna.gz'}
        ]
    }

    object_version_pattern = re.compile(r'^[0-9]+\/1$')
    results = impl.save_assemblies_from_fastas(context, params)[0]['results']
    ws = _get_workspace(config)
    for idx, res in enumerate(results):
        # check the workspace is returning UPAs
        assert res['filtered_input'] is None
        assert object_version_pattern.match("/".join(res['upa'].split("/")[-2:]))

        # check relevant object info fields
        obj = ws.get_objects2({'objects': [{'ref': res['upa']}]})['data'][0]
        info = obj['info']
        assert info[1] == file_names[idx]
        assert info[2].split('-')[0] == 'KBaseGenomeAnnotations.Assembly'
        assert info[6] == config['ws_id']
        assert info[10] == object_metas[idx]


def test_invalid_thread_param(config, context, scratch):
    tmp_dir = scratch / ("test_invalid_thread_param" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    file_names = [
        'GCA_002506415.1_ASM250641v1_genomic.fna.gz',
        'GCF_000007065.1_ASM706v1_genomic.fna.gz'
    ]

    # copy 2 assembly files into the data dir
    for file_name in file_names:
        shutil.copy(data / file_name, tmp_dir)

    params = {
        'workspace_id': config['ws_id'],
        'inputs': [
            {'file': tmp_dir / 'GCA_002506415.1_ASM250641v1_genomic.fna.gz', 'assembly_name': 'GCA_002506415.1_ASM250641v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000007065.1_ASM706v1_genomic.fna.gz', 'assembly_name': 'GCF_000007065.1_ASM706v1_genomic.fna.gz'}
        ]
    }

    var_name = "MAX_THREADS"
    config['file_config']["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "10e"
    impl = AssemblyUtil(config['file_config'])
    with raises(Exception) as got:
        impl.save_assemblies_from_fastas(context, params)
    assert_exception_correct(got.value, ValueError(f"Cannot evaluate the string {var_name}"))

    config['file_config']["KBASE_SECURE_CONFIG_PARAM_MAX_THREADS"] = "-1"
    impl = AssemblyUtil(config['file_config'])
    with raises(Exception) as got:
        impl.save_assemblies_from_fastas(context, params)
    assert_exception_correct(got.value, ValueError(f"{var_name} cannot be negative"))


def test_generator_overflow(config, scratch):
    tmp_dir = scratch / ("test_generator_overflow" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    file_names = [
        'GCA_002506415.1_ASM250641v1_genomic.fna.gz',
        'GCF_000007065.1_ASM706v1_genomic.fna.gz',
        'GCF_000970165.1_ASM97016v1_genomic.fna.gz',
        'GCF_000970185.1_ASM97018v1_genomic.fna.gz',
        'GCF_000970205.1_ASM97020v1_genomic.fna.gz',
        'GCF_000970245.1_ASM97024v1_genomic.fna.gz'
    ]

    # copy 6 assembly files into the data dir
    for file_name in file_names:
        shutil.copy(data / file_name, tmp_dir)

    params = {
        'workspace_id': config['ws_id'],
        'inputs': [
            {'file': tmp_dir / 'GCA_002506415.1_ASM250641v1_genomic.fna.gz', 'assembly_name': 'GCA_002506415.1_ASM250641v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000007065.1_ASM706v1_genomic.fna.gz', 'assembly_name': 'GCF_000007065.1_ASM706v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970165.1_ASM97016v1_genomic.fna.gz', 'assembly_name': 'GCF_000970165.1_ASM97016v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970185.1_ASM97018v1_genomic.fna.gz', 'assembly_name': 'GCF_000970185.1_ASM97018v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970205.1_ASM97020v1_genomic.fna.gz', 'assembly_name': 'GCF_000970205.1_ASM97020v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970245.1_ASM97024v1_genomic.fna.gz', 'assembly_name': 'GCF_000970245.1_ASM97024v1_genomic.fna.gz'}
        ]
    }

    dfu = DataFileUtil(config['callback_url'], token=config['token'])
    fta = FastaToAssembly(dfu, scratch)
    fta.import_fasta_mass(params, 1, 10, 38300)


def test_invalid_max_cumsize(config, scratch):
    tmp_dir = scratch / ("test_invalid_max_cumsize" + str(uuid.uuid4()))
    os.makedirs(tmp_dir)
    data = Path('data')
    file_names = [
        'GCA_002506415.1_ASM250641v1_genomic.fna.gz',
        'GCF_000007065.1_ASM706v1_genomic.fna.gz',
        'GCF_000970165.1_ASM97016v1_genomic.fna.gz',
        'GCF_000970185.1_ASM97018v1_genomic.fna.gz',
        'GCF_000970205.1_ASM97020v1_genomic.fna.gz',
        'GCF_000970245.1_ASM97024v1_genomic.fna.gz'
    ]

    # copy 6 assembly files into the data dir
    for file_name in file_names:
        shutil.copy(data / file_name, tmp_dir)

    params = {
        'workspace_id': config['ws_id'],
        'inputs': [
            {'file': tmp_dir / 'GCA_002506415.1_ASM250641v1_genomic.fna.gz', 'assembly_name': 'GCA_002506415.1_ASM250641v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000007065.1_ASM706v1_genomic.fna.gz', 'assembly_name': 'GCF_000007065.1_ASM706v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970165.1_ASM97016v1_genomic.fna.gz', 'assembly_name': 'GCF_000970165.1_ASM97016v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970185.1_ASM97018v1_genomic.fna.gz', 'assembly_name': 'GCF_000970185.1_ASM97018v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970205.1_ASM97020v1_genomic.fna.gz', 'assembly_name': 'GCF_000970205.1_ASM97020v1_genomic.fna.gz'},
            {'file': tmp_dir / 'GCF_000970245.1_ASM97024v1_genomic.fna.gz', 'assembly_name': 'GCF_000970245.1_ASM97024v1_genomic.fna.gz'}
        ]
    }
    dfu = DataFileUtil(config['callback_url'], token=config['token'])
    fta = FastaToAssembly(dfu, scratch)
    with raises(Exception) as got:
        fta.import_fasta_mass(params, 1, 10, "38300")
    assert_exception_correct(got.value, ValueError("max_cumsize must be an integer or decimal"))

    with raises(Exception) as got:
        fta.import_fasta_mass(params, 1, 10, -1)
    assert_exception_correct(got.value, ValueError("max_cumsize must be > 0"))

    with raises(Exception) as got:
        fta.import_fasta_mass(params, 1, 10, 1024 * 1024 * 1024 * 1024)
    assert_exception_correct(got.value, ValueError(f"max_cumsize must be <= {1024 * 1024 * 1024 * 0.95}"))