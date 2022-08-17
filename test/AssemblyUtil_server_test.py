import json
import os
import shutil
import time
import unittest
import uuid
from os import environ
from pathlib import Path
from pprint import pprint
from configparser import ConfigParser

from AssemblyUtil.authclient import KBaseAuth as _KBaseAuth
from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
from AssemblyUtil.AssemblyUtilServer import MethodContext
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as workspaceService


class AssemblyUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('AssemblyUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        authServiceUrl = cls.cfg.get('auth-service-url',
                                     'https://kbase.us/services/authorization/Sessions/Login')
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'AssemblyUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = AssemblyUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_AssemblyUtil_" + str(suffix)
        self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def check_fasta_file(self, ws_obj_name, expected_file):
        assemblyUtil = self.getImpl()

        print('attempting download')
        fasta = assemblyUtil.get_assembly_as_fasta(
            self.getContext(),
            {'ref': self.getWsName() + "/" + ws_obj_name}
        )[0]
        pprint(fasta)
        # let's compare files pointed from fasta['path'] and expected_file
        with open(expected_file, 'r') as f:
            expected_data = f.read()
        with open(fasta['path'], 'r') as f:
            actual_data = f.read()
        self.assertEqual(actual_data, expected_data)

    # @unittest.skip('x')
    def test_basic_upload_and_download(self):
        assemblyUtil = self.getImpl()

        # use a temp dir for input files since when the files are downloaded to test them
        # in _check_fasta_file they always go in the root directory, creating the potential
        # for clobberation.
        # Get assemblies should be updated to allow specifying a target dir
        tmp_dir = Path(self.cfg['scratch']) / (
            "test_basic_upload_and_download_input" + str(uuid.uuid4()))
        os.makedirs(tmp_dir)
        file_name = "test.fna"
        shutil.copy(os.path.join("data", file_name), tmp_dir)
        fasta_path = os.path.join(tmp_dir, file_name)
        print('attempting upload')
        ws_obj_name = 'MyNewAssembly'
        result = assemblyUtil.save_assembly_from_fasta(self.getContext(),
                                                       {'file': {'path': fasta_path},
                                                        'workspace_name': self.getWsName(),
                                                        'assembly_name': ws_obj_name,
                                                        'taxon_ref': 'ReferenceTaxons/unknown_taxon',
                                                        })
        pprint(result)
        self.check_fasta_file(ws_obj_name, fasta_path)

        print('attempting upload of gzipped file')
        # DFU can't see the file unless it's in scratch
        shutil.copy(Path("data") / "test2.fna.gz", tmp_dir)
        result = assemblyUtil.save_assembly_from_fasta(
            self.getContext(),
            {
                'file': {'path': str(Path(tmp_dir) / "test2.fna.gz")},
                'workspace_name': self.getWsName(),
                'assembly_name': 'MyNewAssembly_gzip',
                'taxon_ref': 'ReferenceTaxons/unknown_taxon',
            })
        pprint(result)
        self.check_fasta_file('MyNewAssembly_gzip', "data/test2.fna")

        print('attempting upload through shock')
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        shock_id = data_file_cli.file_to_shock({'file_path': fasta_path})['shock_id']
        ws_obj_name2 = 'MyNewAssembly.2'
        result2 = assemblyUtil.save_assembly_from_fasta(self.getContext(),
                                                        {'shock_id': shock_id,
                                                         'workspace_name': self.getWsName(),
                                                         'assembly_name': ws_obj_name2
                                                         })
        pprint(result2)
        self.check_fasta_file(ws_obj_name2, fasta_path)

    def test_empty_file_error_message(self):
        assemblyUtil = self.getImpl()
        tmp_dir = self.__class__.cfg['scratch']
        file_name = "empty.FASTA"
        shutil.copy(os.path.join("data", file_name), tmp_dir)
        empty_file_path = os.path.join(tmp_dir, file_name)
        try:
            assemblyUtil.save_assembly_from_fasta(self.getContext(),
                                                    {'file': {'path': empty_file_path},
                                                     'workspace_name': self.getWsName(),
                                                     'assembly_name': 'empty',
                                                     'taxon_ref': 'ReferenceTaxons/unknown_taxon',
                                                     "min_contig_length": 500
                                                     })
            raise Exception("Wrong error message")
        except ValueError as e:
            self.assertEqual(
                'Either the original FASTA file contained no sequences or they were all '
                + 'filtered out based on the min_contig_length parameter',
                str(e))

    def test_legacy_contigset_download(self):
        ws_obj_name4 = 'LegacyContigs'
        contigset_data = {'id': 'contigset',
                          'md5': '5217cb4e8684817ba925ebe125caf54a',
                          'source': 'file',
                          'source_id': 'file',
                          'contigs': [{
                                        'id': 's1',
                                        'description': 'something',
                                        'sequence': 'agtggggg'
                                       }, {
                                        'id': 's2',
                                        'description': 'something2',
                                        'sequence': 'gacgattt'
                                       }, {
                                        'id': 's3',
                                        'sequence': 'ctgtgtttgtgtgtgtgt'
                                       }]
                          }

        obj_info = self.getWsClient().save_objects({'workspace': self.getWsName(),
                                                    'objects': [{'type': 'KBaseGenomes.ContigSet',
                                                                 'name': ws_obj_name4,
                                                                 'data': contigset_data
                                                                 }]
                                                    })[0]

        assemblyUtil = self.getImpl()
        fasta = assemblyUtil.get_assembly_as_fasta(self.getContext(),
                                                   {'ref': self.getWsName() + '/' + obj_info[1],
                                                    'filename': 'legacy.fa'})[0]

        file_name = "legacy_test.fna"
        scratch_dir = self.__class__.cfg['scratch']
        shutil.copy(os.path.join("data", file_name), scratch_dir)
        expected_fasta_path = os.path.join(scratch_dir, file_name)
        with open(expected_fasta_path, 'r') as f:
            expected_data = f.read()
        with open(fasta['path'], 'r') as f:
            actual_data = f.read()
        self.assertEqual(actual_data, expected_data)

    def test_load_with_filter_and_options(self):
        assemblyUtil = self.getImpl()
        tmp_dir = self.__class__.cfg['scratch']
        file_name = "legacy_test.fna"
        shutil.copy(os.path.join("data", file_name), tmp_dir)
        fasta_path = os.path.join(tmp_dir, file_name)
        print('attempting upload')
        ws_obj_name = 'FilteredAssembly'
        result = assemblyUtil.save_assembly_from_fasta(self.getContext(),
                                                       {'file': {'path': fasta_path},
                                                        'workspace_name': self.getWsName(),
                                                        'assembly_name': ws_obj_name,
                                                        'min_contig_length': 9,
                                                        'external_source': 'someplace',
                                                        'external_source_id': 'id',
                                                        'external_source_origination_date': 'sunday',
                                                        'type': 'metagenome',
                                                        'contig_info': {'s3': {'is_circ': 0, 'description': 'somethin'}}
                                                        })

        dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])
        assembly = dfu.get_objects({'object_refs': [result[0]]})['data'][0]['data']

        self.assertEqual(len(assembly['contigs']), 1)
        self.assertEqual(assembly['contigs']['s3']['md5'], '4f339bd56e5f43ecb52e8682a790a111')
        self.assertEqual(assembly['contigs']['s3']['contig_id'], 's3')
        self.assertEqual(assembly['contigs']['s3']['length'], 18)
        self.assertEqual(assembly['contigs']['s3']['is_circ'], 0)
        self.assertEqual(assembly['contigs']['s3']['description'], 'somethin')

        self.assertEqual(assembly['dna_size'], 18)
        self.assertEqual(assembly['gc_content'], 0.44444)
        self.assertEqual(assembly['md5'], 'eba4d1771060e19671a56832d159526e')
        self.assertEqual(assembly['num_contigs'], 1)
        self.assertEqual(assembly['type'], 'metagenome')
        self.assertEqual(assembly['external_source'], 'someplace')
        self.assertEqual(assembly['external_source_id'], 'id')
        self.assertEqual(assembly['external_source_origination_date'], 'sunday')

    def test_assembly_does_not_exist(self):
        assemblyUtil = self.getImpl()
        tmp_dir = self.__class__.cfg['scratch']
        file_name = "not_a_real_file.fna"
        fasta_path = tmp_dir + "/" + file_name
        print('attempting upload')
        ws_obj_name = 'FilteredAssembly'
        with self.assertRaisesRegex(ValueError, 'KBase Assembly Utils tried to save an assembly'):
            assemblyUtil.save_assembly_from_fasta(self.getContext(),
                                                       {'file': {'path': fasta_path},
                                                        'workspace_name': self.getWsName(),
                                                        'assembly_name': ws_obj_name,
                                                        'min_contig_length': 500
                                                        })
