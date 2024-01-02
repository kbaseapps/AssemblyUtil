# import json
# import os
# import shutil
# import time
# import unittest
# import uuid
# from os import environ
# from pathlib import Path
# from pprint import pprint
# from configparser import ConfigParser

# from AssemblyUtil.authclient import KBaseAuth as _KBaseAuth
# from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
# from AssemblyUtil.AssemblyUtilServer import MethodContext
# from installed_clients.DataFileUtilClient import DataFileUtil
# from installed_clients.WorkspaceClient import Workspace as workspaceService


# class AssemblyUtilTest(unittest.TestCase):

#     @classmethod
#     def setUpClass(cls):
#         token = environ.get('KB_AUTH_TOKEN', None)
#         config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
#         cls.cfg = {}
#         config = ConfigParser()
#         config.read(config_file)
#         for nameval in config.items('AssemblyUtil'):
#             cls.cfg[nameval[0]] = nameval[1]
#         authServiceUrl = cls.cfg['auth-service-url']
#         auth_client = _KBaseAuth(authServiceUrl)
#         user_id = auth_client.get_user(token)
#         # WARNING: don't call any logging methods on the context object,
#         # it'll result in a NoneType error
#         cls.ctx = MethodContext(None)
#         cls.ctx.update({'token': token,
#                         'user_id': user_id,
#                         'provenance': [
#                             {'service': 'AssemblyUtil',
#                              'method': 'please_never_use_it_in_production',
#                              'method_params': []
#                              }],
#                         'authenticated': 1})
#         cls.wsURL = cls.cfg['workspace-url']
#         cls.wsClient = workspaceService(cls.wsURL, token=token)
#         ws = cls.wsClient.create_workspace(
#             {'workspace': 'test_AssemblyUtil_' + str(int(time.time() * 1000))})
#         cls.ws_name = ws[1]
#         cls.ws_id = ws[0]
#         cls.serviceImpl = AssemblyUtil(cls.cfg)
#         cls.scratch = cls.cfg['scratch']
#         cls.callback_url = os.environ['SDK_CALLBACK_URL']

#     @classmethod
#     def tearDownClass(cls):
#         if hasattr(cls, 'ws_name'):
#             cls.wsClient.delete_workspace({'workspace': cls.ws_name})
#             print('Test workspace was deleted')

#     def getWsClient(self):
#         return self.wsClient

#     def getImpl(self):
#         return self.serviceImpl

#     def getContext(self):
#         return self.ctx

#     def check_fasta_file(self, ref, expected_file):
#         assemblyUtil = self.getImpl()

#         print('attempting download')
#         fasta = assemblyUtil.get_assembly_as_fasta(self.getContext(), {'ref': ref})[0]
#         pprint(fasta)
#         # let's compare files pointed from fasta['path'] and expected_file
#         with open(expected_file, 'r') as f:
#             expected_data = f.read()
#         with open(fasta['path'], 'r') as f:
#             actual_data = f.read()
#         self.assertEqual(actual_data, expected_data)

#     def check_ws_name(self, upa, name):
#         res = self.wsClient.get_object_info3({'objects': [{'ref': upa}]})
#         assert res['infos'][0][1] == name

#     def get_save_funcs(self, assemblyUtil):
#         return [assemblyUtil.save_assembly_from_fasta, assemblyUtil.save_assembly_from_fasta2]

#     # @unittest.skip('x')
#     def test_basic_upload_and_download(self):
#         assemblyUtil = self.getImpl()

#         # use a temp dir for input files since when the files are downloaded to test them
#         # in _check_fasta_file they always go in the root directory, creating the potential
#         # for clobberation.
#         # Get assemblies should be updated to allow specifying a target dir
#         tmp_dir = Path(self.cfg['scratch']) / (
#             "test_basic_upload_and_download_input" + str(uuid.uuid4()))
#         os.makedirs(tmp_dir)
#         file_name = "test.fna"
#         shutil.copy(os.path.join("data", file_name), tmp_dir)
#         fasta_path = os.path.join(tmp_dir, file_name)
#         print('attempting upload')
#         params = {'file': {'path': fasta_path},
#                   'workspace_name': self.ws_name,
#                   'taxon_ref': 'ReferenceTaxons/unknown_taxon',
#                   }
#         self._basic_upload_and_download_check(
#             assemblyUtil, params, "MyNewAssembly", "MyNewAssembly2", fasta_path)

#         print('attempting upload of gzipped file using workspace id vs name')
#         # DFU can't see the file unless it's in scratch
#         shutil.copy(Path("data") / "test2.fna.gz", tmp_dir)
#         params = {
#             'file': {'path': str(Path(tmp_dir) / "test2.fna.gz")},
#             'workspace_id': self.ws_id,
#             'taxon_ref': 'ReferenceTaxons/unknown_taxon',
#         }
#         self._basic_upload_and_download_check(
#             assemblyUtil, params, "MyNewAssembly_gzip", "MyNewAssembly_gzip2", "data/test2.fna")

#         print('attempting upload through shock')
#         data_file_cli = DataFileUtil(self.callback_url)
#         shock_id = data_file_cli.file_to_shock({'file_path': fasta_path})['shock_id']
#         params = {
#             'shock_id': shock_id,
#             'workspace_name': self.ws_name,
#         }
#         self._basic_upload_and_download_check(
#             assemblyUtil, params, "MyNewAssembly_shock", "MyNewAssembly_shock2", fasta_path)

#         # todo: add checks here on ws object

#         result4 = assemblyUtil.export_assembly_as_fasta(
#             self.getContext(),
#             {'input_ref': self.ws_name + '/' + 'MyNewAssembly_shock'})
#         pprint(result4)

#     def _basic_upload_and_download_check(self, assemblyUtil, params, name1, name2, fasta_path):
#         params['assembly_name'] = name1
#         result = assemblyUtil.save_assembly_from_fasta(self.getContext(), params)
#         pprint(result)
#         self.check_fasta_file(result[0], fasta_path)
#         self.check_ws_name(result[0], name1)
#         params['assembly_name'] = name2
#         result = assemblyUtil.save_assembly_from_fasta2(self.getContext(), params)
#         pprint(result)
#         self.check_fasta_file(result[0]['upa'], fasta_path)
#         self.check_ws_name(result[0]['upa'], name2)
#         assert result[0]['filtered_input'] is None

#     def test_empty_file_error_message(self):
#         assemblyUtil = self.getImpl()
#         tmp_dir = self.cfg['scratch']
#         file_name = "empty.FASTA"
#         shutil.copy(os.path.join("data", file_name), tmp_dir)
#         empty_file_path = os.path.join(tmp_dir, file_name)
#         for func in self.get_save_funcs(assemblyUtil):
#             try:
#                 func(self.getContext(),
#                     {'file': {'path': empty_file_path},
#                     'workspace_name': self.ws_name,
#                     'assembly_name': 'empty',
#                     'taxon_ref': 'ReferenceTaxons/unknown_taxon',
#                     "min_contig_length": 500
#                     })
#                 raise Exception("Wrong error message")
#             except ValueError as e:
#                 self.assertIn(
#                     'Either the original FASTA file contained no sequences or they were all '
#                     + 'filtered out based on the min_contig_length parameter for file '
#                     + '/kb/module/work/tmp/import_fasta_',
#                     str(e))
#                 self.assertIn('/empty.FASTA.filtered.fa', str(e))

#     def test_legacy_contigset_download(self):
#         ws_obj_name4 = 'LegacyContigs'
#         contigset_data = {'id': 'contigset',
#                           'md5': '5217cb4e8684817ba925ebe125caf54a',
#                           'source': 'file',
#                           'source_id': 'file',
#                           'contigs': [{
#                                         'id': 's1',
#                                         'description': 'something',
#                                         'sequence': 'agtggggg'
#                                        }, {
#                                         'id': 's2',
#                                         'description': 'something2',
#                                         'sequence': 'gacgattt'
#                                        }, {
#                                         'id': 's3',
#                                         'sequence': 'ctgtgtttgtgtgtgtgt'
#                                        }]
#                           }

#         obj_info = self.getWsClient().save_objects({'workspace': self.ws_name,
#                                                     'objects': [{'type': 'KBaseGenomes.ContigSet',
#                                                                  'name': ws_obj_name4,
#                                                                  'data': contigset_data
#                                                                  }]
#                                                     })[0]

#         assemblyUtil = self.getImpl()
#         fasta = assemblyUtil.get_assembly_as_fasta(self.getContext(),
#                                                    {'ref': self.ws_name + '/' + obj_info[1],
#                                                     'filename': 'legacy.fa'})[0]

#         file_name = "legacy_test.fna"
#         scratch_dir = self.__class__.cfg['scratch']
#         shutil.copy(os.path.join("data", file_name), scratch_dir)
#         expected_fasta_path = os.path.join(scratch_dir, file_name)
#         with open(expected_fasta_path, 'r') as f:
#             expected_data = f.read()
#         with open(fasta['path'], 'r') as f:
#             actual_data = f.read()
#         self.assertEqual(actual_data, expected_data)

#     def test_load_with_filter_and_options(self):
#         assemblyUtil = self.getImpl()
#         tmp_dir = self.__class__.cfg['scratch']
#         file_name = "legacy_test.fna"
#         shutil.copy(os.path.join("data", file_name), tmp_dir)
#         fasta_path = os.path.join(tmp_dir, file_name)
#         print('attempting upload')
#         ws_obj_name = 'FilteredAssembly'
#         params = {'file': {'path': fasta_path},
#                   'workspace_name': self.ws_name,
#                   'assembly_name': ws_obj_name,
#                   'min_contig_length': 9,
#                   'external_source': 'someplace',
#                   'external_source_id': 'id',
#                   'external_source_origination_date': 'sunday',
#                   'type': 'metagenome',
#                   'contig_info': {'s3': {'is_circ': 0, 'description': 'somethin'}}
#                   }
#         result = assemblyUtil.save_assembly_from_fasta(self.getContext(), params)
#         self._load_with_filter_and_options_check(result[0])
#         params['assembly_name'] = 'FilteredAssembly2'
#         result = assemblyUtil.save_assembly_from_fasta2(self.getContext(), params)
#         self._load_with_filter_and_options_check(result[0]['upa'])
#         self.assertTrue(result[0]['filtered_input'].endswith('/legacy_test.fna.filtered.fa'))

#     def _load_with_filter_and_options_check(self, upa):
#         dfu = DataFileUtil(self.callback_url)
#         assembly = dfu.get_objects({'object_refs': [upa]})['data'][0]['data']

#         self.assertEqual(len(assembly['contigs']), 1)
#         self.assertEqual(assembly['contigs']['s3']['md5'], '4f339bd56e5f43ecb52e8682a790a111')
#         self.assertEqual(assembly['contigs']['s3']['contig_id'], 's3')
#         self.assertEqual(assembly['contigs']['s3']['length'], 18)
#         self.assertEqual(assembly['contigs']['s3']['is_circ'], 0)
#         self.assertEqual(assembly['contigs']['s3']['description'], 'somethin')

#         self.assertEqual(assembly['dna_size'], 18)
#         self.assertEqual(assembly['gc_content'], 0.44444)
#         self.assertEqual(assembly['md5'], 'eba4d1771060e19671a56832d159526e')
#         self.assertEqual(assembly['num_contigs'], 1)
#         self.assertEqual(assembly['type'], 'metagenome')
#         self.assertEqual(assembly['external_source'], 'someplace')
#         self.assertEqual(assembly['external_source_id'], 'id')
#         self.assertEqual(assembly['external_source_origination_date'], 'sunday')

#     def test_assembly_does_not_exist(self):
#         assemblyUtil = self.getImpl()
#         tmp_dir = self.__class__.cfg['scratch']
#         file_name = "not_a_real_file.fna"
#         fasta_path = tmp_dir + "/" + file_name
#         print('attempting upload')
#         ws_obj_name = 'FilteredAssembly'
#         err = ("KBase Assembly Utils tried to save an assembly, but the calling "
#                + f"application specified a file \\('{tmp_dir}/not_a_real_file.fna'\\) "
#                + "that is missing. Please check the application logs for details.")
#         for func in self.get_save_funcs(assemblyUtil):
#             with self.assertRaisesRegex(ValueError, err):
#                 func(self.getContext(),
#                     {'file': {'path': fasta_path},
#                         'workspace_name': self.ws_name,
#                         'assembly_name': ws_obj_name,
#                         'min_contig_length': 500
#                     })
