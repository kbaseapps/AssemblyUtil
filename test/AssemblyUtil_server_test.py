import unittest
import os
import json
import time
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
from AssemblyUtil.AssemblyUtilServer import MethodContext
from DataFileUtil.DataFileUtilClient import DataFileUtil


class AssemblyUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'provenance': [
                            {'service': 'AssemblyUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('AssemblyUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = AssemblyUtil(cls.cfg)

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
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    def check_fasta_file(self, ws_obj_name, expected_file):
        assemblyUtil = self.getImpl()

        print('attempting download')
        fasta = assemblyUtil.get_assembly_as_fasta(self.getContext(), 
            {'ref':self.getWsName() + "/MyNewAssembly"})[0]
        pprint(fasta)
        # let's compare files pointed from fasta['path'] and expected_file
        expected_data = None
        with open(expected_file, 'r') as f:
            expected_data = f.read()
        actual_data = None
        with open(fasta['path'], 'r') as f:
            actual_data = f.read()
        self.assertEqual(actual_data, expected_data)

    def test_stuff(self):
        assemblyUtil = self.getImpl()

        tmp_dir = self.__class__.cfg['scratch']
        file_name = "test.fna"
        shutil.copy(os.path.join("data", file_name), tmp_dir)
        fasta_path = os.path.join(tmp_dir, file_name)
        print('attempting upload')
        ws_obj_name = 'MyNewAssembly'
        result = assemblyUtil.save_assembly_from_fasta(self.getContext(), 
            {
                'file':{'path':fasta_path},
                'workspace_name':self.getWsName(),
                'assembly_name':ws_obj_name
            });
        pprint(result)
        self.check_fasta_file(ws_obj_name, fasta_path)


        print('attempting upload through shock')
        data_file_cli = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=
                               self.__class__.ctx['token'],
                               service_ver='dev')
        shock_id = data_file_cli.file_to_shock({'file_path': fasta_path})['shock_id']
        ws_obj_name2 = 'MyNewAssembly.2'
        result2 = assemblyUtil.save_assembly_from_fasta(self.getContext(), 
            {
                'shock_id':shock_id,
                'workspace_name':self.getWsName(),
                'assembly_name':ws_obj_name2
            });
        pprint(result2)
        self.check_fasta_file(ws_obj_name2, fasta_path)

        print('attempting upload via ftp url')
        ftp_url = 'ftp://ftp.ensemblgenomes.org/pub/release-29/bacteria//fasta/bacteria_8_collection/acaryochloris_marina_mbic11017/dna/Acaryochloris_marina_mbic11017.GCA_000018105.1.29.dna.genome.fa.gz'
        ws_obj_name3 = 'MyNewAssembly.3'
        result3 = assemblyUtil.save_assembly_from_fasta(self.getContext(), 
            {
                'ftp_url':ftp_url,
                'workspace_name':self.getWsName(),
                'assembly_name':ws_obj_name3
            });
        pprint(result3)
        #todo: add checks here on ws object

