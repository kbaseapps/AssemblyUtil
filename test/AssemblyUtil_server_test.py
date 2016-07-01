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
from shock_util.shock_utilClient import shock_util


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



    def test_stuff(self):
        assemblyUtil = self.getImpl()

        tmp_dir = self.__class__.cfg['scratch']
        file_name = "test.fna"
        shutil.copy(os.path.join("data", file_name), tmp_dir)
        fasta_path = os.path.join(tmp_dir, file_name)
        print('attempting upload')
        result = assemblyUtil.save_assembly_from_fasta(self.getContext(), 
            {
                'file':{'path':fasta_path},
                'workspace_name':self.getWsName(),
                'assembly_name':'MyNewAssembly'
            });
        pprint(result)

        print('attempting download')
        fasta = assemblyUtil.get_assembly_as_fasta(self.getContext(), 
            {'ref':self.getWsName() + "/MyNewAssembly"})[0]
        pprint(fasta)

        print('attempting upload through shock')
        shock_cli = shock_util(os.environ['SDK_CALLBACK_URL'], token=
                               self.__class__.ctx['token'],
                               service_ver='dev')
        shock_id = shock_cli.file_to_node({'file_path': fasta_path})['shock_id']
        result2 = assemblyUtil.save_assembly_from_fasta(self.getContext(), 
            {
                'shock_id':shock_id,
                'workspace_name':self.getWsName(),
                'assembly_name':'MyNewAssembly.2'
            });
        pprint(result2)


