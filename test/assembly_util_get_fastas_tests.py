# -*- coding: utf-8 -*-
import os
import time
import unittest
import json
import shutil
from configparser import ConfigParser

from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
from AssemblyUtil.FastaToAssembly import FastaToAssembly
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from AssemblyUtil.AssemblyUtilServer import MethodContext
from AssemblyUtil.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace

class AssemblyUtil_FastaTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('AssemblyUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
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
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = AssemblyUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        shutil.copy2('data/NC_021490.fasta', cls.scratch)

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
        wsName = "test_kb_cmash_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def test_metagenome_binned_input(self):

        # Setup
        path = "data/binnedContigs.json"
        ws_path = '/kb/module/work/tmp'
        assembly_path = "data/CCESR16_SPAdes.assembly.fa"
        shutil.copy2(path, ws_path)
        shutil.copy2(assembly_path, ws_path)
        dfu = DataFileUtil(self.callback_url)
        wsName = self.getWsName()
        ws_id = dfu.ws_name_to_id(wsName)

        # FASTA to assembly object
        Fasta_assembly_dict = {"path": '/kb/module/work/tmp/CCESR16_SPAdes.assembly.fa', "assembly_name": "meta_assembly"}
        assembly_params = {"file": Fasta_assembly_dict, "workspace_name": wsName, "assembly_name": "test_assembly"}
        meta_assembly_ref = self.getImpl().save_assembly_from_fasta(self.ctx, assembly_params)[0]

        # Upload genome, copy genome to workspace folder, & genome data dictionary input
        meta_data = json.load(open(path))
        meta_data['assembly_ref'] = meta_assembly_ref
        meta_dict = [{'name': 'Meta_test',
                  'type': 'KBaseMetagenomes.BinnedContigs',
                  'data': meta_data}]

        # Create .Genome object in workspace with save_objects
        binned_obj = dfu.save_objects({'id': ws_id, 'objects': meta_dict})

        binned_obj_info = binned_obj[0]
        binned_obj_ref = str(binned_obj_info[6]) + '/' + str(binned_obj_info[0]) + '/' + str(binned_obj_info[4])

        # Get FASTA
        ret = self.getImpl().get_fastas(self.callback_url, [binned_obj_ref])

    def test_genome_set_input(self):

        # Setup: copy data file to workspace and get workspace id
        path = "data/TestGenome.json"
        ws_path = '/kb/module/work/tmp'
        shutil.copy2(path, ws_path)
        dfu = DataFileUtil(self.callback_url)
        wsName = self.getWsName()
        ws_id = dfu.ws_name_to_id(wsName)

        # Initiate Dictionaries
        genome_dict, genome_set_dict, dfu_genomeset_dict, dfu_genomeset_dict_2, dfu_genome_search_dict, dfu_genome_search_dict_2 = {}, {}, {}, {}, {}, {}

        # Upload genome & genome data dictionary input
        data = json.load(open(path))
        objs1 = [{'name': 'genome_test',
                  'type': 'KBaseGenomes.Genome',
                  'data': data }]
        # Create .Genome object in workspace with save_objects
        genome_obj = dfu.save_objects({'id': ws_id, 'objects': objs1})

        # Get .Genome object reference
        genome_info = genome_obj[0]
        genome_ref = str(genome_info[6]) + '/' + str(genome_info[0]) + '/' + str(genome_info[4])

        # Create genome object info dictionary
        genome_dict.update({"label": "GenomeSetTest", "ref": genome_ref})

        # Create genome set object dictionary
        genome_set_dict.update({"description": " ", "items": [genome_dict]})

        # Create DataFileUtil dictionaries for genome set data
        dfu_genomeset_dict.update({"type": "KBaseSets.GenomeSet", "data": genome_set_dict,
                                "name": "Genome_Set_Test"})
        dfu_genomeset_dict_2.update({'id': ws_id, 'objects': [dfu_genomeset_dict]})

        # Lastly, create .GenomeSet object with save_objects and get GenomeSet object reference
        genome_set_obj = dfu.save_objects(dfu_genomeset_dict_2)
        genome_set_info = genome_set_obj[0]
        genome_set_ref = str(genome_set_info[6]) + '/' + str(genome_set_info[0]) + '/' + str(genome_set_info[4])

       # Test KBaseSearch.GenomeSet
        genome_set_dict.pop('items', None)
        genome_set_dict['elements'] = {"Set1" : genome_dict}

        # Create DataFileUtil dictionaries for KBaseSearch.GenomeSet data
        dfu_genome_search_dict.update({"type": "KBaseSearch.GenomeSet", "data": genome_set_dict,
                                "name": "Genome_Set_Test_2"})
        dfu_genome_search_dict_2.update({'id': ws_id, 'objects': [dfu_genome_search_dict]})

        # Lastly, create .GenomeSet object with save_objects and get GenomeSet object reference
        search_genome_obj = dfu.save_objects(dfu_genome_search_dict_2)
        search_genome_info = search_genome_obj[0]
        search_set_ref = str(search_genome_info[6]) + '/' + str(search_genome_info[0]) + '/' + str(search_genome_info[4])

        # Get FASTAS
        ret = self.getImpl().get_fastas(self.callback_url, [genome_set_ref, search_set_ref])

    def test_AssemblySet_input(self):

        # Initiate empty data dictionaries and get data_util
        dfu = DataFileUtil(self.callback_url)
        assembly_dict = dict()
        assembly_set_dict = dict()
        dfu_dict = dict()
        dfu_dict_2 = dict()
        # Get workspace id and name
        wsName = self.getWsName()
        ws_id = dfu.ws_name_to_id(wsName)

        # FASTA to assembly object
        Fasta_assembly_dict = {"path": "/kb/module/work/tmp/NC_021490.fasta", "assembly_name": "test_assembly" }
        params = {"file": Fasta_assembly_dict, "workspace_name": wsName, "assembly_name":"test_assembly"}
        ref = self.getImpl().save_assembly_from_fasta(self.ctx, params)

        # Create assembly data dictionaries
        assembly_dict.update({"label":"assemblySetTest", "ref" : ref[0]})
        assembly_set_dict.update({"description": " ", "items" :[assembly_dict]})
        # Create DataFileUtil dictionaries
        dfu_dict.update({"type" : "KBaseSets.AssemblySet", "data": assembly_set_dict,
                         "name": "Assembly_Test"})
        dfu_dict_2.update({'id': ws_id, 'objects': [dfu_dict]})

        # Create assembly set object
        assembly_set_obj = dfu.save_objects(dfu_dict_2)
        assembly_set_ref = [str(assembly_set_obj[0][6]) + '/' + str(assembly_set_obj[0][0]) + '/' + str(assembly_set_obj[0][4])]

        # Get FASTA
        ret = self.getImpl().get_fastas(self.callback_url, assembly_set_ref)

    def test_genome_input(self):
        ref_list = ["27079/16/1"]
        ret = self.getImpl().get_fastas(self.callback_url, ref_list)

    def test_annotations_assembly_input(self):
        ref_list = ["27079/3/1", "23594/10/1"]
        ret = self.getImpl().get_fastas(self.callback_url, ref_list)

    @unittest.skip
    def test_assorted_refs(self):
        ref_list = ["27079/3/1", "44333/4/1", "27079/2/1"]
        ret = self.getImpl().get_fastas(self.callback_url, ref_list)
