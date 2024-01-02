# -*- coding: utf-8 -*-
import os
import time
import unittest
import json
import shutil
from configparser import ConfigParser

from AssemblyUtil.AssemblyUtilImpl import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
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
        cls.testWS = 'KBaseTestData'
        cls.serviceImpl = AssemblyUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'])

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

    def getContext(self):
        return self.__class__.ctx

    def create_Genome_Obj(self, ws_id, genome_obj):

        genome_dict = {}

        # Create .Genome object in workspace with save_objects
        genome_obj = self.dfu.save_objects({'id': ws_id, 'objects': genome_obj})

        # Get .Genome object reference
        genome_info = genome_obj[0]
        genome_ref = str(genome_info[6]) + '/' + str(genome_info[0]) + '/' + str(genome_info[4])

        # GenomeSet info dictionary
        genome_dict.update({"label": "GenomeSetTest", "ref": genome_ref})

        return genome_dict

    def create_KBaseSet_GenomeSetObj(self, ws_id, genome_dict):

        genome_set_dict, dfu_genomeset_dict, dfu_genomeset_dict_2 = {}, {}, {}

        # GenomeSet object dictionary
        genome_set_dict.update({"description": " ", "items": [genome_dict]})

        # DataFileUtil dictionaries to create GenomeSet object
        dfu_genomeset_dict.update({"type": "KBaseSets.GenomeSet", "data": genome_set_dict,
                                   "name": "Genome_Set_Test"})
        dfu_genomeset_dict_2.update({'id': ws_id, 'objects': [dfu_genomeset_dict]})

        # Create .GenomeSet object with save_objects and get GenomeSet object reference
        genome_set_obj = self.dfu.save_objects(dfu_genomeset_dict_2)
        genome_set_info = genome_set_obj[0]
        genome_set_ref = str(genome_set_info[6]) + '/' \
                         + str(genome_set_info[0]) + '/' \
                         + str(genome_set_info[4])

        return (genome_set_ref, genome_set_dict)

    def create_KBaseSearch_GenomeSetObj(self, ws_id, genome_dict, genome_set_dict):

        dfu_genome_search_dict, dfu_genome_search_dict_2 = {}, {}

        # Change GenomeSet dictionary to KBaseSearch.GenomeSet specifications
        genome_set_dict.pop('items', None)
        genome_set_dict['elements'] = {"Set1": genome_dict}

        # DataFileUtil dictionaries to create KBaseSearch.GenomeSet objects
        dfu_genome_search_dict.update({"type": "KBaseSearch.GenomeSet", "data": genome_set_dict,
                                       "name": "Genome_Set_Test_2"})
        dfu_genome_search_dict_2.update({'id': ws_id, 'objects': [dfu_genome_search_dict]})

        # Create .GenomeSet object with save_objects and get GenomeSet object reference
        search_genome_obj = self.dfu.save_objects(dfu_genome_search_dict_2)
        search_genome_info = search_genome_obj[0]
        search_set_ref = str(search_genome_info[6]) + '/' \
                         + str(search_genome_info[0]) + '/' \
                         + str(search_genome_info[4])

        return search_set_ref

    def _assert_inputs(self, ret, ref_lst):
        """assert that all the input references appear in the output"""
        refs = []
        for ref, value in ret.items():
            if value.get('parent_refs'):
                refs += value['parent_refs']
            else:
                refs.append(ref)
        self.assertCountEqual(refs, ref_lst)

    def _assert_outputs(self, ret):
        """assert """
        for ref, value in ret.items():
            for path in value['paths']:
                self.assertTrue(os.path.isfile(path))
            self.assertTrue('type' in value)

    # @unittest.skip('skip')
    def test_genome_set_input(self):
        """ test_genome_set_input tests get_fasta for KBaseSets.GenomeSet and KBaseSearch.GenomeSet objects.
            A genome is saved to the workspace with save_objects and then used to create a GenomeSet.
            The GenomeSet object consists of a single genome: [genome]
            First a KBaseSet.GenomeSet object is made and the reference saved. Then the original genome dictionary
            is modified to fit the type spec of KBaseSearch.GenomeSet; "items" -> "elements"
            With the genome object correct dictionary, an KBaseSearch.Genome object is made an its reference is saved.
            Both object references are placed in an array and fed in get_fasta for testing. """

        # Setup: copy data file to workspace and get workspace id
        wsName = self.getWsName()
        ws_id = self.dfu.ws_name_to_id(wsName)
        ctx = self.getContext()

        path = "data/GenomeSet_TestData/ListeriaMonocytogenes.json"
        ws_path = '/kb/module/work/tmp'
        shutil.copy2(path, ws_path)

        # Upload genome & genome data dictionary input
        data = json.load(open(path))
        genome_obj = [{'name': 'genome_test',
                  'type': 'KBaseGenomes.Genome',
                  'data': data}]

        # Create Genome object from json file
        genome_data_dict = self.create_Genome_Obj(ws_id, genome_obj)
        # Create KBaseSets.GenomeSet object and ref
        GenomeSet_tuple = self.create_KBaseSet_GenomeSetObj(ws_id, genome_data_dict)
        # Seperate ref and dictionary from genomeset data tuple
        GenomeSet_ref = GenomeSet_tuple[0]
        GenomeSet_dict = GenomeSet_tuple[1]
        # Create KBaseSearch.GenomeSet object
        KB_Search_GenomeSet_ref = self.create_KBaseSearch_GenomeSetObj(ws_id, genome_data_dict, GenomeSet_dict)

        # Get FASTAS for KBaseSets.GenomeSet and KBaseSearch.GenomeSet references

        ret = self.getImpl().get_fastas(ctx, {"ref_lst": [GenomeSet_ref, KB_Search_GenomeSet_ref]})[0]
        self._assert_inputs(ret, [GenomeSet_ref, KB_Search_GenomeSet_ref])
        self._assert_outputs(ret)

    # @unittest.skip('skip')
    def test_AssemblySet_input(self):
        """test_AssemblySet_input tests get_fasta for KBaseSets.AssemblySet.
            From theNC_021490.fasta file, a assembly object is saved to the workspace
            and KBaseSets.AssemblySet is made.
            It's reference is then input into get_fasta."""

        # Initiate empty data dictionaries and get workspace ID
        assembly_dict, assembly_set_dict, dfu_dict, dfu_dict_2 = {}, {}, {}, {}
        wsName = self.getWsName()
        ctx = self.getContext()
        ws_id = self.dfu.ws_name_to_id(wsName)
        # Copy data file to workspace
        path = "data/AssemblySet_TestData/NC_021490.fasta"
        ws_path = '/kb/module/work/tmp'
        shutil.copy2(path, ws_path)

        # FASTA to assembly object
        Fasta_assembly_dict = {"path": "/kb/module/work/tmp/NC_021490.fasta", "assembly_name": "test_assembly" }
        params = {"file": Fasta_assembly_dict, "workspace_name": wsName, "assembly_name":"test_assembly"}
        ref = self.getImpl().save_assembly_from_fasta(ctx, params)

        # Create assembly info dictionarary
        assembly_dict.update({"label":"assemblySetTest", "ref" : ref[0]})
        # Create assembly object dictionarary
        assembly_set_dict.update({"description": " ", "items" :[assembly_dict]})

        # Create DataFileUtil dictionaries to build AssemblySet object
        dfu_dict.update({"type" : "KBaseSets.AssemblySet", "data": assembly_set_dict,
                         "name": "Assembly_Test"})
        dfu_dict_2.update({'id': ws_id, 'objects': [dfu_dict]})

        # Create AssemblySet object and get reference
        assembly_set_obj = self.dfu.save_objects(dfu_dict_2)
        assembly_set_ref = [str(assembly_set_obj[0][6]) + '/' \
                            + str(assembly_set_obj[0][0]) + '/' \
                            + str(assembly_set_obj[0][4])]

        # Get FASTA for AssemblySet object
        ret = self.getImpl().get_fastas(ctx, {"ref_lst": assembly_set_ref})[0]
        self._assert_inputs(ret, assembly_set_ref)
        self._assert_outputs(ret)

    # @unittest.skip('skip')
    def test_annotated_metagenome_input(self):
        """"""
        wsName = self.getWsName()
        ret = self.wsClient.copy_object({'from': {'workspace': self.testWS, 'name': 'metagenome.gff_metagenome.assembly.fa_assembly'},
                                         'to': {'workspace': wsName, 'name': 'metagenome.gff_metagenome.assembly.fa_assembly'}})
        upa = '{}/{}/{}'.format(ret[6], ret[0], ret[4])
        ret = self.getImpl().get_fastas(self.getContext(), {'ref_lst':  [upa]})[0]
        self._assert_outputs(ret)

    # @unittest.skip('skip')
    def test_metagenome_binned_input(self):
        """test_metagenome_binned_input tests get_fasta for KBaseMetagenomes.BinnedContigs.
            From the CCESR16 Contig object and assembly, a KBaseMetagenomes.BinnedContigs is made
            and it's reference is input into get_fasta."""

        # Setup & copy data file to workspace
        path = "data/Metagenome_TestData/binnedContigs.json"
        ws_path = '/kb/module/work/tmp'
        assembly_path = "data/Metagenome_TestData/SPAdes_Test.assembly.fa"
        shutil.copy2(path, ws_path)
        shutil.copy2(assembly_path, ws_path)
        wsName = self.getWsName()
        ctx = self.getContext()
        ws_id = self.dfu.ws_name_to_id(wsName)

        # FASTA to assembly object
        Fasta_assembly_dict = {"path": '/kb/module/work/tmp/SPAdes_Test.assembly.fa', "assembly_name": "meta_assembly"}
        assembly_params = {"file": Fasta_assembly_dict, "workspace_name": wsName, "assembly_name": "test_assembly"}
        meta_assembly_ref = self.getImpl().save_assembly_from_fasta(ctx, assembly_params)[0]

        # Upload genome, copy genome to workspace folder, & genome data dictionary input
        meta_data = json.load(open(path))
        meta_data['assembly_ref'] = meta_assembly_ref
        meta_dict = [{'name': 'Meta_test',
                      'type': 'KBaseMetagenomes.BinnedContigs',
                      'data': meta_data}]

        # Create .Genome object in workspace with save_objects
        binned_obj = self.dfu.save_objects({'id': ws_id, 'objects': meta_dict})

        binned_obj_info = binned_obj[0]
        binned_obj_ref = str(binned_obj_info[6]) + '/' + str(binned_obj_info[0]) + '/' + str(binned_obj_info[4])

        # Get FASTA
        ret = self.getImpl().get_fastas(ctx, {"ref_lst": [binned_obj_ref]})[0]
        self._assert_inputs(ret, [binned_obj_ref])
        self._assert_outputs(ret)

    # @unittest.skip('skip')
    def test_genome_input(self):
        wsName = self.getWsName()
        ret = self.wsClient.copy_object(
            {'from': {'workspace': self.testWS, 'name': 'KBase_derived_16_paired_trim_MEGAHIT.contigs.fa_genome.gff_genome'},
             'to': {'workspace': wsName, 'name': 'KBase_derived_16_paired_trim_MEGAHIT.contigs.fa_genome.gff_genome'}})
        upa = '{}/{}/{}'.format(ret[6], ret[0], ret[4])
        ref_list = {"ref_lst": [upa]}
        ret = self.getImpl().get_fastas(self.getContext(), ref_list)[0]
        self._assert_inputs(ret, ref_list['ref_lst'])
        self._assert_outputs(ret)

    # @unittest.skip('skip')
    def test_annotations_assembly_input(self):
        wsName = self.getWsName()
        ret1 = self.wsClient.copy_object(
            {'from': {'workspace': self.testWS,
                      'name': 'archaea_test.fa_assembly.fa_assembly'},
             'to': {'workspace': wsName, 'name': 'archaea_test.fa_assembly.fa_assembly'}})
        upa1 = '{}/{}/{}'.format(ret1[6], ret1[0], ret1[4])
        ret2 = self.wsClient.copy_object(
            {'from': {'workspace': self.testWS,
                      'name': '3300011599_1.fa_assembly.fa_assembly'},
             'to': {'workspace': wsName, 'name': '3300011599_1.fa_assembly.fa_assembly'}})
        upa2 = '{}/{}/{}'.format(ret2[6], ret2[0], ret2[4])
        ref_list = {"ref_lst": [upa1, upa2]}
        ret = self.getImpl().get_fastas(self.getContext(), ref_list)[0]
        self._assert_inputs(ret, ref_list['ref_lst'])
        self._assert_outputs(ret)
