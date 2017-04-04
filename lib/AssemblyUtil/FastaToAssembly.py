import os
import sys
import uuid

from pprint import pprint
from hashlib import md5
from collections import Counter

from Bio import SeqIO

from DataFileUtil.DataFileUtilClient import DataFileUtil


class FastaToAssembly:

    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)

        # Note added X due to kb|g.1886.fasta
        self.valid_chars = "-ACGTUWSMKRYBDHVNX"
        self.amino_acid_specific_characters = "PLIFQE"


    def import_fasta(self, ctx, params):
        print('validating parameters')
        self.validate_params(params)

        print('staging input files')
        fasta_file_path = self.stage_input(params)

        if 'min_contig_length' in params:
            min_contig_length = int(params['min_contig_length'])
            print('filtering fasta file by contig length (min len=' + str(min_contig_length) + 'bp)')
            fasta_file_path = self.filter_contigs_by_length(fasta_file_path, min_contig_length)

        assembly_data = self.parse_fasta(fasta_file_path, params)
        print('parsing FASTA file: ' + str(fasta_file_path))
        pprint(assembly_data)


        print('saving assembly to KBase')

        # save file to shock and build handle
        fasta_file_handle_info = self.save_fasta_file_to_shock(fasta_file_path)
        pprint(fasta_file_handle_info)

        # construct the output object
        assembly_object_to_save = self.build_assembly_object(assembly_data,
                                                             fasta_file_handle_info,
                                                             params)

        # save to WS and return
        if 'workspace_id' in params:
            workspace_id = int(params['workspace_id'])
        else:
            workspace_id = self.dfu.ws_name_to_id(params['workspace_name'])
        assembly_info = self.save_assembly_object(workspace_id,
                                                  params['assembly_name'],
                                                  assembly_object_to_save)

        return assembly_info


    def build_assembly_object(self, assembly_data, fasta_file_handle_info, params):
        ''' construct the WS object data to save based on the parsed info and params '''
        assembly_data['assembly_id'] = params['assembly_name']
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
            assembly_data['external_source_origination_date'] = params['external_source_origination_date']

        return assembly_data


    def parse_fasta(self, fasta_file_path, params):
        ''' Do the actual work of inspecting each contig '''

        # variables to store running counts of things
        total_length = 0
        base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        md5_list = []

        # map from contig_id to contig_info
        all_contig_data = {}

        for record in SeqIO.parse(fasta_file_path, "fasta"):
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
                'length': len(record.seq),
                'is_circular': 'Unknown',


                'num_bytes': 0,
                'start_position': 0,
                'end_position': 0
            }

            # 1) compute sequence character statistics running total
            total_length += contig_info['length']
            sequence_count_table = dict(Counter(sequence))
            for character in sequence_count_table:
                if character in base_counts:
                    base_counts[character] = base_counts[character] + sequence_count_table[character]
                else:
                    base_counts[character] = sequence_count_table[character]
                if character not in self.valid_chars:
                    if character in self.amino_acid_specific_characters:
                        raise ValueError('This fasta file may have amino acids in it instead ' +
                                         'of the required nucleotides.')
                    raise ValueError("This FASTA file has non nucleic acid characters : {0}".format(character))

            # 2) record number of 'N' characters
            Ncount = 0
            if 'N' in sequence_count_table:
                Ncount = sequence_count_table['N']
            contig_info['Ncount'] = Ncount

            # 3) record md5 checksum
            contig_md5 = md5(sequence).hexdigest()
            contig_info['md5'] = contig_md5
            md5_list.append(contig_md5)

            # 4) record the all important GC to ~3 significant digits
            GC_count = 0
            for base in ['G', 'C']:
                if base in sequence_count_table:
                    GC_count += sequence_count_table[base]
            contig_info['gc_content'] = round(float(GC_count) / float(contig_info['length']), 3)

            # 5) add to contig list
            if contig_info['contig_id'] in all_contig_data:
                raise ValueError('The fasta header key ' + contig_info['contig_id'] +
                                 'appears more than once in the file')
            all_contig_data[contig_info['contig_id']] = contig_info

        # Aggregate stats for the data
        assembly_data = {
            'md5': md5(",".join(sorted(md5_list))).hexdigest(),
            'base_counts': base_counts,
            'dna_size': total_length,
            'gc_content': round(float(GC_count) / float(contig_info['length']), 3),
            'contigs': all_contig_data,
            'num_contigs': len(all_contig_data)
        }
        return assembly_data



    def filter_contigs_by_length(self, fasta_file_path, min_contig_length):
        ''' removes all contigs less than the min_contig_length provided '''
        filtered_fasta_file_path = fasta_file_path + '.filtered.fa'

        fasta_record_iter = SeqIO.parse(fasta_file_path, 'fasta')
        short_seq_iterator = (record for record in fasta_record_iter if len(record.seq) >= min_contig_length)
        SeqIO.write(short_seq_iterator, filtered_fasta_file_path, 'fasta')

        return filtered_fasta_file_path



#     contig_set_dict["external_source"] = source
#     contig_set_dict["external_source_id"] = os.path.basename(fasta_file_name) 
# #    contig_set_dict["external_source_origination_date"] = str(os.stat(fasta_file_name).st_ctime)

#     if date_string is not None:
#         contig_set_dict["external_source_origination_date"] = date_string
#     contig_set_dict["contigs"] = fasta_dict
#     contig_set_dict["dna_size"] = total_length
#     contig_set_dict["gc_content"] = float(gc_length)/float(total_length)
# #    print "Fasta dict Keys :"+",".join(fasta_dict.keys())+":" 
#     contig_set_dict["num_contigs"] = len(fasta_dict.keys())
#     contig_set_dict["type"] = "Unknown"
#     contig_set_dict["notes"] = "Note MD5s are generated from uppercasing the sequences" 
#     contig_set_dict["base_counts"] = base_counts 
#     if taxon_reference is not None:
#         contig_set_dict["taxon_ref"] = taxon_reference





#     input_file_handle = TextFileDecoder.open_textdecoder(fasta_file_name, 'ISO-8859-1')    
#     fasta_header = None
#     fasta_description = None
#     sequence_list = []
#     fasta_dict = dict()
#     first_header_found = False
#     contig_set_md5_list = []
#     # Pattern for replacing white space
#     pattern = re.compile(r'\s+')
#     sequence_exists = False
    
#     total_length = 0
#     gc_length = 0
#     #Note added X and x due to kb|g.1886.fasta
#     valid_chars = "-AaCcGgTtUuWwSsMmKkRrYyBbDdHhVvNnXx"
#     amino_acid_specific_characters = "PpLlIiFfQqEe" 

#     #Base_counts - is dict of base characters and their counts.
#     base_counts = dict()

#     sequence_start = 0
#     sequence_stop = 0

#     current_line = input_file_handle.readline()
#     while current_line != None and len(current_line) > 0:
# #        print "CURRENT LINE: " + current_line
#         if (current_line[0] == ">"):
#             # found a header line
#             # Wrap up previous fasta sequence
#             if (not sequence_exists) and first_header_found:
#                 logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header))        
#                 raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header))
#             if not first_header_found:
#                 first_header_found = True
#                 sequence_start = 0
#             else:
#                 sequence_stop = input_file_handle.tell() - len(current_line)
#                 # build up sequence and remove all white space
#                 total_sequence = ''.join(sequence_list)
#                 total_sequence = re.sub(pattern, '', total_sequence)
#                 if not total_sequence :
#                     logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header)) 
#                     raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header))
# #                for character in total_sequence:
# #                    if character not in valid_chars:
# #                        if character in amino_acid_specific_characters:
# #                            raise Exception("This fasta file may have amino acids in it instead of the required nucleotides.")
# #                        raise Exception("This FASTA file has non nucleic acid characters : {0}".format(character))
#                 seq_count = collections.Counter(total_sequence.upper())
#                 seq_dict = dict(seq_count)
#                 for character in seq_dict:
#                     if character in base_counts:
#                         base_counts[character] =  base_counts[character] + seq_dict[character]
#                     else:
#                         base_counts[character] =  seq_dict[character]
#                     if character not in valid_chars:
#                         if character in amino_acid_specific_characters:
#                             raise Exception("This fasta file may have amino acids in it instead of the required nucleotides.")
#                         raise Exception("This FASTA file has non nucleic acid characters : {0}".format(character))

#                 contig_dict = dict() 
#                 Ncount = 0
#                 if "N" in seq_dict:
#                     Ncount = seq_dict["N"]
#                 contig_dict["Ncount"] = Ncount 
#                 length = len(total_sequence)
#                 total_length = total_length + length
#                 contig_gc_length = len(re.findall('G|g|C|c',total_sequence))
#                 contig_dict["gc_content"] = float(contig_gc_length)/float(length) 
#                 gc_length = gc_length + contig_gc_length
#                 fasta_key = fasta_header.strip()
#                 contig_dict["contig_id"] = fasta_key 
#                 contig_dict["length"] = length 
#                 contig_dict["name"] = fasta_key 
#                 contig_md5 = hashlib.md5(total_sequence.upper()).hexdigest() 
#                 contig_dict["md5"] = contig_md5 
#                 contig_set_md5_list.append(contig_md5)

#                 contig_dict["is_circular"] = "Unknown"
#                 if fasta_description is not None: 
#                     contig_dict["description"] = fasta_description
#                 if contig_information_dict is not None:
#                     if contig_information_dict[fasta_key] is not None:
#                         if contig_information_dict[fasta_key]["definition"] is not None:
#                             contig_dict["description"] = contig_information_dict[fasta_key]["definition"]
#                         if contig_information_dict[fasta_key]["is_circular"] is not None:
#                             contig_dict["is_circular"] = contig_information_dict[fasta_key]["is_circular"]
#                 contig_dict["start_position"] = sequence_start
#                 contig_dict["num_bytes"] = sequence_stop - sequence_start

# #                print "Sequence Start: " + str(sequence_start) + "Fasta: " + fasta_key
# #                print "Sequence Stop: " + str(sequence_stop) + "Fasta: " + fasta_key

#                 if fasta_key in fasta_dict:
#                     raise Exception("The fasta header {0} appears more than once in the file ".format(fasta_key))
#                 else: 
#                     fasta_dict[fasta_key] = contig_dict
               
#                 # get set up for next fasta sequence
#                 sequence_list = []
#                 sequence_exists = False
                
# #               sequence_start = input_file_handle.tell()               
#             sequence_start = 0            

#             fasta_header_line = current_line.strip().replace('>','')
#             try:
#                 fasta_header , fasta_description = fasta_header_line.split(' ',1)
#             except:
#                 fasta_header = fasta_header_line
#                 fasta_description = None
#         else:
#             if sequence_start == 0:
#                 sequence_start = input_file_handle.tell() - len(current_line) 
#             sequence_list.append(current_line)
#             sequence_exists = True
#         current_line = input_file_handle.readline()
# #        print "ENDING CURRENT LINE: " + current_line

#     # wrap up last fasta sequence
#     if (not sequence_exists) and first_header_found: 
#         logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header))        
#         raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header)) 
#     elif not first_header_found :
#         logger.error("There are no contigs in this file") 
#         raise Exception("There are no contigs in this file") 
#     else: 
#         sequence_stop = input_file_handle.tell()
#         # build up sequence and remove all white space      
#         total_sequence = ''.join(sequence_list)
#         total_sequence = re.sub(pattern, '', total_sequence)
#         if not total_sequence :
#             logger.error("There is no sequence related to FASTA record : {0}".format(fasta_header)) 
#             raise Exception("There is no sequence related to FASTA record : {0}".format(fasta_header)) 

# #        for character in total_sequence: 
#         seq_count = collections.Counter(total_sequence.upper()) 
#         seq_dict = dict(seq_count) 
#         for character in seq_dict:
#             if character in base_counts:
#                 base_counts[character] =  base_counts[character] + seq_dict[character]
#             else:
#                 base_counts[character] =  seq_dict[character]
#             if character not in valid_chars: 
#                 if character in amino_acid_specific_characters:
#                     raise Exception("This fasta file may have amino acids in it instead of the required nucleotides.")
#                 raise Exception("This FASTA file has non nucleic acid characters : {0}".format(character))

#         contig_dict = dict() 
#         Ncount = 0
#         if "N" in seq_dict:
#             Ncount = seq_dict["N"]
#         contig_dict["Ncount"] = Ncount 
#         length = len(total_sequence)
#         total_length = total_length + length
#         contig_gc_length = len(re.findall('G|g|C|c',total_sequence))
#         contig_dict["gc_content"] = float(contig_gc_length)/float(length) 
#         gc_length = gc_length + contig_gc_length
#         fasta_key = fasta_header.strip()
#         contig_dict["contig_id"] = fasta_key 
#         contig_dict["length"] = length
#         contig_dict["name"] = fasta_key

#         contig_dict["is_circular"] = "Unknown"
#         if fasta_description is not None:
#             contig_dict["description"] = fasta_description
#         if contig_information_dict is not None: 
#             if contig_information_dict[fasta_key] is not None:
#                 if contig_information_dict[fasta_key]["definition"] is not None:
#                     contig_dict["description"] = contig_information_dict[fasta_key]["definition"]
#                 if contig_information_dict[fasta_key]["is_circular"] is not None:
#                     contig_dict["is_circular"] = contig_information_dict[fasta_key]["is_circular"]
#         contig_md5 = hashlib.md5(total_sequence.upper()).hexdigest()
#         contig_dict["md5"]= contig_md5
#         contig_set_md5_list.append(contig_md5)
#         contig_dict["start_position"] = sequence_start
#         contig_dict["num_bytes"] = sequence_stop - sequence_start
        
#         if fasta_key in fasta_dict:
#             raise Exception("The fasta header {0} appears more than once in the file ".format(fasta_key))
#         else: 
#             fasta_dict[fasta_key] = contig_dict
#         input_file_handle.close()

#     contig_set_dict = dict()
#     contig_set_dict["md5"] = hashlib.md5(",".join(sorted(contig_set_md5_list))).hexdigest()
#     contig_set_dict["assembly_id"] = assembly_name
#     contig_set_dict["name"] = assembly_name
#     contig_set_dict["external_source"] = source
#     contig_set_dict["external_source_id"] = os.path.basename(fasta_file_name) 
# #    contig_set_dict["external_source_origination_date"] = str(os.stat(fasta_file_name).st_ctime)

#     if date_string is not None:
#         contig_set_dict["external_source_origination_date"] = date_string
#     contig_set_dict["contigs"] = fasta_dict
#     contig_set_dict["dna_size"] = total_length
#     contig_set_dict["gc_content"] = float(gc_length)/float(total_length)
# #    print "Fasta dict Keys :"+",".join(fasta_dict.keys())+":" 
#     contig_set_dict["num_contigs"] = len(fasta_dict.keys())
#     contig_set_dict["type"] = "Unknown"
#     contig_set_dict["notes"] = "Note MD5s are generated from uppercasing the sequences" 
#     contig_set_dict["base_counts"] = base_counts 
#     if taxon_reference is not None:
#         contig_set_dict["taxon_ref"] = taxon_reference









    def save_assembly_object(self, workspace_id, assembly_name, obj_data):
        print('Saving Assembly to Workspace')
        sys.stdout.flush()
        obj_info = self.dfu.save_objects({'id': workspace_id,
                                          'objects': [{'type': 'KBaseGenomeAnnotations.Assembly',
                                                       'data': obj_data,
                                                       'name': assembly_name
                                                       }]
                                          })[0]
        return obj_info


    def save_fasta_file_to_shock(self, fasta_file_path):
        ''' Given the path to the file, upload to shock and return Handle information
            returns:
                typedef structure {
                    string shock_id;
                    Handle handle;
                    string node_file_name;
                    string size;
                } FileToShockOutput;

        '''
        print('Uploading fasta file (' + str(fasta_file_path) + ') to SHOCK')
        sys.stdout.flush()
        return self.dfu.file_to_shock({'file_path': fasta_file_path, 'make_handle': 1})


    def stage_input(self, params):
        ''' Setup the input_directory by fetching the files and returning the path to the file'''
        file_path = None
        if 'file' in params:
            file_path = os.path.abspath(params['file']['path'])
        if 'shock_id' in params:
            print('Downloading file from SHOCK node: ' + str(params['shock_id']))
            sys.stdout.flush()
            input_directory = os.path.join(self.scratch, 'assembly-upload-staging-' + str(uuid.uuid4()))
            os.makedirs(input_directory)
            file_name = self.dfu.shock_to_file({'file_path': input_directory,
                                                'shock_id': params['shock_id']
                                                })['node_file_name']
            file_path = os.path.join(input_directory, file_name)
        if 'ftp_url' in params:
            print('Downloading file from: ' + str(params['ftp_url']))
            sys.stdout.flush()
            file_path = self.dfu.download_web_file({'file_url': params['ftp_url'],
                                                    'download_type': 'FTP'
                                                    })['copy_file_path']

        # extract the file if it is compressed
        if file_path is not None:
            unpacked_file = self.dfu.unpack_file({'file_path': file_path})
            return unpacked_file['file_path']

        raise ValueError('No valid fasta could be extracted based on the input parameters')





    def validate_params(self, params):
        for key in ('workspace_name', 'assembly_name'):
            if key not in params:
                raise ValueError('required "' + key + '" field was not defined')

        # one and only one of either 'file', 'shock_id', or ftp_url is required
        input_count = 0
        for key in ('file', 'shock_id', 'ftp_url'):
            if key in params and params[key] is not None:
                input_count = input_count + 1
                if key == 'file':
                    if not isinstance(params[key], dict) or 'path' not in params[key]:
                        raise ValueError('when specifying a fasta file input, "path" field was not defined in "file"')

        if input_count == 0:
            raise ValueError('required fasta file as input, set as either "file", "shock_id", or "ftp_url"')
        if input_count > 1:
            raise ValueError('required exactly one fasta file as input source, you set more than one of ' +
                             'these fields: "file", "shock_id", or "ftp_url"')















#         # transform scripts assume the file gets dumped in its own unique directory, so create that special
#         # directory and move the file there temporarily
#         input_directory =  os.path.join(self.sharedFolder, 'assembly-upload-staging-'+str(uuid.uuid4()))
#         os.makedirs(input_directory)

#         fasta_file_path = None
#         if 'file' not in params:
#             if 'shock_id' not in params:
#                 if 'ftp_url' not in params:
#                     raise ValueError('No input file (either file.path, shock_id, or ftp_url) provided')
#                 else:
#                     # TODO handle ftp - this creates a directory for us, so update the input directory
#                     print('calling Transform download utility: script_utils.download');
#                     print('URL provided = '+params['ftp_url']);
#                     script_utils.download_from_urls(
#                             working_directory = input_directory,
#                             token = ctx['token'], # not sure why this requires a token to download from a url...
#                             urls  = {
#                                         'ftpfiles': params['ftp_url']
#                                     }
#                         );
#                     input_directory = os.path.join(input_directory,'ftpfiles')
#                     # unpack everything in input directory
#                     dir_contents = os.listdir(input_directory)
#                     print('downloaded directory listing:')
#                     pprint(dir_contents)
#                     dir_files = []
#                     for f in dir_contents:
#                         if os.path.isfile(os.path.join(input_directory, f)):
#                             dir_files.append(f)

#                     print('processing files in directory...')
#                     for f in dir_files:
#                         # unpack if needed using the standard transform utility
#                         print('unpacking '+f)
#                         script_utils.extract_data(filePath=os.path.join(input_directory,f))

#             else:
#                 # handle shock file
#                 dfUtil = DataFileUtil(self.callback_url)
#                 file_name = dfUtil.shock_to_file({
#                                     'file_path': input_directory,
#                                     'shock_id': params['shock_id']
#                                 })['node_file_name']
#                 fasta_file_path = os.path.join(input_directory, file_name)
#         else:
#             # copy the local file to the input staging directory
#             # (NOTE: could just move it, but then this method would have the side effect of moving your
#             # file which another SDK module might have an open handle on - could also get around this
#             # by having the script take a file list directly rather working in a directory)
#             if 'path' not in params['file']:
#                 raise ValueError('file.path field was not defined')
#             local_file_path = params['file']['path']
#             fasta_file_path = os.path.join(input_directory, os.path.basename(local_file_path))
#             shutil.copy2(local_file_path, fasta_file_path)

#         if fasta_file_path is not None:
#             print("input fasta file =" + fasta_file_path)

#             # unpack if needed using the standard transform utility
#             script_utils.extract_data(filePath=fasta_file_path)

#         print("Handle URL: " + self.handleURL)
#         # do the upload
#         result = uploader.upload_assembly(
#                 logger=None,

#                 shock_service_url = self.shockURL,
#                 handle_service_url = self.handleURL,
#                 workspace_service_url = self.workspaceURL,

#                 input_directory=input_directory,
#                 workspace_name=params['workspace_name'],
#                 assembly_name=params['assembly_name'],

#                 provenance=ctx.provenance()
#             )

# #                    taxon_reference = None, 
# #                    source = None, 
# #                    date_string = None,
# #                    contig_information_dict = None,
# #                    logger = None):

#         # clear the temp directory
#         shutil.rmtree(input_directory)

#         # get WS metadata to return the reference to the object
#         ws = Workspace(url=self.workspaceURL)
#         info = ws.get_object_info_new({'objects':[{'ref':params['workspace_name'] + '/' + result}],'includeMetadata':0, 'ignoreErrors':0})[0]

#         ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
