from Bio import AlignIO
from Bio.Align import Alignment
import time
from config import Config
from utils import save_to_json, load_from_json
import os
import logging
from Bio import SeqIO


def reduce_sequences_in_alignments(alignments:dict):
    
    alignments_with_one_sequence = {}

    for number, alignment in alignments.items():
        sequences = alignment["sequences"]
        new_sequences = []
        for seq in sequences:
            if '-' not in seq:
                new_sequences.append(seq)
                break
        
        if len(new_sequences)!=0:
            d={}
            d["id"] = alignment["id"]
            d["name"] = alignment['name']
            d["strand"] = alignment['strand']
            d["genome_filename"] = alignment.get("genome_filename","")
            d["genome_length"] = alignment.get("genome_length","")
            d["genome_header"] = alignment.get("genome_header","")
            d["sequences"]= new_sequences
            alignments_with_one_sequence[number] = d
        #else:
        #    print("all alignments have '-' on it")

    return alignments_with_one_sequence    
                  

    


#-----------------------------------------------------------------------------------------------------------------------
'''Save filtered alignments to a location as a json file'''
def save_alignments(alignments:dict, config_args: Config, alignment_filename:str):
    
    alignment_filepath = os.path.join(config_args.results_path,alignment_filename)
    success = save_to_json(dictio=alignments, filename=alignment_filepath)
    if success:
        config_args.filtered_xmfa_path = alignment_filepath

    return success


#------------------------------------------------------------------------------------------------------------------------
'''Returns the sequence lenght of an alignment''' 
def compute_alignment_length(alignment:Alignment):
    l=0
    if len(alignment)>0:
        l =len(alignment[0].seq)
    
    return l

#------------------------------------------------------------------------------------------------------------------------
'''Computes how much coverage an alignment has on a set of genomes (i.e., how many genomes are in the alignment)'''
def compute_alignment_coverage(alignment:Alignment, max_records: int):

    coverage_percentage = len(alignment) / max_records

    return coverage_percentage


#------------------------------------------------------------------------------------------------------------------------
'''Computes the percentage of identity of an alignment'''
def compute_alignment_percentage_of_identity(alignment: Alignment):
    total_positions = alignment.get_alignment_length()
    total_identity = 0
    for i in range(total_positions):
        column = alignment[:, i]
        column_set = set(column)
        if len(column_set) == 1 and "-" not in column:
            total_identity = total_identity + 1

    percentage_identity = (total_identity / total_positions) 

    return percentage_identity


#------------------------------------------------------------------------------------------------------------------------
'''Computes the percentage of identity of two sequences'''
def compute_pairwise_percentage_of_identity(sequence1: str, sequence2:str):
    total_positions = min(len(sequence1), len(sequence2))
    total_identity = 0
    for i in range(total_positions):
        column = [sequence1[i], sequence2[i]]
        column_set = set(column)
        if len(column_set) == 1 and "-" not in column:
            total_identity = total_identity + 1

    percentage_identity = (total_identity / total_positions) 

    return percentage_identity

#------------------------------------------------------------------------------------------------------------------------
def compute_average_alignment_percentage_of_identity(alignment:Alignment):

    sequences = []
    for s in alignment:
        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    
    i=0
    sum_percentages = 0 
    total_pairs = 0
    for i in range(0, len(sequences)-1):
        sequence_1 = sequences[i]
        j= i + 1
        while (j<len(sequences)):
            #print("Pairs {},{} \n".format(i,j))
            sequence_2 = sequences[j]
            percentage_of_identity = compute_pairwise_percentage_of_identity(sequence1=sequence_1, sequence2=sequence_2)
            total_pairs = total_pairs + 1 
            sum_percentages = sum_percentages + percentage_of_identity
            j = j + 1
    
    avg_percentage = sum_percentages / total_pairs
    return avg_percentage
        



#------------------------------------------------------------------------------------------------------------------------
'''Returns all the sequences (as a list of strings) contained in an alignment'''
'''If all the sequences are identical, it will return only one sequence'''
def get_alignment_sequences(alignment: Alignment, percentage_identity: None):
    
    #if percentage_identity==None:
    #    percentage_identity = 0

    sequences = []

    #If the percentage of identity is 1 (100%) all the sequences are the same. Just return one sequence. 
    #if percentage_identity == 1:
     #   sequences.append(str(alignment[0].seq))
    #else:
    #    for s in alignment:
    #        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object

    for s in alignment:
        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    


    return sequences



#------------------------------------------------------------------------------------------------------------------------
def load_from_file(filename: str, format:str, type=None):
        success = True
        if os.path.isfile(filename):
            try:
                name = os.path.basename(filename)
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
            except Exception as e:
                logging.info(e)
                success = False
        else:
            logging.info ("Filename %s does not exists.", filename)
            success = False
        return seq_records


#------------------------------------------------------------------------------------------------------------------------
def load_file_1(filename: str):
    success = True
    if os.path.isfile(filename):
        try:
            file = open(filename)
            line = file.readline()
            while (line[0]!='>'):
                try:
                    line = file.readline()
                    print(line)
                except Exception as e:
                    print(e)        
        except Exception as e:
            print(e)
    else:
        logging.info ("Filename %s does not exists.", filename)
        success = False
        
    file.close()
    
    return success


#------------------------------------------------------------------------------------------------------------------------
def fill_dict(sequence_file: str, sequence_header:str, sequence_length:str):
    d = {}
    d = {"sequence_file" :sequence_file,
            "sequence_header" :sequence_header,
            "sequence_length" : sequence_length}
    return d          

#------------------------------------------------------------------------------------------------------------------------
def load_file(filename: str):
    
    sequence_info = {}
 
    if os.path.isfile(filename):
        try:
            file = open(filename)
            sequence_index = 0
            sequence_file = ""
            sequence_header = ""
            sequence_length= ""

            for line in file:
                #print(line)
                #print(sequence_index)
                if (line[0]=='>'):
                    #We finished reading all the sequences info
                    #Fill the info of the last sequence
                    d = fill_dict(sequence_file=sequence_file,
                                sequence_header=sequence_header,
                                sequence_length=sequence_length)
                    sequence_info[sequence_index] = d
                          
                    break
                else:
                    if "##SequenceIndex" in line:
                        d = {}
                        if sequence_index >= 1:
                            d = fill_dict(sequence_file=sequence_file,
                                          sequence_header=sequence_header,
                                          sequence_length=sequence_length)
                            sequence_info[sequence_index] = d
                            sequence_file = ""
                            sequence_header = ""
                            sequence_length= ""
                        
                        sequence_index = sequence_index + 1

                    if '##SequenceFile' in line:
                        sequence_file = line[len('##SequenceFile'):len(line)-1]
                    if '##SequenceHeader' in line:
                        sequence_header = line[len('##SequenceHeader'):len(line)-1]
                    if '##SequenceLength' in line:
                        sequence_length = line[len("##SequenceLength"):len(line)-1]
        except Exception as e:
            print(e)
    else:
        logging.info ("Filename %s does not exists.", filename)
        
    file.close()
    
    return sequence_info



#------------------------------------------------------------------------------------------------------------------------
'''Filters all alignments on a xmfa file based on three parameters: elenght, coverage, and identity'''
'''Returns the filtered alignments as a dictionary'''
'''Saves the alignments as a json file into a location'''
def filter_alignments(alignments_file: str, min_alignment_length: int, min_alignment_coverage: int, min_alignment_identity: float, ingroup_size:int):

    filtered_alignments = {}
    
    try:
        align = AlignIO.parse(alignments_file, "mauve")
        alignments = list(align)
    except Exception as e:
        logging.info (e)
        return
    
    try:
        alignment_info = load_file(alignments_file)
    except Exception as e:
        logging.info (e)
    
    
    
    count = 1

    #logging.info("Filtering using pairwise \n")

    for alignment in alignments:
        #logging.info("Filtering cluster {} of {}".format(count,len(alignments)))
        alignment_length = compute_alignment_length(alignment=alignment)
        
        if alignment_length >= min_alignment_length:
            alignment_presence = compute_alignment_coverage(alignment=alignment, max_records=ingroup_size)

            if alignment_presence >= min_alignment_coverage:
                alignment_percentage_identity = compute_alignment_percentage_of_identity(alignment)

                #alignment_percentage_identity = compute_average_alignment_percentage_of_identity(alignment=alignment)

                #alignment_percentage_identity = 1
                
                if alignment_percentage_identity >= min_alignment_identity:
                    id = alignment[0].id
                    alignment_id = id.split()[0]
                    alignment_number = alignment_id[7:len(alignment_id)]
                    d = {}
                    d["id"] = id
                    d["name"] = alignment[0].name
                    d["strand"] = alignment[0].annotations['strand']
                    d["sequences"]= get_alignment_sequences(alignment=alignment, percentage_identity=alignment_percentage_identity)
                    
                    try:
                        sequence_index = int(d.get("name",0))
                        if sequence_index in alignment_info:
                            a = alignment_info[sequence_index]
                            d["genome_filename"] = a.get("sequence_file","")
                            d["genome_length"] = a.get("sequence_length","")
                            d["genome_header"] = a.get("sequence_header","")
                    except:
                            d["genome_filename"] = ""
                            d["genome_length"] = ""
                            d["genome_header"] = ""
                        
                        
                    filtered_alignments[alignment_number] = d
    
        count = count + 1

    #---------------------------------------------------------
    #This could be useful for a class
    logging.info ("Total alignments: {}".format(len(alignments)))
    logging.info ("Filtered Alignments: {}".format(len(filtered_alignments)))
    logging.info ("Percentage kept: {}".format (len(filtered_alignments)/len(alignments)*100))

    return filtered_alignments


#----------------------------------------
'''Runs the filter alignments'''
def run_filter(config_args: Config):

    filtered_alignments = {}
    start = time.time()
    
    
    filtered_alignments = filter_alignments(alignments_file=config_args.xmfa_file_path, 
                                            min_alignment_length = config_args.minimum_alignment_length, 
                                            min_alignment_coverage = config_args.minimum_alignment_coverage,
                                            min_alignment_identity =  config_args.minimum_alignment_percentage_identity,
                                            ingroup_size=config_args.ingroup_size)
    
    success = save_alignments(alignments=filtered_alignments, config_args=config_args, alignment_filename = config_args.filtered_xmfa_name)
    print_stats(filtered_alignments=filtered_alignments, config=config_args)
    

    #d =  load_from_json(filename = config_args.filtered_xmfa_path)

    filtered_alignments_reduced = reduce_sequences_in_alignments(alignments=filtered_alignments)
    success = save_alignments(alignments=filtered_alignments_reduced, config_args=config_args, alignment_filename = config_args.reduced_filtered_xmfa_name)
    
    end = time.time()
    mins = (end-start)/60
    
    logging.info("Filtering Runtime mins {}".format(mins))
    return success


#----------------------------------------
def print_stats(filtered_alignments:dict, config):
    
    max = 0
    file = open(config.results_path+'/lenghts.txt','w+') 
    for alignment in filtered_alignments.items():
        for sequence in alignment[1]["sequences"]:
            if len(sequence)>max:
                max = len(sequence)
        file.write(str(max))
        file.write('\n')
        max = 0
    file.close()




