from Bio import AlignIO
from Bio.Align import Alignment
import time
from config import Config
from utils import save_to_json, load_from_json
import os
import logging




#-----------------------------------------------------------------------------------------------------------------------
'''Save filtered alignments to a location as a json file'''
def save_alignments(alignments:dict, config_args: Config):
    
    alignment_filename = config_args.filtered_xmfa_name
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
'''Returns all the sequences (as a list of strings) contained in an alignment'''
'''If all the sequences are identical, it will return only one sequence'''
def get_alignment_sequences(alignment: Alignment, percentage_identity: None):
    
    if percentage_identity==None:
        percentage_identity = 0

    sequences = []

    #If the percentage of identity is 1 (100%) all the sequences are the same. Just return one sequence. 
    if percentage_identity == 1:
        sequences.append(str(alignment[0].seq))
    else:
        for s in alignment:
            sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    
    return sequences


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
    
    logging.info("Filtering alignments...")
    for alignment in alignments:

        alignment_length = compute_alignment_length(alignment=alignment)
        
        if alignment_length >= min_alignment_length:
            alignment_presence = compute_alignment_coverage(alignment=alignment, max_records=ingroup_size)

            if alignment_presence >= min_alignment_coverage:
                alignment_percentage_identity = compute_alignment_percentage_of_identity(alignment)
                #alignment_percentage_identity = 1
                
                if alignment_percentage_identity >= min_alignment_identity:
                    id = alignment[0].id
                    alignment_id = id.split()[0]
                    alignment_number = alignment_id[7:len(alignment_id)]
                    d = {}
                    d["id"] = id
                    d["sequences"]= get_alignment_sequences(alignment=alignment, percentage_identity=alignment_percentage_identity)
                    
                    filtered_alignments[alignment_number] = d

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
    
    success = save_alignments(alignments=filtered_alignments , config_args=config_args)
    

    d =  load_from_json(filename = config_args.filtered_xmfa_path)
    

    end = time.time()
    mins = (end-start)/60
    
    logging.info("Filtering Runtime mins {}".format(mins))


