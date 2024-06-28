from Bio import AlignIO
from Bio.Align import Alignment
import time
from config import Config
from utils import save_to_json, load_from_json
import os
import logging
from Bio import SeqIO
import concurrent.futures
import math


def reduce_sequences_in_alignments(alignments:dict):
    
    alignments_with_one_sequence = {}

    for number, alignment in alignments.items():
        sequences = alignment["sequences"]
        new_sequences = []
        for seq in sequences:
            if '-' not in seq:
                new_sequences.append(seq)
                break
        
        if len(new_sequences)==0:
            seq = sequences[0]
            new_seq = seq.replace("-","")
            new_sequences.append(new_seq)
        
        if len(new_sequences)!=0:
            d={}
            d["id"] = alignment["id"]
            d["name"] = alignment['name']
            d["strand"] = alignment['strand']
            d["genome_filename"] = alignment.get("genome_filename","")
            d["genome_length"] = alignment.get("genome_length","")
            d["genome_header"] = alignment.get("genome_header","")
            d["percentage_identity"] = alignment.get("percentage_identity","")
            d["sequences"]= new_sequences
            alignments_with_one_sequence[number] = d
       
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
'''Computes the percentage of identity of two sequences (alternative)'''
def compute_pairwise_percentage_of_identity(sequence1: str, sequence2:str):
    
    length = len(sequence1)

    sum = 0
    for i in range(0,length):
        if sequence1[i]==sequence2[i] and sequence1[i]!="-":
            sum = sum + 1 

    percentage_identity  = sum / length

    return percentage_identity


#------------------------------------------------------------------------------------------------------------------------
def compute_average_alignment_percentage_of_identity_sequential(alignment:Alignment):
    
    sequences = []
    for s in alignment:
        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    
    n = len(sequences)
    total_identity = 0
    count = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            identity = compute_pairwise_percentage_of_identity(sequence1=sequences[i], sequence2=sequences[j])
            total_identity = total_identity + identity
            count = count + 1
    
    average_identity = total_identity / count if count > 0 else 0


    
    return average_identity



#------------------------------------------------------------------------------------------------------------
def __compute_average_percentage_of_identity(sequences:list, processor_pairs:list, config_args:Config):
    
    total_identity = 0
    for pair in processor_pairs:
        i = pair[0]
        j = pair[1]
        identity = compute_pairwise_percentage_of_identity(sequence1=sequences[i], sequence2=sequences[j])
        total_identity = total_identity + identity
        
    return total_identity




#------------------------------------------------------------------------------------------------------------
def compute_average_alignment_percentage_of_identity_parallel(alignment:Alignment, config_args:Config, executor):
    
    sequences = []
    for s in alignment:
        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    
    n = len(sequences) 
    total_pairs = math.comb(n,2)
    pairs_per_processor = math.floor(total_pairs / config_args.processors_number)
    
    if pairs_per_processor<1:
        pairs_per_processor = 1

    processors_jobs = []
    pair_list = []
    count = 0
    processor_index = 1
    for i in range(n):
        for j in range(i + 1, n):
            pair = (i,j)
            pair_list.append(pair)
            count = count + 1
            if count == pairs_per_processor and processor_index<config_args.processors_number:
                processors_jobs.append(pair_list)
                count = 0
                pair_list = []
                processor_index = processor_index + 1
    

    processors_jobs.append(pair_list)
    processor_results = []

    future_results = [executor.submit(__compute_average_percentage_of_identity, sequences, processor_pairs, config_args) for processor_pairs in processors_jobs]   
              
    for finished in concurrent.futures.as_completed(future_results, timeout=600):
        try:
            processor_results.append(finished.result())
        except concurrent.futures._base.TimeoutError:
            logging.error("Process took to long to complete")
        except Exception as exc:
            logging.error("Exception occurred")
            logging.error(exc)

    sum = 0
    for result in processor_results:
        sum = sum + result

    average_percentage_identity = sum / total_pairs

    return average_percentage_identity
    
#------------------------------------------------------------------------------------------------------------------------
'''Returns all the sequences (as a list of strings) contained in an alignment'''
'''If all the sequences are identical, it will return only one sequence'''
def get_alignment_sequences(alignment: Alignment, percentage_identity: None):
    

    sequences = []
    
    for s in alignment:
        sequences.append(str(s.seq))  #If more information is needed, store the whole sequence object
    


    return sequences



#------------------------------------------------------------------------------------------------------------------------
def load_from_file(filename: str, format:str, type=None):
        seq_records=None
        if os.path.isfile(filename):
            try:
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
            except Exception as e:
                logging.info(e)
        else:
            logging.info ("Filename %s does not exists.", filename)
        return seq_records



#------------------------------------------------------------------------------------------------------------------------
def fill_dict(sequence_file: str, sequence_header:str, sequence_length:str):
    d = {}
    d = {"sequence_file" :sequence_file,
            "sequence_header" :sequence_header,
            "sequence_length" : sequence_length}
    return d          

#------------------------------------------------------------------------------------------------------------------------
def parse_metadata_xmfa(filename: str):
    
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
def filter_alignments_parallel(alignments_file: str, min_alignment_length: int, min_alignment_coverage: int, min_alignment_identity: float, ingroup_size:int, config_args:Config):

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=config_args.processors_number)

    filtered_alignments = {}
    
    #Parse the alignments file (xmfa) returned by parsnp to 
    #obtain alignments 
    try:
        align = AlignIO.parse(alignments_file, "mauve")
        alignments = list(align)
        config_args.stats.alignments_found_by_parsnp = len(alignments)
    except Exception as e:
        logging.info (e)
        return
    
    #Parse the alignments file (xmfa) returned by parsnp to
    #obtain metadata (number of genomes, filename, strand, etc)
    try:
        alignment_info = parse_metadata_xmfa(alignments_file)
    except Exception as e:
        logging.info (e)
    
    
    
    count = 1
    alignments_kept_by_length = 0
    alignments_kept_by_coverage = 0
    alignments_kept_by_identity = 0
    

    alignments_discarded_by_length = 0
    alignments_discarded_by_coverage = 0
    alignments_discarded_by_identity = 0
    
   

    for alignment in alignments:
        logging.info("Filtering cluster  {} of {}".format(count,len(alignments)))
        
        alignment_length = compute_alignment_length(alignment=alignment)
        
        if alignment_length >= min_alignment_length:
            alignments_kept_by_length = alignments_kept_by_length + 1

            alignment_presence = compute_alignment_coverage(alignment=alignment, max_records=ingroup_size)

            if alignment_presence >= min_alignment_coverage:
                alignments_kept_by_coverage = alignments_kept_by_coverage + 1
                
                alignment_percentage_identity = compute_average_alignment_percentage_of_identity_parallel(alignment=alignment, config_args=config_args, executor=executor)
                
                if alignment_percentage_identity >= min_alignment_identity:
                    alignments_kept_by_identity = alignments_kept_by_identity + 1

                    id = alignment[0].id
                    alignment_id = id.split()[0]
                    alignment_number = alignment_id[7:len(alignment_id)]
                    d = {}  
                    d["id"] = id
                    d["name"] = alignment[0].name
                    d["strand"] = alignment[0].annotations['strand']
                    d["sequences"]= get_alignment_sequences(alignment=alignment, percentage_identity=alignment_percentage_identity)
                    d["percentage_identity"] = alignment_percentage_identity

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
            
                else:
                    alignments_discarded_by_identity = alignments_discarded_by_identity + 1
       
            else:
                alignments_discarded_by_coverage = alignments_discarded_by_coverage + 1
        else:
            alignments_discarded_by_length = alignments_discarded_by_length + 1
    
        count = count + 1

    #---------------------------------------------------------
    #This could be useful for a class
    logging.info ("Total alignments: {}".format(len(alignments)))
    logging.info ("Filtered Alignments: {}".format(len(filtered_alignments)))
    logging.info ("Percentage kept: {}".format (len(filtered_alignments)/len(alignments)*100))
    logging.info("Discarded by identity {}".format(alignments_discarded_by_identity))


    config_args.stats.set_alignments_discarded_by_length(count=alignments_discarded_by_length)
    config_args.stats.set_alignments_discarded_by_coverage(count=alignments_discarded_by_coverage)
    config_args.stats.set_alignments_discared_by_percentage_of_identity(count=alignments_discarded_by_identity)

    return filtered_alignments

#----------------------------------------
'''Runs the filter alignments'''
def run_filter(config_args: Config):

    filtered_alignments = {}
    start = time.time()
    
    filtered_alignments = filter_alignments_parallel(alignments_file=config_args.xmfa_file_path, 
                                            min_alignment_length = config_args.minimum_alignment_length, 
                                            min_alignment_coverage = config_args.minimum_alignment_coverage,
                                            min_alignment_identity =  config_args.minimum_alignment_percentage_identity,
                                            ingroup_size=config_args.ingroup_size,
                                            config_args=config_args)
            

    success = save_alignments(alignments=filtered_alignments, config_args=config_args, alignment_filename = config_args.filtered_xmfa_name)


    filtered_alignments_reduced = reduce_sequences_in_alignments(alignments=filtered_alignments)
        
    success = save_alignments(alignments=filtered_alignments_reduced, config_args=config_args, alignment_filename = config_args.reduced_filtered_xmfa_name)
    
    config_args.stats.compute_total_alignnments_kept()
    
    end = time.time()
    mins = (end-start)/60
    
    logging.info("Filtering Runtime mins {}".format(mins))
    config_args.stats.filtering_runtime = mins

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




