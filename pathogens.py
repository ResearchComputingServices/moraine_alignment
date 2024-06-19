from genome import Genome
from genome import Subsequence
from genome import Alignment
import logging
import os
import shutil
from config import Config
from utils import get_files_in_dir_recursive, create_dir, get_today_datetime, has_subfolders
from Bio import SeqIO
from utils import load_from_json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ncbi_blast import blast_subsequences_against_genomes, blast_subsequences_against_genomes__parallel
import math



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''Removes subsequences that have hits with more than max_percentage_genomes% in the outgroup'''
def filter_subsequences_with_maximum_hits(genomes:list, outgroup_size:int, max_percentage_genomes=None):
    
    if max_percentage_genomes==None:
        return
    

    max_allowed_genomes = round(outgroup_size*(max_percentage_genomes/100))

    for genome in genomes:
        subsequences_to_keep = []
        for subsequence in genome.subsequences:
            counter_genomes_hit=0
            for gen_id, match in subsequence.outgroup_hits.items():
                if (match!=None) and len(match)>0:
                    if "hit_score" in match:
                        counter_genomes_hit = counter_genomes_hit + 1
               
            if counter_genomes_hit <= max_allowed_genomes:
                subsequences_to_keep.append(subsequence)

        genome.subsequences = subsequences_to_keep



#-------------------------------------------------------------------------------------------------
'''Write a (one) sequence to a fasta file'''
def write_sequence_to_file(genomes:dict, file_location:str)->list:
    
    genome_files = {}
    for genome_id, genome in genomes.items():
        file_name = create_fasta_file(genome,file_location)    
        if len(file_name)>0:
                genome_files[genome.id]= file_name
        
        
    return genome_files


#--------------------------------------------------------------------------------------------------------------------
'''Creates a fasta file with multiple sequences'''
def create_multifasta_file(genome: Genome, start:int, end:int, file_location:str)->str:
    
    list = []    
    index = start
    fasta_file_path = ""

    for sequence in genome.subsequences[start:end+1]:
        str_id = "alignment_{} subseqindex_{} subseqid_{}".format(genome.id, index, sequence.id)
        record = SeqRecord(Seq(sequence.value), id=str_id)
        index = index + 1
        list.append(record)


    # Write the records to a FASTA file
    if len(list)>0:
        try:
            fasta_file = "Alignment_{}_S_{}_to_{}.fasta".format(genome.id,start,end)
            fasta_file_path = os.path.join(file_location,fasta_file)
            os.makedirs(file_location, exist_ok=True)
            SeqIO.write(list, fasta_file_path, "fasta")  
        except Exception as e:
            logging.error(e)


    return fasta_file_path



#--------------------------------------------------------------------------------------------------------------------
'''Creates a fasta file with a sequence'''
def create_fasta_file(genome:Genome, file_location:str)->str:
      
    fasta_file_path = ""

    try:
        str_id = "alignment_{}_{}".format(genome.id,genome.description)
        record = SeqRecord(Seq(genome.genome_sequence), id=str_id)
        
        fasta_file = "Alignment_{}.fasta".format(genome.id)
        fasta_file_path = os.path.join(file_location,fasta_file)
        os.makedirs(file_location, exist_ok=True)
        SeqIO.write(record, fasta_file_path, "fasta")  
    except Exception as e:
            logging.error(e)

    return fasta_file_path


#--------------------------------------------------------------------------------------------------------------------
'''Write multiple sequences to a fasta file'''
def write_subsequences_to_files(genomes:list, file_location:str, number_processors:int)->list:
    
    genome_subseq_files ={}
    for genome in genomes:
        logging.info ("Genome {} -> # Sequences to write: {}".format(genome.id, len(genome.subsequences)))
        #print ("Genome {} -> # Sequences to write: {}".format(genome.filename, len(genome.subsequences)))
        number_subsequences_per_processor = math.floor(len(genome.subsequences) / number_processors)
        
        if number_subsequences_per_processor<=0:
            #There are very few sequences.
            #Less than number of processors
            #All subsequences will go in one file only
            number_processors = 1
            number_subsequences_per_processor = len(genome.subsequences)

            
        
        logging.info ("Creating multifasta files for genome {}...".format(genome.id))
        #print("Creating multifasta files for genome {}...".format(genome.filename))
        processor_work={}
        ini = 0
        end = number_subsequences_per_processor-1
        for i in range(0,number_processors):
            processor_work[i]={"ini": ini, "end":end}
            ini = end + 1
            end = ini + number_subsequences_per_processor - 1
            
        #Adjust the last processor's work for the number of sequences
        processor_work[number_processors-1]["end"]=len(genome.subsequences)-1
                
        #print ("CPU WORK")
        #pprint.pprint(processor_work)

        file_list=[]
        for i in range(0,number_processors):
            start = processor_work[i]["ini"]
            end = processor_work[i]["end"]
            file_name = create_multifasta_file(genome,start,end,file_location)    
            if len(file_name)>0:
                file_list.append(file_name)
        
        if len(file_list)>0:
            genome_subseq_files[genome.id] = file_list
    
    return genome_subseq_files

#----------------------------------------------------------------------------------------------------
'''Write multiple sequences to a fasta file'''
def write_subsequences_to_files_fixed_files(genomes:list, file_location:str, number_processors:int)->list:
    
    genome_subseq_files ={}
    for genome in genomes:
        logging.info ("Genome {} -> # Sequences to write: {}".format(genome.id, len(genome.subsequences)))
        
        number_subsequences_per_processor = 500 
        number_processors = math.floor(len(genome.subsequences) / number_subsequences_per_processor)
        

        #print ("Genome {} -> # Sequences to write: {}".format(genome.filename, len(genome.subsequences)))
        #number_subsequences_per_processor = math.floor(len(genome.subsequences) / number_processors)
        
        if number_subsequences_per_processor<=0 or number_processors==0:
            #There are very few sequences.
            #Less than number of processors
            #All subsequences will go in one file only
            number_processors = 1
            number_subsequences_per_processor = len(genome.subsequences)

            
        
        logging.info ("Creating multifasta files for genome {}...".format(genome.id))
        #print("Creating multifasta files for genome {}...".format(genome.filename))
        processor_work={}
        ini = 0
        end = number_subsequences_per_processor-1
        for i in range(0,number_processors):
            processor_work[i]={"ini": ini, "end":end}
            ini = end + 1
            end = ini + number_subsequences_per_processor - 1
            
        #Adjust the last processor's work for the number of sequences
        processor_work[number_processors-1]["end"]=len(genome.subsequences)-1
                
        #print ("CPU WORK")
        #pprint.pprint(processor_work)

        file_list=[]
        for i in range(0,number_processors):
            start = processor_work[i]["ini"]
            end = processor_work[i]["end"]
            file_name = create_multifasta_file(genome,start,end,file_location)    
            if len(file_name)>0:
                file_list.append(file_name)
        
        if len(file_list)>0:
            genome_subseq_files[genome.id] = file_list
    
    return genome_subseq_files

    

#----------------------------------------------------------------------------------------------------
'''Check if all the genomes in the outgroup are in a single folder wihtout subfolders'''
'''Otherwise, copy the genomes into a single directory'''
def verify_outgroup_folder(config_args:Config):

    proceed = False

    genome_filepaths = get_files_in_dir_recursive(folder_path=config_args.outgroup_path)

    #There aren't files to process
    if len(genome_filepaths)<=0:
        return proceed
    
    #The ingroup size was not specified, set it to the count of files in the directory
    if config_args.outgroup_size==None:
        config_args.ingroup_size = len(genome_filepaths)

    #We will copy the files to a different location if the ingroup location has directories
    #or if the ingroup size specified is less than the number of files in the ingroup directory
    if has_subfolders(config_args.outgroup_path) or (config_args.outgroup_size<len(genome_filepaths)):

        logging.info("The outgroup location has subfolders")
        
        #Copy files into a single folder
        dir_path = config_args.data_folder
        dir_name = "Temp_outgroup_" + get_today_datetime() 
        new_outgroup_path = create_dir(dir_path, dir_name)
        
        if os.path.isdir(new_outgroup_path):

            logging.info("Copying files into temporal directory: {}".format(new_outgroup_path))
            #Copy files to temporary directory
            count = 0
            while (count < config_args.ingroup_size):
                file = genome_filepaths[count]
                try:
                    shutil.copy(file, new_outgroup_path)
                except Exception as e:
                    logging.error("Error when copying files.")
                    logging.error(e)
                count = count + 1
    
            #Set the ingroup location to this new folder
            config_args.outgroup_path = new_outgroup_path
            config_args.copied_outgroup_folder = True

    proceed = True
    return proceed

#--------------------------------------------------------------------------------------------------
'''For stress testing only'''
'''Read the genomes located in location'''
def read_genomes_duplicate(config_args: Config, larger_size=int, type=None):

    proceed = False 
    genomes = {}
    

    #Loads genomes into a folder
    logging.info("Checking outgroup folder")
    proceed = verify_outgroup_folder(config_args=config_args)

    if not proceed:
        return genomes


    genome_filepaths = get_files_in_dir_recursive(folder_path=config_args.outgroup_path)
    repeated_times = round(larger_size / len(genome_filepaths))

    for i in range(0, repeated_times):
        for file_path in genome_filepaths:
            try:
                genome = Genome(filename=file_path, format=config_args.format, type=type)
                if genome!=None:
                    genomes[genome.id]=genome
            except Exception as e:
                logging.error(e)

    return genomes



#--------------------------------------------------------------------------------------------------
'''Read the genomes located in location'''
def read_genomes(config_args: Config, type=None):

    proceed = False 
    genomes = {}
    

    #Loads genomes into a folder
    logging.info("Checking outgroup folder")
    proceed = verify_outgroup_folder(config_args=config_args)

    if not proceed:
        return genomes

    genome_filepaths = get_files_in_dir_recursive(folder_path=config_args.outgroup_path)
    for file_path in genome_filepaths:
        try:
            genome = Genome(filename=file_path, format=config_args.format, type=type)
            if genome!=None:
                genomes[genome.id]=genome
        except Exception as e:
            logging.error(e)

    return genomes





#-------------------------------------------------------------------------------------------------
def get_alignment_subsequences(config_args:Config, alignment:dict)->list:
    #Get all the subsequences if the alignmnent sequence is larger than a threshold
    sequence = alignment["sequences"][0]
    alignment =Alignment(sequence=sequence, 
                         description = alignment["id"], 
                         strand = alignment["strand"], 
                         type="Alignment", 
                         length = len(sequence),
                         genome_header = alignment["genome_header"] ,
                         genome_filename = alignment["genome_filename"], 
                         genome_length = alignment["genome_length"])
    
    
    if len(alignment.sequence) >= config_args.sequence_length:
        alignment.compute_all_subsequences(subseq_len=config_args.sequence_length, shift=config_args.shift)
    else:
        #Save sequence into a single fasta file
        subsequence = Subsequence(value=alignment.genome_sequence, 
                                  alignment_id=alignment.id, position=0, 
                                  alignment_position_in_genome=alignment.position_in_genome)
        alignment.subsequences.append(subsequence)
        alignment.original_subsequences_count=len(alignment.subsequences)
    
    #subsequence_files = write_subsequences_to_files(genomes=[alignment], file_location=config_args.in_process_path, number_processors=config_args.processors_number)
    subsequence_files = write_subsequences_to_files_fixed_files(genomes=[alignment], file_location=config_args.in_process_path, number_processors=config_args.processors_number)
  
    return alignment, subsequence_files


#-------------------------------------------------------------------------------------------------
def get_pathogens_from_alignments(config_args:Config):
    
    success = False

    #Read outgroup files and create a genome object
    logging.info("Reding Outgroup ------------------------------------------------ \n")
    outgroup_genomes = read_genomes(config_args=config_args, type="Outgroup")
    
    if len(outgroup_genomes)==0:
        return success 


    #Create files for the genomes 
    outgroup_sequence_files = write_sequence_to_file(genomes=outgroup_genomes, file_location=config_args.multifasta_path)
    if len(outgroup_sequence_files)==0:
        return success
    
    #pprint.pprint(outgroup)
    
    logging.info("\n")

    count = 1
    alignment_genomes = {}
    #There exists a filtered alignment file
    if config_args.filtered_xmfa_path!="":
        success, filtered_alignments =  load_from_json(filename = config_args.filtered_xmfa_path)
        if success:
            for alignment in filtered_alignments.items():
                logging.info("\n")
                logging.info("Procesing genome {} of {} --------------------------------".format(count, len(filtered_alignments)))
                #Get all the subsequences if the alignmnent sequence is larger than a threshold
                alignment_genome, subsequence_files = get_alignment_subsequences(config_args=config_args, alignment=alignment[1])
                
                #Blast the alignment subsequences (or sequence) against the outgroup
                if len(subsequence_files)>0:
                    blast_subsequences_against_genomes(genomes_query=[alignment_genome], genomes_subject=outgroup_genomes, query_subseq_fasta_paths=subsequence_files,
                                                   subject_sequence_fasta_paths=outgroup_sequence_files,config_args=config_args)


                #Replace for file saving to avoid memory explosion again.
                alignment_genomes[alignment_genome.id]=alignment_genome
                count = count + 1    
    
    return alignment_genomes



#-------------------------------------------------------------------------------------------------
def get_pathogens_from_alignments(config_args:Config):
    
    success = False

    #Read outgroup files and create a genome object
    logging.info("Reding Outgroup ------------------------------------------------ \n")
    outgroup_genomes = read_genomes(config_args=config_args, type="Outgroup")
    
    if len(outgroup_genomes)==0:
        return success 


    #Create files for the genomes 
    outgroup_sequence_files = write_sequence_to_file(genomes=outgroup_genomes, file_location=config_args.multifasta_path)
    if len(outgroup_sequence_files)==0:
        return success
    
    #pprint.pprint(outgroup)
    
    logging.info("\n")

    count = 1
    alignment_genomes = {}
    #There exists a filtered alignment file
    if config_args.filtered_xmfa_path!="":
        success, filtered_alignments =  load_from_json(filename = config_args.filtered_xmfa_path)
        if success:
            for alignment in filtered_alignments.items():
                logging.info("\n")
                logging.info("Procesing genome {} of {} --------------------------------".format(count, len(filtered_alignments)))
                #Get all the subsequences if the alignmnent sequence is larger than a threshold
                alignment, subsequence_files = get_alignment_subsequences(config_args=config_args, alignment=alignment[1])
                
                #Blast the alignment subsequences (or sequence) against the outgroup
                if len(subsequence_files)>0:
                    blast_subsequences_against_genomes(genomes_query=[alignment], genomes_subject=outgroup_genomes, query_subseq_fasta_paths=subsequence_files,
                                                   subject_sequence_fasta_paths=outgroup_sequence_files,config_args=config_args)


                #Replace for file saving to avoid memory explosion again.
                alignment_genomes[alignment.id]=alignment
                count = count + 1    
    
    return alignment_genomes


#-------------------------------------------------------------------------------------------------
def get_pathogens_from_alignments__parallel(config_args:Config):
    
    success = False

    #Read outgroup files and create a genome object
    logging.info("Reding Outgroup ------------------------------------------------ \n")
    #outgroup_genomes = read_genomes(config_args=config_args, type="Outgroup")
    outgroup_genomes = read_genomes_duplicate(config_args=config_args, larger_size=100, type="Outgroup")


    if len(outgroup_genomes)==0:
        return success 


    #Create files for the genomes 
    outgroup_sequence_files = write_sequence_to_file(genomes=outgroup_genomes, file_location=config_args.multifasta_path)
    if len(outgroup_sequence_files)==0:
        return success
    
    #pprint.pprint(outgroup)
    
    logging.info("\n")

    count = 1
    alignment_genomes = {}
    total_subsequences_from_alignments = 0
    
    #There exists a filtered alignment file
    if config_args.filtered_xmfa_path!="":
        success, filtered_alignments =  load_from_json(filename = config_args.filtered_xmfa_path)
        
        #The pipeline didn't run the filtering part
        if config_args.stats.total_alignments_kept==0:
            config_args.stats.total_alignments_kept = len(filtered_alignments)
            
        if success:
            for filtered_alignment in filtered_alignments.items():
                logging.info("\n")
                logging.info("Procesing genome {} of {} --------------------------------".format(count, len(filtered_alignments)))
                #Get all the subsequences if the alignmnent sequence is larger than a threshold
                alignment_object, subsequence_files = get_alignment_subsequences(config_args=config_args, alignment=filtered_alignment[1])
                total_subsequences_from_alignments = total_subsequences_from_alignments + alignment_object.original_subsequences_count
                
                #Blast the alignment subsequences (or sequence) against the outgroup
                if len(subsequence_files)>0:
                    #logging.info("{} multifasta files to blast against {} outgroup".format(len(subsequence_files),len(outgroup_genomes)))
                    blast_subsequences_against_genomes__parallel(genomes_query=[alignment_object], genomes_subject=outgroup_genomes, query_subseq_fasta_paths=subsequence_files,
                                                   subject_sequence_fasta_paths=outgroup_sequence_files,config_args=config_args)
                    

                #Filter subsequences
                #filter_subsequences_with_maximum_hits(genomes=[alignment_genome], outgroup_size=len(outgroup_genomes), max_percentage_genomes=config_args.outgroup_filter_percentage)
                #Replace for file saving to avoid memory explosion again.
                
                alignment_genomes[alignment_object.id]=alignment_object
                count = count + 1    

                if count==3:
                    break
    
    config_args.stats.total_subsequences_from_alignments = total_subsequences_from_alignments
    return alignment_genomes


 #-------------------------------------------------------------------------------------------------   
