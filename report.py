import pandas as pd
import logging
from genome import Genome
from genome import Alignment
from config import Config
from utils import *
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''Removes subsequences that have hits with more than max_percentage_genomes% in the outgroup'''
def filter_subsequences_with_maximum_hits(alignments:list, outgroup_size:int, max_percentage_genomes=None):
    
    if max_percentage_genomes==None:
        return
    
    max_allowed_genomes = round(outgroup_size*(max_percentage_genomes/100))

    for alignment in alignments:
        subsequences_to_keep = []
        for subsequence in alignment.subsequences:
            counter_genomes_hit=0
            for gen_id, match in subsequence.outgroup_hits.items():
                if (match!=None) and len(match)>0:
                    if "hit_score" in match:
                        counter_genomes_hit = counter_genomes_hit + 1
               
            if counter_genomes_hit <= max_allowed_genomes:
                subsequences_to_keep.append(subsequence)

        alignment.subsequences = subsequences_to_keep


'------------------------------------------------------------------------------------------------------------'
'''Creates a fasta file with multiple sequences'''
def create_multifasta_file_for_list(alignments: List[Alignment], config_args:Config)->str:
    
    list = []    
    fasta_file_path = ""

    filter_subsequences_with_maximum_hits(alignments=alignments, outgroup_size=config_args.outgroup_size) 

    for id, alignment in alignments.items():
        for sequence in alignment.subsequences:
            str_id = "genome_{} position_{}".format(alignment.genome_filename, sequence.position_in_genome)
            record = SeqRecord(Seq(sequence.value), id=str_id)
            list.append(record)
       
    if len(list)>0:
        fasta_file = "all_candidates.fasta"
        fasta_file_path = os.path.join(config_args.results_path, fasta_file)
        os.makedirs(config_args.results_path, exist_ok=True)
        SeqIO.write(list, fasta_file_path, "fasta")  
    
    return fasta_file_path




def general_header(config_args:Config, total_time:float):
    
    header_dict = {}
    header_dict["Ingroup count: "]=config_args.ingroup_size
    header_dict["Outgroup count: "]=config_args.outgroup_size
    header_dict["Subsequence Length: "] = config_args.sequence_length
    header_dict["Shift: "] = config_args.shift
    header_dict["E value - Outgroup: "]= config_args.e_cutoff_outgroup
    header_dict["Perc Identity - Outgroup: "] = config_args.perc_identity_outgroup
    header_dict["Maximum Percentage Hits in Outgroup: "] = config_args.outgroup_filter_percentage
    header_dict["Processing time in mins: "]=total_time
    header_dict["Report date: "] = get_today_datetime()

    header_fr = pd.DataFrame(header_dict, index=[0])
    header_fr_t = header_fr.T
    header_fr_t.columns =[""]
    

    return header_fr_t

def stats_page(config_args:Config):

    header_dict = {}
    header_dict["Total Alignments Found by Parsnp: "] = config_args.stats.alignments_found_by_parsnp
    
    header_dict["Alignments Discarded by Length: "] = config_args.stats.alignments_discarded_by_length
    header_dict["Percentage Alignments Discarded by Length: "] = str(round(config_args.stats.percentage_discarded_by_length*100)) + "%"
   
    header_dict["Alignments Discarded by Coverage: "] = config_args.stats.alignments_discarded_by_coverage
    header_dict["Percentage Alignments Discarded by Coverage: "] = str(round(config_args.stats.percentage_discarded_by_coverage*100))+"%"
    
    header_dict["Alignments Discarded by Percentage of Identity: "] = config_args.stats.alignments_discared_by_percentage_of_identity
    header_dict["Percentage Alignments Discarded by Percentage of Identity: "] = str(round(config_args.stats.percentage_discared_by_percentage_of_identity*100))+"%"


    header_dict["Total Alignment Discarded (%): "] =  config_args.stats.total_alignments_taken
    header_dict["Percentage Total Alignment Discarded: "] = str(round(config_args.stats.percentage_alignment_taken*100))+"%"
    
    header_dict["Total Alignment Kept: "] = config_args.stats.total_alignments_kept
    header_dict["Percentage Total Alignment Kept: "] = str(round(config_args.stats.percentage_alignment_kept*100))+"%"
    

    header_dict["Total Subsequences From Alignments: "] = config_args.stats.total_subsequences_from_alignments
    header_dict["Total Candidates: "] = config_args.stats.total_candidates
    header_dict["Candidates Percentage: "] = config_args.stats.percentage_candidates

    header_dict["Parsnp runtime: "] =  config_args.stats.parsnp_runtime
    header_dict["Filtering runtime: "] = config_args.stats.filtering_runtime
    header_dict["Blast runtime: "] = config_args.stats.blast_runtime
    header_dict["Total runtime: "] =  config_args.stats.total_runtime

    header_fr = pd.DataFrame(header_dict, index=[0])
    header_fr_t = header_fr.T
    header_fr_t.columns =[""]
    

    return header_fr_t

   

    

    
    
def candidates_main(pathogen_candidates:List[Alignment]):             

    main_dict = {}
    index = 0
    for id, genome in pathogen_candidates.items():
        for subsequence in genome.subsequences:
            d={}
            #d["Subsequence Id"]=subsequence.id
            d["Subsequence"] = subsequence.value
            d["Outgroup Hits"] = len(subsequence.outgroup_hits)
            d["Subsequence position in Genome"] = subsequence.position_in_genome
            #d["Outgroup Score"] = "TBC"
            d["Percentage of Identity"] = genome.percentage_of_identity
            d["Alignment Info"] = genome.description
            d["Alignment Strand"] = genome.strand
            #d["Alignment Position in Genome"] = genome.position_in_genome 
            d["Alignment Length"] = genome.length
            #d["Subsequence Position"] = subsequence.position
            d["Genome Filename"] = genome.genome_filename
            d["Genome Description"] = genome.genome_header
            d["Genome Length"] = genome.genome_length
            main_dict[index] = d
            index = index + 1

    
    main_fr = pd.DataFrame(main_dict)
    
    if main_fr.empty:
        return main_fr
    
    main_fr_t = main_fr.T
    try:
        main_fr_t_s = main_fr_t.sort_values(['Outgroup Hits', 'Percentage of Identity'], ascending=[True, False])
        #main_fr_t_s = main_fr_t.sort_values(['Outgroup Hits'], ascending=[True])
    except Exception as e:
        logging.info(e)
        #print(e)
   
    return main_fr_t_s



#----------------------------------------------------------------------------------------------------
def create_general_report(config_args:Config,
                          pathogen_candidates: List[Genome],
                          mins:float):
 
    logging.info("\nCreating general report ------------------------------------------ \n")
    filename_path = generate_filename(config_args.results_path, "report","xlsx")  
    writer = pd.ExcelWriter(filename_path, engine='xlsxwriter')
   
    #General information
    #Create header dataframe
    header_fr = general_header(config_args=config_args, total_time=mins)
    header_fr.to_excel(writer, sheet_name="Pathogens Candidates", index=True)
  
    #All candidates in one sheet
    sheet_name = "Pathogens Candidates"
    main_fr = candidates_main(pathogen_candidates=pathogen_candidates)
    if not main_fr.empty:
        main_fr.to_excel(writer, startrow=15, sheet_name=sheet_name, index=False, header=True)

    #Stats page
    stats_fr = stats_page(config_args=config_args)
    sheet_name = "Stats"
    if not stats_fr.empty:
        stats_fr.to_excel(writer, startrow=4, sheet_name=sheet_name, index=True)


    workbook = writer.book
    workbook.formats[0].set_font_size(14)

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    create_multifasta_file_for_list(alignments=pathogen_candidates, config_args=config_args)
    #zip_directory(output_dir=results_dir_name, zip_dir=results_dir_name)

    return filename_path


