import pandas as pd
import logging
from genome import Genome
from genome import Alignment
from config import Config
from utils import *
from typing import List

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
        #main_fr_t_s = main_fr_t.sort_values(['Outgroup Hits', 'Outgroup Score'], ascending=[True, True])
        main_fr_t_s = main_fr_t.sort_values(['Outgroup Hits'], ascending=[True])
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

    workbook = writer.book
    workbook.formats[0].set_font_size(14)

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()


    #create_multifasta_file_for_list(genomes=ingroup, file_location=results_dir_name)
    #zip_directory(output_dir=results_dir_name, zip_dir=results_dir_name)

    return filename_path


