import os
from config import Config
from parsnp_alignment import run_parsnp
from alignments_filter import run_filter
from utils import delete_files, get_today_datetime,create_dir
from ncbi_blast import run_blast


def main():
    #Load input parameters
    config_args = Config()
    config_args._print()


    #Verify output paths exit; otherwise create them.
    for path in [config_args.temp_path, config_args.output_path, config_args.output_parsnp_path]:
        os.makedirs(path, exist_ok=True)


    #Create results folder
    dir_name = get_today_datetime()
    results_path = create_dir(config_args.output_path,dir_name)
    config_args.results_path = results_path
 

    #Run parsnp 
    success = True
    if config_args.xmfa_file_path=="":
        success = run_parsnp(config_args=config_args, reference_genome=None)
        
    
    #If success, filter alignments
    if success:
        run_filter(config_args=config_args)


    #Delete at some point the temp directory where the ingroup was copied (if one was created)
    if config_args.copied_ingroup_folder:
        delete_files(config_args.ingroup_folder)

if __name__ == '__main__':
    #main()
    run_blast()    
