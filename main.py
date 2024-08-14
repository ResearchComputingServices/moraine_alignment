import os
from config import Config
from parsnp_alignment import run_parsnp
from alignments_filter import run_filter
from utils import delete_files, get_today_datetime,create_dir, unzip_folder, cleanup_folder, remove_folder_older_than
from ncbi_blast import run_blast
import logging
import time
import pathogens
from report import create_general_report
    
def config_logging(log_dir:str):
    output_filename_prefix = "LOG"
    log_file_name = f"{output_filename_prefix}_{time.strftime('%y_%m_%d_%H_%M')}.log"
    log_file_name_path = os.path.join(log_dir, log_file_name)
    logging.basicConfig(level=logging.INFO,
                        filename=log_file_name_path, 
                        filemode='w', 
                        format='%(asctime)s %(levelname)s %(filename)s %(funcName)s(%(lineno)d) %(message)s')
    logging.info("Logs have been setup!")


'''This routine setup the config file based on the input parameters'''
'''and the log file.'''
def setup_run():

    config_args = None

    #Load input parameters
    config_args = Config()
   
    #Verify output paths exit; otherwise create them.
    for path in [config_args.temp_path, config_args.output_path, config_args.output_parsnp_path, config_args.multifasta_path, config_args.in_process_path]:
        os.makedirs(path, exist_ok=True)

    #Create results folder
    dir_name = get_today_datetime()
    results_path = create_dir(config_args.output_path,dir_name)
    config_args.results_path = results_path
    # Write the config file to the results folder
    config_args.to_file()

    #Setup the logs
    config_logging(log_dir=results_path)

    return config_args


def main():
    
    start = time.time()


    #Setup all the running parameters
    config_args = setup_run()
    config_args.stats.start_time = start

    if config_args==None:
        return

    #Log input parameters and config
    config_args._log()
    
    #Run parsnp 
    success = True
    if config_args.xmfa_file_path=="" and config_args.filtered_xmfa_path=="":
        success = run_parsnp(config_args=config_args, reference_genome=None)
        
    
    #If success, filter alignments
    if success and config_args.filtered_xmfa_path=="":
         success = run_filter(config_args=config_args)

    
    if success:
        alignment_genomes = pathogens.get_pathogens_from_alignments(config_args=config_args)
        end = time.time()
        mins = (end-start)/60
        config_args.stats.end_time =end
        config_args.stats.total_runtime = mins
        create_general_report(config_args=config_args, pathogen_candidates=alignment_genomes, mins=mins)

    # #Delete at some point the temp directory where the ingroup was copied (if one was created)
    if config_args.copied_ingroup_folder:
         delete_files(config_args.ingroup_path)
         os.rmdir(config_args.ingroup_path)
 
    if config_args.copied_outgroup_folder:
        delete_files(config_args.outgroup_path)
        os.rmdir(config_args.outgroup_path)

    logging.info("Deleting temporal files -------------------------------------\n")
    delete_files(config_args.multifasta_path)
    delete_files(config_args.in_process_path)
    
    end = time.time()
    mins = (end-start)/60
    
    logging.info("\n\n Total Runtime mins {}".format(mins))


if __name__ == '__main__':
    main()
    