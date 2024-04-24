from config import Config
from utils import get_files_in_dir_recursive, create_dir, get_today_datetime, has_subfolders
import shutil
import os
import shlex, subprocess

'''File with routines that run and execute parsnp'''


#------------------------------------------------------------------------------------------------------------------------
'''Calls parsnp'''
'''parsnp should be already installed'''
def parsnp_wrapper(input_directory:str, reference_genome_with_path:str,  output_directory: str, number_threads: int):   
    result = None
    try:
        #command = "./parsnp -r " + reference_genome_with_path +  " -d " + input_directory +  " -o " + output_directory 
        #the -c forces all the genomes in the input_directory to be included in the alignment
        command = "parsnp  -p" + str(number_threads) + " -r " + reference_genome_with_path +  " -d " + input_directory +  " -o " + output_directory + " -c"
        args = shlex.split(command)
        #result = subprocess.run(args, capture_output=True, text=True)
        result = subprocess.call(args)
    except Exception as e:
        print(e)
    
    return result






#----------------------------------------------------------------------------------------------------
'''Check if all the genomes are in a single folder wihtout subfolders'''
'''Otherwise, copy the genomes into a single directory'''
def verify_ingroup_folder(config_args:Config):

    proceed = False

    genome_filepaths = get_files_in_dir_recursive(folder_path=config_args.ingroup_path)

    #There aren't files to process
    if len(genome_filepaths)<=0:
        return proceed
    
    #The ingroup size was not specified, set it to the count of files in the directory
    if config_args.ingroup_size==None:
        config_args.ingroup_size = len(genome_filepaths)

    #We will copy the files to a different location if the ingroup location has directories
    #or if the ingroup size specified is less than the number of files in the ingroup directory
    if has_subfolders(config_args.ingroup_path) or (config_args.ingroup_size<len(genome_filepaths)):

        print ("The ingroup location has subfolders")
        
        #Copy files into a single folder
        dir_path = config_args.data_folder
        dir_name = "Temp_parsnp_" + get_today_datetime() 
        new_ingroup_path = create_dir(dir_path, dir_name)
        
        if os.path.isdir(new_ingroup_path):

            print ("Copying files into directory: {}".format(new_ingroup_path))
            #Copy files to temporary directory
            count = 0
            while (count < config_args.ingroup_size):
                file = genome_filepaths[count]
                try:
                    shutil.copy(file, new_ingroup_path)
                except Exception as e:
                    print("Error when copying files.")
                    print(e)
                count = count + 1
    
            #Set the ingroup location to this new folder
            config_args.ingroup_path = new_ingroup_path
            config_args.copied_ingroup_folder = True

    proceed = True
    return proceed

    



#--------------------------------------------------------------------------------------------------
'''Loads the genomes to run parsnp into a single directory (in case subfolders exists)'''
'''Runs parsnp using a reference_genome (optional)'''
'''Saves the output into a designated output folder'''
def run_parsnp(config_args: Config, reference_genome:str=None):
    
    #Loads genomes into a folder
    proceed = verify_ingroup_folder(config_args=config_args)

    if proceed:
        if reference_genome==None or (not os.path.isfile(reference_genome)):
            reference_genome="!"   #Choose randomly a reference genome
        

    dir_name = get_today_datetime()
    xmfa_filepath = create_dir(config_args.output_parsnp_path, dir_name)
    
    success = False
    if len(xmfa_filepath)>0:
        output = parsnp_wrapper(input_directory=config_args.ingroup_path, 
                                reference_genome_with_path=reference_genome,  
                                output_directory=xmfa_filepath, 
                                number_threads=config_args.parsnp_number_threads)   
        xmfa_filelocation = os.path.join(xmfa_filepath,config_args.xmfa_file_name)
        if output==0 or os.path.isfile(xmfa_filelocation):
            success = True
            config_args.xmfa_file_path = xmfa_filelocation
    
    return success