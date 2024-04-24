import os
from datetime import datetime
from zipfile import ZipFile
import json


####################################################################################################
def save_to_json(dictio:dict, filename:str):
    success = True
    try:
        with open(filename, "w") as outfile:
            json.dump(dictio, outfile)
    except Exception as e:
        success = False
        print(e)
    
    return success


####################################################################################################
def load_from_json(filename:str):
    success = True
    dictio = {}
    try:
        with open(filename, 'r') as openfile:
            # Reading from json file
            dictio = json.load(openfile)
    except Exception as e:
        success = False
        print(e)
    
    return success, dictio
    
####################################################################################################
def create_dir(dir_path, dir_name):
    
    fullpath_dir=''
    try:
        fullpath_dir = os.path.join(dir_path,dir_name)
        if not os.path.isdir(fullpath_dir):
            os.mkdir(fullpath_dir)
    except Exception as e:
        fullpath_dir=''
        print(e)
    
    return fullpath_dir

####################################################################################################
def get_today_datetime():
    now = datetime.now()
    #dt_string = now.strftime("%Y_%m_%d")
    dt_string = now.strftime("%Y_%m_%d_%H_%M")
    return dt_string

####################################################################################################
def generate_filename(path:str, name:str, extension:str):
    today_date = get_today_datetime()
    filename = today_date + '_' + name + '.' + extension
    filename_path = os.path.join(path, filename)
    return filename_path
 

####################################################################################################
def delete_files(directory_path):
    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            os.remove(os.path.join(directory_path, filename))

####################################################################################################
def get_filename(file_path):
    try:
        full_name = os.path.basename(file_path)
        #file_name = os.path.splitext(full_name)
    except:
         full_name = file_path
    return full_name


####################################################################################################
def is_directory(dir: str):
    if os.path.isdir(dir):
         return True
    else:
         return False


####################################################################################################
def get_files_in_dir(folder_path: str) -> list:
    try:
        file_path_list = []
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isfile(file_path):
                file_path_list.append(file_path)
    except Exception as e:
            print ("Directory does not exist.")
            print (e)
            
            #logging.info ("Directory does not exist.")
            #logging.info(e)
        
    return file_path_list




def get_files_in_dir_recursive(folder_path: str) -> list:
    try:
        file_path_list = []
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isfile(file_path):
                file_path_list.append(file_path)
            else:
                if os.path.isdir(file_path):
                    file_path_list.extend(get_files_in_dir_recursive(file_path))
    except Exception as e:
            print("Directory does not exist.")
            print(e)
    return file_path_list




def clear_file_paths(file_paths: list) -> list:
    clear_list = []
    for path in file_paths:
        if os.path.isfile(path):
            supported_formats = ['.fasta', '.fna', '.fsa_nt']
            extension = os.path.splitext(path)[1]
            if extension in supported_formats:
                clear_list.append(path)
    return clear_list



def get_sublist(ori_list: list, size: int) -> list:
    
    if size == None:
        size = len(ori_list)

    if size>len(ori_list):
        size = len(ori_list)
    
    sublist = ori_list[0:size]
    
    return sublist

####################################################################################################
def zip_directory(output_dir, zip_dir):
    try:

        name = "results"
        extension = "zip"
        today_date = get_today_datetime()
        zip_filename = today_date + '_' + name + '.' + extension
        zip_filename_path = os.path.join(output_dir, zip_filename)
   

        with ZipFile(zip_filename_path, 'w') as zip_object:
            # Traverse all files in directory
            for folder_name, sub_folders, file_names in os.walk(zip_dir):
                for filename in file_names:
                    if filename != zip_filename:
                        # Create filepath of files in directory
                        file_path = os.path.join(folder_name, filename)
                        # Add files to zip file
                        zip_object.write(file_path, os.path.basename(file_path))

        if os.path.exists(zip_filename_path):
            print("ZIP file created")
        else:
            print("ZIP file not created")
    except Exception as e:
        print ("An exception has occurred")
        print(e)


'''Get files in directory visiting all subdirectories'''
def get_files_in_dir_recursive(folder_path: str) -> list:
    try:
        file_path_list = []
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isfile(file_path):
                file_path_list.append(file_path)
            else:
                if os.path.isdir(file_path):
                    recursive_list = get_files_in_dir_recursive(folder_path=file_path)
                    file_path_list = file_path_list +recursive_list
                    #file_path_list.extend(get_files_in_dir_recursive(folder_path=file_path))
    except Exception as e:
            #logging.info ("Directory does not exist.")
            #logging.info(e)
            print("Directory does not exist.")
            print(e)

    return file_path_list



'''Returns true if a folder has subfolders'''
def has_subfolders(folder_path: str) -> bool:
    subfolders = False
    try:
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isdir(file_path):
                subfolders = True
                break
    except Exception as e:
            #logging.info ("Directory does not exist.")
            #logging.info(e)
            print("Directory does not exist.")
            print(e)

    return subfolders







  