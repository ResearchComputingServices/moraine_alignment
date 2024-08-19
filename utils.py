import logging
import os
from datetime import datetime
from zipfile import ZipFile
import shutil
import json
from sys import platform


def save_to_json(dictio: dict, filename: str) -> bool:
    """
    Save a dictionary to a JSON file.

    Keyword arguments:
        dictio (dict): The dictionary to be saved.
        filename (str): The path to the JSON file.

    Returns:
        bool: True if the dictionary was successfully saved, False otherwise.
    """

    success = True
    try:
        with open(filename, "w") as outfile:
            json.dump(dictio, outfile)
    except Exception as e:
        success = False
        logging.error(e)

    return success


def load_from_json(filename: str):
    """
    Load data from a JSON file.

    Keyword arguments:
        filename (str): The path to the JSON file.

    Returns:
        tuple: A tuple containing a boolean indicating the success of loading the data and a dictionary containing the loaded data.
    """

    success = True
    dictio = {}
    try:
        with open(filename, 'r') as openfile:
            # Reading from json file
            dictio = json.load(openfile)
    except Exception as e:
        success = False
        logging.error(e)

    return success, dictio


def create_dir(dir_path, dir_name):
    """
    Create a directory with the given directory path and directory name.

    Keyword arguments:
        dir_path (str): The path of the parent directory.
        dir_name (str): The name of the directory to be created.

    Returns:
        str: The full path of the created directory.

    Raises:
        OSError: If an error occurs while creating the directory.
    """

    fullpath_dir = ''
    try:
        fullpath_dir = os.path.join(dir_path, dir_name)
        if not os.path.isdir(fullpath_dir):
            os.mkdir(fullpath_dir)
    except Exception as e:
        fullpath_dir = ''
        logging.error(e)

    return fullpath_dir


def get_today_datetime():
    """
    Returns the current date and time in the format "YYYY_MM_DD_HH_MM".

    Returns:
        str: The current date and time in the format "YYYY_MM_DD_HH_MM".
    """

    now = datetime.now()
    dt_string = now.strftime("%Y_%m_%d_%H_%M")
    return dt_string


def generate_filename(path: str, name: str, extension: str):
    """
    Generate a filename with the given path, name, and extension.

    Parameters:
        path (str): The path where the file will be saved.
        name (str): The name of the file.
        extension (str): The extension of the file.

    Returns:
        str: The full path of the generated filename.
    """

    today_date = get_today_datetime()
    filename = today_date + '_' + name + '.' + extension
    filename_path = os.path.join(path, filename)
    return filename_path


def delete_files(directory_path):
    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            os.remove(os.path.join(directory_path, filename))


def get_filename(file_path):
    """
    Returns the filename from the given file path.

    Keyword arguments:
        file_path (str): The path of the file.

    Returns:
        str: The filename extracted from the file path.
    """

    try:
        full_name = os.path.basename(file_path)
    except:
        full_name = file_path
    return full_name


def is_directory(dir: str):
    """
    Check if the given path is a directory.

    Keyword arguments:
        dir (str): The path to check.

    Returns:
        bool: True if the path is a directory, False otherwise.
    """

    if os.path.isdir(dir):
        return True
    else:
        return False


def get_files_in_dir(folder_path: str) -> list:
    """
    Retrieves a list of file paths in the specified directory.

    Keyword arguments:
        folder_path (str): The path of the directory.

    Returns:
        list: A list of file paths in the directory.

    Raises:
        Exception: If the directory does not exist.
    """

    try:
        file_path_list = []
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isfile(file_path):
                file_path_list.append(file_path)
    except Exception as e:
        logging.error("Directory does not exist.")
        logging.error(e)

    return file_path_list


def clear_file_paths(file_paths: list) -> list:
    """
    Filters a list of file paths and returns only the paths that have supported formats.

    Keyword arguments:
        file_paths (list): A list of file paths.

    Returns:
        list: A list of file paths with supported formats.
    """

    clear_list = []
    for path in file_paths:
        if os.path.isfile(path):
            supported_formats = ['.fasta', '.fna', '.fsa_nt']
            extension = os.path.splitext(path)[1]
            if extension in supported_formats:
                clear_list.append(path)
    return clear_list


def get_sublist(ori_list: list, size: int) -> list:
    """
    Returns a sublist of the original list with the specified size.

    Keyword arguments:
        ori_list (list): The original list.
        size (int): The size of the sublist.

    Returns:
        list: The sublist of the original list.
    """

    if size == None:
        size = len(ori_list)

    if size > len(ori_list):
        size = len(ori_list)

    sublist = ori_list[0:size]

    return sublist


def zip_directory(output_dir, zip_dir):
    """
    Compresses the contents of a directory into a zip file.

    Keyword arguments:
        output_dir (str): The path to the directory where the zip file will be created.
        zip_dir (str): The path to the directory whose contents will be compressed.

    Raises:
        Exception: If an error occurs during the compression process.

    Returns:
        None
    """

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
                        zip_object.write(
                            file_path, os.path.basename(file_path))

        if os.path.exists(zip_filename_path):
            logging.info("ZIP file created")
        else:
            logging.info("ZIP file not created")
    except Exception as e:
        logging.info("An exception has occurred")
        logging.info(e)


def get_files_in_dir_recursive(folder_path: str) -> list:
    """
    Recursively retrieves all file paths within a given directory.

    Keyword arguments:
        folder_path (str): The path to the directory.

    Returns:
        list: A list of file paths within the directory and its subdirectories.
    """

    try:
        file_path_list = []
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isfile(file_path):
                file_path_list.append(file_path)
            else:
                if os.path.isdir(file_path):
                    recursive_list = get_files_in_dir_recursive(
                        folder_path=file_path)
                    file_path_list = file_path_list + recursive_list
    except Exception as e:
        logging.info("Directory does not exist.")
        logging.info(e)

    return file_path_list


def has_subfolders(folder_path: str) -> bool:
    """
    Check if a given folder path contains subfolders.

    Keyword arguments:
        folder_path (str): The path of the folder to check.

    Returns:
        bool: True if the folder contains subfolders, False otherwise.
    """

    subfolders = False
    try:
        for f in os.listdir(folder_path):
            file_path = os.path.join(folder_path, f)
            if os.path.isdir(file_path):
                subfolders = True
                break
    except Exception as e:
        logging.info("Directory does not exist.")
        logging.info(e)

    return subfolders


def unzip_folder(zip_path: str, dest_path: str):
    """
    Unzips a zip file to the specified destination path.

    Keyword arguments:
        zip_path (str): The path to the zip file.
        dest_path (str): The destination path where the zip file will be extracted.

    Returns:
        bool: True if the unzip operation is successful, False otherwise.
    """

    success = True
    if os.path.isfile(zip_path):

        supported_formats = ['.zip']

        if os.path.splitext(zip_path)[1] not in supported_formats:
            logging.error("The file is not a zip file.")
            success = False
            return success

        try:
            _format = os.path.splitext(zip_path)[1][1:]
            dir_path = os.path.splitext(zip_path)[0]
            if not os.path.isdir(dest_path):
                dest_path = dir_path

            shutil.unpack_archive(zip_path, dest_path, _format)

        except Exception as e:
            logging.info("Error when unzip file")
            logging.info(e)
            success = False
    else:
        logging.info("The path does not exists.")
        success = False

    return success


def remove_file_older_than(filepath, days):
    """
    Removes a file if it is older than the specified number of days.

    Keyword arguments:
        filepath (str): The path to the file.
        days (int): The number of days.

    Returns:
        float: The number of days the file is older than, or -1 if the file does not exist.

    Raises:
        OSError: If there is an error while removing the file.
    """

    old_days = -1
    if os.path.isfile(filepath):
        if days >= 0:
            if platform == "linux" or platform == "linux2":
                file_created_datetime = os.stat(filepath).st_ctime
            elif platform == "darwin":
                file_created_datetime = os.stat(filepath).st_birthtime
            current_datetime = datetime.timestamp(datetime.now())
            dif = current_datetime - file_created_datetime  # Dif in seconds between two days
            # 86400 = Number of secs in 1 day
            old_days = dif / 86400
            if old_days >= days:
                try:
                    os.remove(filepath)
                except OSError as e:
                    logging.error("Error: %s : %s" % (filepath, e.strerror))
    return old_days


def remove_folder_older_than(folder, days=None):
    """
    Remove a folder if it is older than a specified number of days.

    Keyword arguments:
        folder (str): The path to the folder.
        days (int, optional): The number of days. If not provided or negative, defaults to 0.

    Returns:
        None

    Raises:
        OSError: If an error occurs while removing the folder.
    """

    if not os.path.exists(folder):
        return

    if not days or days < 0:
        days = 0
    try:
        if platform == "linux" or platform == "linux2":
            file_created_datetime = os.stat(folder).st_ctime
        elif platform == "darwin":
            file_created_datetime = os.stat(folder).st_birthtime

        current_datetime = datetime.timestamp(datetime.now())
        dif = current_datetime - file_created_datetime  # Dif in seconds between two days
        # 86400 = Number of secs in 1 day
        old_days = dif / 86400
        if old_days >= days:
            try:
                shutil.rmtree(folder)
            except Exception as e:
                logging.error(e)
    except Exception as e:
        logging.error(e)


def cleanup_folder(folder, days=None):
    """
    Cleans up the specified folder by removing files and subfolders that are older than the specified number of days.

    Keyword arguments:
        folder (str): The path to the folder to be cleaned up.
        days (int, optional): The number of days. If not provided or set to a negative value, all files and subfolders will be removed. Defaults to None.

    Raises:
        OSError: If an error occurs while accessing or removing files.
    """

    try:
        if not os.path.exists(folder):
            return

        if not days or days < 0:
            days = 0

        for path in os.listdir(folder):
            if not path.startswith('.'):
                fullpath = os.path.join(folder, path)
                if os.path.isdir(fullpath):
                    remove_folder_older_than(folder=fullpath, days=days)
                else:
                    if os.path.isfile(fullpath):
                        remove_file_older_than(filepath=fullpath, days=days)

    except Exception as e:
        logging.error("Error: %s : %s" % (folder, e.strerror))
