import logging
from pathlib import Path
import os
import json
import copy
import input_parameters as ip
from stats import Stats
import argparse

config_parser = argparse.ArgumentParser(description='Configuration parameters for the pipeline')
config_parser.add_argument('--data_folder', type=str, help='Path to the data folder')
config_parser.add_argument('--multifasta_folder', type=str, help='Path to the multifasta folder')
config_parser.add_argument('--temp_folder', type=str, help='Path to the temp folder')
config_parser.add_argument('--ingroup_folder', type=str, help='Path to the ingroup folder')
config_parser.add_argument('--outgroup_folder', type=str, help='Path to the outgroup folder')
config_parser.add_argument('--output_parsnp_folder', type=str, help='Path to the output parsnp folder')
config_parser.add_argument('--output_folder', type=str, help='Path to the output folder')
config_parser.add_argument('--ingroup_size', type=int, help='Size of the ingroup')
config_parser.add_argument('--outgroup_size', type=int, help='Size of the outgroup')
config_parser.add_argument('--minimum_alignment_length', type=int, help='Minimum alignment length')
config_parser.add_argument('--minimum_alignment_coverage', type=float, help='Minimum alignment coverage')
config_parser.add_argument('--minimum_alignment_percentage_identity', type=float, help='Minimum alignment percentage identity')
config_parser.add_argument('--format', type=str, help='Alignment format')
config_parser.add_argument('--sequence_length', type=int, help='Length of the sequence')
config_parser.add_argument('--shift', type=int, help='Shift value')
config_parser.add_argument('--outgroup_filter_percentage', type=float, help='Outgroup filter percentage')
config_parser.add_argument('--xmfa_file_path', type=str, help='Path to the xmfa file')
config_parser.add_argument('--xmfa_file_name', type=str, help='Name of the xmfa file')
config_parser.add_argument('--copied_ingroup_folder', type=bool, help='Flag indicating if the ingroup folder is copied')
config_parser.add_argument('--copied_outgroup_folder', type=bool, help='Flag indicating if the outgroup folder is copied')
config_parser.add_argument('--filtered_xmfa_name', type=str, help='Name of the filtered xmfa file')
config_parser.add_argument('--reduced_filtered_xmfa_name', type=str, help='Name of the reduced filtered xmfa file')
config_parser.add_argument('--filtered_xmfa_path', type=str, help='Path to the filtered xmfa file')
config_parser.add_argument('--results_path', type=str, help='Path to the results folder')
config_parser.add_argument('--parsnp_number_threads', type=int, help='Number of threads for parsnp')
config_parser.add_argument('--processors_number', type=int, help='Number of processors')
config_parser.add_argument('--e_cutoff_outgroup', type=float, help='E cutoff for outgroup')
config_parser.add_argument('--perc_identity_outgroup', type=float, help='Percentage identity for outgroup')
config_parser.add_argument('--max_hsps', type=int, help='Maximum number of HSPs')
config_parser.add_argument('--cleanup_days', type=int, help='Number of days for cleanup')



class Config:
    def __init__(self, **kwargs):
        """
        Initializes the Config object with the given parameters.
        Keyword arguments:
            - data_folder (str): The path to the data folder.
            - multifasta_folder (str): The name of the multifasta folder.
            - temp_folder (str): The name of the temporary folder.
            - ingroup_folder (str): The name of the ingroup folder.
            - outgroup_folder (str): The name of the outgroup folder.
            - output_parsnp_folder (str): The name of the output parsnp folder.
            - output_folder (str): The name of the output folder.
            - ingroup_size (int): The size of the ingroup.
            - outgroup_size (int): The size of the outgroup.
            - minimum_alignment_length (int): The minimum alignment length.
            - minimum_alignment_coverage (float): The minimum alignment coverage.
            - minimum_alignment_percentage_identity (float): The minimum alignment percentage identity.
            - format (str): The format of the data.
            - sequence_length (int): The length of the sequence.
            - shift (int): The shift value.
            - outgroup_filter_percentage (float): The outgroup filter percentage.
            - xmfa_file_path (str): The path to the XMFA file.
            - xmfa_file_name (str): The name of the XMFA file.
            - copied_ingroup_folder (bool): Indicates if the ingroup folder has been copied.
            - copied_outgroup_folder (bool): Indicates if the outgroup folder has been copied.
            - filtered_xmfa_name (str): The name of the filtered XMFA file.
            - reduced_filtered_xmfa_name (str): The name of the reduced filtered XMFA file.
            - filtered_xmfa_path (str): The path to the filtered XMFA file.
            - results_path (str): The path to the results.
            - parsnp_number_threads (int): The number of threads for Parsnp.
            - processors_number (int): The number of processors.
            - e_cutoff_outgroup (float): The E cutoff for the outgroup.
            - perc_identity_outgroup (float): The percentage identity for the outgroup.
            - max_hsps (int): The maximum number of HSPs.
            - stats (Stats): The statistics object.
            - cleanup_days (int): The number of days for cleanup.
        
        Additional parameters can be provided as keyword arguments and will overwrite the default values.
        
        Returns:
            None
        """
        
        # Input parameters
        
        self.data_folder = ip.DATA_FOLDER
        self.multifasta_folder = ip.MULTIFASTA_FOLDER
        self.temp_folder = ip.TEMP_FOLDER
        self.ingroup_folder = ip.INGROUP_FOLDER
        self.outgroup_folder = ip.OUTGROUP_FOLDER
        self.output_parsnp_folder = ip.OUTPUT_PARSNP_FOLDER
        self.output_folder = ip.OUTPUT_FOLDER

        self.ingroup_size = ip.INGROUP_SIZE
        self.outgroup_size = ip.OUTGROUP_SIZE

        self.multifasta_path = os.path.join(
            self.data_folder, self.multifasta_folder)
        self.in_process_path = os.path.join(
            self.multifasta_path, ip.IN_PROCESS_FOLDER)
        self.ingroup_path = os.path.join(self.data_folder, self.ingroup_folder)
        self.outgroup_path = os.path.join(
            self.data_folder, self.outgroup_folder)

        base_dir = Path(__file__).resolve().parent.as_posix()
        self.temp_path = os.path.join(base_dir, self.temp_folder)
        self.output_path = os.path.join(base_dir, self.output_folder)
        self.output_parsnp_path = os.path.join(
            base_dir, self.output_parsnp_folder)

        self.ingroup_size = ip.INGROUP_SIZE
        self.outgroup_size = ip.OUTGROUP_SIZE

        self.minimum_alignment_length = ip.MINIMUM_ALIGNMENT_LENGTH
        self.minimum_alignment_coverage = ip.MINIMUM_ALIGNMENT_COVERAGE/100
        self.minimum_alignment_percentage_identity = ip.MINIMUM_ALIGNMENT_PERCENTAGE_IDENTITY/100

        self.format = ip.FORMAT
        self.sequence_length = ip.SEQUENCE_LENGTH
        self.shift = ip.SHIFT
        self.outgroup_filter_percentage = ip.OUTGROUP_FILTER_PERCENTAGE

        # Other variables
        self.xmfa_file_path = ""
        self.xmfa_file_name = "parsnp.xmfa"
        self.copied_ingroup_folder = False
        self.copied_outgroup_folder = False
        self.filtered_xmfa_name = "filtered_parsnp.json"
        self.reduced_filtered_xmfa_name = "reduced_filtered_parsnp.json"
        self.filtered_xmfa_path = ""
        self.results_path = ""

        if ip.PARSNP_XMFA_FILE_LOCATION != "" and os.path.isfile(ip.PARSNP_XMFA_FILE_LOCATION):
            self.xmfa_file_path = ip.PARSNP_XMFA_FILE_LOCATION

        if ip.FILTERED_XMFA_FILE_LOCATION != "" and os.path.isfile(ip.FILTERED_XMFA_FILE_LOCATION):
            self.filtered_xmfa_path = ip.FILTERED_XMFA_FILE_LOCATION

        self.parsnp_number_threads = 8

        self.processors_number = os.cpu_count()

        self.e_cutoff_outgroup = ip.E_CUTOFF_OUTGROUP
        self.perc_identity_outgroup = ip.PERC_IDENTITY_OUTGROUP
        self.max_hsps = 1

        self.stats = Stats()

        self.cleanup_days = ip.CLEANUP_DAYS
        
        # Overwrite default values with the ones provided in the command line
        for key in vars(self):
            if key in kwargs and kwargs[key] is not None:
                setattr(self, key, kwargs[key])

    def to_file(self):
        """
        Save the configuration parameters to a file.
        
        Raises:
            ValueError: If `results_path` is not defined.
        """
        
        variables = copy.deepcopy(vars(self))
        to_remove = ['stats']
        for var in to_remove:
            variables.pop(var)
        if self.results_path == "":
            raise ValueError("results_path is not defined")
        else:
            os.makedirs(self.results_path, exist_ok=True)
            with open(os.path.join(self.results_path, "config_parameters.txt"), "w", encoding='utf-8') as f:
                json.dump(variables, f, indent=2)

    def _print(self):
        """
        Print the configuration parameters.
        This method prints all the configuration parameters of the object.
        
        Keyword arguments:
            - self: The object instance.
        
        Returns:
            - None
        """
        
        print("Configuration parameters:")
        for key in vars(self):
            print(f"{key}: {getattr(self, key)}")

    def _log(self):
        '''Logs the configuration parameters'''
        info_line_1 = f"Ingroup location: {self.ingroup_path} - Outgroup location: {self.outgroup_path}"
        info_line_2 = f"Ingroup size: {self.ingroup_size} "\
            f"- Outgroup size: {self.outgroup_size}"
        info_line_3 = f"Alignment filtering parameters:  Min. Length: {self.minimum_alignment_length}, Min. % Coverage {self.minimum_alignment_coverage}, Min. % Identity {self.minimum_alignment_percentage_identity}"
        info_line_4 = f"xmfa location: {self.output_parsnp_path} - Filtered xmfa location: {self.results_path}"

        logging.info(info_line_1)
        logging.info(info_line_2)
        logging.info(info_line_3)
        logging.info(info_line_4)
