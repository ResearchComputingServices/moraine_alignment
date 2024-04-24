from pathlib import Path
import os
import input_parameters as ip


class Config:
    def __init__(self):
        #Input parameters
        self.data_folder = ip.DATA_FOLDER
        self.multifasta_folder = ip.MULTIFASTA_FOLDER
        self.temp_folder = ip.TEMP_FOLDER
        self.ingroup_folder = ip.INGROUP_FOLDER
        self.outgroup_folder = ip.OUTGROUP_FOLDER
        self.output_parsnp_folder = ip.OUTPUT_PARSNP_FOLDER
        self.output_folder = ip.OUTPUT_FOLDER
        
        self.ingroup_size = ip.INGROUP_SIZE
        self.outgroup_size = ip.OUTGROUP_SIZE

        self.multifasta_path = os.path.join(self.data_folder,self.multifasta_folder)   
        self.ingroup_path = os.path.join(self.data_folder,self.ingroup_folder)
        self.outgroup_path = os.path.join(self.data_folder,self.outgroup_folder)

        base_dir = Path(__file__).resolve().parent.as_posix()
        self.temp_path = os.path.join(base_dir,self.temp_folder)
        self.output_path = os.path.join(base_dir,self.output_folder)
        self.output_parsnp_path = os.path.join(base_dir,self.output_parsnp_folder)
        
        self.ingroup_size = ip.INGROUP_SIZE
        self.outgroup_size =ip.OUTGROUP_SIZE

        self.minimum_alignment_length = ip.MINIMUM_ALIGNMENT_LENGTH                  
        self.minimum_alignment_coverage = ip.MINIMUM_ALIGNMENT_COVERAGE                
        self.minimum_alignment_percentage_identity = ip.MINIMUM_ALIGNMENT_PERCENTAGE_IDENTITY    


        #Other variables
        self.xmfa_file_path=""
        self.xmfa_file_name = "parsnp.xmfa"
        self.copied_ingroup_folder = False
        self.filtered_xmfa_name = "filtered_parsnp.json"
        self.filtered_xmfa_path=""
        self.results_path = ""

        if ip.PARSNP_XMFA_FILE_LOCATION!="" and os.path.isfile(ip.PARSNP_XMFA_FILE_LOCATION):
            self.xmfa_file_path = ip.PARSNP_XMFA_FILE_LOCATION

        if os.cpu_count()>1:
           self.parsnp_number_threads = 8 
        else:
           self.parsnp_number_threads = 1 

        
    def _print(self):
        print(self.data_folder)
        print(self.multifasta_folder)
        print(self.temp_folder) 
        print(self.ingroup_folder) 
        print(self.outgroup_folder) 
        print(self.output_parsnp_folder) 
        print(self.output_folder)
        
        print(self.ingroup_size) 
        print(self.outgroup_size) 

        print(self.multifasta_path)
        print(self.ingroup_path)
        print(self.outgroup_path)

        print(self.temp_path)
        print(self.output_path) 
        print(self.output_parsnp_path) 

        print(self.ingroup_size)
        print(self.outgroup_size) 

        print(self.minimum_alignment_length)                  
        print(self.minimum_alignment_coverage)                
        print(self.minimum_alignment_percentage_identity) 


