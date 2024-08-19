# Do NOT change these parameters (lines 2-7).
DATA_FOLDER = "/data"
MULTIFASTA_FOLDER = "Multifasta/"
IN_PROCESS_FOLDER = "in_process"
TEMP_FOLDER = "Temp/"
OUTPUT_FOLDER = "Output/"
OUTPUT_PARSNP_FOLDER = "Output_parsnp/"

# Ingroup and Outgroup folders WITHOUT path
INGROUP_FOLDER = "Ingroup_Demo"
OUTGROUP_FOLDER = "Outgroup_Demo"

INGROUP_SIZE = None         #Leave None to use all genomes in the ingroup folder or specify the number of genomes to use
OUTGROUP_SIZE = None        #Leave None to use all genomes in the outgroup  folder or specify the number of genomes to use

# Filtering parameters for the parsnp alignment
MINIMUM_ALIGNMENT_LENGTH = 300                   #Minimum sequence length in an alignment
MINIMUM_ALIGNMENT_COVERAGE = 99                  #Minimum percentage of genomes that there should be an alignment values from [0,100]
MINIMUM_ALIGNMENT_PERCENTAGE_IDENTITY = 99  

# Add the location of the parsnp file to avoid re-run parsnp
PARSNP_XMFA_FILE_LOCATION = ""
FILTERED_XMFA_FILE_LOCATION = ""

# Maximum sequence length for a sequence alignment to be blast
SEQUENCE_LENGTH = 300
SHIFT = 250
E_CUTOFF_OUTGROUP = 1
PERC_IDENTITY_OUTGROUP = 60

# Format of the input files
FORMAT = "fasta"

# Removing older output files
CLEANUP_DAYS = 60

# Maximum percentage of matches (0-100) 0 = There should be zero hits, and 100 we allow all the hits
OUTGROUP_FILTER_PERCENTAGE = 50
