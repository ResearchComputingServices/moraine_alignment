DATA_FOLDER = "/data"
MULTIFASTA_FOLDER = "Multifasta/"
IN_PROCESS_FOLDER = "in_process"
TEMP_FOLDER = "Temp/"
INGROUP_FOLDER = "Ingroup_Echoli_1000"
OUTGROUP_FOLDER = "Outgroup_Wong/"
OUTPUT_FOLDER = "Output/"
OUTPUT_PARSNP_FOLDER = "Output_parsnp/"

INGROUP_SIZE = 500
OUTGROUP_SIZE = 4

#Filtering parameters for the parsnp alignment
MINIMUM_ALIGNMENT_LENGTH = 300                   #Minimum sequence length in an alignment
MINIMUM_ALIGNMENT_COVERAGE = 0.99                #Minimum percentage of genomes that there should be an alignment values from [0,1]
MINIMUM_ALIGNMENT_PERCENTAGE_IDENTITY = 0.99  


#Add the location of the parsnp file to avoid re-run parsnp
#PARSNP_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output_parsnp/Ingroup_3_Outgroup_4/parsnp.xmfa"
#FILTERED_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output/filtered_test_2/filtered_parsnp.json"
#PARSNP_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output_parsnp/Ingroup_500_Outgroup_500/parsnp.xmfa"
PARSNP_XMFA_FILE_LOCATION = ""
FILTERED_XMFA_FILE_LOCATION = ""
#FILTERED_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output/2024_06_19_16_25/reduced_filtered_parsnp.json"
#FILTERED_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output/2024_06_10_11_16/reduced_filtered_parsnp.json"
#FILTERED_XMFA_FILE_LOCATION = "/home/jazminromero/alignment/Output/2024_06_10_10_38/reduced_filtered_parsnp.json"


#Maximum sequence length for a sequence alignment to be blast
SEQUENCE_LENGTH = 300
SHIFT = 250
E_CUTOFF_OUTGROUP = 1
PERC_IDENTITY_OUTGROUP = 60
#Maximum percentage of matches (0-100) 0 = There should be zero hits, and 100 we allow all the hits
OUTGROUP_FILTER_PERCENTAGE = 25

#Format
FORMAT = "fasta"


#Removing older output files
CLEANUP_DAYS = 60