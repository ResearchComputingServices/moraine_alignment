DATA_FOLDER = "/data"
MULTIFASTA_FOLDER = "Multifasta/"
TEMP_FOLDER = "Temp/"
INGROUP_FOLDER = "Ingroup_Echoli_1000"
OUTGROUP_FOLDER = "Outgroup_Wong/"
OUTPUT_FOLDER = "Output/"
OUTPUT_PARSNP_FOLDER = "Output_parsnp/"

INGROUP_SIZE = 500
OUTGROUP_SIZE = 4

#Filtering parameters for the parsnp alignment
MINIMUM_ALIGNMENT_LENGTH = 150                   #Minimum sequence length in an alignment
MINIMUM_ALIGNMENT_COVERAGE = 0.95                #Minimum percentage of genomes that there should be an alignment values from [0,1]
MINIMUM_ALIGNMENT_PERCENTAGE_IDENTITY = 0.99   


#Add the location of the parsnp file to avoid re-run parsnp
PARSNP_XMFA_FILE_LOCATION = ""
