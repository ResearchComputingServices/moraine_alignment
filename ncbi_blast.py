import shlex, subprocess
import os
from utils import get_files_in_dir_recursive, clear_file_paths, get_sublist, create_dir, delete_files, get_today_datetime, get_files_in_dir
import shutil
import time
from Bio.Blast.Applications import NcbiblastnCommandline



'''Calls ncbi_blast using biopython wrapper'''
def ncbi_blast_bio(
    query_file:str,
    subject_file:str, 
    e_cutoff: float, 
    identity_perc_cutoff: float, 
    max_hsps=None):
    
    results = None 
   
    delete_file = False
    try:
        if not os.path.isfile(subject_file) or not os.path.isfile(query_file):
            print("File %s does not exists OR", subject_file)
            print("File %s does not exists ", query_file)
            return results
    except Exception as exc:
        print("Error accessing file %s", subject_file)
        print(exc)

    #Execute blast
    try:
        if max_hsps==None:
            output = NcbiblastnCommandline(query=query_file, subject=subject_file, evalue=e_cutoff, perc_identity=identity_perc_cutoff, outfmt=5, task="blastn-short")()[0]
        else:
            output = NcbiblastnCommandline(query=query_file, subject=subject_file, evalue=e_cutoff, perc_identity=identity_perc_cutoff, max_hsps=max_hsps, outfmt=5, task="blastn-short")()[0]
    
        print(output)
        #Uncomment to try any of these alternatives to call when Ncbiblast is not longer available
        #output2 = blastn_wrapper_cl(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)
        #output3 = blastn_wrapper_sub(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)

    
        #results = __parse_blast_xml_output(output, subject_id)
        #Uncomment these to parse the results obtained by the alternatives wrappers
        #results2 = __parse_blast_xml_output(output2, subject_id)
        #results3 = __parse_blast_xml_output(output3, subject_id)

    
    except Exception as exc:
        print("Error when executing blast")
        print(exc)

    #Delete subject file if created
    try:
         if delete_file and os.path.isfile(subject_file):
             os.remove(subject_file)
    except Exception as exc:
         print("Could not delete file %s", subject_file)
       
    return results



'''Alternative Blast wrapper using os.popen'''
def ncbi_blast_wrapper_cl(query_file:str, subject_file:str,  e_cutoff: float, identity_perc_cutoff: float, max_hsps: int):
    
    try:
        if max_hsps==None:
            command = "blastn -outfmt 5 -query " + query_file +  " -evalue " + str(e_cutoff) + " -subject " + subject_file + " -task=blastn-short " + " -perc_identity " + str(identity_perc_cutoff) 
        else:
            command = "blastn -outfmt 5 -query " + query_file +  " -evalue " + str(e_cutoff) + " -subject " + subject_file + " -task=blastn-short " + " -perc_identity " + str(identity_perc_cutoff) + " -max_hsps " + str(max_hsps) 
            
        output = os.popen(command).read()
    except Exception as e:
        #logging.error(e)
        print(e)
    return output


'''Alternative Blast wrapper using subprocess: to be used when biopython wrapper is removed'''
def ncbi_blast_wrapper_sub(query_file:str, subject_file:str,  e_cutoff: float, identity_perc_cutoff: float, max_hsps:int):
    
    try:
        if max_hsps==None:
            command = "blastn -outfmt 5 -query " + query_file +  " -evalue " + str(e_cutoff) + " -subject " + subject_file + " -task=blastn-short " + " -perc_identity " + str(identity_perc_cutoff)
        else:
            command = "blastn -outfmt 5 -query " + query_file +  " -evalue " + str(e_cutoff) + " -subject " + subject_file + " -task=blastn-short " + " -perc_identity " + str(identity_perc_cutoff) + " -max_hsps " + str(max_hsps) 
            
            
        args = shlex.split(command)
        result = subprocess.run(args, capture_output=True, text=True)
        output = result.stdout.strip()
    except Exception as e:
        #logging.error(e)
        print(e)
    return output



def run_blast():
    query_file = "/data/Test/my_query.fasta"
    subject_file = "/data/Test/my_subject.fasta"

    ncbi_blast_bio(query_file=query_file, subject_file=subject_file, e_cutoff=1E-40, identity_perc_cutoff=90, max_hsps=1)
    
        