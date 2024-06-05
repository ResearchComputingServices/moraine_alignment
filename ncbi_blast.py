import shlex, subprocess
import os
import time
from Bio.Blast.Applications import NcbiblastnCommandline
import logging
import concurrent.futures
from typing import List
from genome import Genome
from config import Config
import uuid
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
import io


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





'''Parser Blast output and creates a dictionary'''
'''The output is a dictionary'''
def __parse_blast_xml_output(output:str, subject_id:str):

    results={}
    index = 0

    for blast_result_record in NCBIXML.parse(io.StringIO(output)):
        try:
            list_hits = []

            for alignment in blast_result_record.alignments:
                for hsp in alignment.hsps:
                        list_hits.append(hsp)


            query_id = blast_result_record.query.split()
            genome_id = query_id[0].split("_")[1]
            subseq_index = query_id[1].split("_")[1]
            subseq_id = query_id[2].split("_")[1]

            #We add the subsequence results only it got hits.
            d={}
            hits = {}
                
            if len(list_hits)>0:
                hits["query_id"]=subseq_id
                hits["hit_genome_id"]=subject_id
                hits["hit_score"]=list_hits[0].score
                hits["hit_str"]=str(list_hits[0])

            d["genome_id"]=genome_id
            d["subseq_index"]=subseq_index
            d["subseq_id"]=subseq_id
            d["hits"] = hits

            
            results[subseq_id] = d
        
        except Exception as e:
            logging.error(e)
            print(e)

        index = index + 1
    return results



'''Call blast wrapper from biopython to blast a multifasta query_file with a subject'''
'''Returns a blast result per sequence in the query_file'''
def ncbi_blast_multifasta(
    query_file:str,
    subject_file:str, 
    subject_sequence: str, 
    subject_id:str, 
    e_cutoff: float, 
    identity_perc_cutoff: float, 
    max_hsps=None,
    config_args: Config = None):
    
    results = None 



    delete_file = False
    
    #If the subject file is not provided 
    #Write the subject sequence to a file
    if subject_file=='' or not os.path.isfile(subject_file):
        if subject_sequence!="":
            #Create subject fasta file
            subject_file_name = str(uuid.uuid4()) + '.fasta'
            os.makedirs(config_args.TEMP_PATH, exist_ok=True)
            subject_file = os.path.join(config_args.TEMP_PATH,subject_file_name) 
            subject_sequence_record = SeqRecord(Seq(subject_sequence), id="subject")
            SeqIO.write(subject_sequence_record, subject_file, "fasta")
            delete_file = True
        else:
            return results

    #Check if the subject file exists
    if not os.path.isfile(subject_file) or not os.path.isfile(query_file):
        logging.error("File %s does not exists", subject_file)
        logging.error("File %s does not exists", query_file)
        print("File %s does not exists OR", subject_file)
        print("File %s does not exists ", query_file)
        return results

    #Execute blast
    if max_hsps==None:
        output = NcbiblastnCommandline(query=query_file, subject=subject_file, evalue=e_cutoff, perc_identity=identity_perc_cutoff, outfmt=5, task="blastn-short")()[0]
    else:
        output = NcbiblastnCommandline(query=query_file, subject=subject_file, evalue=e_cutoff, perc_identity=identity_perc_cutoff, max_hsps=max_hsps, outfmt=5, task="blastn-short")()[0]
    

    #Uncomment to try any of these alternatives to call when Ncbiblast is not longer available
    #output2 = blastn_wrapper_cl(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)
    #output3 = blastn_wrapper_sub(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)

    
    results = __parse_blast_xml_output(output, subject_id)
    #Uncomment these to parse the results obtained by the alternatives wrappers
    #results2 = __parse_blast_xml_output(output2, subject_id)
    #results3 = __parse_blast_xml_output(output3, subject_id)
    
    #Delete subject file if created
    if delete_file and os.path.isfile(subject_file):
        os.remove(subject_file)
       
    return results


'''Maps the results obtained by blast to the object in the class '''
def map_blast_results_to_genome(genome: Genome, blast_results:dict):

    if blast_results!=None:
        for subseq_id, result in blast_results.items():
            if result["hits"]!=None and len(result["hits"])>0:
                try:
                    index = int(result["subseq_index"])
                    genome_subject_id = result["hits"]["hit_genome_id"]
                    genome.subsequences[index].outgroup_hits[genome_subject_id] = result["hits"]    
                    if genome.subsequences[index].id!=subseq_id:
                        logging.error("Something went wrong")
                except Exception as e:
                    logging.info(e)






'''Blast the set of subsequences of a genome against a genome sequence'''
def blast_subsequences_against_genomes(genomes_query: List[Genome], 
                                       genomes_subject: List[Genome], 
                                       query_subseq_fasta_paths: List[str], 
                                       subject_sequence_fasta_paths: List[str],          
                                       config_args: Config = None):
    start = time.time()
    
    #cluster_run = False


    executor = concurrent.futures.ProcessPoolExecutor(max_workers=config_args.processors_number)

    for query_index  in range(0,len(genomes_query)):
        
        #print ("Genome {}".format(query_index))
        query_genome_id = genomes_query[query_index].id
        query_genome_description = genomes_query[query_index].description
        
        #Get the file with the subsequences filepaths corresponding to blast
        if query_genome_id not in query_subseq_fasta_paths:
            continue
        
        query_fasta_paths = query_subseq_fasta_paths[query_genome_id]
        
        subject_index = 0
        for subject_genome_id, subject_genome in genomes_subject.items():

            subject_fasta = subject_sequence_fasta_paths[subject_genome_id]
            subject_genome_filename = subject_genome.filename
            subject_genome_sequence = ""
            
            #We don't want to blast a genome against itself
            if (query_genome_id!=subject_genome_id):
                
                logging.info("Blast -> subs in Alignment-{}_{} - Outgroup-{}_{}".format(query_index,query_genome_description,subject_index,subject_genome_filename))
                #logging.info("Blast -> Query {} - Subject {} ".format(query_index,subject_genome_fasta_file))
                
                #for fp in file_paths:
                #    ncbi_blast_multifasta(query_file=fp, subject_file=subject_genome_fasta_file, subject_id=subject_genome_id, subject_sequence=subject_genome_sequence, e_cutoff=config_args.e_cutoff_outgroup,
                #                          identity_perc_cutoff=config_args.perc_identity_outgroup,max_hsps=config_args.max_hsps,config_args=config_args)


                processors_results = []
                future_results = [executor.submit(ncbi_blast_multifasta, fp, subject_fasta, subject_genome_sequence, subject_genome_id, config_args.e_cutoff_outgroup, config_args.perc_identity_outgroup, config_args.max_hsps) for fp in query_fasta_paths]
    
                for finished in concurrent.futures.as_completed(future_results, timeout=600):
                    try:
                        processors_results.append(finished.result())
                        #logging.info(finished.result())
                    except concurrent.futures._base.TimeoutError:
                        logging.error("Process took to long to complete")
                        #print("Process took to long to compete")
                    except Exception as exc:
                        logging.error("Exception occurred")
                        logging.error(exc)
                        #print ("Exception occurred: ")
                        #print (exc)

                for blast_results in processors_results:
                    map_blast_results_to_genome(genome=genomes_query[query_index], blast_results=blast_results)
                    #save_blast_result(blast_results=blast_results)


            subject_index  = subject_index + 1   
                        
    #logging.info("total json files: {}.".format(total_json_files))            
    end = time.time()
    mins = (end-start)/60
    
    logging.info("Blast mins {}".format(mins))
















def run_blast():
    query_file = "/data/Test/my_query.fasta"
    subject_file = "/data/Test/my_subject.fasta"

    ncbi_blast_bio(query_file=query_file, subject_file=subject_file, e_cutoff=1E-40, identity_perc_cutoff=90, max_hsps=1)
    
        