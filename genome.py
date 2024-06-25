import uuid
import os
from Bio import SeqIO
import logging



class Subsequence:
    def __init__(self, value: str, alignment_id, position:int, alignment_position_in_genome:int):
        self.id = str(uuid.uuid4())
        self.value = value
        self.alignment_id = alignment_id
        self.position = position
        self.blast_score = 0
        self.outgroup_hits={}
        if alignment_position_in_genome!=-1:
            self.position_in_genome = alignment_position_in_genome + position
        else:
            self.position_in_genome = -1



class Genome:
    def __init__(self, filename=None, format=None, type=None):
        self.id = str(uuid.uuid4())
        self.filename = ""
        self.genome_sequence = ""
        self.type=""
        self.description=""
        if type!=None:
            self.type=type

        if filename!=None and format!=None:
            self.load_from_file(filename=filename, format=format)



    def load_from_file(self, filename: str, format:str, type=None):
        success = True
        if os.path.isfile(filename):
            try:
                name = os.path.basename(filename)
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
                self.compute_genome_sequence(seq_records=seq_records)
                self.filename= name
                if type!=None:
                    self.type=type

            except Exception as e:
                logging.info(e)
                success = False
        else:
            logging.info ("Filename %s does not exists.", filename)
            success = False
        return success
    
    '''A genome sequence is returned by records by biopython'''
    '''Here, the whole sequence in stored in one variable.'''
    def compute_genome_sequence(self, seq_records):
        sequence=""
        for rec in seq_records:
            sequence = sequence + str(rec.seq)
        self.genome_sequence = sequence




class Alignment:
    def __init__(self, sequence=None, 
                 description=None, 
                 strand=None, 
                 type=None, 
                 length=None, 
                 format=None,
                 genome_filename=None,
                 genome_header=None,
                 genome_length=None,
                 percentage_identity = None):
        
        self.id = str(uuid.uuid4())
        self.sequence = sequence
        self.subsequences = []
        self.original_subsequences_count = 0
        self.description = description
        self.strand = strand
        self.type=""
        self.length=length
        self.position_in_genome = 0
        self.percentage_of_identity = percentage_identity
    
        self.genome_filename = genome_filename
        self.genome_header = genome_header
        self.genome_length = genome_length
        
        if type!=None:
            self.type=type
        
        if self.description!=None:
            self.get_position_in_genome()


    def load_from_file(self, filename: str, format:str, type=None):
        success = True
        if os.path.isfile(filename):
            try:
                name = os.path.basename(filename)
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
                self.compute_sequence(seq_records=seq_records)
                self.filename= name
                if type!=None:
                    self.type=type

            except Exception as e:
                logging.info(e)
                success = False
        else:
            logging.info ("Filename %s does not exists.", filename)
            success = False
        return success


    '''A genome sequence is returned by records by biopython'''
    '''Here, the whole sequence in stored in one variable.'''
    def compute_sequence(self, seq_records):
        sequence=""
        for rec in seq_records:
            sequence = sequence + str(rec.seq)
        self.sequence = sequence


    '''Computes all subsequences of length subseq_len with a shift displacement'''
    '''Each subsequence is a string'''
    def compute_all_subsequences(self, subseq_len: int, shift: int):
        
        subseq_list = []
        seq_len = len(self.sequence)
        index = 0
        
        while (index <= seq_len - subseq_len):
            subseq_value = self.sequence[index:index + subseq_len]
            subseq = Subsequence(value=subseq_value, alignment_id=self.id, position=index, alignment_position_in_genome=self.position_in_genome) 
            subseq_list.append(subseq)
            index = index + shift

        #Getting the last ending subsequence 
        if index < seq_len:
            subseq_value = self.sequence[seq_len-subseq_len:seq_len]
            subseq = Subsequence(value=subseq_value, alignment_id=self.id, position=seq_len-subseq_len, alignment_position_in_genome=self.position_in_genome) 
            subseq_list.append(subseq)
            
            
        self.subsequences = subseq_list
        self.original_subsequences_count = len(self.subsequences)
    

    def get_position_in_genome(self):
        try:
            p = self.description.find("p")
            e = self.description.find("/")
            position = self.description[p+1:e]
            self.position_in_genome= int(position)
        except:
            self.position_in_genome = -1