import uuid
import os
from Bio import SeqIO
import logging


class Subsequence:
    """
    Represents a subsequence in a genome alignment.
    
    Attributes:
        id (str): The unique identifier of the subsequence.
        value (str): The value of the subsequence.
        alignment_id: The identifier of the alignment the subsequence belongs to.
        position (int): The position of the subsequence.
        blast_score (int): The blast score of the subsequence.
        outgroup_hits (dict): A dictionary containing the outgroup hits of the subsequence.
        position_in_genome (int): The position of the subsequence in the genome alignment.
    """

    def __init__(
        self, value: str, alignment_id, position: int, alignment_position_in_genome: int
    ):
        self.id = str(uuid.uuid4())
        self.value = value
        self.alignment_id = alignment_id
        self.position = position
        self.blast_score = 0
        self.outgroup_hits = {}
        if alignment_position_in_genome != -1:
            self.position_in_genome = alignment_position_in_genome + position
        else:
            self.position_in_genome = -1


class Genome:
    """
    Represents a genome.
    """

    def __init__(self, filename=None, format=None, type=None):
        """
        Initializes a Genome object.
        
        Keyword arguments:
            filename (str, optional): The filename of the genome file. Defaults to None.
            format (str, optional): The format of the genome file. Defaults to None.
            type (str, optional): The type of the genome. Defaults to None.
        """

        self.id = str(uuid.uuid4())
        self.filename = ""
        self.genome_sequence = ""
        self.type = ""
        self.description = ""
        if type != None:
            self.type = type

        if filename != None and format != None:
            self.load_from_file(filename=filename, format=format)

    def load_from_file(self, filename: str, format: str, type=None) -> bool:
        """
        Load genome data from a file.
        
        Keyword arguments:
            filename (str): The path to the file containing the genome data.
            format (str): The format of the file.
            type (optional): The type of the genome data.
        
        Returns:
            bool: True if the genome data is successfully loaded, False otherwise.
        """

        success = True
        if os.path.isfile(filename):
            try:
                name = os.path.basename(filename)
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
                self.compute_genome_sequence(seq_records=seq_records)
                self.filename = name
                if type != None:
                    self.type = type

            except Exception as e:
                logging.info(e)
                success = False
        else:
            logging.info("Filename %s does not exists.", filename)
            success = False
        return success

    """A genome sequence is returned by records by biopython"""
    """Here, the whole sequence in stored in one variable."""

    def compute_genome_sequence(self, seq_records):
        sequence = ""
        for rec in seq_records:
            sequence = sequence + str(rec.seq)
        self.genome_sequence = sequence


class Alignment:
    """
    Represents an alignment object.
    """

    def __init__(
        self,
        sequence=None,
        description=None,
        strand=None,
        type=None,
        length=None,
        format=None,
        genome_filename=None,
        genome_header=None,
        genome_length=None,
        percentage_identity=None,
    ):
        """
        Initialize a Genome object.
        
        Keyword arguments:
            sequence (str, optional): The sequence of the genome. Defaults to None.
            description (str, optional): The description of the genome. Defaults to None.
            strand (str, optional): The strand of the genome. Defaults to None.
            type (str, optional): The type of the genome. Defaults to None.
            length (int, optional): The length of the genome. Defaults to None.
            format (str, optional): The format of the genome. Defaults to None.
            genome_filename (str, optional): The filename of the genome. Defaults to None.
            genome_header (str, optional): The header of the genome. Defaults to None.
            genome_length (int, optional): The length of the genome. Defaults to None.
            percentage_identity (float, optional): The percentage identity of the genome. Defaults to None.
        """

        self.id = str(uuid.uuid4())
        self.sequence = sequence
        self.subsequences = []
        self.original_subsequences_count = 0
        self.description = description
        self.strand = strand
        self.type = ""
        self.length = length
        self.position_in_genome = 0
        self.percentage_of_identity = percentage_identity

        self.genome_filename = genome_filename
        self.genome_header = genome_header
        self.genome_length = genome_length

        if type != None:
            self.type = type

        if self.description != None:
            self.get_position_in_genome()

    def load_from_file(self, filename: str, format: str, type=None) -> bool:
        """
        Load genome data from a file.
        
        Keyword arguments:
            filename (str): The path to the file containing the genome data.
            format (str): The format of the file.
            type (Optional): The type of the genome data.
        
        Returns:
            bool: True if the genome data is successfully loaded, False otherwise.
        """

        success = True
        if os.path.isfile(filename):
            try:
                name = os.path.basename(filename)
                logging.info("Reading genome from file: " + filename)
                seq_iterator = SeqIO.parse(filename, format)
                seq_records = list(seq_iterator)
                self.compute_sequence(seq_records=seq_records)
                self.filename = name
                if type != None:
                    self.type = type

            except Exception as e:
                logging.info(e)
                success = False
        else:
            logging.info("Filename %s does not exists.", filename)
            success = False
        return success

    def compute_sequence(self, seq_records):
        """
        A genome sequence is returned by records by biopython
        Here, the whole sequence in stored in one variable.
        
        Keyword arguments:
            seq_records (list): A list of sequence records.
        
        Returns:
            None
        
        Raises:
            None
        """

        sequence = ""
        for rec in seq_records:
            sequence = sequence + str(rec.seq)
        self.sequence = sequence

    def compute_all_subsequences(self, subseq_len: int, shift: int):
        """
        Computes all subsequences of length subseq_len with a shift displacement
        Each subsequence is a string
        
        Keyword arguments:
            subseq_len (int): The length of the subsequences to be computed.
            shift (int): The amount of shift between each subsequence.
        
        Returns:
            None
        
        Raises:
            None
        """

        subseq_list = []
        seq_len = len(self.sequence)
        index = 0

        while index <= seq_len - subseq_len:
            subseq_value = self.sequence[index : index + subseq_len]
            subseq = Subsequence(
                value=subseq_value,
                alignment_id=self.id,
                position=index,
                alignment_position_in_genome=self.position_in_genome,
            )
            subseq_list.append(subseq)
            index = index + shift

        # Getting the last ending subsequence
        if index < seq_len:
            subseq_value = self.sequence[seq_len - subseq_len : seq_len]
            subseq = Subsequence(
                value=subseq_value,
                alignment_id=self.id,
                position=seq_len - subseq_len,
                alignment_position_in_genome=self.position_in_genome,
            )
            subseq_list.append(subseq)

        self.subsequences = subseq_list
        self.original_subsequences_count = len(self.subsequences)

    def get_position_in_genome(self):
        """
        Retrieves the position of the genome.
        
        Returns:
            int: The position of the genome in the description.
        """

        try:
            p = self.description.find("p")
            e = self.description.find("/")
            position = self.description[p + 1 : e]
            self.position_in_genome = int(position)
        except:
            self.position_in_genome = -1
