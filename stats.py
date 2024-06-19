class Stats:
    def __init__(self):
        self.alignments_found_by_parsnp = 0
        self.alignments_discarded_by_length = 0
        self.alignments_discarded_by_coverage = 0 
        self.alignments_discared_by_percentage_of_identity = 0
        self.total_alignments_kept = 0
        self.total_alignments_taken = 0

        self.percentage_discarded_by_length = 0
        self.percentage_discarded_by_coverage = 0 
        self.percentage_discared_by_percentage_of_identity = 0
        self.percentage_alignment_kept = 0
        self.percentage_alignment_taken = 0
        

        self.total_subsequences_from_alignments = 0
        self.total_candidates = 0
        self.percentage_candidates = 0
        

        self.top_five_hitted_outgroups_index = []
        self.top_five_hitted_outgroups_filenames = []

        self.start_time = 0
        self.end_time = 0
        self.parsnp_runtime = 0
        self.filtering_runtime = 0
        self.blast_runtime = 0
        self.total_runtime = 0


    def set_alignments_found_by_parsnp(self, total_alignments):
            self.alignments_found_by_parsnp = total_alignments

    def set_alignments_discarded_by_length(self, count):
        self.alignments_discarded_by_length = count
        if self.alignments_found_by_parsnp>0: 
            self.percentage_discarded_by_length = self.alignments_discarded_by_length/self.alignments_found_by_parsnp
    
    def set_alignments_discarded_by_coverage(self, count):
        self.alignments_discarded_by_coverage = count
        if self.alignments_found_by_parsnp>0: 
            self.percentage_discarded_by_coverage = self.alignments_discarded_by_coverage/self.alignments_found_by_parsnp
    
    def set_alignments_discared_by_percentage_of_identity(self, count):
        self.alignments_discared_by_percentage_of_identity = count
        if self.alignments_found_by_parsnp>0: 
            self.percentage_discared_by_percentage_of_identity = self.alignments_discared_by_percentage_of_identity / self.alignments_found_by_parsnp
        

    def compute_total_alignnments_kept(self):
        self.total_alignments_taken = self.alignments_discarded_by_length+self.alignments_discarded_by_coverage+self.alignments_discared_by_percentage_of_identity
        self.total_alignments_kept = self.alignments_found_by_parsnp - self.total_alignments_taken
        if self.alignments_found_by_parsnp>0:
            self.percentage_alignment_taken = self.total_alignments_taken / self.alignments_found_by_parsnp
            self.percentage_alignment_kept = self.total_alignments_kept / self.alignments_found_by_parsnp

    
    def compute_candidates_percentage(self, candidates_count):
        self.total_candidates = candidates_count
        if self.total_subsequences_from_alignments>0:
            self.percentage_candidates = 0
        




