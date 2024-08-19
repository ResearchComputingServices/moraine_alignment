class Stats:
    """
    Class representing statistics for alignment analysis.
    """

    def __init__(self):
        """
        Initializes the statistics object.

        Attributes:
        - alignments_found_by_parsnp: Number of alignments found by Parsnp.
        - alignments_discarded_by_length: Number of alignments discarded due to length.
        - alignments_discarded_by_coverage: Number of alignments discarded due to coverage.
        - alignments_discarded_by_percentage_of_identity: Number of alignments discarded due to percentage of identity.
        - total_alignments_kept: Total number of alignments kept.
        - total_alignments_taken: Total number of alignments taken.
        - percentage_discarded_by_length: Percentage of alignments discarded due to length.
        - percentage_discarded_by_coverage: Percentage of alignments discarded due to coverage.
        - percentage_discarded_by_percentage_of_identity: Percentage of alignments discarded due to percentage of identity.
        - percentage_alignment_kept: Percentage of alignments kept.
        - percentage_alignment_taken: Percentage of alignments taken.
        - total_subsequences_from_alignments: Total number of subsequences from alignments.
        - total_candidates: Total number of candidates.
        - percentage_candidates: Percentage of candidates.
        - top_five_hitted_outgroups_index: List of indices of the top five hit outgroups.
        - top_five_hitted_outgroups_filenames: List of filenames of the top five hit outgroups.
        - start_time: Start time of the process.
        - end_time: End time of the process.
        - parsnp_runtime: Runtime of Parsnp.
        - filtering_runtime: Runtime of filtering.
        - blast_runtime: Runtime of BLAST.
        - total_runtime: Total runtime.
        """

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
        """
        Sets the number of alignments found by Parsnp.

        Keyword arguments:
            total_alignments: The total number of alignments found by Parsnp.

        Returns:
            None
        """
        self.alignments_found_by_parsnp = total_alignments

    def set_alignments_discarded_by_length(self, count):
        """
        Sets the number of alignments discarded by length.

        Keyword arguments:
            count (int): The number of alignments discarded by length.

        Returns:
            None
        """

        self.alignments_discarded_by_length = count
        if self.alignments_found_by_parsnp > 0:
            self.percentage_discarded_by_length = (
                self.alignments_discarded_by_length / self.alignments_found_by_parsnp
            )

    def set_alignments_discarded_by_coverage(self, count):
        """
        Sets the number of alignments discarded by coverage.

        Keyword arguments:
            count (int): The number of alignments discarded by coverage.

        Returns:
            None
        """

        if self.alignments_found_by_parsnp > 0:
            self.percentage_discarded_by_coverage = (
                self.alignments_discarded_by_coverage / self.alignments_found_by_parsnp
            )

        self.alignments_discarded_by_coverage = count
        if self.alignments_found_by_parsnp > 0:
            self.percentage_discarded_by_coverage = (
                self.alignments_discarded_by_coverage / self.alignments_found_by_parsnp
            )

    def set_alignments_discared_by_percentage_of_identity(self, count):
        """
        Sets the number of alignments discarded by percentage of identity.

        Keyword arguments:
            count (int): The count of alignments discarded.

        Returns:
            None
        """

        self.alignments_discared_by_percentage_of_identity = count
        if self.alignments_found_by_parsnp > 0:
            self.percentage_discared_by_percentage_of_identity = (
                self.alignments_discared_by_percentage_of_identity
                / self.alignments_found_by_parsnp
            )

    def compute_total_alignnments_kept(self):
        """
        Computes the total number of alignments kept based on the number of alignments discarded by length,
        coverage, and percentage of identity. It also calculates the percentage of alignments taken and
        the percentage of alignments kept.

        Returns:
            None
        """

        self.total_alignments_taken = (
            self.alignments_discarded_by_length
            + self.alignments_discarded_by_coverage
            + self.alignments_discared_by_percentage_of_identity
        )
        self.total_alignments_kept = (
            self.alignments_found_by_parsnp - self.total_alignments_taken
        )
        if self.alignments_found_by_parsnp > 0:
            self.percentage_alignment_taken = (
                self.total_alignments_taken / self.alignments_found_by_parsnp
            )
            self.percentage_alignment_kept = (
                self.total_alignments_kept / self.alignments_found_by_parsnp
            )

    def compute_candidates_percentage(self, candidates_count):
        """
        Computes the percentage of candidates based on the given candidates count.

        Keyword arguments:
            candidates_count (int): The total count of candidates.

        Returns:
            None
        """
        self.total_candidates = candidates_count
        if self.total_subsequences_from_alignments > 0:
            self.percentage_candidates = 0
