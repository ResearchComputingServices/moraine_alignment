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


def __parse_blast_xml_output(output: str, subject_id: str) -> dict:
    """
    Parser Blast output and creates a dictionary.

    Keyword arguments:
        output (str): The XML output of the BLAST search.
        subject_id (str): The ID of the subject sequence.

    Returns:
        dict: A dictionary containing the parsed results. The keys are subsequence IDs and the values are dictionaries containing the following information:
            - genome_id (str): The ID of the genome.
            - subseq_index (str): The index of the subsequence.
            - subseq_id (str): The ID of the subsequence.
            - hits (dict): A dictionary containing information about the hits. The keys are as follows:
                - query_id (str): The ID of the query subsequence.
                - hit_genome_id (str): The ID of the subject genome.
                - hit_score (float): The score of the hit.
                - hit_str (str): A string representation of the hit.
    """
    results = {}
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

            # We add the subsequence results only it got hits.
            d = {}
            hits = {}

            if len(list_hits) > 0:
                hits["query_id"] = subseq_id
                hits["hit_genome_id"] = subject_id
                hits["hit_score"] = list_hits[0].score
                hits["hit_str"] = str(list_hits[0])

            d["genome_id"] = genome_id
            d["subseq_index"] = subseq_index
            d["subseq_id"] = subseq_id
            d["hits"] = hits

            results[subseq_id] = d

        except Exception as e:
            logging.error(e)
            print(e)

        index = index + 1
    return results


def ncbi_blast_multifasta(
    query_file: str,
    subject_file: str,
    subject_sequence: str,
    subject_id: str,
    e_cutoff: float,
    identity_perc_cutoff: float,
    max_hsps: int | None,
    config_args: Config = None,
):
    """
    Call blast wrapper from biopython to blast a multifasta query_file with a subject
    Returns a blast result per sequence in the query_file

    Keyword arguments:
        query_file (str): Path to the query file.
        subject_file (str): Path to the subject file.
        subject_sequence (str): Subject sequence.
        subject_id (str): Subject ID.
        e_cutoff (float): E-value cutoff.
        identity_perc_cutoff (float): Identity percentage cutoff.
        max_hsps (int, optional): Maximum number of HSPs to report. Defaults to None.
        config_args (Config, optional): Configuration arguments. Defaults to None.

    Returns:
        results (object): BLAST search results.

    Raises:
        FileNotFoundError: If the subject file or query file does not exist.
    """
    results = None

    delete_file = False

    # If the subject file is not provided
    # Write the subject sequence to a file
    if subject_file == "" or not os.path.isfile(subject_file):
        if subject_sequence != "":
            # Create subject fasta file
            subject_file_name = str(uuid.uuid4()) + ".fasta"
            os.makedirs(config_args.TEMP_PATH, exist_ok=True)
            subject_file = os.path.join(
                config_args.TEMP_PATH, subject_file_name)
            subject_sequence_record = SeqRecord(
                Seq(subject_sequence), id="subject")
            SeqIO.write(subject_sequence_record, subject_file, "fasta")
            delete_file = True
        else:
            return results

    # Check if the subject file exists
    if not os.path.isfile(subject_file) or not os.path.isfile(query_file):
        logging.error("File %s does not exists", subject_file)
        logging.error("File %s does not exists", query_file)
        print("File %s does not exists OR", subject_file)
        print("File %s does not exists ", query_file)
        return results

    # Execute blast
    if max_hsps is None:
        output = NcbiblastnCommandline(
            query=query_file,
            subject=subject_file,
            evalue=e_cutoff,
            perc_identity=identity_perc_cutoff,
            outfmt=5,
            task="blastn-short",
        )()[0]
    else:
        output = NcbiblastnCommandline(
            query=query_file,
            subject=subject_file,
            evalue=e_cutoff,
            perc_identity=identity_perc_cutoff,
            max_hsps=max_hsps,
            outfmt=5,
            task="blastn-short",
        )()[0]

    # Uncomment to try any of these alternatives to call when Ncbiblast is not longer available
    # output2 = blastn_wrapper_cl(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)
    # output3 = blastn_wrapper_sub(query_file, subject_file,  e_cutoff, identity_perc_cutoff, max_hsps)

    results = __parse_blast_xml_output(output, subject_id)
    # Uncomment these to parse the results obtained by the alternatives wrappers
    # results2 = __parse_blast_xml_output(output2, subject_id)
    # results3 = __parse_blast_xml_output(output3, subject_id)

    # Delete subject file if created
    if delete_file and os.path.isfile(subject_file):
        os.remove(subject_file)

    return results


def map_blast_results_to_genome(genome: Genome, blast_results: dict):
    """
    Maps the results obtained by blast to the object in the class.

    Keyword arguments:
        genome (Genome): The genome object to map the blast results to.
        blast_results (dict): The blast results dictionary.

    Returns:
        None

    Raises:
        None
    """
    if blast_results is not None:
        for subseq_id, result in blast_results.items():
            if result["hits"] is not None and len(result["hits"]) > 0:
                try:
                    index = int(result["subseq_index"])
                    genome_subject_id = result["hits"]["hit_genome_id"]
                    genome.subsequences[index].outgroup_hits[genome_subject_id] = (
                        result["hits"]
                    )
                    if genome.subsequences[index].id != subseq_id:
                        logging.error("Something went wrong")
                except Exception as e:
                    logging.info(e)


def blast_multifasta_files_vs_subject(
    subject_info: dict, query_fasta_paths: List[str], config_args
):
    """
    Parallelization by outgroup
    Each processor will blast the subsequences of an alignment to the outgroup' sequence

    Keyword arguments:
        subject_info (dict): Information about the subject sequence, including the subject FASTA file and genome ID.
        query_fasta_paths (List[str]): List of file paths to the query FASTA files.
        config_args: Configuration arguments for the BLAST search.

    Returns:
        List: List of BLAST results for each query FASTA file.
    """

    subject_fasta = subject_info["subject_fasta"]
    subject_genome_sequence = ""
    subject_genome_id = subject_info["subject_genome_id"]

    query_blast_results = []
    for fp in query_fasta_paths:
        blast_result = ncbi_blast_multifasta(
            query_file=fp,
            subject_file=subject_fasta,
            subject_sequence=subject_genome_sequence,
            subject_id=subject_genome_id,
            e_cutoff=config_args.e_cutoff_outgroup,
            identity_perc_cutoff=config_args.perc_identity_outgroup,
            max_hsps=config_args.max_hsps,
        )
        query_blast_results.append(blast_result)

    return query_blast_results


def blast_query_subsequences_vs_outgroup(
    config_args: Config,
    query_fasta_paths: List[str],
    genomes_subject: dict,
    subject_sequence_fasta_paths: List[str],
):
    """
    Blast query subsequences against outgroup genomes.

    Keyword arguments:
        config_args (Config): Configuration arguments.
        query_fasta_paths (List[str]): List of file paths for query FASTA sequences.
        genomes_subject (dict): Dictionary of outgroup genomes.
        subject_sequence_fasta_paths (List[str]): List of file paths for subject sequence FASTA files.

    Returns:
        List: List of results from the blast process.
    """

    executor = concurrent.futures.ProcessPoolExecutor(
        max_workers=config_args.processors_number
    )

    processors_results = []
    processor_index = 0
    subject_index = 0
    processors_tasks = []
    total_subjects = len(genomes_subject)

    for subject_genome_id, subject_genome in genomes_subject.items():
        if (
            processor_index < config_args.processors_number
            and subject_index < total_subjects
        ):

            subject_fasta = subject_sequence_fasta_paths[subject_genome_id]
            subject_info = {}
            subject_info["subject_genome_id"] = subject_genome_id
            subject_info["subject_fasta"] = subject_fasta
            processors_tasks.append(subject_info)
            processor_index = processor_index + 1
            subject_index = subject_index + 1

            if (
                processor_index >= config_args.processors_number
                or subject_index >= total_subjects
            ):

                future_results = [
                    executor.submit(
                        blast_multifasta_files_vs_subject,
                        subject_info,
                        query_fasta_paths,
                        config_args,
                    )
                    for subject_info in processors_tasks
                ]

                for finished in concurrent.futures.as_completed(
                    future_results, timeout=600
                ):
                    try:
                        processors_results.append(finished.result())
                    except concurrent.futures._base.TimeoutError:
                        logging.error("Process took to long to complete")
                    except Exception as exc:
                        logging.error("Exception occurred")
                        logging.error(exc)

                processors_tasks = []
                processor_index = 0

    return processors_results


def blast_subsequences_against_genomes(
    genomes_query: List[Genome],
    genomes_subject: List[Genome],
    query_subseq_fasta_paths: List[str],
    subject_sequence_fasta_paths: List[str],
    config_args: Config = None,
):
    """
    Blast the set of subsequences of a genome against a genome sequence.
    Parallelized by Outgroup.

    Keyword arguments:
        genomes_query (List[Genome]): List of query genomes.
        genomes_subject (List[Genome]): List of subject genomes.
        query_subseq_fasta_paths (List[str]): List of file paths for query subsequences in FASTA format.
        subject_sequence_fasta_paths (List[str]): List of file paths for subject sequences in FASTA format.
        config_args (Config, optional): Configuration arguments. Defaults to None.

    Returns:
        None
    """
    start = time.time()

    logging.info(
        "Blast {} multifasta files vs {} outgroups \n".format(
            len(query_subseq_fasta_paths), len(subject_sequence_fasta_paths)
        )
    )
    for query_index in range(0, len(genomes_query)):

        query_genome_id = genomes_query[query_index].id

        # Get the file with the subsequences filepaths corresponding to blast
        if query_genome_id not in query_subseq_fasta_paths:
            continue

        query_fasta_paths = query_subseq_fasta_paths[query_genome_id]
        blast_results = blast_query_subsequences_vs_outgroup(
            config_args=config_args,
            query_fasta_paths=query_fasta_paths,
            genomes_subject=genomes_subject,
            subject_sequence_fasta_paths=subject_sequence_fasta_paths,
        )

        for blast_result_list in blast_results:
            for blast_result in blast_result_list:
                map_blast_results_to_genome(
                    genome=genomes_query[query_index], blast_results=blast_result
                )

    end = time.time()
    mins = (end - start) / 60
    config_args.stats.blast_runtime = config_args.stats.blast_runtime + mins

    logging.info("Blast mins {}".format(mins))
