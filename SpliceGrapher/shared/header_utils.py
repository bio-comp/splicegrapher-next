"""Header parsing helpers extracted from shared.utils."""


def process_fastq_header(header):
    """Find the chromosome identifier within a FASTQ header."""
    return header.split()[0]


def process_fasta_header(header):
    """Find the chromosome identifier within a FASTA header."""
    return header.split()[0]


def process_labeled_fasta_header(header):
    """
    Extract a sequence ID and its label from a fasta file.
    Used with training data FASTA files, so it assumes
    headers have the form "ID label=#".
    """
    id, labelToken = header.split()
    return id.strip(), labelToken.split("=")[1]
