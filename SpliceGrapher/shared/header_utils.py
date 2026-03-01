"""Header parsing helpers extracted from shared.utils."""


def _first_header_token(header: str) -> str:
    """Return the first token in a header string."""
    tokens = header.split()
    if not tokens:
        raise ValueError("Header is empty")
    return tokens[0]


def process_fastq_header(header: str) -> str:
    """Find the chromosome identifier within a FASTQ header."""
    return _first_header_token(header)


def process_fasta_header(header: str) -> str:
    """Find the chromosome identifier within a FASTA header."""
    return _first_header_token(header)


def process_labeled_fasta_header(header: str) -> tuple[str, str]:
    """
    Extract a sequence ID and its label from a fasta file.
    Used with training data FASTA files, so it assumes
    headers have the form "ID label=#".
    """
    sequence_id = _first_header_token(header).strip()
    label_token = next((token for token in header.split()[1:] if token.startswith("label=")), None)
    if label_token is None:
        raise ValueError("Header missing label= token")

    _, label = label_token.split("=", 1)
    if not label:
        raise ValueError("Header has empty label value")
    return sequence_id, label
