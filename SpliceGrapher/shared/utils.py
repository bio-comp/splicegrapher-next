"""Compatibility shim for legacy imports from :mod:`SpliceGrapher.shared.utils`.

Use focused modules under ``SpliceGrapher.shared`` for new code.
"""

from SpliceGrapher.shared.collection_utils import as_list, as_set, binary_search
from SpliceGrapher.shared.config_utils import configMap, getEnvironmentValue
from SpliceGrapher.shared.file_utils import (
    ez_open,
    file_len,
    file_prefix,
    find_file,
    make_graph_list_file,
    validate_dir,
    validate_file,
)
from SpliceGrapher.shared.format_utils import (
    comma_format,
    dict_string,
    list_string,
    substring_after,
    substring_before,
    substring_between,
    time_string,
    timestamp,
    to_numeric,
)
from SpliceGrapher.shared.header_utils import (
    process_fasta_header,
    process_fastq_header,
    process_labeled_fasta_header,
)
from SpliceGrapher.shared.process_utils import (
    getAttribute,
    idFactory,
    logMessage,
    runCommand,
    writeStartupMessage,
)
from SpliceGrapher.shared.progress import ProgressIndicator, RandomListIterator

__all__ = [
    "as_list",
    "as_set",
    "binary_search",
    "comma_format",
    "configMap",
    "dict_string",
    "ez_open",
    "file_len",
    "file_prefix",
    "find_file",
    "getAttribute",
    "getEnvironmentValue",
    "idFactory",
    "list_string",
    "logMessage",
    "make_graph_list_file",
    "process_fastq_header",
    "process_fasta_header",
    "process_labeled_fasta_header",
    "runCommand",
    "substring_after",
    "substring_before",
    "substring_between",
    "timestamp",
    "time_string",
    "to_numeric",
    "validate_dir",
    "validate_file",
    "writeStartupMessage",
    "ProgressIndicator",
    "RandomListIterator",
]
