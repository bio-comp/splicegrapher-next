"""Compatibility shim for legacy imports from :mod:`SpliceGrapher.shared.utils`.

Use focused modules under ``SpliceGrapher.shared`` for new code.
"""

from SpliceGrapher.shared.collection_utils import asList, asSet, bsearch
from SpliceGrapher.shared.config_utils import configMap, getEnvironmentValue
from SpliceGrapher.shared.file_utils import (
    ezopen,
    fileLen,
    filePrefix,
    findFile,
    makeGraphListFile,
    validateDir,
    validateFile,
)
from SpliceGrapher.shared.format_utils import (
    commaFormat,
    dictString,
    listString,
    substringAfter,
    substringBefore,
    substringBetween,
    timestamp,
    timeString,
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
    "asList",
    "asSet",
    "bsearch",
    "commaFormat",
    "configMap",
    "dictString",
    "ezopen",
    "fileLen",
    "filePrefix",
    "findFile",
    "getAttribute",
    "getEnvironmentValue",
    "idFactory",
    "listString",
    "logMessage",
    "makeGraphListFile",
    "process_fastq_header",
    "process_fasta_header",
    "process_labeled_fasta_header",
    "runCommand",
    "substringAfter",
    "substringBefore",
    "substringBetween",
    "timestamp",
    "timeString",
    "to_numeric",
    "validateDir",
    "validateFile",
    "writeStartupMessage",
    "ProgressIndicator",
    "RandomListIterator",
]
