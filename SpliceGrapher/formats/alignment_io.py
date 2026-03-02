"""Alignment I/O helpers for SAM/BAM/CRAM inputs."""

import io
import os
import re
import sys

import pysam
import structlog

from SpliceGrapher.core.enums import SamHeaderLine, SamHeaderTag
from SpliceGrapher.formats.shortread_compat import (
    SpliceJunction,
    is_depths_file,
    read_depths,
)
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.process_utils import getAttribute
from SpliceGrapher.shared.progress import ProgressIndicator

LOGGER = structlog.get_logger(__name__)

# Header tags
HEADER_HD_TAG = SamHeaderTag.HD
HEADER_LN_TAG = SamHeaderTag.LN
HEADER_SN_TAG = SamHeaderTag.SN
HEADER_SO_TAG = SamHeaderTag.SO
HEADER_SQ_TAG = SamHeaderTag.SQ
HEADER_VN_TAG = SamHeaderTag.VN
HEADER_SQ_LINE = SamHeaderLine.SQ

# SAM files have the following tab-delimited columns:
# Column  Value
# 0       QNAME -- Query pair NAME if paired; or Query NAME if unpaired
# 1       FLAG  -- bitwise FLAG (Sam spec., section 2.2.2)
# 2       RNAME -- Reference sequence NAME 3
# 3       POS   -- 1-based leftmost POSition/coordinate of the clipped sequence
# 4       MAPQ  -- MAPping Quality (phred-scaled posterior probability that the
#                  mapping position of this read is incorrect)
# 5       CIGAR -- extended CIGAR string
# 6       MRNM  -- Mate Reference sequence NaMe
# 7       MPOS  -- 1-based leftmost Mate POSition of the clipped sequence
# 8       ISIZE -- inferred Insert SIZE
# 9       SEQ   -- query SEQuence; = for match to the reference; n/N/. for
#                  ambiguity; cases are not maintained
# 10      QUAL  -- query QUALity; ASCII-33 gives the Phred base quality
# 11-??   Optional TAG:VTYPE:VALUE triplets
#
# A -- Printable character
# i -- Signed 32-bit integer
# f -- Single-precision floating number
# Z -- Printable string, including space
# H -- Byte array in the Hex format
# B -- Integer or numeric array
VALID_VTYPES = ["A", "i", "f", "Z", "H", "B"]

QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE, SEQ, QUAL, TAGS = range(12)
ALL_COLUMNS = [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, MRNM, MPOS, ISIZE, SEQ, QUAL, TAGS]
REQUIRED_COLUMNS = ALL_COLUMNS[:-1]
INT_COLUMNS = [FLAG, POS, MPOS, MAPQ]

# Tag used by Cufflinks and SpliceGrapher to identify alignment strand
STRAND_TAG = "XS"

# Tag used by SpliceGrapher SAM files to identify junction type
JCT_CODE_TAG = "YC"
KNOWN_JCT_TAG = "YC:A:K"
RECOMBINED_JCT_TAG = "YC:A:U"
PREDICTED_JCT_TAG = "YC:A:P"

READ_LEN_TAG = "NM"

# Match patterns of valid CIGAR symbol sequences
MATCH_RE = re.compile("^M(NM)*$")
CIGAR_MATCH = re.compile("^([0-9]+[SH])?[0-9]+[MSH]([0-9]+[NID][0-9]+[MSH])*([0-9]+[SH])?$")
CIGAR_TOKEN = re.compile("[0-9]+[A-Z]")

PYSAM_CIGAR = "MIDNSHP=X"

# Defined in pysam documentation, but evidently
# not in the Python code.
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8

# Our own lists for convenience
PYSAM_MERGE_OPS = [BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP]
PYSAM_IGNORE_OPS = [BAM_CDEL, BAM_CPAD]
PYSAM_SAVE_OPS = [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF, BAM_CREF_SKIP]

NULL_CIGAR = "*"
NULL_CHROMOSOME = "*"
EXACT_CIGAR = "="

# Flag values used to assess/interpret records
UNMAPPED_FLAG = 4
REVERSE_FLAG = 16


class ChromosomeTracker(object):
    def __init__(self):
        self.chromosome = None
        self.finished = {}

    def allFinished(self, cSet):
        """Return True if all chromosomes in the given set/list are finished."""
        for c in cSet:
            if not self.isFinished(c):
                return False
        return True

    def isFinished(self, c):
        """Return True if the given chromosome has already been seen."""
        try:
            return self.finished[c.lower()]
        except KeyError:
            return False

    def getFinished(self):
        """Returns a list of all chromosomes that have been seen."""
        return [c for c in self.finished if self.finished[c]]

    def finish(self):
        """Tells the tracker to finish off the last chromosome."""
        if self.chromosome:
            self.finished[self.chromosome] = True

    def update(self, c):
        """Main feature of this class: when it detects a change from one
        chromosome to another, marks the old one as finished and starts a new one."""
        chrom = c.lower()
        if self.chromosome == chrom:
            return False
        # Mismatch:
        if self.chromosome:
            self.finished[self.chromosome] = True
        self.finished[chrom] = False
        self.chromosome = chrom
        return True


def _is_alignment_path(source):
    """Return ``True`` when ``source`` is a filesystem path to an alignment file."""
    return isinstance(source, str) and os.path.isfile(source)


def isCramFile(filePath):
    """Simple heuristic returns True if path is a CRAM file; false otherwise."""
    return filePath.lower().endswith(".cram")


def _open_alignment_file(path, **args):
    """Open SAM/BAM/CRAM input with pysam and CRAM reference safeguards."""
    reference_fasta = getAttribute("reference_fasta", None, **args)
    if isBamFile(path):
        return pysam.AlignmentFile(path, "rb")
    if isCramFile(path):
        try:
            if reference_fasta:
                return pysam.AlignmentFile(path, "rc", reference_filename=str(reference_fasta))
            return pysam.AlignmentFile(path, "rc")
        except ValueError as exc:
            if "reference" in str(exc).lower():
                raise ValueError(
                    "Unable to decode CRAM without reference. "
                    "Re-run with reference_fasta=<path-to-fasta>."
                ) from exc
            raise
    return pysam.AlignmentFile(path, "r")


def _tokenize_cigar_tuples(cigar_tuples):
    """Convert pysam CIGAR tuples into normalized depth tokens.

    The token stream uses:
    - ``M`` for match/equal/diff (reference-consuming aligned blocks),
    - ``N`` for reference skips (splice junction introns),
    - ``D`` for deletions (reference-consuming but not aligned sequence).

    Query-only operations (``I/S/H``) and pads (``P``) are ignored for
    reference-coordinate depth stepping.
    """
    if not cigar_tuples:
        return []

    tokens = []
    for op, size in cigar_tuples:
        if op in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]:
            code = "M"
        elif op == BAM_CREF_SKIP:
            code = "N"
        elif op == BAM_CDEL:
            code = "D"
        elif op in [BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD]:
            continue
        else:
            continue

        if tokens and tokens[-1][-1] == code:
            prev_size = int(tokens[-1][:-1])
            tokens[-1] = f"{prev_size + size}{code}"
        else:
            tokens.append(f"{size}{code}")
    return tokens


def _pysam_record_to_tokens(rec):
    """Return normalized CIGAR tokens for a pysam aligned segment."""
    if rec.cigartuples is None:
        return []
    tokens = _tokenize_cigar_tuples(rec.cigartuples)
    if not tokens:
        read_len = rec.query_length or 0
        if read_len > 0:
            return [f"{read_len}M"]
    return tokens


def _next_match_anchor(cigar_tuples, start_idx):
    """Return the contiguous downstream match length used as junction anchor."""
    anchor = 0
    for op, size in cigar_tuples[start_idx:]:
        if op in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]:
            anchor += size
            continue
        if op in [BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD]:
            continue
        break
    return anchor


def _build_splice_junctions(chromosome, start_pos, cigar_tuples, strand, jct_code=""):
    """Build splice junctions directly from pysam CIGAR tuples.

    Junction coordinates are 1-based closed intervals:
    - donor is the last reference base before an ``N`` skip
    - acceptor is the first reference base after the skip
    """
    if not cigar_tuples:
        return []

    junctions = []
    reference_pos = start_pos
    left_anchor = 0

    for idx, (op, size) in enumerate(cigar_tuples):
        if op in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]:
            left_anchor += size
            reference_pos += size
            continue

        if op in [BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD]:
            continue

        if op == BAM_CDEL:
            reference_pos += size
            left_anchor = 0
            continue

        if op != BAM_CREF_SKIP:
            continue

        donor = reference_pos - 1
        reference_pos += size
        acceptor = reference_pos
        right_anchor = _next_match_anchor(cigar_tuples, idx + 1)
        if left_anchor > 0 and right_anchor > 0:
            junctions.append(
                SpliceJunction(
                    chromosome,
                    donor,
                    acceptor,
                    [left_anchor, right_anchor],
                    jct_code,
                    strand,
                )
            )
        left_anchor = 0

    return junctions


def _collect_pysam_data(path, **args):
    """Collect depths/junctions/alignment spans using pysam-backed iteration."""
    alignments = getAttribute("alignments", False, **args)
    chromosomes = getAttribute("chromosomes", None, **args)
    chrom_set = makeChromosomeSet(chromosomes)
    get_junctions = getAttribute("junctions", True, **args)
    maxpos = getAttribute("maxpos", sys.maxsize, **args)
    minanchor = getAttribute("minanchor", 0, **args)
    minjct = getAttribute("minjct", 1, **args)

    align = {}
    depths = {}
    junction_tmp = {}
    limit = {}

    with _open_alignment_file(path, **args) as alignment:
        for rec in alignment.fetch(until_eof=True):
            if rec.is_unmapped:
                continue

            chrom = alignment.get_reference_name(rec.reference_id).lower()
            if chrom_set and chrom not in chrom_set:
                continue

            start_pos = rec.reference_start + 1
            if start_pos > maxpos:
                continue

            if chrom not in depths:
                depths[chrom] = [0] * (maxpos + 1) if maxpos < sys.maxsize else [0]
            if chrom not in limit:
                limit[chrom] = 0
            if alignments and chrom not in align:
                align[chrom] = []

            tokens = _pysam_record_to_tokens(rec)
            if not tokens:
                continue

            cur_pos = start_pos
            for token in tokens:
                code = token[-1]
                delta = int(token[:-1])
                nxt_pos = cur_pos + delta
                while nxt_pos + 1 > len(depths[chrom]):
                    depths[chrom] += [0] * len(depths[chrom])

                if code == "M":
                    for i in range(cur_pos, nxt_pos):
                        depths[chrom][i] += 1
                    if alignments:
                        align[chrom].append((cur_pos, delta))
                cur_pos = nxt_pos

            limit[chrom] = max(limit[chrom], cur_pos + 1)

            if get_junctions and rec.cigartuples:
                if rec.has_tag(STRAND_TAG):
                    strand = rec.get_tag(STRAND_TAG)
                else:
                    strand = "-" if rec.is_reverse else "+"
                if rec.has_tag(JCT_CODE_TAG):
                    jct_code = rec.get_tag(JCT_CODE_TAG)
                else:
                    jct_code = ""
                jct_list = _build_splice_junctions(
                    chrom,
                    start_pos,
                    rec.cigartuples,
                    strand,
                    jct_code,
                )
                if not jct_list:
                    continue

                junction_tmp.setdefault(chrom, {})
                junction_tmp[chrom].setdefault(strand, {})
                for new_jct in jct_list:
                    junction_tmp[chrom][strand].setdefault(new_jct.p1, {})
                    junction_tmp[chrom][strand][new_jct.p1].setdefault(new_jct.p2, None)
                    existing = junction_tmp[chrom][strand][new_jct.p1][new_jct.p2]
                    if existing is None:
                        junction_tmp[chrom][strand][new_jct.p1][new_jct.p2] = new_jct
                    else:
                        existing.update(new_jct)

    for chrom in depths:
        max_depth = max(maxpos, limit[chrom]) if maxpos < sys.maxsize else limit[chrom]
        if max_depth < len(depths[chrom]):
            depths[chrom] = depths[chrom][:max_depth]

    junctions = {}
    if get_junctions:
        for chrom in sorted(junction_tmp.keys()):
            junctions[chrom] = []
            for strand in sorted(junction_tmp[chrom].keys()):
                for anchor in sorted(junction_tmp[chrom][strand].keys()):
                    for acceptor in sorted(junction_tmp[chrom][strand][anchor].keys()):
                        jct = junction_tmp[chrom][strand][anchor][acceptor]
                        if jct.count >= minjct and jct.minAnchor() >= minanchor:
                            junctions[chrom].append(jct)

    if alignments:
        return depths, junctions, align
    return depths, junctions


def acceptSAMRecord(s, counter, **args):
    """Converts a string to a SAM record and returns the record if it's valid;
    returns None for unmapped SAM records; throws an exception for illegal records.
    NB: assumes record string has been strip()-ed."""
    if not s:
        return None, None  # blank lines
    if s[0] == "@":
        return None, None  # SAM comments

    try:
        rec = AlignmentRecord(s)
    except ValueError as ve:
        LOGGER.error("invalid_sam_record", line_number=counter, line=s)
        raise ve

    if rec.attrs[RNAME] == NULL_CHROMOSOME:
        return None, None
    if bool(rec.attrs[FLAG] & UNMAPPED_FLAG):
        return None, None
    if not validCigarString(rec.attrs[CIGAR]):
        return None, None

    ungap = "%dM" % len(rec.attrs[SEQ])
    if rec.attrs[CIGAR] in [ungap, EXACT_CIGAR]:
        tokens = [ungap]
    else:
        tokens = cigarTokens(rec.attrs[CIGAR])
        # Convert S,H,I and D tokens into M or N:
        try:
            newCigar = convertCigarString(tokens, len(rec.attrs[SEQ]))
        except ValueError as ve:
            raise ve
        rec.attrs[CIGAR] = newCigar
        tokens = cigarTokens(newCigar)
        newLen = sum([int(s[:-1]) for s in tokens if s[-1] == "M"])
        if newLen < len(rec.attrs[SEQ]):
            rec.attrs[SEQ] = rec.attrs[SEQ][:newLen]
            rec.attrs[QUAL] = rec.attrs[QUAL][:newLen]
        elif newLen > len(rec.attrs[SEQ]):
            delta = newLen - len(rec.attrs[SEQ])
            rec.attrs[SEQ] += rec.attrs[SEQ][0] * delta
            rec.attrs[QUAL] += rec.attrs[QUAL][0] * delta

    return rec, tokens


def bamIterator(path, **args):
    """Return an iterator over alignment records rendered as SAM strings."""
    with _open_alignment_file(path, **args) as bam_stream:
        chr_map = pysamChromosomeMap(bam_stream)
        headers = pysamHeaders(bam_stream)
        for h in headers:
            yield h
        for rec in bam_stream:
            yield pysamReadToString(rec, chr_map)


def cigarTokens(cigar):
    """Takes a CIGAR string and returns a list of cigar tokens.
    For example, '12M304N20M' will return ['12M','304N','20M']"""
    return [cigar[m.start() : m.end()] for m in CIGAR_TOKEN.finditer(cigar)]


def convertCigarString(tokens, requiredSize):
    """Converts a CIGAR string that contains S, H, I and D tokens
    into one that consists only of M and N tokens.  For read-mapping
    purposes, these are the only aspects that really matter.  Assumes
    that the string is a valid CIGAR string."""
    if not tokens:
        raise ValueError('Unrecognized CIGAR string (no tokens): "%s"' % tokens)
    elif len(tokens) == 1:
        if (tokens[0][-1] == "M") or (tokens[0] == EXACT_CIGAR):
            return tokens[0]
        else:
            raise ValueError('Unrecognized CIGAR string (single token): "%s"' % tokens)

    # NB: to get here, len(tokens) > 1
    result = []
    prev = tokens[0]
    if prev[-1] == "M":
        result.append(prev)

    for curr in tokens[1:]:
        symbol = curr[-1]
        try:
            # Tokens must be an integer followed by a single character
            size = int(curr[:-1])
        except ValueError:
            raise ValueError('Unrecognized CIGAR string: "%s"' % curr)

        if symbol == "M":
            # If a match follows a hard/soft clip, just
            # extend the match to include clipped region.
            if prev[-1] in ["H", "S"]:
                size += int(prev[:-1])
                result.append("%dM" % size)
            elif prev[-1] == "M":
                # After an insertion/deletion, the
                # previous token might also be a match
                size += int(prev[:-1])
                result[-1] = "%dM" % size
            else:
                # Ignore an I or N
                result.append("%dM" % size)
        elif symbol == "I":
            # Extra characters not in reference, so skip them
            prev = result[-1]
            continue  # do not update prev
        elif symbol == "D":
            # A deletion must follow a match and
            # the lengths should be added together.
            # assert(prev[-1] == 'M')
            size += int(prev[:-1])
            result[-1] = "%dM" % size
            # Pretend this token was a match
            # instead of a delete
            curr = result[-1]
        elif symbol == "N":
            # Splice junction added as is
            result.append(curr)
        elif symbol in ["H", "S"]:
            # Soft and hard clips should only
            # come at the end, after a match
            # assert(prev[-1] == 'M')
            size += int(prev[:-1])
            result[-1] = "%dM" % size
        prev = curr

    return "".join(result)


# A -- Printable character
# i -- Signed 32-bit integer
# f -- Single-precision floating number
# Z -- Printable string, including space
# H -- Byte array in the Hex format
# B -- Integer or numeric array
def pysamAttributesToStrings(pysamAttributes):
    """Converts a pysam list of attribute values (as duples) into a
    canonical SAM attribute string."""
    result = []
    for name, value in pysamAttributes:
        if type(value) is int:
            new_string = "%s:i:%d" % (name, value)
        elif type(value) is float:
            new_string = "%s:i:%f" % (name, value)
        elif type(value) is list:
            new_string = "%s:B:%s" % (name, value)
        elif type(value) is chr:
            new_string = "%s:A:%c" % (name, value)
        elif type(value) is str:
            new_string = "%s:Z:%s" % (name, value)
        else:
            continue
        result.append(new_string)
    return result


def pysamChromosomeMap(pysamFile):
    """Builds a map of pysam indexes to their corresponding chromosome names."""
    result = {}
    if HEADER_SQ_TAG in pysamFile.header:
        chromList = pysamFile.header[HEADER_SQ_TAG]
        for i in range(len(chromList)):
            result[str(i)] = chromList[i][HEADER_SN_TAG]
    return result


def pysamCigarToString(pysamCigar):
    """Converts a pysam list of CIGAR values (as duples) into a
    canonical CIGAR string.  (The pysam cigarstring attribute does
    not use the canonical form.)"""
    return "".join(["%d%s" % (d[1], PYSAM_CIGAR[d[0]]) for d in pysamCigar])


def pysamHeaders(pysamFile):
    """Returns a list of header strings associated with a pysam file."""
    result = []
    hdict = pysamFile.header
    if HEADER_HD_TAG in hdict:
        valdict = hdict[HEADER_HD_TAG]
        hstring = "@%s" % HEADER_HD_TAG
        for k in [HEADER_VN_TAG, HEADER_SO_TAG]:
            if k in valdict:
                hstring += "\t%s:%s" % (k, valdict[k])
        result.append(hstring + "\n")

    if HEADER_SQ_TAG in hdict:
        seqList = hdict[HEADER_SQ_TAG]
        for seq in seqList:
            # @SQ SN:Chr4 LN:18585056
            result.append(
                "@%s\t%s:%s\t%s:%s\n"
                % (
                    HEADER_SQ_TAG,
                    HEADER_SN_TAG,
                    seq[HEADER_SN_TAG],
                    HEADER_LN_TAG,
                    seq[HEADER_LN_TAG],
                )
            )
    return result


def pysamMergeCigar(pysamCigar):
    """Accepts pysam CIGAR duples that include adjustments (soft/hard
    clipping, inserts, deletes) and returns a set of just match/skip duples
    for corresponding locations in the reference sequence."""
    if len(pysamCigar) == 1:
        return pysamCigar
    i = 0
    # skip over initial ignored operations
    while i < len(pysamCigar) and pysamCigar[i][0] in PYSAM_IGNORE_OPS:
        i += 1
    # starts with an intron; something's wrong
    prev = pysamCigar[i]
    if prev[0] == BAM_CREF_SKIP:
        raise ValueError("CIGAR error: starts with splice junction")

    result = [prev] if prev[0] in PYSAM_SAVE_OPS else []

    for curr in pysamCigar[i + 1 :]:
        oper = curr[0]
        length = curr[1]
        if prev[0] in PYSAM_MERGE_OPS:
            if oper == BAM_CREF_SKIP:
                raise ValueError("CIGAR error: splice junction preceded by clip/insert/delete")
            curr = (oper, length + prev[1])
        elif oper == prev[0]:
            result[-1] = (prev[0], length + prev[1])
            continue
        elif oper in PYSAM_MERGE_OPS and prev[0] == BAM_CMATCH:
            result[-1] = (prev[0], length + prev[1])

        if oper in PYSAM_SAVE_OPS:
            if result and oper == result[-1][0]:
                result[-1] = (oper, length + result[-1][1])
            else:
                result.append(curr)
        prev = curr
    return result


def pysamReadDepths(bamFile, chromosome, gene, **args):
    """Returns a relative start position and an array of read depths
    for the given gene based on input from a BAM file."""
    margin = getAttribute("margin", 0, **args)
    verbose = getAttribute("verbose", False, **args)

    loBound = max(0, gene.minpos - margin)
    upBound = gene.maxpos + margin
    allreads = set([r for r in bamFile.fetch(chromosome, loBound, upBound)])
    nSpliced = 0
    nUngapped = 0
    result = [0] * (upBound - loBound + 1)
    for r in allreads:
        if len(r.cigar) > 1:  # Spliced alignment
            if pysamStrand(r) != gene.strand:
                continue
            nSpliced += 1
            cigar = pysamMergeCigar(r.cigar)
            pos = r.pos + 1
            for tok in cigar:
                if tok[0] == BAM_CMATCH:
                    start = max(pos, loBound) - loBound
                    end = min(upBound, pos + tok[1]) - loBound
                    for i in range(start, end):
                        result[i] += 1
                pos += tok[1]
        else:  # ungapped alignment
            # assert(len(r.cigar) == 1)
            nUngapped += 1
            start = max(r.pos + 1, loBound) - loBound
            end = min(r.pos + r.qlen + 1, upBound) - loBound
            for i in range(start, end + 1):
                result[i] += 1

    if verbose:
        LOGGER.info(
            "pysam_read_depths_loaded",
            ungapped_reads=nUngapped,
            spliced_reads=nSpliced,
            gene_id=gene.id,
        )
    return loBound, result


def pysamReadToString(pysamRecord, chromMap):
    """Converts a pysam record into a normal SAM record.  The chromosome
    map links pysam chromosome indexes to their corresponding names
    (see method pysamChromosomeMap)."""
    parts = str(pysamRecord).split("\t")
    try:
        parts[2] = chromMap[parts[2]]
    except KeyError:
        raise ValueError("Chromosome map is missing pysam index %s:\n%s" % (parts[2], chromMap))

    # pysam records adjust the position by -1
    parts[3] = str(int(parts[3]) + 1)
    parts[5] = pysamCigarToString(pysamRecord.cigar)
    attrs = pysamAttributesToStrings(pysamRecord.tags)
    parts = parts[:-1] + attrs
    return "\t".join(parts)


def pysamSpliceJunctions(pysamRecord, chrMap):
    """Converts a pysam record to a list of splice junction records.
    Returns a list of SpliceJunction objects if the SAM record
    represents a spliced alignment.  Otherwise returns None.
    Note: assumes the CIGAR string has already been validated
    and that there are at least 2 matches in the list."""
    # Example matches: ['16M', '204N', '20M', '355N', '40M']
    cigarDuples = pysamMergeCigar(pysamRecord.cigar)
    sizes = [d[1] for d in cigarDuples]
    pos = [pysamRecord.pos + 1 + sizes[0] - 1]
    for b in sizes[1:-1]:
        nextPos = pos[-1] + b
        pos.append(nextPos)

    # Look for junction type: known='YC:A:K' recombined='YC:A:U' predicted='YC:A:P'
    # (Given by SpliceGrapher, but not other programs.)
    tagDict = dict(pysamRecord.tags)
    jctCode = ""
    try:
        jctCode = tagDict[JCT_CODE_TAG]
    except KeyError:
        pass

    result = []
    exons = [d[1] for d in cigarDuples if d[0] == BAM_CMATCH]
    k = 0
    chrom = chrMap[str(pysamRecord.tid)]
    strand = pysamStrand(pysamRecord, tagDict)
    for i in range(0, len(pos), 2):
        jct = SpliceJunction(
            chrom, pos[i], pos[i + 1] + 1, [exons[k], exons[k + 1]], jctCode, strand
        )
        result.append(jct)
        k += 1
    return result


def pysamStrand(pysamRecord, tagDict=None):
    """Convenience method returns the strand given by the pysam record."""
    if not tagDict:
        tagDict = dict(pysamRecord.tags)
    try:
        return tagDict[STRAND_TAG]
    except KeyError:
        return "-" if pysamRecord.is_reverse else "+"


def getNextSamChromosome(samStream, seed=None, **args):
    """Fetches all records for the next chromosome in the
    stream, along with the first record from the next chromosome."""
    verbose = getAttribute("verbose", False, **args)
    result = []
    if seed:
        result.append(seed)
        targetChrom = seed.split("\t")[2]

    seedLine = None
    targetChrom = None
    indicator = ProgressIndicator(10000000, verbose=verbose)
    for line in samStream:
        indicator.update()
        if line.startswith("@"):
            continue
        s = line.strip()
        chrom = s.split("\t")[2]
        if not targetChrom:
            targetChrom = chrom

        if chrom == targetChrom:
            result.append(s)
        else:
            seedLine = s
            break
    indicator.finish()
    return result, seedLine


def getSamAlignments(samRecords, **args):
    """Reads a SAM file and fills a dictionary with the start
    position and length of every read in the file.  Spliced
    alignments are stored as separate short reads."""
    verbose = getAttribute("verbose", False, **args)
    chromosomes = getAttribute("chromosomes", None, **args)
    chromSet = makeChromosomeSet(chromosomes)
    indicator = ProgressIndicator(1000000, verbose=verbose)
    tracker = ChromosomeTracker()
    result = {}
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith("@"):
            continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec:
            continue

        c = rec.chromosome()
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet):
                break
            if c not in chromSet:
                continue
        if c not in result:
            result[c] = []

        pos = rec.attrs[POS]
        for m in matches:
            try:
                code = m[-1]
                delta = int(m[:-1])
            except Exception as e:
                LOGGER.error(
                    "sam_parse_error",
                    line_number=indicator.ctr,
                    line=line,
                )
                raise e

            # Update reads for matches only
            if code == "M":
                result[c].append((pos, delta))
            pos += delta

    indicator.finish()
    return result


def getSamDepths(samRecords, **args):
    """Reads a SAM file and fills a read-depth array for each
    chromosome found in the file.  Returns a dictionary of
    read-depth lists indexed by chromosome."""
    verbose = getAttribute("verbose", False, **args)
    maxpos = getAttribute("maxpos", sys.maxsize, **args)
    chromosomes = getAttribute("chromosomes", None, **args)
    chromSet = makeChromosomeSet(chromosomes)

    is_native_alignment = (
        isBamFile(samRecords) or isCramFile(samRecords) or samRecords.lower().endswith(".sam")
    )
    if _is_alignment_path(samRecords) and is_native_alignment:
        depths, _ = _collect_pysam_data(samRecords, junctions=False, **args)
        return depths

    if is_depths_file(samRecords):
        depths, jcts = read_depths(samRecords, junctions=False, **args)
        return depths

    depths = {}
    limit = {}
    if verbose and maxpos < sys.maxsize:
        LOGGER.info("loading_sam_records", max_position=maxpos)
    indicator = ProgressIndicator(1000000, verbose=verbose)

    tracker = ChromosomeTracker()
    omitted = 0
    total = 0
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith("@"):
            continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec:
            omitted += 1
            continue

        c = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet):
                break
            if c not in chromSet:
                continue

        total += 1

        if c not in limit:
            limit[c] = 0
        if c not in depths:
            depths[c] = [0] * (maxpos + 1) if maxpos < sys.maxsize else [0]

        prvPos = rec.attrs[POS]
        if prvPos > maxpos:
            break

        for m in matches:
            try:
                code = m[-1]
                delta = int(m[:-1])
            except Exception as e:
                LOGGER.error(
                    "sam_parse_error",
                    line_number=indicator.ctr,
                    line=line,
                )
                raise e
            curPos = prvPos + delta

            # Expand lists as necessary:
            while curPos + 1 > len(depths[c]):
                depths[c] += [0] * len(depths[c])

            # Update depth for matches only
            if code == "M":
                for i in range(prvPos, curPos):
                    depths[c][i] += 1
            prvPos = curPos

        # SAM records should be sorted, but a spliced read may
        # have a lower start position than subsequent reads
        # yet still have a higher end position:
        limit[c] = max(limit[c], curPos + 1)

    indicator.finish()

    for c in depths:
        maxDepth = max(maxpos, limit[c]) if maxpos < sys.maxsize else limit[c]
        if maxDepth < len(depths[c]):
            depths[c] = depths[c][:maxDepth]

    if verbose:
        loaded = total - omitted
        LOGGER.info(
            "sam_depths_loaded",
            loaded_records=loaded,
            omitted_records=omitted,
            total_records=total,
        )

    return depths


def getSamHeaders(samRecords, **args):
    """Reads a SAM file and returns just the header strings as a list."""
    if _is_alignment_path(samRecords):
        with _open_alignment_file(samRecords, **args) as sam_stream:
            return [h.strip() for h in pysamHeaders(sam_stream)]

    result = []
    instream = samInput(samRecords)
    for line in instream:
        if not line.startswith("@"):
            break
        result.append(line.strip())
    if isinstance(instream, io.IOBase):
        instream.close()
    return result


def getSamHeaderInfo(samStream, **args):
    """Parses header records in a SAM file and returns a list of
    possible chromosomes along with the first alignment record."""
    if _is_alignment_path(samStream):
        with _open_alignment_file(samStream, **args) as stream:
            chroms = set([name for name in stream.references])
            return chroms, None

    verbose = getAttribute("verbose", False, **args)
    seedLine = None
    result = set()
    for line in samStream:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                parts = line.strip().split("\t")
                # assert(len(parts) == 3)
                # assert(parts[1].startswith('SN'))
                result.add(parts[1].split(":")[1])
        else:
            seedLine = line
            break
    if verbose:
        LOGGER.info("sam_header_chromosome_count", chromosome_count=len(result))
        LOGGER.info("sam_header_chromosomes", chromosomes=",".join(sorted(result)))
    return result, seedLine


def getSamJunctions(samRecords, **args):
    """Reads a SAM file and builds lists of splice-junctions
    for each chromosome found in the file.  Returns a dictionary of
    junction lists indexed by chromosome."""
    verbose = getAttribute("verbose", False, **args)
    maxpos = getAttribute("maxpos", sys.maxsize, **args)
    minanchor = getAttribute("minanchor", 0, **args)
    minjct = getAttribute("minjct", 1, **args)
    chromosomes = getAttribute("chromosomes", None, **args)

    is_native_alignment = (
        isBamFile(samRecords) or isCramFile(samRecords) or samRecords.lower().endswith(".sam")
    )
    if _is_alignment_path(samRecords) and is_native_alignment:
        _, junctions = _collect_pysam_data(samRecords, **args)
        return junctions

    if is_depths_file(samRecords):
        depths, jcts = read_depths(samRecords, depths=False, **args)
        return jcts

    # Restrict records to the given chromosomes
    chromSet = makeChromosomeSet(chromosomes)

    if verbose:
        LOGGER.info(
            "loading_splice_junctions",
            min_anchor=minanchor,
            min_junction_support=minjct,
        )

    indicator = ProgressIndicator(1000000, verbose=verbose)
    junctions = {}
    tracker = ChromosomeTracker()
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith("@"):
            continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec or len(matches) < 2:
            continue

        if rec.attrs[POS] > maxpos:
            break
        c = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet):
                break
            if c not in chromSet:
                continue

        jctList = recordToSpliceJunction(rec, matches)
        if not jctList:
            continue

        strand = rec.attrs[STRAND_TAG]
        if c not in junctions:
            junctions[c] = {}
        if strand not in junctions[c]:
            junctions[c][strand] = {}
        for newJct in jctList:
            junctions[c][strand].setdefault(newJct.p1, {})
            junctions[c][strand][newJct.p1].setdefault(newJct.p2, None)

            if junctions[c][strand][newJct.p1][newJct.p2] is None:
                junctions[c][strand][newJct.p1][newJct.p2] = newJct
            else:
                junctions[c][strand][newJct.p1][newJct.p2].update(newJct)

    result = {}
    for c in sorted(junctions.keys()):
        result[c] = []
        for s in sorted(junctions[c].keys()):
            for a in sorted(junctions[c][s].keys()):
                for b in sorted(junctions[c][s][a].keys()):
                    jct = junctions[c][s][a][b]
                    if jct.count >= minjct and jct.minAnchor() >= minanchor:
                        result[c].append(jct)

    indicator.finish()
    return result


def getSamReadData(samRecords, **args):
    """Reads a SAM file and fills a read-depth array for each
    chromosome found in the file along with a junction array
    for each chromosome.  By default, returns two dictionaries:
    one with read-depth lists indexed by chromosome and another
    with junction lists indexed by chromosome.  If the 'alignments'
    option is set, returns a third dictionary of (pos,len) tuples
    representing every read alignment."""
    alignments = getAttribute("alignments", False, **args)
    chromosomes = getAttribute("chromosomes", None, **args)
    chromSet = makeChromosomeSet(chromosomes)
    getJct = getAttribute("junctions", True, **args)
    maxpos = getAttribute("maxpos", sys.maxsize, **args)
    minanchor = getAttribute("minanchor", 0, **args)
    minjct = getAttribute("minjct", 1, **args)
    verbose = getAttribute("verbose", False, **args)

    is_native_alignment = (
        isBamFile(samRecords) or isCramFile(samRecords) or samRecords.lower().endswith(".sam")
    )
    if _is_alignment_path(samRecords) and is_native_alignment:
        return _collect_pysam_data(samRecords, **args)

    if is_depths_file(samRecords):
        depths, jcts = read_depths(samRecords, **args)
        return depths, jcts

    align = {}
    depths = {}
    jctTmp = {}
    limit = {}
    if verbose:
        if maxpos < sys.maxsize:
            LOGGER.info("loading_sam_records", max_position=maxpos)
        LOGGER.info(
            "loading_splice_junctions",
            min_anchor=minanchor,
            min_junction_support=minjct,
        )

    indicator = ProgressIndicator(1000000, verbose=verbose)
    tracker = ChromosomeTracker()
    curPos = 0
    for line in samInput(samRecords):
        indicator.update()
        if line.startswith("@"):
            continue
        rec, matches = acceptSAMRecord(line.strip(), indicator.ctr)
        if not rec:
            continue

        c = rec.chromosome()

        # If chromosome changes, check to see if we have loaded all of them
        if chromSet:
            if tracker.update(c) and tracker.allFinished(chromSet):
                break
            if c not in chromSet:
                continue

        if c not in depths:
            depths[c] = [0] * (maxpos + 1) if maxpos < sys.maxsize else [0]

        if c not in limit:
            limit[c] = 0

        if alignments and c not in align:
            align[c] = []

        prvPos = rec.attrs[POS]
        if prvPos > maxpos:
            break

        for m in matches:
            try:
                code = m[-1]
                delta = int(m[:-1])
            except Exception as e:
                LOGGER.error(
                    "invalid_cigar_string",
                    line_number=indicator.ctr,
                    record=line,
                    cigar=rec.cigar(),
                )
                raise e

            curPos = prvPos + delta

            # Expand lists as necessary:
            while curPos + 1 > len(depths[c]):
                depths[c] += [0] * len(depths[c])

            # Update depth for matches only
            if code == "M":
                for i in range(prvPos, curPos):
                    depths[c][i] += 1
                if alignments:
                    align[c].append((prvPos, delta))
            prvPos = curPos

        # SAM records should be sorted, but a spliced read may
        # have a lower start position than subsequent reads
        # yet still have a higher end position:
        limit[c] = max(limit[c], curPos + 1)

        # Now look for junctions
        if getJct and len(matches) > 1:
            jctList = recordToSpliceJunction(rec, matches)
            if not jctList:
                continue

            strand = rec.attrs[STRAND_TAG]
            if c not in jctTmp:
                jctTmp[c] = {}
            if strand not in jctTmp[c]:
                jctTmp[c][strand] = {}
            for newJct in jctList:
                jctTmp[c][strand].setdefault(newJct.p1, {})
                jctTmp[c][strand][newJct.p1].setdefault(newJct.p2, None)

                if jctTmp[c][strand][newJct.p1][newJct.p2] is None:
                    jctTmp[c][strand][newJct.p1][newJct.p2] = newJct
                else:
                    jctTmp[c][strand][newJct.p1][newJct.p2].update(newJct)

    indicator.finish()

    for c in depths:
        maxDepth = max(maxpos, limit[c]) if maxpos < sys.maxsize else limit[c]
        if maxDepth < len(depths[c]):
            depths[c] = depths[c][:maxDepth]

    junctions = {}
    if getJct:
        for c in sorted(jctTmp.keys()):
            junctions[c] = []
            for s in sorted(jctTmp[c].keys()):
                for a in sorted(jctTmp[c][s].keys()):
                    for b in sorted(jctTmp[c][s][a].keys()):
                        jct = jctTmp[c][s][a][b]
                        if jct.count >= minjct and jct.minAnchor() >= minanchor:
                            junctions[c].append(jctTmp[c][s][a][b])

    # Really need a better way to do this
    if alignments:
        return depths, junctions, align
    else:
        return depths, junctions


def getSamSequences(samRecords, **args):
    """Reads a SAM file and returns a dictionary of sequence names and lengths."""
    result = {}
    headers = getSamHeaders(samRecords, **args)
    for line in headers:
        if not line.startswith(HEADER_SQ_LINE):
            continue
        parts = line.strip().split("\t")
        sdict = dict([tuple(s.split(":")) for s in parts[1:]])
        try:
            result[sdict[HEADER_SN_TAG]] = sdict[HEADER_LN_TAG]
        except KeyError:
            continue
    return result


def isBamFile(filePath):
    """Simple heuristic returns True if path is a BAM file; false otherwise."""
    return filePath.lower().endswith(".bam")


def loadSAMRecords(f):
    """Loads SAM records from a file and returns them in a list."""
    result = []
    for s in ez_open(f):
        if not s.startswith("@"):
            result.append(AlignmentRecord(s))
    return result


def makeChromosomeSet(chromList):
    """Converts a string or a list of chromosomes into a set of unique values."""
    chromSet = None
    if chromList:
        if isinstance(chromList, str):
            chromList = [chromList]
        chromSet = set([c.lower() for c in chromList])
    return chromSet


def recordToSpliceJunction(samRec, matches):
    """Converts a SAM record to a list of splice junction records.
    Returns a list of SpliceJunction objects if the SAM record
    represents a spliced alignment.  Otherwise returns None.
    Note: assumes the CIGAR string has already been validated
    and that there are at least 2 matches in the list."""

    # Example matches: ['16M', '204N', '20M', '355N', '40M']
    sizes = [int(m[:-1]) for m in matches]
    pos = [samRec.attrs[POS] + sizes[0] - 1]
    for b in sizes[1:-1]:
        nextPos = pos[-1] + b
        pos.append(nextPos)

    # Look for junction type: known='YC:A:K' recombined='YC:A:U' predicted='YC:A:P'
    # (Given by SpliceGrapher, but not other programs.)
    try:
        jctCode = samRec[JCT_CODE_TAG]
    except KeyError:
        jctCode = ""

    result = []
    exons = [int(m[:-1]) for m in matches if m[-1] == "M"]
    k = 0
    for i in range(0, len(pos), 2):
        jct = SpliceJunction(
            samRec.chromosome(),
            pos[i],
            pos[i + 1] + 1,
            [exons[k], exons[k + 1]],
            jctCode,
            samRec.attrs[STRAND_TAG],
        )
        result.append(jct)
        k += 1
    return result


def samInput(source):
    """Convenience method that returns an iterator over a set of SAM records.
    If source is a list, set or file stream it returns the object.  A string
    is interpreted as a file path to be opened for input and the stream returned."""
    if isinstance(source, (list, tuple, set, io.IOBase)):
        return source
    elif isinstance(source, str):
        return samIterator(source)
    else:
        raise ValueError("Unrecognized type %s for SAM input" % type(source))


def samIterator(path, isBam=False, **args):
    """Return an iterator over SAM/BAM/CRAM records."""
    if isBam or isBamFile(path) or isCramFile(path):
        return bamIterator(path, **args)
    else:
        return ez_open(path)


def validCigarString(s):
    """Return True when the given CIGAR string is usable."""
    if s == NULL_CIGAR:
        return False
    elif s == EXACT_CIGAR:
        return True
    else:
        return CIGAR_MATCH.match(s) is not None


class AlignmentRecord(object):
    """Encapsulates all the information in a SAM record."""

    def __init__(self, s):
        self.attrs = {}
        self.vtypes = {}
        parts = s.split("\t")
        col = 0
        while col in REQUIRED_COLUMNS:
            try:
                if col in INT_COLUMNS:
                    self.attrs[col] = int(parts[col])
                    self.vtypes[col] = "i"
                else:
                    self.attrs[col] = parts[col]
                    self.vtypes[col] = "A"
            except IndexError:
                raise ValueError(
                    "Too few columns (%d<%d) in input record:\n%s"
                    % (len(parts), len(REQUIRED_COLUMNS), s)
                )
            col += 1

        while col < len(parts):
            triplets = parts[col].split()
            for triplet in triplets:
                try:
                    (tag, vtype, val) = triplet.split(":")
                    if vtype not in VALID_VTYPES:
                        raise ValueError()
                    self.vtypes[tag] = vtype
                    self.attrs[tag] = val
                except ValueError:
                    raise ValueError(
                        'Illegal tag:type:value SAM attribute "%s"; '
                        "record contains %d columns" % (triplet, len(parts))
                    )
            col += 1

        # Only do these once:
        if STRAND_TAG not in self.attrs:
            self.attrs[STRAND_TAG] = "-" if bool(self.attrs[FLAG] & REVERSE_FLAG) else "+"
            self.vtypes[STRAND_TAG] = "A"

    def chromosome(self):
        """Convenience method returns the chromosome."""
        return self.attrs[RNAME].lower()

    def cigar(self):
        """Convenience method returns the cigar string."""
        return self.attrs[CIGAR]

    def __cmp__(self, other):
        return self.attrs[POS] - other.attrs[POS]

    def flag(self):
        """Convenience method returns the bitwise flag value as an int."""
        return self.attrs[FLAG]

    def __getitem__(self, key):
        if key not in self.vtypes:
            return self.attrs[key]
        elif self.vtypes[key] in ["A", "Z"]:
            return self.attrs[key]
        elif self.vtypes[key] == "i":
            return int(self.attrs[key])
        elif self.vtypes[key] == "f":
            return float(self.attrs[key])
        elif self.vtypes[key] == "H":
            return int(self.attrs[key], 16)
        else:
            return self.attrs[key]

    # Required for creating dicts/sets
    def __hash__(self):
        return str(self).__hash__()

    def matchpos(self):
        """Convenience method returns the mate position value as an int."""
        return self.attrs[MPOS]

    def pos(self):
        """Convenience method returns the position value as an int."""
        return self.attrs[POS]

    def quality(self):
        """Convenience method returns the phred quality value as an int."""
        return self.attrs[MAPQ]

    def query(self):
        """Convenience method returns the original query sequence."""
        return self.attrs[SEQ]

    def read(self):
        """Convenience method returns the read id."""
        return self.attrs[QNAME]

    def setAttributeString(self, valueString):
        """Sets an attribute in SAM key:vtype:value format.
        Accordingly, the valueString must have the form 'key:vtype:value'."""
        self.setAttribute(valueString.split(":"))

    def setAttribute(self, values):
        """Set an attribute from [key, type, value] SAM key:type:value tokens."""
        if type(values) not in [set, list]:
            raise ValueError("Value must be a set or a list; received %s" % type(values))
        if not (2 <= len(values) <= 3):
            raise ValueError("Value list must have 2 or 3 elements; received %d" % len(values))
        key, vtype, value = values
        if vtype not in VALID_VTYPES:
            raise ValueError("Invalid vtype: %s" % vtype)
        self.vtypes[key] = vtype
        self.attrs[key] = value

    def __str__(self):
        reqd = [str(self.attrs[x]) for x in REQUIRED_COLUMNS]
        tags = [
            "%s:%s:%s" % (x, self.vtypes[x], self.attrs[x])
            for x in self.attrs
            if x not in REQUIRED_COLUMNS
        ]
        vals = reqd + tags
        return "\t".join([x for x in vals if x is not None])

    def strand(self):
        """Convenience method returns the strand given by the bitwise flag."""
        return self.attrs[STRAND_TAG]
