#!/usr/bin/env python3
"""Generate expected Bio.SeqIO results from BioPython as JSON.

Each test case is run through BioPython's Bio.SeqIO and the result is
emitted as a JSON value between case markers. The MoonBit companion
implementation can produce the same JSON structure so the two can be
diffed.

Output format:
    === CASE: test_name ===
    <json>
    === END ===
"""
import json
import sys
import warnings
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Suppress BioPython warnings so output stays clean.
warnings.filterwarnings("ignore")


def record_to_dict(rec, include_qual=False):
    """Convert a SeqRecord to a JSON-serializable dict."""
    d = {
        "id": rec.id,
        "name": rec.name,
        "description": rec.description,
        "seq": str(rec.seq),
        "length": len(rec.seq),
    }
    if include_qual:
        d["qual"] = list(rec.letter_annotations.get("phred_quality", []))
    return d


def emit(case_name, payload):
    """Print a labeled test case as JSON between markers."""
    print(f"=== CASE: {case_name} ===")
    print(json.dumps(payload))
    print("=== END ===")


# ---------------------------------------------------------------------------
# Test FASTA content
# ---------------------------------------------------------------------------

FASTA_MULTI = """>seq1 Description for seq1
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>seq2 Another description
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
>seq3
ACGTACGTACGTACGT
"""

FASTA_LOWERCASE = """>lower1 lowercase sequence
atggccattgtaatgggccgctgaaagggtgcccgatag
>lower2 MixedCase Sequence
atgGCCattGTAatgGGCcgctGAAagggtgccCGATAG
"""

FASTA_MULTILINE = """>multi1 Multi-line sequence
ATGGCCATTGTA
ATGGGCCGCTGAA
AGGGTGCCCGATAG
>multi2 Another multi-line
ACGT
ACGT
ACGT
ACGT
"""

FASTA_SINGLE = """>only One record only
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
"""

# ---------------------------------------------------------------------------
# Test FASTQ content (Sanger format, Phred+33 offset)
# ---------------------------------------------------------------------------

FASTQ_MULTI = """@seq1 Description for seq1
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@seq2 Another description
ACGTACGT
+
IIIIIIII
@seq3 No quality worth mentioning
atggcc
+
;;;;;;
"""


def case_parse_fasta_multiple():
    records = list(SeqIO.parse(StringIO(FASTA_MULTI), "fasta"))
    emit("parse_fasta_multiple", [record_to_dict(r) for r in records])


def case_parse_fasta_lowercase():
    records = list(SeqIO.parse(StringIO(FASTA_LOWERCASE), "fasta"))
    emit("parse_fasta_lowercase", [record_to_dict(r) for r in records])


def case_parse_fasta_multiline():
    records = list(SeqIO.parse(StringIO(FASTA_MULTILINE), "fasta"))
    emit("parse_fasta_multiline", [record_to_dict(r) for r in records])


def case_read_single_record():
    rec = SeqIO.read(StringIO(FASTA_SINGLE), "fasta")
    emit("read_single_record", record_to_dict(rec))


def case_read_multiple_error():
    """SeqIO.read() on a multi-record file should raise ValueError."""
    try:
        SeqIO.read(StringIO(FASTA_MULTI), "fasta")
        emit("read_multiple_error", {"error": None})
    except Exception as e:
        emit("read_multiple_error",
             {"error": type(e).__name__, "message": str(e)})


def case_to_dict():
    d = SeqIO.to_dict(SeqIO.parse(StringIO(FASTA_MULTI), "fasta"))
    payload = {k: record_to_dict(v) for k, v in d.items()}
    emit("to_dict", payload)


def case_write_fasta():
    recs = [
        SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
                  id="written1", name="written1",
                  description="written1 First written record"),
        SeqRecord(Seq("acgtacgtacgt"),
                  id="written2", name="written2",
                  description="written2 lowercase input"),
        SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
                  id="written3", name="written3",
                  description="written3"),
    ]
    out = StringIO()
    SeqIO.write(recs, out, "fasta")
    emit("write_fasta", out.getvalue())


def case_parse_fastq():
    records = list(SeqIO.parse(StringIO(FASTQ_MULTI), "fastq"))
    emit("parse_fastq",
         [record_to_dict(r, include_qual=True) for r in records])


def case_write_fastq():
    recs = [
        SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
                  id="fq1", name="fq1",
                  description="fq1 written fastq record",
                  letter_annotations={"phred_quality": [40] * 39}),
        SeqRecord(Seq("acgt"),
                  id="fq2", name="fq2",
                  description="fq2 lowercase",
                  letter_annotations={"phred_quality": [10, 20, 30, 40]}),
    ]
    out = StringIO()
    SeqIO.write(recs, out, "fastq")
    emit("write_fastq", out.getvalue())


# ---------------------------------------------------------------------------
# Test GenBank content
# ---------------------------------------------------------------------------

GENBANK_SAMPLE = """LOCUS       NM_001              48 bp    DNA     linear   PRI 01-JAN-2020
DEFINITION  Test gene for comparison.
ACCESSION   NM_001
VERSION     NM_001.1
SOURCE      Test organism
  ORGANISM  Test organism
            Eukaryota; Test.
FEATURES             Location/Qualifiers
     source          1..48
                     /organism="Test organism"
     gene            1..48
                     /gene="test_gene"
     CDS             1..48
                     /translation="MKRHPGSTH"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgc
//
"""


def case_parse_genbank():
    records = list(SeqIO.parse(StringIO(GENBANK_SAMPLE), "genbank"))
    if records:
        emit("parse_genbank", record_to_dict(records[0]))
    else:
        emit("parse_genbank", {})


def main():
    case_parse_fasta_multiple()
    case_parse_fasta_lowercase()
    case_parse_fasta_multiline()
    case_read_single_record()
    case_read_multiple_error()
    case_to_dict()
    case_write_fasta()
    case_parse_fastq()
    case_write_fastq()
    case_parse_genbank()


if __name__ == "__main__":
    main()
