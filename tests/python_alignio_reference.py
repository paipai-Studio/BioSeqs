#!/usr/bin/env python3
"""Generate expected Bio.AlignIO results from BioPython as JSON.

Each test case is run through BioPython's Bio.AlignIO and the result is
emitted as a JSON value between case markers. The MoonBit companion
implementation can produce the same JSON structure so the two can be
diffed.

Output format:
    === CASE: test_name ===
    <json>
    === END ===
"""
import json
import warnings
from io import StringIO

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

warnings.filterwarnings("ignore")


def alignment_to_dict(align):
    """Convert a MultipleSeqAlignment to a JSON-serializable dict."""
    return {
        "num_records": len(align),
        "alignment_length": align.get_alignment_length(),
        "records": [
            {
                "id": rec.id,
                "name": rec.name,
                "description": rec.description,
                "seq": str(rec.seq),
                "length": len(rec.seq),
            }
            for rec in align
        ],
    }


def emit(case_name, payload):
    """Print a labeled test case as JSON between markers."""
    print(f"=== CASE: {case_name} ===")
    print(json.dumps(payload))
    print("=== END ===")


# ---------------------------------------------------------------------------
# Test ClustalW content
# ---------------------------------------------------------------------------

CLUSTAL_SAMPLE = """CLUSTAL X (1.81) multiple sequence alignment

AF191659.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
AF191658.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
AF191661.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
AF191660.1      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG

AF191659.1      AGA
AF191658.1      AGA
AF191661.1      AGA
AF191660.1      AGA
"""

# ---------------------------------------------------------------------------
# Test PHYLIP interlaced content
# ---------------------------------------------------------------------------

PHYLIP_SAMPLE = """3 20
seq1        ATGCCGATCGATCGATCGAT
seq2        ATGCCGATCGATCGATCGAT
seq3        ATGCCGATCGATCGATCGAT
"""

# ---------------------------------------------------------------------------
# Test PHYLIP sequential content
# ---------------------------------------------------------------------------

PHYLIP_SEQ_SAMPLE = """3 20
seq1        ATGCCGATCGATCGATCGAT
seq2        ATGCCGATCGATCGATCGAT
seq3        ATGCCGATCGATCGATCGAT
"""

# ---------------------------------------------------------------------------
# Test FASTA alignment content
# ---------------------------------------------------------------------------

FASTA_ALIGN_SAMPLE = """>seq1 Description 1
ATGCCGATCGATCGATCGATCGATCGATCGATCGATC
>seq2 Description 2
ATGCCGATCGATCGATCGATCGATCGATCGATCGATC
>seq3 Description 3
ATGCCGATCGATCGATCGATCGATCGATCGATCGATC
"""


def case_parse_clustal():
    align = AlignIO.read(StringIO(CLUSTAL_SAMPLE), "clustal")
    emit("parse_clustal", alignment_to_dict(align))


def case_parse_phylip():
    align = AlignIO.read(StringIO(PHYLIP_SAMPLE), "phylip")
    emit("parse_phylip", alignment_to_dict(align))


def case_parse_phylip_sequential():
    align = AlignIO.read(StringIO(PHYLIP_SEQ_SAMPLE), "phylip-sequential")
    emit("parse_phylip_sequential", alignment_to_dict(align))


def case_parse_fasta_align():
    align = AlignIO.read(StringIO(FASTA_ALIGN_SAMPLE), "fasta")
    emit("parse_fasta_align", alignment_to_dict(align))


def case_write_clustal():
    recs = [
        SeqRecord(Seq("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
                  id="AF191659.1", name="AF191659.1", description="AF191659.1"),
        SeqRecord(Seq("TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
                  id="AF191658.1", name="AF191658.1", description="AF191658.1"),
        SeqRecord(Seq("TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAGA"),
                  id="AF191660.1", name="AF191660.1", description="AF191660.1"),
    ]
    align = MultipleSeqAlignment(recs)
    out = StringIO()
    AlignIO.write([align], out, "clustal")
    emit("write_clustal", out.getvalue())


def case_write_phylip():
    recs = [
        SeqRecord(Seq("ATGCCGATCGATCGATCGAT"),
                  id="seq1", name="seq1", description="seq1"),
        SeqRecord(Seq("ATGCCGATCGATCGATCGAT"),
                  id="seq2", name="seq2", description="seq2"),
        SeqRecord(Seq("ATGCCGATCGATCGATCGAT"),
                  id="seq3", name="seq3", description="seq3"),
    ]
    align = MultipleSeqAlignment(recs)
    out = StringIO()
    AlignIO.write([align], out, "phylip")
    emit("write_phylip", out.getvalue())


def case_write_fasta_align():
    recs = [
        SeqRecord(Seq("ATGCCGATCGATCGATCGATCGATCGATCGATCGATC"),
                  id="seq1", name="seq1", description="seq1 desc"),
        SeqRecord(Seq("ATGCCGATCGATCGATCGATCGATCGATCGATCGATC"),
                  id="seq2", name="seq2", description="seq2 desc"),
    ]
    align = MultipleSeqAlignment(recs)
    out = StringIO()
    AlignIO.write([align], out, "fasta")
    emit("write_fasta_align", out.getvalue())


def case_alignment_length():
    align = AlignIO.read(StringIO(CLUSTAL_SAMPLE), "clustal")
    emit("alignment_length", {
        "num_records": len(align),
        "alignment_length": align.get_alignment_length(),
    })


def case_alignment_column():
    align = AlignIO.read(StringIO(CLUSTAL_SAMPLE), "clustal")
    column_0 = "".join([str(rec.seq[0]) for rec in align])
    column_5 = "".join([str(rec.seq[5]) for rec in align])
    emit("alignment_column", {
        "column_0": column_0,
        "column_5": column_5,
    })


def case_read_multiple_error():
    """AlignIO.read() on multiple alignments should raise ValueError."""
    try:
        AlignIO.read(StringIO(CLUSTAL_SAMPLE + "\n" + CLUSTAL_SAMPLE), "clustal")
        emit("read_multiple_error", {"error": None})
    except Exception as e:
        emit("read_multiple_error",
             {"error": type(e).__name__, "message": str(e)})


def main():
    case_parse_clustal()
    case_parse_phylip()
    case_parse_phylip_sequential()
    case_parse_fasta_align()
    case_write_clustal()
    case_write_phylip()
    case_write_fasta_align()
    case_alignment_length()
    case_alignment_column()
    case_read_multiple_error()


if __name__ == "__main__":
    main()
