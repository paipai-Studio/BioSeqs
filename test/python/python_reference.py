#!/usr/bin/env python3
"""Generate expected Seq results from BioPython as JSON.

Each test case is run through BioPython's Bio.Seq.Seq and the result is
emitted as a JSON object. The MoonBit companion program (cmd/main/main.mbt)
produces the same JSON structure so the two can be diffed.
"""
import json
import sys
import warnings

from Bio.Seq import Seq

# Suppress BioPython's partial-codon / dual-coding warnings so the output
# stays clean (we only compare return values, not warnings).
warnings.filterwarnings("ignore")


def run_cases():
    """Return a list of (name, result) pairs for every test case."""
    cases = []

    def record(name, result):
        cases.append({"name": name, "result": str(result)})

    # --- Basic construction ---
    s = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    record("to_string", s)
    record("length", len(s))

    # --- complement / complement_rna ---
    record("complement_CGA", Seq("CGA").complement())
    record("complement_CGAUT", Seq("CGAUT").complement())
    record("complement_rna_CGAUT", Seq("CGAUT").complement_rna())
    record("complement_lower", Seq("cgaut").complement())
    record("complement_rna_lower", Seq("cgaut").complement_rna())

    # --- reverse_complement ---
    record("reverse_complement", s.reverse_complement())
    record("reverse_complement_rna", s.reverse_complement_rna())
    # Note: BioPython removed Seq.reverse(); use slicing [::-1] instead.
    record("reverse_ACGT", Seq("ACGT")[::-1])

    # --- transcribe / back_transcribe ---
    record("transcribe", s.transcribe())
    rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
    record("back_transcribe", rna.back_transcribe())

    # --- translate ---
    record("translate_default", s.translate())
    record("translate_to_stop", s.translate(to_stop=True))
    record("translate_stop_at", s.translate(stop_symbol="@"))
    record("translate_partial", Seq("ATGGCCAT").translate())
    record("translate_TAN", Seq("TAN").translate())
    record("translate_NNN", Seq("NNN").translate())

    # Valid CDS
    cds = Seq("ATGGCCATTGGCCTGAAAGGGTGCCCCCCGTAG")
    record("translate_cds", cds.translate(cds=True))

    # translate cds error cases (expect error string)
    try:
        Seq("AAACCCTAG").translate(cds=True)
        record("translate_cds_bad_start", "NO_ERROR")
    except Exception as e:
        record("translate_cds_bad_start", type(e).__name__)

    try:
        s.translate(cds=True)
        record("translate_cds_internal_stop", "NO_ERROR")
    except Exception as e:
        record("translate_cds_internal_stop", type(e).__name__)

    try:
        Seq("TA?").translate()
        record("translate_invalid_codon", "NO_ERROR")
    except Exception as e:
        record("translate_invalid_codon", type(e).__name__)

    # --- upper / lower ---
    record("upper", Seq("acgtACGT").upper())
    record("lower", Seq("acgtACGT").lower())

    # --- count / count_overlap ---
    c = Seq("ATATGAAATTTGAAAA")
    record("count_AAA", c.count("AAA"))
    record("count_overlap_AAA", c.count_overlap("AAA"))
    record("count_ATA_ATATA", Seq("ATATA").count("ATA"))
    record("count_overlap_ATA_ATATA", Seq("ATATA").count_overlap("ATA"))

    # --- find / rfind ---
    record("find_GCC", s.find("GCC"))
    record("rfind_GCC", s.rfind("GCC"))
    record("find_XYZ", s.find("XYZ"))
    record("rfind_XYZ", s.rfind("XYZ"))

    # --- startswith / endswith ---
    record("startswith_ATG", s.startswith("ATG"))
    record("endswith_TAG", s.endswith("TAG"))
    record("startswith_GTG", s.startswith("GTG"))

    # --- contains ---
    record("contains_AAA", "AAA" in Seq("ATATGAAATTTGAAAA"))
    record("contains_XYZ", "XYZ" in Seq("ATATGAAATTTGAAAA"))

    # --- add / mul ---
    record("add", Seq("AC") + Seq("GT"))
    record("mul3", Seq("AC") * 3)

    # --- slice ---
    sl = Seq("ATGGCCATTGTA")
    record("slice_1_4", sl[1:4])
    record("slice_0_3", sl[0:3])

    # --- strip / lstrip / rstrip ---
    record("strip_dash", Seq("---ACGT---").strip("-"))
    record("lstrip_dash", Seq("---ACGT---").lstrip("-"))
    record("rstrip_dash", Seq("---ACGT---").rstrip("-"))

    # --- replace ---
    record("replace_TU", Seq("ATGGCC").replace("T", "U"))

    # --- split ---
    parts = Seq("ACGTACGTACGT").split("AC")
    record("split_count", len(parts))
    record("split_1", parts[1])

    # --- join ---
    record("join", Seq("-").join([Seq("AC"), Seq("GT"), Seq("TT")]))

    # --- equality ---
    record("eq_same", Seq("ACGT") == Seq("ACGT"))
    record("eq_diff", Seq("ACGT") == Seq("TTTT"))

    return cases


def main():
    cases = run_cases()
    json.dump(cases, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
