// Blackbox tests for the MoonBit Seq reimplementation, cross-checked against
// BioPython 1.87's Bio.Seq.Seq. Each assertion's expected value was produced
// by running the equivalent BioPython code (see tests/python_reference.py).
// Run `moon test` to execute.

///|
test "seq_construction_and_to_string" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_eq(s.to_string(), "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_eq(s.repr(), "Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')")
  assert_eq(s.length(), 39)
}

///|
test "seq_complement_and_complement_rna" {
  // Seq("CGA").complement() -> Seq('GCT')
  assert_eq(@src.Seq::new("CGA").complement().to_string(), "GCT")
  // Seq("CGAUT").complement() -> Seq('GCTAA')  (U treated as T)
  assert_eq(@src.Seq::new("CGAUT").complement().to_string(), "GCTAA")
  // Seq("CGAUT").complement_rna() -> Seq('GCUAA')  (T treated as U)
  assert_eq(@src.Seq::new("CGAUT").complement_rna().to_string(), "GCUAA")
  // lowercase preserved
  assert_eq(@src.Seq::new("cgaut").complement().to_string(), "gctaa")
  assert_eq(
    @src.Seq::new("cgaut").complement_rna().to_string(),
    "gctaa".replace(old="t", new="u"),
  )
}

///|
test "reverse_complement" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.reverse_complement() -> 'CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT'
  assert_eq(
    s.reverse_complement().to_string(),
    "CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT",
  )
  // complement then reverse == reverse then complement (round-trip check)
  assert_eq(
    s.reverse().complement().to_string(),
    s.complement().reverse().to_string(),
  )
}

///|
test "seq_transcribe_and_back_transcribe" {
  let dna = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // dna.transcribe() -> 'AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG'
  assert_eq(
    dna.transcribe().to_string(),
    "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG",
  )
  let rna = @src.Seq::new("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGTGCCCGAUAG")
  assert_eq(
    rna.back_transcribe().to_string(),
    "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
  )
  // transcribe then back_transcribe is identity
  assert_eq(dna.transcribe().back_transcribe().to_string(), dna.to_string())
}

///|
test "translate default" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate() -> 'MAIVMGR*KGAR*'
  assert_eq(s.translate().to_string(), "MAIVMGR*KGAR*")
}

///|
test "seq_translate_to_stop" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate(to_stop=True) -> 'MAIVMGR'
  assert_eq(s.translate(to_stop=true).to_string(), "MAIVMGR")
}

///|
test "seq_translate_cds_valid" {
  // A valid CDS: ATG ... TAG with no in-frame internal stop codons.
  // Codons: ATG GCC ATT GGC CTG AAA GGG TGC CCC CCG TAG (33 bp)
  // Translation (cds=True): M A I G L K G C P P
  let cds = @src.Seq::new("ATGGCCATTGGCCTGAAAGGGTGCCCCCCGTAG")
  assert_eq(cds.translate(cds=true).to_string(), "MAIGLKGCPP")
}

///|
test "seq_translate_cds_invalid_start" {
  // First codon is AAA (not a start codon).
  let bad = @src.Seq::new("AAACCCTAG")
  try {
    ignore(bad.translate(cds=true))
    fail("expected SeqError for invalid start codon")
  } catch {
    _ => ()
  }
}

///|
test "seq_translate_cds_internal_stop" {
  // Has an internal TGA stop codon.
  let bad = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  try {
    ignore(bad.translate(cds=true))
    fail("expected SeqError for internal stop codon")
  } catch {
    _ => ()
  }
}

///|
test "seq_translate_partial_codon" {
  // Length not a multiple of 3; trailing bases are dropped (BioPython warns).
  // BioPython: Seq("ATGGCCAT").translate() -> 'MA'
  assert_eq(@src.Seq::new("ATGGCCAT").translate().to_string(), "MA")
}

///|
test "seq_translate_ambiguous_codon" {
  // TAN -> X (possible stop, all valid letters)
  assert_eq(@src.Seq::new("TAN").translate().to_string(), "X")
  // NNN -> X
  assert_eq(@src.Seq::new("NNN").translate().to_string(), "X")
}

///|
test "seq_translate_invalid_codon_raises" {
  // TA? contains ? which is not a valid nucleotide letter.
  try {
    ignore(@src.Seq::new("TA?").translate())
    fail("expected SeqError for invalid codon")
  } catch {
    _ => ()
  }
}

///|
test "seq_translate_stop_symbol" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate(stop_symbol="@") -> 'MAIVMGR@KGAR@'
  assert_eq(s.translate(stop_symbol='@').to_string(), "MAIVMGR@KGAR@")
}

///|
test "upper and lower" {
  assert_eq(@src.Seq::new("acgtACGT").upper().to_string(), "ACGTACGT")
  assert_eq(@src.Seq::new("acgtACGT").lower().to_string(), "acgtacgt")
}

///|
test "count" {
  let s = @src.Seq::new("ATATGAAATTTGAAAA")
  // BioPython: s.count("AAA") -> 2 (non-overlapping, at indices 5 and 12)
  assert_eq(s.count("AAA"), 2)
  // count_overlap: "AAA" at indices 5, 12, 13 -> 3 overlapping
  assert_eq(s.count_overlap("AAA"), 3)
  // overlapping pattern: "ATA" in "ATATA" -> 2 overlapping
  assert_eq(@src.Seq::new("ATATA").count_overlap("ATA"), 2)
  assert_eq(@src.Seq::new("ATATA").count("ATA"), 1)
}

///|
test "seq_find_and_rfind" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.find("GCC") -> 3, s.rfind("GCC") -> 30
  assert_eq(s.find("GCC"), 3)
  assert_eq(s.rfind("GCC"), 30)
  assert_eq(s.find("XYZ"), -1)
  assert_eq(s.rfind("XYZ"), -1)
}

///|
test "seq_startswith_and_endswith" {
  let s = @src.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_true(s.startswith("ATG"))
  assert_true(s.endswith("TAG"))
  assert_false(s.startswith("GTG"))
  assert_false(s.endswith("ATG"))
}

///|
test "contains" {
  let s = @src.Seq::new("ATATGAAATTTGAAAA")
  assert_true(s.contains("AAA"))
  assert_false(s.contains("XYZ"))
}

///|
test "seq_add_and_mul" {
  assert_eq((@src.Seq::new("AC") + @src.Seq::new("GT")).to_string(), "ACGT")
  assert_eq(@src.Seq::new("AC").add_str("GT").to_string(), "ACGT")
  assert_eq(@src.Seq::new("AC").mul(3).to_string(), "ACACAC")
  assert_eq(@src.Seq::new("AC").mul(0).to_string(), "")
}

///|
test "slice" {
  let s = @src.Seq::new("ATGGCCATTGTA")
  // Python s[1:4] -> 'TGG'
  assert_eq(s.slice(start=1, end=4).to_string(), "TGG")
  // Python s[0:3] -> 'ATG'
  assert_eq(s.slice(start=0, end=3).to_string(), "ATG")
}

///|
test "get_char" {
  let s = @src.Seq::new("ACGT")
  assert_eq(s.get_char(0), Some('A'))
  assert_eq(s.get_char(3), Some('T'))
  match s.get_char(4) {
    None => ()
    Some(_) => fail("expected None for out-of-bounds index")
  }
}

///|
test "seq_strip_lstrip_rstrip" {
  assert_eq(@src.Seq::new("---ACGT---").strip(chars="-").to_string(), "ACGT")
  assert_eq(
    @src.Seq::new("---ACGT---").lstrip(chars="-").to_string(),
    "ACGT---",
  )
  assert_eq(
    @src.Seq::new("---ACGT---").rstrip(chars="-").to_string(),
    "---ACGT",
  )
  // default strips whitespace
  assert_eq(@src.Seq::new("  ACGT  ").strip().to_string(), "ACGT")
}

///|
test "replace" {
  assert_eq(
    @src.Seq::new("ATGGCC").replace(old="T", new="U").to_string(),
    "AUGGCC",
  )
}

///|
test "split" {
  let parts = @src.Seq::new("ACGTACGTACGT").split("AC")
  // "ACGTACGTACGT".split("AC") -> ['', 'GT', 'GT', 'GT']
  assert_eq(parts.length(), 4)
  assert_eq(parts[0].to_string(), "")
  assert_eq(parts[1].to_string(), "GT")
  assert_eq(parts[2].to_string(), "GT")
  assert_eq(parts[3].to_string(), "GT")
}

///|
test "join" {
  let spacer = @src.Seq::new("-")
  let parts = [@src.Seq::new("AC"), @src.Seq::new("GT"), @src.Seq::new("TT")]
  // "-".join(["AC", "GT", "TT"]) -> 'AC-GT-TT'
  assert_eq(spacer.join(parts).to_string(), "AC-GT-TT")
}

///|
test "equality" {
  let a = @src.Seq::new("ACGT")
  let b = @src.Seq::new("ACGT")
  let c = @src.Seq::new("TTTT")
  assert_true(a == b)
  assert_false(a == c)
}

///|
test "reverse" {
  assert_eq(@src.Seq::new("ACGT").reverse().to_string(), "TGCA")
}

///|
test "seq_complement_rna_round_trip" {
  let dna = @src.Seq::new("ATCG")
  // reverse_complement of DNA, then back should equal complement of reverse
  let rc = dna.reverse_complement()
  assert_eq(rc.to_string(), "CGAT")
}
