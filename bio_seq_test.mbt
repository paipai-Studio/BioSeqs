// Blackbox tests for the MoonBit Seq reimplementation, cross-checked against
// BioPython 1.87's Bio.Seq.Seq. Each assertion's expected value was produced
// by running the equivalent BioPython code (see tests/python_reference.py).
// Run `moon test` to execute.

///|
test "construction and to_string" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_eq(s.to_string(), "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_eq(s.repr(), "Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')")
  assert_eq(s.length(), 39)
}

///|
test "complement and complement_rna" {
  // Seq("CGA").complement() -> Seq('GCT')
  assert_eq(@bio_seq.Seq::new("CGA").complement().to_string(), "GCT")
  // Seq("CGAUT").complement() -> Seq('GCTAA')  (U treated as T)
  assert_eq(@bio_seq.Seq::new("CGAUT").complement().to_string(), "GCTAA")
  // Seq("CGAUT").complement_rna() -> Seq('GCUAA')  (T treated as U)
  assert_eq(@bio_seq.Seq::new("CGAUT").complement_rna().to_string(), "GCUAA")
  // lowercase preserved
  assert_eq(@bio_seq.Seq::new("cgaut").complement().to_string(), "gctaa")
  assert_eq(@bio_seq.Seq::new("cgaut").complement_rna().to_string(), "gctaa".replace(old="t", new="u"))
}

///|
test "reverse_complement" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
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
test "transcribe and back_transcribe" {
  let dna = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // dna.transcribe() -> 'AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG'
  assert_eq(
    dna.transcribe().to_string(),
    "AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG",
  )
  let rna = @bio_seq.Seq::new("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGTGCCCGAUAG")
  assert_eq(
    rna.back_transcribe().to_string(),
    "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
  )
  // transcribe then back_transcribe is identity
  assert_eq(dna.transcribe().back_transcribe().to_string(), dna.to_string())
}

///|
test "translate default" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate() -> 'MAIVMGR*KGAR*'
  assert_eq(s.translate().to_string(), "MAIVMGR*KGAR*")
}

///|
test "translate to_stop" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate(to_stop=True) -> 'MAIVMGR'
  assert_eq(s.translate(to_stop=true).to_string(), "MAIVMGR")
}

///|
test "translate cds valid" {
  // A valid CDS: ATG ... TAG with no in-frame internal stop codons.
  // Codons: ATG GCC ATT GGC CTG AAA GGG TGC CCC CCG TAG (33 bp)
  // Translation (cds=True): M A I G L K G C P P
  let cds = @bio_seq.Seq::new("ATGGCCATTGGCCTGAAAGGGTGCCCCCCGTAG")
  assert_eq(cds.translate(cds=true).to_string(), "MAIGLKGCPP")
}

///|
test "translate cds invalid start" {
  // First codon is AAA (not a start codon).
  let bad = @bio_seq.Seq::new("AAACCCTAG")
  try {
    ignore(bad.translate(cds=true))
    fail("expected SeqError for invalid start codon")
  } catch {
    _ => ()
  }
}

///|
test "translate cds internal stop" {
  // Has an internal TGA stop codon.
  let bad = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  try {
    ignore(bad.translate(cds=true))
    fail("expected SeqError for internal stop codon")
  } catch {
    _ => ()
  }
}

///|
test "translate partial codon" {
  // Length not a multiple of 3; trailing bases are dropped (BioPython warns).
  // BioPython: Seq("ATGGCCAT").translate() -> 'MA'
  assert_eq(@bio_seq.Seq::new("ATGGCCAT").translate().to_string(), "MA")
}

///|
test "translate ambiguous codon" {
  // TAN -> X (possible stop, all valid letters)
  assert_eq(@bio_seq.Seq::new("TAN").translate().to_string(), "X")
  // NNN -> X
  assert_eq(@bio_seq.Seq::new("NNN").translate().to_string(), "X")
}

///|
test "translate invalid codon raises" {
  // TA? contains ? which is not a valid nucleotide letter.
  try {
    ignore(@bio_seq.Seq::new("TA?").translate())
    fail("expected SeqError for invalid codon")
  } catch {
    _ => ()
  }
}

///|
test "translate stop_symbol" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.translate(stop_symbol="@") -> 'MAIVMGR@KGAR@'
  assert_eq(s.translate(stop_symbol='@').to_string(), "MAIVMGR@KGAR@")
}

///|
test "upper and lower" {
  assert_eq(@bio_seq.Seq::new("acgtACGT").upper().to_string(), "ACGTACGT")
  assert_eq(@bio_seq.Seq::new("acgtACGT").lower().to_string(), "acgtacgt")
}

///|
test "count" {
  let s = @bio_seq.Seq::new("ATATGAAATTTGAAAA")
  // BioPython: s.count("AAA") -> 2 (non-overlapping, at indices 5 and 12)
  assert_eq(s.count("AAA"), 2)
  // count_overlap: "AAA" at indices 5, 12, 13 -> 3 overlapping
  assert_eq(s.count_overlap("AAA"), 3)
  // overlapping pattern: "ATA" in "ATATA" -> 2 overlapping
  assert_eq(@bio_seq.Seq::new("ATATA").count_overlap("ATA"), 2)
  assert_eq(@bio_seq.Seq::new("ATATA").count("ATA"), 1)
}

///|
test "find and rfind" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  // BioPython: s.find("GCC") -> 3, s.rfind("GCC") -> 30
  assert_eq(s.find("GCC"), 3)
  assert_eq(s.rfind("GCC"), 30)
  assert_eq(s.find("XYZ"), -1)
  assert_eq(s.rfind("XYZ"), -1)
}

///|
test "startswith and endswith" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  assert_true(s.startswith("ATG"))
  assert_true(s.endswith("TAG"))
  assert_false(s.startswith("GTG"))
  assert_false(s.endswith("ATG"))
}

///|
test "contains" {
  let s = @bio_seq.Seq::new("ATATGAAATTTGAAAA")
  assert_true(s.contains("AAA"))
  assert_false(s.contains("XYZ"))
}

///|
test "add and mul" {
  assert_eq((@bio_seq.Seq::new("AC") + @bio_seq.Seq::new("GT")).to_string(), "ACGT")
  assert_eq(@bio_seq.Seq::new("AC").add_str("GT").to_string(), "ACGT")
  assert_eq(@bio_seq.Seq::new("AC").mul(3).to_string(), "ACACAC")
  assert_eq(@bio_seq.Seq::new("AC").mul(0).to_string(), "")
}

///|
test "slice" {
  let s = @bio_seq.Seq::new("ATGGCCATTGTA")
  // Python s[1:4] -> 'TGG'
  assert_eq(s.slice(start=1, end=4).to_string(), "TGG")
  // Python s[0:3] -> 'ATG'
  assert_eq(s.slice(start=0, end=3).to_string(), "ATG")
}

///|
test "get_char" {
  let s = @bio_seq.Seq::new("ACGT")
  assert_eq(s.get_char(0), Some('A'))
  assert_eq(s.get_char(3), Some('T'))
  match s.get_char(4) {
    None => ()
    Some(_) => fail("expected None for out-of-bounds index")
  }
}

///|
test "strip lstrip rstrip" {
  assert_eq(@bio_seq.Seq::new("---ACGT---").strip(chars="-").to_string(), "ACGT")
  assert_eq(@bio_seq.Seq::new("---ACGT---").lstrip(chars="-").to_string(), "ACGT---")
  assert_eq(@bio_seq.Seq::new("---ACGT---").rstrip(chars="-").to_string(), "---ACGT")
  // default strips whitespace
  assert_eq(@bio_seq.Seq::new("  ACGT  ").strip().to_string(), "ACGT")
}

///|
test "replace" {
  assert_eq(
    @bio_seq.Seq::new("ATGGCC").replace(old="T", new="U").to_string(),
    "AUGGCC",
  )
}

///|
test "split" {
  let parts = @bio_seq.Seq::new("ACGTACGTACGT").split("AC")
  // "ACGTACGTACGT".split("AC") -> ['', 'GT', 'GT', 'GT']
  assert_eq(parts.length(), 4)
  assert_eq(parts[0].to_string(), "")
  assert_eq(parts[1].to_string(), "GT")
  assert_eq(parts[2].to_string(), "GT")
  assert_eq(parts[3].to_string(), "GT")
}

///|
test "join" {
  let spacer = @bio_seq.Seq::new("-")
  let parts = [
    @bio_seq.Seq::new("AC"), @bio_seq.Seq::new("GT"), @bio_seq.Seq::new("TT"),
  ]
  // "-".join(["AC", "GT", "TT"]) -> 'AC-GT-TT'
  assert_eq(spacer.join(parts).to_string(), "AC-GT-TT")
}

///|
test "equality" {
  let a = @bio_seq.Seq::new("ACGT")
  let b = @bio_seq.Seq::new("ACGT")
  let c = @bio_seq.Seq::new("TTTT")
  assert_true(a == b)
  assert_false(a == c)
}

///|
test "reverse" {
  assert_eq(@bio_seq.Seq::new("ACGT").reverse().to_string(), "TGCA")
}

///|
test "complement_rna round trip" {
  let dna = @bio_seq.Seq::new("ATCG")
  // reverse_complement of DNA, then back should equal complement of reverse
  let rc = dna.reverse_complement()
  assert_eq(rc.to_string(), "CGAT")
}
