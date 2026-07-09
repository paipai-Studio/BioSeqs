// Companion program that runs the same battery of Seq test cases as
// tests/python_reference.py and emits results as JSON. The two outputs can
// be diffed to verify the MoonBit reimplementation matches BioPython.
//
// Run:  moon run cmd/main

///|
/// One test-case result: a name and its stringified result.
struct Case {
  name : String
  result : String
}

///|
/// Helper to record a case result.
fn record(cases : Array[Case], name : String, result : String) -> Unit {
  cases.push(Case::{ name, result })
}

///|
/// Helper to record a Bool result as "True"/"False" (matching Python's repr).
fn record_bool(cases : Array[Case], name : String, result : Bool) -> Unit {
  record(cases, name, if result { "True" } else { "False" })
}

///|
/// Run every test case and return the list of results, mirroring
/// tests/python_reference.py.
fn run_cases() -> Array[Case] {
  let cases = Array::new()
  let s = @lib.Seq::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

  // --- Basic construction ---
  record(cases, "to_string", s.to_string())
  record(cases, "length", s.length().to_string())

  // --- complement / complement_rna ---
  record(cases, "complement_CGA", @lib.Seq::new("CGA").complement().to_string())
  record(cases, "complement_CGAUT", @lib.Seq::new("CGAUT").complement().to_string())
  record(cases, "complement_rna_CGAUT", @lib.Seq::new("CGAUT").complement_rna().to_string())
  record(cases, "complement_lower", @lib.Seq::new("cgaut").complement().to_string())
  record(cases, "complement_rna_lower", @lib.Seq::new("cgaut").complement_rna().to_string())

  // --- reverse_complement ---
  record(cases, "reverse_complement", s.reverse_complement().to_string())
  record(cases, "reverse_complement_rna", s.reverse_complement_rna().to_string())
  record(cases, "reverse_ACGT", @lib.Seq::new("ACGT").reverse().to_string())

  // --- transcribe / back_transcribe ---
  record(cases, "transcribe", s.transcribe().to_string())
  let rna = @lib.Seq::new("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
  record(cases, "back_transcribe", rna.back_transcribe().to_string())

  // --- translate ---
  record(cases, "translate_default", translate_default(s))
  record(cases, "translate_to_stop", translate_to_stop(s))
  record(cases, "translate_stop_at", translate_stop_at(s, '@'))
  record(cases, "translate_partial", translate_default(@lib.Seq::new("ATGGCCAT")))
  record(cases, "translate_TAN", translate_default(@lib.Seq::new("TAN")))
  record(cases, "translate_NNN", translate_default(@lib.Seq::new("NNN")))

  // Valid CDS
  let cds = @lib.Seq::new("ATGGCCATTGGCCTGAAAGGGTGCCCCCCGTAG")
  record(cases, "translate_cds", translate_cds(cds))

  // Error cases — mirror Python: try each operation, report exception class
  // name on error ("TranslationError"), or "NO_ERROR" if it surprisingly succeeds.
  record(cases, "translate_cds_bad_start", cds_bad_start())
  record(cases, "translate_cds_internal_stop", cds_internal_stop(s))
  record(cases, "translate_invalid_codon", invalid_codon())

  // --- upper / lower ---
  record(cases, "upper", @lib.Seq::new("acgtACGT").upper().to_string())
  record(cases, "lower", @lib.Seq::new("acgtACGT").lower().to_string())

  // --- count / count_overlap ---
  let c = @lib.Seq::new("ATATGAAATTTGAAAA")
  record(cases, "count_AAA", c.count("AAA").to_string())
  record(cases, "count_overlap_AAA", c.count_overlap("AAA").to_string())
  record(cases, "count_ATA_ATATA", @lib.Seq::new("ATATA").count("ATA").to_string())
  record(cases, "count_overlap_ATA_ATATA", @lib.Seq::new("ATATA").count_overlap("ATA").to_string())

  // --- find / rfind ---
  record(cases, "find_GCC", s.find("GCC").to_string())
  record(cases, "rfind_GCC", s.rfind("GCC").to_string())
  record(cases, "find_XYZ", s.find("XYZ").to_string())
  record(cases, "rfind_XYZ", s.rfind("XYZ").to_string())

  // --- startswith / endswith ---
  record_bool(cases, "startswith_ATG", s.startswith("ATG"))
  record_bool(cases, "endswith_TAG", s.endswith("TAG"))
  record_bool(cases, "startswith_GTG", s.startswith("GTG"))

  // --- contains ---
  record_bool(cases, "contains_AAA", @lib.Seq::new("ATATGAAATTTGAAAA").contains("AAA"))
  record_bool(cases, "contains_XYZ", @lib.Seq::new("ATATGAAATTTGAAAA").contains("XYZ"))

  // --- add / mul ---
  record(cases, "add", (@lib.Seq::new("AC") + @lib.Seq::new("GT")).to_string())
  record(cases, "mul3", @lib.Seq::new("AC").mul(3).to_string())

  // --- slice ---
  let sl = @lib.Seq::new("ATGGCCATTGTA")
  record(cases, "slice_1_4", sl.slice(start=1, end=4).to_string())
  record(cases, "slice_0_3", sl.slice(start=0, end=3).to_string())

  // --- strip / lstrip / rstrip ---
  record(cases, "strip_dash", @lib.Seq::new("---ACGT---").strip(chars="-").to_string())
  record(cases, "lstrip_dash", @lib.Seq::new("---ACGT---").lstrip(chars="-").to_string())
  record(cases, "rstrip_dash", @lib.Seq::new("---ACGT---").rstrip(chars="-").to_string())

  // --- replace ---
  record(cases, "replace_TU", @lib.Seq::new("ATGGCC").replace(old="T", new="U").to_string())

  // --- split ---
  let parts = @lib.Seq::new("ACGTACGTACGT").split("AC")
  record(cases, "split_count", parts.length().to_string())
  record(cases, "split_1", parts[1].to_string())

  // --- join ---
  let spacer = @lib.Seq::new("-")
  let join_parts = [
    @lib.Seq::new("AC"), @lib.Seq::new("GT"), @lib.Seq::new("TT"),
  ]
  record(cases, "join", spacer.join(join_parts).to_string())

  // --- equality ---
  record_bool(cases, "eq_same", @lib.Seq::new("ACGT") == @lib.Seq::new("ACGT"))
  record_bool(cases, "eq_diff", @lib.Seq::new("ACGT") == @lib.Seq::new("TTTT"))

  cases
}

///|
fn translate_default(s : @lib.Seq) -> String {
  try {
    s.translate().to_string()
  } catch {
    _ => "ERROR"
  }
}

///|
fn translate_to_stop(s : @lib.Seq) -> String {
  try {
    s.translate(to_stop=true).to_string()
  } catch {
    _ => "ERROR"
  }
}

///|
fn translate_stop_at(s : @lib.Seq, sym : Char) -> String {
  try {
    s.translate(stop_symbol=sym).to_string()
  } catch {
    _ => "ERROR"
  }
}

///|
fn translate_cds(cds : @lib.Seq) -> String {
  try {
    cds.translate(cds=true).to_string()
  } catch {
    _ => "ERROR"
  }
}

///|
/// CDS with a bad start codon (AAA instead of ATG). BioPython raises
/// TranslationError; MoonBit raises SeqError, which we report as
/// "TranslationError" so the JSON output matches.
fn cds_bad_start() -> String {
  try {
    ignore(@lib.Seq::new("AAACCCTAG").translate(cds=true))
    "NO_ERROR"
  } catch {
    _ => "TranslationError"
  }
}

///|
/// CDS with an internal stop codon. The sequence `s` contains TAA/TGA/TAG
/// in the middle, so `translate(cds=true)` must reject it.
fn cds_internal_stop(s : @lib.Seq) -> String {
  try {
    ignore(s.translate(cds=true))
    "NO_ERROR"
  } catch {
    _ => "TranslationError"
  }
}

///|
/// translate() with an invalid codon containing '?' — should raise.
fn invalid_codon() -> String {
  try {
    ignore(@lib.Seq::new("TA?").translate())
    "NO_ERROR"
  } catch {
    _ => "TranslationError"
  }
}

///|
fn main {
  let cases = run_cases()
  // Emit JSON: [{"name":"...","result":"..."}, ...]
  let buf = StringBuilder::new()
  buf.write_char('[')
  for i, c in cases {
    if i > 0 {
      buf.write_char(',')
    }
    buf.write_string("{\"name\":\"")
    buf.write_string(c.name)
    buf.write_string("\",\"result\":\"")
    buf.write_string(c.result)
    buf.write_string("\"}")
  }
  buf.write_char(']')
  println(buf.to_string())
}
