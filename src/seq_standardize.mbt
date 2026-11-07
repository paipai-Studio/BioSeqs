/// seq_standardize.mbt — sequence standardization pipeline
/// (MolStandardize analogue, item B3)
///
/// The "U+T" of the chemistry world: uppercase, strip terminal/unknown
/// markers and gaps, validate against a target alphabet with an explicit
/// ambiguity-code policy, and report every action taken so the caller can
/// audit what happened to the sequence.

///|
/// Policy for IUPAC ambiguity codes (R/Y/N/...): keep them as-is or
/// reject the record.
pub(all) enum SeqAmbiguityPolicy {
  Keep
  Reject
}

///|
/// Standardization options.
pub struct SeqStandardizeOptions {
  uppercase : Bool
  strip_stars : Bool
  strip_gaps : Bool
  strip_spaces : Bool
  ambiguity_policy : SeqAmbiguityPolicy
  target_alphabet : String // "dna" | "rna" | "protein" | "" (auto-detect)
}

///|
/// Audit report: actions applied in order, invalid characters found
/// (position, char) and the number of removed characters.
pub struct SeqStandardizeReport {
  actions : Array[String]
  invalid_chars : Array[(Int32, String)]
  removed_count : Int32
}

///|
/// Sensible defaults: uppercase + strip stars/gaps/spaces, keep
/// ambiguity codes, auto-detect the alphabet.
pub fn seq_standardize_defaults() -> SeqStandardizeOptions {
  SeqStandardizeOptions::{
    uppercase: true,
    strip_stars: true,
    strip_gaps: true,
    strip_spaces: true,
    ambiguity_policy: Keep,
    target_alphabet: "",
  }
}

///|
/// Public constructor for external callers (the struct itself is
/// read-only outside this package).
pub fn SeqStandardizeOptions::new(
  uppercase~ : Bool = true,
  strip_stars~ : Bool = true,
  strip_gaps~ : Bool = true,
  strip_spaces~ : Bool = true,
  ambiguity_policy~ : SeqAmbiguityPolicy = Keep,
  target_alphabet~ : String = "",
) -> SeqStandardizeOptions {
  SeqStandardizeOptions::{
    uppercase,
    strip_stars,
    strip_gaps,
    strip_spaces,
    ambiguity_policy,
    target_alphabet,
  }
}

///|
fn seq_std_is_iupac_ambiguity(c : Char) -> Bool {
  match c {
    'R' | 'Y' | 'K' | 'M' | 'S' | 'W' | 'B' | 'D' | 'H' | 'V' | 'N' => true
    _ => false
  }
}

///|
fn seq_std_is_core_dna(c : Char) -> Bool {
  match c {
    'A' | 'C' | 'G' | 'T' => true
    _ => false
  }
}

///|
fn seq_std_is_core_rna(c : Char) -> Bool {
  match c {
    'A' | 'C' | 'G' | 'U' => true
    _ => false
  }
}

///|
fn seq_std_is_core_protein(c : Char) -> Bool {
  match c {
    'A' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'K' | 'L' | 'M' | 'N' | 'P' | 'Q' | 'R' | 'S' | 'T' | 'V' | 'W' | 'Y' => true
    'B' | 'X' | 'Z' | 'J' | 'U' | 'O' => true
    _ => false
  }
}

///|
/// Auto-detect the alphabet: RNA when U is present without T, protein
/// when protein-specific residues appear, DNA otherwise.
pub fn seq_std_detect_alphabet(s : String) -> String {
  let n = s.length()
  let mut has_u = false
  let mut has_t = false
  let mut protein_specific = false
  let mut i = 0
  while i < n {
    let c = s.unsafe_get(i).unsafe_to_char()
    match c {
      'U' => has_u = true
      'T' => has_t = true
      'D' | 'E' | 'F' | 'I' | 'L' | 'K' | 'M' | 'P' | 'Q' | 'W' | 'Y' => protein_specific = true
      _ => ()
    }
    i = i + 1
  }
  if has_u && !has_t {
    "rna"
  } else if protein_specific {
    "protein"
  } else {
    "dna"
  }
}

///|
/// Standardize one record following the options. Returns Ok of the
/// cleaned record (with an "standardize" annotation) plus the audit
/// report, or Err of the report listing the invalid characters that
/// blocked acceptance.
pub fn seq_standardize(
  record : SeqRecord,
  opts : SeqStandardizeOptions,
) -> Result[(SeqRecord, SeqStandardizeReport), SeqStandardizeReport] {
  let actions : Array[String] = Array(capacity=0)
  let invalid_chars : Array[(Int32, String)] = Array(capacity=0)
  let s = record.seq.to_string()
  let n = s.length()
  let sb = StringBuilder()
  let mut removed = 0
  let mut i = 0
  // Pass 1: case + strip
  while i < n {
    let c = s.unsafe_get(i).unsafe_to_char()
    let ci = c.to_int()
    let is_lower = ci >= 97 && ci <= 122
    let mut skip = false
    if opts.strip_stars && c == '*' {
      skip = true
    }
    if opts.strip_gaps && (c == '-' || c == '.') {
      skip = true
    }
    if opts.strip_spaces && (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
      skip = true
    }
    if skip {
      removed = removed + 1
    } else if is_lower && opts.uppercase {
      let _ = sb.write_char((ci - 32).unsafe_to_char())
    } else {
      let _ = sb.write_char(c)
    }
    i = i + 1
  }
  if opts.uppercase {
    actions.push("uppercase")
  }
  if removed > 0 {
    actions.push("stripped " + removed.to_string() + " chars")
  }
  let cleaned = sb.to_string()
  // Pass 2: alphabet validation
  let alphabet = if opts.target_alphabet.length() > 0 {
    opts.target_alphabet
  } else {
    seq_std_detect_alphabet(cleaned)
  }
  actions.push("alphabet=" + alphabet)
  let m = cleaned.length()
  let mut j = 0
  while j < m {
    let c = cleaned.unsafe_get(j).unsafe_to_char()
    let is_core = match alphabet {
      "rna" => seq_std_is_core_rna(c)
      "protein" => seq_std_is_core_protein(c)
      _ => seq_std_is_core_dna(c)
    }
    if seq_std_is_iupac_ambiguity(c) {
      match opts.ambiguity_policy {
        Keep => ()
        Reject => invalid_chars.push((j, c.to_string()))
      }
    } else if !is_core {
      invalid_chars.push((j, c.to_string()))
    }
    j = j + 1
  }
  let report : SeqStandardizeReport = {
    actions,
    invalid_chars,
    removed_count: removed,
  }
  if report.invalid_chars.length() > 0 {
    Err(report)
  } else {
    let new_seq = Seq::new(cleaned)
    let annotations : Map[String, String] = Map([], capacity=2)
    for k, v in record.annotations {
      annotations[k] = v
    }
    annotations["alphabet"] = alphabet
    let out : SeqRecord = {
      seq: new_seq,
      id: record.id,
      name: record.name,
      description: record.description,
      letter_annotations: record.letter_annotations,
      annotations,
      features: record.features,
    }
    Ok((out, report))
  }
}

///|
/// Standardize a batch: returns (accepted records, rejected records with
/// their reports).
pub fn seq_standardize_batch(
  records : Array[SeqRecord],
  opts : SeqStandardizeOptions,
) -> (Array[SeqRecord], Array[(SeqRecord, SeqStandardizeReport)]) {
  let accepted : Array[SeqRecord] = Array(capacity=0)
  let rejected : Array[(SeqRecord, SeqStandardizeReport)] = Array(capacity=0)
  for r in records {
    match seq_standardize(r, opts) {
      Ok(pair) => {
        let (clean, _) = pair
        accepted.push(clean)
      }
      Err(report) => rejected.push((r, report))
    }
  }
  (accepted, rejected)
}
