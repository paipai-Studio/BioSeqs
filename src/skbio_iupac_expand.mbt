/// IUPAC ambiguity code expansion and regex conversion.
///
/// scikit-bio equivalent: IUPAC ambiguity code handling
///
/// IUPAC nucleotide codes:
/// A=adenine, C=cytosine, G=guanine, T=thymine, U=uracil
/// R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C
/// B=C/G/T, A=A/C/G, D=A/G/T, H=A/C/T, V=A/C/G, N=A/C/G/T

///|
/// Expand an IUPAC ambiguity code to all possible nucleotides.
/// Returns array of single-character strings.
pub fn skbio_expand_iupac(code : Char) -> Array[Char] {
  match code {
    'A' => ['A']
    'C' => ['C']
    'G' => ['G']
    'T' => ['T']
    'U' => ['U']
    'R' => ['A', 'G']
    'Y' => ['C', 'T']
    'S' => ['G', 'C']
    'W' => ['A', 'T']
    'K' => ['G', 'T']
    'M' => ['A', 'C']
    'B' => ['C', 'G', 'T']
    'D' => ['A', 'G', 'T']
    'H' => ['A', 'C', 'T']
    'V' => ['A', 'C', 'G']
    'N' => ['A', 'C', 'G', 'T']
    _ => [code]
  }
}

///|
/// Expand an IUPAC sequence to all possible concrete sequences.
/// Warning: exponential blowup. Only practical for short sequences (<= ~20 ambiguous positions).
pub fn skbio_expand_iupac_seq(seq : String) -> Array[String] {
  let results : Array[String] = Array(capacity=0)
  results.push("")
  for c in seq {
    let bases = skbio_expand_iupac(c)
    let new_results : Array[String] = Array(capacity=0)
    for r in results {
      for b in bases {
        let sb = StringBuilder::new()
        sb.write_string(r)
        sb.write_char(b)
        new_results.push(sb.to_string())
      }
    }
    results.clear()
    for r in new_results {
      results.push(r)
    }
  }
  results
}

///|
/// Convert an IUPAC sequence to a regex pattern.
/// Each IUPAC code becomes a character class.
/// Example: "ARY" → "[A][AG][CT]"
pub fn skbio_iupac_to_regex(seq : String) -> String {
  let sb = StringBuilder::new()
  for c in seq {
    let bases = skbio_expand_iupac(c)
    if bases.length() == 1 {
      sb.write_char(bases[0])
    } else {
      sb.write_char('[')
      for b in bases {
        sb.write_char(b)
      }
      sb.write_char(']')
    }
  }
  sb.to_string()
}

///|
/// Count the number of ambiguous (non-ACGTU) positions in a sequence.
pub fn skbio_count_ambiguous(seq : String) -> Int32 {
  let mut count : Int32 = 0
  for c in seq {
    match c {
      'A' => ()
      'C' => ()
      'G' => ()
      'T' => ()
      'U' => ()
      'a' => ()
      'c' => ()
      'g' => ()
      't' => ()
      'u' => ()
      _ => count = count + 1
    }
  }
  count
}

///|
/// Check if a sequence contains any IUPAC ambiguity codes.
pub fn skbio_has_ambiguous(seq : String) -> Bool {
  for c in seq {
    match c {
      'A' => ()
      'C' => ()
      'G' => ()
      'T' => ()
      'U' => ()
      'a' => ()
      'c' => ()
      'g' => ()
      't' => ()
      'u' => ()
      _ => return true
    }
  }
  false
}

///|
/// Get the IUPAC code for a set of nucleotides.
/// Inverse of expand_iupac: {A,G} → R, {C,T} → Y, etc.
pub fn skbio_iupac_code(bases : Array[Char]) -> Char {
  // Sort bases for consistent lookup
  let sorted = bases.copy()
  sorted.sort()
  let key = {
    let sb = StringBuilder::new()
    for b in sorted {
      sb.write_char(b)
    }
    sb.to_string()
  }
  match key {
    "A" => 'A'
    "C" => 'C'
    "G" => 'G'
    "T" => 'T'
    "U" => 'U'
    "AG" => 'R'
    "CT" => 'Y'
    "CG" => 'S'
    "AT" => 'W'
    "GT" => 'K'
    "AC" => 'M'
    "CGT" => 'B'
    "AGT" => 'D'
    "ACT" => 'H'
    "ACG" => 'V'
    "ACGT" => 'N'
    _ => 'N'
  }
}

///|
/// Merge two IUPAC codes (union of their expansions).
pub fn skbio_iupac_merge(a : Char, b : Char) -> Char {
  let bases_a = skbio_expand_iupac(a)
  let bases_b = skbio_expand_iupac(b)
  let merged : Set[Char] = Set::from_array(bases_a)
  for c in bases_b {
    merged.add(c)
  }
  skbio_iupac_code(merged.to_array())
}
