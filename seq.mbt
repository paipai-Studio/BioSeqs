// A reimplementation of BioPython's Bio.Seq.Seq in MoonBit.
//
// Seq is an immutable biological sequence object (essentially a string with
// biological methods). It mirrors the public API of Bio.Seq.Seq:
//   - complement / complement_rna / reverse_complement / reverse_complement_rna
//   - transcribe / back_transcribe
//   - translate (Standard NCBI table 1, with to_stop / cds / gap options)
//   - string-like methods: upper, lower, count, count_overlap, find, rfind,
//     startswith, endswith, contains, strip, lstrip, rstrip, split, join,
//     replace, add (+), mul (*), slicing, length, equality.

///|
/// Error raised by `translate` for invalid codons or invalid CDS sequences,
/// mirroring Bio.Data.CodonTable.TranslationError.
pub suberror SeqError {
  SeqError(String)
}

///|
/// Find the first index of `sub` at or after `from` (character index), or
/// None. Implemented manually because MoonBit's String::find has no start
/// offset parameter.
fn find_from(s : String, sub : String, from : Int) -> Int? {
  let n = s.char_length()
  let m = sub.char_length()
  if m == 0 {
    return Some(from)
  }
  if from < 0 || from + m > n {
    return None
  }
  let mut i = from
  while i + m <= n {
    let mut matched = true
    let mut j = 0
    while j < m {
      match (s.get_char(i + j), sub.get_char(j)) {
        (Some(a), Some(b)) =>
          if a != b {
            matched = false
            break
          }
        _ => {
          matched = false
          break
        }
      }
      j = j + 1
    }
    if matched {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
/// An immutable biological sequence, equivalent to Bio.Seq.Seq.
pub struct Seq {
  data : String
} derive(Eq, Hash)

///|
/// Construct a Seq from a String.
pub fn Seq::new(s : String) -> Seq {
  Seq::{ data: s }
}

///|
/// Construct a Seq from a string literal / StringView.
pub fn Seq::from_str(s : String) -> Seq {
  Seq::{ data: s }
}

///|
/// Return the underlying sequence as a String.
pub fn Seq::to_string(self : Seq) -> String {
  self.data
}

///|
/// Return a debug representation like BioPython's repr: Seq('ACGT').
pub fn Seq::repr(self : Seq) -> String {
  "Seq('" + self.data + "')"
}

///|
/// Number of characters in the sequence.
pub fn Seq::length(self : Seq) -> Int {
  self.data.char_length()
}

///|
/// Concatenate two sequences, like Seq("AC") + Seq("GT"). Implements the
/// `Add` trait so the `+` operator works.
pub impl Add for Seq with fn add(self, other) {
  Seq::{ data: self.data + other.data }
}

///|
/// Concatenate a Seq with a String.
pub fn Seq::add_str(self : Seq, other : String) -> Seq {
  Seq::{ data: self.data + other }
}

///|
/// Repeat the sequence `n` times. (The `*` operator is not used because
/// MoonBit's Mul trait requires both operands to be the same type.)
pub fn Seq::mul(self : Seq, n : Int) -> Seq {
  if n <= 0 {
    Seq::{ data: "" }
  } else {
    Seq::{ data: self.data.repeat(n) }
  }
}

///|
/// Return a subsequence from [start, end) (character indices). Matches
/// Python slicing semantics for non-negative bounds.
pub fn Seq::slice(self : Seq, start~ : Int, end~ : Int) -> Seq {
  Seq::{ data: self.data.substring(start~, end~) }
}

///|
/// Return the character at position `i`, or None if out of bounds.
pub fn Seq::get_char(self : Seq, i : Int) -> Char? {
  self.data.get_char(i)
}

///|
/// Return an uppercased copy.
pub fn Seq::upper(self : Seq) -> Seq {
  Seq::{ data: self.data.to_upper() }
}

///|
/// Return a lowercased copy.
pub fn Seq::lower(self : Seq) -> Seq {
  Seq::{ data: self.data.to_lower() }
}

///|
/// Count non-overlapping occurrences of `sub` in the full sequence, like
/// BioPython's Seq.count.
pub fn Seq::count(self : Seq, sub : String) -> Int {
  if sub.length() == 0 {
    return self.length() + 1
  }
  let mut total = 0
  let data = self.data
  let mut start = 0
  let sub_len = sub.length()
  while start <= data.length() {
    match find_from(data, sub, start) {
      Some(pos) => {
        total = total + 1
        start = pos + sub_len
      }
      None => break
    }
  }
  total
}

///|
/// Count occurrences of `sub`, including overlaps. Mirrors BioPython's
/// Seq.count_overlap.
pub fn Seq::count_overlap(self : Seq, sub : String) -> Int {
  if sub.length() == 0 {
    return self.length() + 1
  }
  let mut total = 0
  let data = self.data
  let mut start = 0
  while start <= data.length() {
    match find_from(data, sub, start) {
      Some(pos) => {
        total = total + 1
        start = pos + 1
      }
      None => break
    }
  }
  total
}

///|
/// Find the first index of `sub`, or -1 if not found (Python str.find
/// semantics).
pub fn Seq::find(self : Seq, sub : String) -> Int {
  match self.data.find(sub) {
    Some(p) => p
    None => -1
  }
}

///|
/// Find the last index of `sub`, or -1 if not found.
pub fn Seq::rfind(self : Seq, sub : String) -> Int {
  match self.data.rev_find(sub) {
    Some(p) => p
    None => -1
  }
}

///|
/// True if the sequence starts with `prefix`.
pub fn Seq::startswith(self : Seq, prefix : String) -> Bool {
  self.data.has_prefix(prefix)
}

///|
/// True if the sequence ends with `suffix`.
pub fn Seq::endswith(self : Seq, suffix : String) -> Bool {
  self.data.has_suffix(suffix)
}

///|
/// True if `sub` is contained in the sequence.
pub fn Seq::contains(self : Seq, sub : String) -> Bool {
  self.data.contains(sub)
}

///|
/// Remove leading and trailing characters in `chars` (default: whitespace),
/// like Seq.strip.
pub fn Seq::strip(self : Seq, chars~ : String = " ") -> Seq {
  Seq::{ data: self.data.trim(chars~).to_owned() }
}

///|
/// Remove leading characters in `chars`.
pub fn Seq::lstrip(self : Seq, chars~ : String = " ") -> Seq {
  Seq::{ data: self.data.trim_start(chars~).to_owned() }
}

///|
/// Remove trailing characters in `chars`.
pub fn Seq::rstrip(self : Seq, chars~ : String = " ") -> Seq {
  Seq::{ data: self.data.trim_end(chars~).to_owned() }
}

///|
/// Replace every occurrence of `old` with `new`.
pub fn Seq::replace(self : Seq, old~ : String, new~ : String) -> Seq {
  Seq::{ data: self.data.replace_all(old~, new~) }
}

///|
/// Split the sequence on `sep` into a list of Seq objects.
pub fn Seq::split(self : Seq, sep : String) -> Array[Seq] {
  let out = Array::new()
  for part in self.data.split(sep) {
    out.push(Seq::{ data: part.to_owned() })
  }
  out
}

///|
/// Join the sequences in `others` with `self` as the spacer, like
/// BioPython's Seq.join.
pub fn Seq::join(self : Seq, others : Array[Seq]) -> Seq {
  let buf = StringBuilder::new()
  for i, s in others {
    if i > 0 {
      buf.write_string(self.data)
    }
    buf.write_string(s.data)
  }
  Seq::{ data: buf.to_string() }
}

// ---------------------------------------------------------------------------
// Biological methods
// ---------------------------------------------------------------------------

///|
/// Return the complement as a DNA sequence. Any U is treated as a T,
/// matching BioPython's Seq.complement().
pub fn Seq::complement(self : Seq) -> Seq {
  Seq::{ data: complement_chars(self.data, dna_complement_map) }
}

///|
/// Return the complement as an RNA sequence. Any T is treated as a U,
/// matching BioPython's Seq.complement_rna().
pub fn Seq::complement_rna(self : Seq) -> Seq {
  Seq::{ data: complement_chars(self.data, rna_complement_map) }
}

///|
/// Reverse the sequence (character-wise).
pub fn Seq::reverse(self : Seq) -> Seq {
  Seq::{ data: self.data.rev() }
}

///|
/// Return the reverse complement as a DNA sequence.
pub fn Seq::reverse_complement(self : Seq) -> Seq {
  Seq::{ data: complement_chars(self.data.rev(), dna_complement_map) }
}

///|
/// Return the reverse complement as an RNA sequence.
pub fn Seq::reverse_complement_rna(self : Seq) -> Seq {
  Seq::{ data: complement_chars(self.data.rev(), rna_complement_map) }
}

///|
/// Transcribe DNA into RNA (T -> U), like BioPython's Seq.transcribe().
pub fn Seq::transcribe(self : Seq) -> Seq {
  Seq::{ data: self.data.replace_all(old="T", new="U").replace_all(old="t", new="u") }
}

///|
/// Back-transcribe RNA into DNA (U -> T), like BioPython's
/// Seq.back_transcribe().
pub fn Seq::back_transcribe(self : Seq) -> Seq {
  Seq::{ data: self.data.replace_all(old="U", new="T").replace_all(old="u", new="t") }
}

///|
/// Translate a nucleotide sequence into a protein sequence using the
/// Standard NCBI table (table 1).
///
/// Arguments mirror BioPython's Seq.translate():
///  - stop_symbol: character used for stop codons (default '*')
///  - to_stop:     if true, translation stops at the first stop codon
///  - cds:         if true, the sequence is treated as a complete CDS;
///                 it must start with a start codon, be a multiple of three,
///                 and end with a single stop codon. Raises SeqError on
///                 failure.
///  - gap:         gap character (default '-'); a codon of three gaps is
///                 translated as a single gap.
pub fn Seq::translate(
  self : Seq,
  stop_symbol~ : Char = '*',
  to_stop~ : Bool = false,
  cds~ : Bool = false,
  gap~ : Char = '-'
) -> Seq raise SeqError {
  let seq = self.data.to_upper()
  let n = seq.char_length()
  let amino = StringBuilder::new()
  // CDS validation (BioPython _translate_str cds branch).
  if cds {
    let first_codon = seq.substring(start=0, end=3)
    if !(standard_start_codons.contains(first_codon)) {
      raise SeqError("First codon '" + first_codon + "' is not a start codon")
    }
    if n % 3 != 0 {
      raise SeqError(
        "Sequence length " + n.to_string() + " is not a multiple of three",
      )
    }
    let last_codon = seq.substring(start=n - 3, end=n)
    if !(standard_stop_codons.contains(last_codon)) {
      raise SeqError("Final codon '" + last_codon + "' is not a stop codon")
    }
    // Skip start (M) and stop codons.
    amino.write_char('M')
    let mut i = 3
    let stop = n - 3
    while i < stop {
      let codon = seq.substring(start=i, end=i + 3)
      if standard_stop_codons.contains(codon) {
        raise SeqError("Extra in frame stop codon '" + codon + "' found.")
      }
      match standard_forward_table.get(codon) {
        Some(aa) => amino.write_char(aa)
        None => {
          if all_valid_letters(codon) {
            amino.write_char('X')
          } else if codon == gap.to_string().repeat(3) {
            amino.write_char(gap)
          } else {
            raise SeqError("Codon '" + codon + "' is invalid")
          }
        }
      }
      i = i + 3
    }
    return Seq::{ data: amino.to_string() }
  }
  // Non-CDS path.
  let mut i = 0
  while i < n - n % 3 {
    let codon = seq.substring(start=i, end=i + 3)
    match standard_forward_table.get(codon) {
      Some(aa) => amino.write_char(aa)
      None => {
        if standard_stop_codons.contains(codon) {
          if to_stop {
            return Seq::{ data: amino.to_string() }
          }
          amino.write_char(stop_symbol)
        } else if all_valid_letters(codon) {
          amino.write_char('X')
        } else if codon == gap.to_string().repeat(3) {
          amino.write_char(gap)
        } else {
          raise SeqError("Codon '" + codon + "' is invalid")
        }
      }
    }
    i = i + 3
  }
  Seq::{ data: amino.to_string() }
}

///|
/// True if every character of `codon` is a valid (ambiguous) nucleotide
/// letter. Used to detect "possible stop" codons such as TAN or NNN.
fn all_valid_letters(codon : String) -> Bool {
  for c in codon.iter() {
    match valid_nucleotide_letters.get(c) {
      Some(_) => continue
      None => return false
    }
  }
  true
}
