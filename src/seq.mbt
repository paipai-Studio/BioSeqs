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
/// Find the first index of `sub` at or after `from` (code-unit index), or
/// None. Uses `unsafe_get` (O(1)) instead of `get_char` (O(n)) for speed.
/// For ASCII biological sequences, code-unit index == character index.
fn find_from(s : String, sub : String, from : Int) -> Int? {
  let n = s.length()
  let m = sub.length()
  if m == 0 {
    return Some(from)
  }
  if from < 0 || from + m > n {
    return None
  }
  let first = sub.unsafe_get(0)
  let mut i = from
  while i + m <= n {
    if s.unsafe_get(i) != first {
      i = i + 1
      continue
    }
    let mut matched = true
    for j in 1..<m {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
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
  Seq::{ data: self.data[start:end].to_owned() }
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
  let sub_len = sub.length()
  if sub_len == 0 {
    return self.data.length() + 1
  }
  let mut total = 0
  let data = self.data
  let mut start = 0
  let n = data.length()
  while start <= n {
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
  let sub_len = sub.length()
  if sub_len == 0 {
    return self.data.length() + 1
  }
  let mut total = 0
  let data = self.data
  let mut start = 0
  let n = data.length()
  while start <= n {
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
/// Implemented manually because MoonBit's `has_prefix` calls `find` (full
/// string search) which is O(n) instead of O(prefix_len).
pub fn Seq::startswith(self : Seq, prefix : String) -> Bool {
  let plen = prefix.length()
  let data = self.data
  if data.length() < plen {
    return false
  }
  for i = 0; i < plen; i = i + 1 {
    if data.unsafe_get(i) != prefix.unsafe_get(i) {
      return false
    }
  }
  true
}

///|
/// True if the sequence ends with `suffix`.
/// Implemented manually because MoonBit's `has_suffix` calls `rev_find`
/// (full string search) which is O(n) instead of O(suffix_len).
pub fn Seq::endswith(self : Seq, suffix : String) -> Bool {
  let slen = suffix.length()
  let data = self.data
  let dlen = data.length()
  if dlen < slen {
    return false
  }
  let offset = dlen - slen
  for i = 0; i < slen; i = i + 1 {
    if data.unsafe_get(offset + i) != suffix.unsafe_get(i) {
      return false
    }
  }
  true
}

///|
/// True if `sub` is contained in the sequence.
pub fn Seq::contains(self : Seq, sub : String) -> Bool {
  self.data.contains(sub)
}

///|
/// Remove leading and trailing characters in `chars` (default: whitespace),
/// like Seq.strip.
pub fn Seq::strip(self : Seq, chars? : String = " ") -> Seq {
  Seq::{ data: self.data.trim(chars~).to_owned() }
}

///|
/// Remove leading characters in `chars`.
pub fn Seq::lstrip(self : Seq, chars? : String = " ") -> Seq {
  Seq::{ data: self.data.trim_start(chars~).to_owned() }
}

///|
/// Remove trailing characters in `chars`.
pub fn Seq::rstrip(self : Seq, chars? : String = " ") -> Seq {
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
  Seq::{ data: complement_chars(self.data, dna_complement_arr) }
}

///|
/// Return the complement as an RNA sequence. Any T is treated as a U,
/// matching BioPython's Seq.complement_rna().
pub fn Seq::complement_rna(self : Seq) -> Seq {
  Seq::{ data: complement_chars(self.data, rna_complement_arr) }
}

///|
/// Reverse the sequence (character-wise). Uses `unsafe_get` from the end for
/// O(1) per-character access and writes directly to a FixedArray[UInt16] for
/// fast string construction.
pub fn Seq::reverse(self : Seq) -> Seq {
  let src = self.data
  let n = src.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    buf[i] = src.unsafe_get(n - 1 - i)
  }
  Seq::{ data: unsafe_from_uint16_fixedarray(buf) }
}

///|
/// Return the reverse complement as a DNA sequence. Uses a single-pass
/// reverse-complement to avoid allocating an intermediate reversed string.
pub fn Seq::reverse_complement(self : Seq) -> Seq {
  Seq::{ data: reverse_complement_chars(self.data, dna_complement_arr) }
}

///|
/// Return the reverse complement as an RNA sequence. Uses a single-pass
/// reverse-complement to avoid allocating an intermediate reversed string.
pub fn Seq::reverse_complement_rna(self : Seq) -> Seq {
  Seq::{ data: reverse_complement_chars(self.data, rna_complement_arr) }
}

///|
/// Transcribe DNA into RNA (T -> U), like BioPython's Seq.transcribe().
/// Single-pass character replacement using FixedArray[UInt16] for fast
/// string construction.
pub fn Seq::transcribe(self : Seq) -> Seq {
  let src = self.data
  let n = src.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  let t : UInt16 = 'T'.to_int().to_uint16()
  let tl : UInt16 = 't'.to_int().to_uint16()
  let u : UInt16 = 'U'.to_int().to_uint16()
  let ul : UInt16 = 'u'.to_int().to_uint16()
  for i = 0; i < n; i = i + 1 {
    let ci = src.unsafe_get(i)
    if ci == t {
      buf[i] = u
    } else if ci == tl {
      buf[i] = ul
    } else {
      buf[i] = ci
    }
  }
  Seq::{ data: unsafe_from_uint16_fixedarray(buf) }
}

///|
/// Back-transcribe RNA into DNA (U -> T), like BioPython's
/// Seq.back_transcribe(). Single-pass character replacement.
pub fn Seq::back_transcribe(self : Seq) -> Seq {
  let src = self.data
  let n = src.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  let u : UInt16 = 'U'.to_int().to_uint16()
  let ul : UInt16 = 'u'.to_int().to_uint16()
  let t : UInt16 = 'T'.to_int().to_uint16()
  let tl : UInt16 = 't'.to_int().to_uint16()
  for i = 0; i < n; i = i + 1 {
    let ci = src.unsafe_get(i)
    if ci == u {
      buf[i] = t
    } else if ci == ul {
      buf[i] = tl
    } else {
      buf[i] = ci
    }
  }
  Seq::{ data: unsafe_from_uint16_fixedarray(buf) }
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
  stop_symbol? : Char = '*',
  to_stop? : Bool = false,
  cds? : Bool = false,
  gap? : Char = '-',
) -> Seq raise SeqError {
  let seq = self.data.to_upper()
  let n = seq.length()
  let max_amino = n / 3
  let amino_buf : FixedArray[UInt16] = FixedArray::make(max_amino, 0)
  let mut amino_len = 0
  let stop_code = stop_symbol.to_int().to_uint16()
  let gap_code = gap.to_int().to_uint16()
  let m_code = 'M'.to_int().to_uint16()
  let x_code = 'X'.to_int().to_uint16()

  if cds {
    if !is_start_codon_at(seq, 0) {
      let first_codon = codon_to_string(seq, 0)
      raise SeqError("First codon '" + first_codon + "' is not a start codon")
    }
    if n % 3 != 0 {
      raise SeqError(
        "Sequence length " + n.to_string() + " is not a multiple of three",
      )
    }
    if !is_stop_codon_at(seq, n - 3) {
      let last_codon = codon_to_string(seq, n - 3)
      raise SeqError("Final codon '" + last_codon + "' is not a stop codon")
    }
    amino_buf[0] = m_code
    amino_len = 1
    let mut i = 3
    let stop = n - 3
    while i < stop {
      if is_stop_codon_at(seq, i) {
        let codon = codon_to_string(seq, i)
        raise SeqError("Extra in frame stop codon '" + codon + "' found.")
      }
      match lookup_codon(seq, i) {
        Some(aa) => {
          amino_buf[amino_len] = aa.to_int().to_uint16()
          amino_len = amino_len + 1
        }
        None =>
          if all_valid_letters_at(seq, i) {
            amino_buf[amino_len] = x_code
            amino_len = amino_len + 1
          } else if is_gap_codon_at(seq, i, gap) {
            amino_buf[amino_len] = gap_code
            amino_len = amino_len + 1
          } else {
            let codon = codon_to_string(seq, i)
            raise SeqError("Codon '" + codon + "' is invalid")
          }
      }
      i = i + 3
    }
    let final_buf : FixedArray[UInt16] = FixedArray::make(amino_len, 0)
    for k = 0; k < amino_len; k = k + 1 {
      final_buf[k] = amino_buf[k]
    }
    return Seq::{ data: unsafe_from_uint16_fixedarray(final_buf) }
  }
  let mut i = 0
  while i < n - n % 3 {
    match lookup_codon(seq, i) {
      Some(aa) => {
        amino_buf[amino_len] = aa.to_int().to_uint16()
        amino_len = amino_len + 1
      }
      None =>
        if is_stop_codon_at(seq, i) {
          if to_stop {
            let final_buf : FixedArray[UInt16] = FixedArray::make(amino_len, 0)
            for k = 0; k < amino_len; k = k + 1 {
              final_buf[k] = amino_buf[k]
            }
            return Seq::{ data: unsafe_from_uint16_fixedarray(final_buf) }
          }
          amino_buf[amino_len] = stop_code
          amino_len = amino_len + 1
        } else if all_valid_letters_at(seq, i) {
          amino_buf[amino_len] = x_code
          amino_len = amino_len + 1
        } else if is_gap_codon_at(seq, i, gap) {
          amino_buf[amino_len] = gap_code
          amino_len = amino_len + 1
        } else {
          let codon = codon_to_string(seq, i)
          raise SeqError("Codon '" + codon + "' is invalid")
        }
    }
    i = i + 3
  }
  let final_buf : FixedArray[UInt16] = FixedArray::make(amino_len, 0)
  for k = 0; k < amino_len; k = k + 1 {
    final_buf[k] = amino_buf[k]
  }
  Seq::{ data: unsafe_from_uint16_fixedarray(final_buf) }
}

///|
/// Read 3 characters starting at index `i` as a String. Used only for error
/// messages.
fn codon_to_string(seq : String, i : Int) -> String {
  let buf = StringBuilder::new(size_hint=3)
  buf.write_char(seq.unsafe_get(i).unsafe_to_char())
  buf.write_char(seq.unsafe_get(i + 1).unsafe_to_char())
  buf.write_char(seq.unsafe_get(i + 2).unsafe_to_char())
  buf.to_string()
}

///|
/// Check if the 3 characters at `i` form a stop codon (TAA/TAG/TGA) without
/// allocating a substring.
fn is_stop_codon_at(seq : String, i : Int) -> Bool {
  let c0 = seq.unsafe_get(i).to_int()
  let c1 = seq.unsafe_get(i + 1).to_int()
  let c2 = seq.unsafe_get(i + 2).to_int()
  let t = 'T'.to_int()
  let a = 'A'.to_int()
  let g = 'G'.to_int()
  // TAA, TAG, TGA
  (c0 == t && c1 == a && (c2 == a || c2 == g)) ||
  (c0 == t && c1 == g && c2 == a)
}

///|
/// Check if the 3 characters at `i` form a start codon (ATG/CTG/TTG) without
/// allocating a substring.
fn is_start_codon_at(seq : String, i : Int) -> Bool {
  let c0 = seq.unsafe_get(i).to_int()
  let c1 = seq.unsafe_get(i + 1).to_int()
  let c2 = seq.unsafe_get(i + 2).to_int()
  let a = 'A'.to_int()
  let c = 'C'.to_int()
  let t = 'T'.to_int()
  let g = 'G'.to_int()
  // ATG, CTG, TTG
  (c0 == a || c0 == c || c0 == t) && c1 == t && c2 == g
}

///|
/// Look up the amino acid for the codon at index `i` without allocating a
/// substring. Uses the precomputed codon-to-Int key for fast Map lookup.
fn lookup_codon(seq : String, i : Int) -> Char? {
  let key = codon_key(seq, i)
  codon_forward_table.get(key)
}

///|
/// Encode 3 characters at index `i` as an Int key for Map lookup.
fn codon_key(seq : String, i : Int) -> Int {
  let c0 = seq.unsafe_get(i).to_int()
  let c1 = seq.unsafe_get(i + 1).to_int()
  let c2 = seq.unsafe_get(i + 2).to_int()
  c0 * 65536 + c1 * 256 + c2
}

///|
/// Check if the 3 characters at `i` are all valid (ambiguous) nucleotide
/// letters, without allocating a substring.
fn all_valid_letters_at(seq : String, i : Int) -> Bool {
  valid_letter(seq.unsafe_get(i).to_int()) &&
  valid_letter(seq.unsafe_get(i + 1).to_int()) &&
  valid_letter(seq.unsafe_get(i + 2).to_int())
}

///|
/// Check if the 3 characters at `i` are all equal to `gap`.
fn is_gap_codon_at(seq : String, i : Int, gap : Char) -> Bool {
  let g = gap.to_int()
  seq.unsafe_get(i).to_int() == g &&
  seq.unsafe_get(i + 1).to_int() == g &&
  seq.unsafe_get(i + 2).to_int() == g
}
