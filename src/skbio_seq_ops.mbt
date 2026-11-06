/// skbio_seq_ops.mbt — Sequence validation, concatenation, hashing, motif search, kmer positions
/// Items 10-14: skbio_seq_validate, skbio_seq_concat, skbio_seq_hash, skbio_motif_search, skbio_seq_kmer_positions

///|
/// Sequence type enum.
pub enum SkbioSeqType {
  SkbioDna
  SkbioRna
  SkbioProtein
}

///|
/// Validate a sequence for the given type.
/// Returns list of error messages (empty = valid).
pub fn skbio_seq_validate(
  seq : String,
  seq_type~ : SkbioSeqType = SkbioDna,
  allow_gaps~ : Bool = true,
  allow_ambiguous~ : Bool = true,
) -> Array[String] {
  let errors : Array[String] = Array(capacity=0)
  let valid_chars = match seq_type {
    SkbioDna => {
      let s = if allow_ambiguous { "ACGTUNRYSWKMBDHV.-" } else { "ACGT.-" }
      if allow_gaps {
        s
      } else {
        "ACGTUNRYSWKMBDHV"
      }
    }
    SkbioRna => {
      let s = if allow_ambiguous { "ACGUUNRYSWKMBDHV.-" } else { "ACGU.-" }
      if allow_gaps {
        s
      } else {
        "ACGUUNRYSWKMBDHV"
      }
    }
    SkbioProtein => {
      let s = if allow_ambiguous {
        "ACDEFGHIKLMNPQRSTVWYBXZJUO.*-"
      } else {
        "ACDEFGHIKLMNPQRSTVWY.-"
      }
      if allow_gaps {
        s
      } else {
        "ACDEFGHIKLMNPQRSTVWYBXZJUO*"
      }
    }
  }
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i].to_int()
    let ch = c.unsafe_to_char().to_string().to_upper()
    let mut found = false
    for j = 0; j < valid_chars.length(); j = j + 1 {
      let vc = valid_chars[j].to_int()
      if vc == c || vc == c + 32 || vc == c - 32 {
        found = true
      }
    }
    if !found {
      errors.push("Invalid character '" + ch + "' at position " + i.to_string())
    }
  }
  errors
}

///|
/// Concatenate two sequences with optional separator.
pub fn skbio_seq_concat(
  seq1 : String,
  seq2 : String,
  separator~ : String = "",
) -> String {
  if separator.length() == 0 {
    seq1 + seq2
  } else {
    seq1 + separator + seq2
  }
}

///|
/// Compute a simple hash/fingerprint of a sequence for deduplication.
/// Uses a rolling polynomial hash (djb2 variant).
pub fn skbio_seq_hash(seq : String, seed~ : Int32 = 5381) -> Int32 {
  let mut hash = seed
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i].to_int()
    hash = (hash * 33 + c) & 0x7FFFFFFF
  }
  hash
}

///|
/// Motif search result.
pub struct SkbioMotifHit {
  position : Int32 // 0-based start position
  matched_seq : String // the matched subsequence
  strand : Int32 // +1 forward, -1 reverse complement
}

///|
/// Search for a motif in a sequence.
/// Supports exact matching and IUPAC ambiguous codes.
/// If search_reverse=true, also searches reverse complement strand.
pub fn skbio_motif_search(
  seq : String,
  motif : String,
  search_reverse~ : Bool = false,
) -> Array[SkbioMotifHit] {
  let results : Array[SkbioMotifHit] = Array(capacity=0)
  let seq_upper = seq.to_upper()
  let motif_upper = motif.to_upper()
  let m_len = motif_upper.length()
  let s_len = seq_upper.length()
  if m_len == 0 || m_len > s_len {
    return results
  }
  // Forward strand search
  for i = 0; i <= s_len - m_len; i = i + 1 {
    let mut is_match = true
    for j = 0; j < m_len; j = j + 1 {
      let s_char = seq_upper[i + j].to_int()
      let m_char = motif_upper[j].to_int()
      if !skbio_iupac_match(m_char, s_char) {
        is_match = false
        break
      }
    }
    if is_match {
      results.push(SkbioMotifHit::{
        position: i,
        matched_seq: seq_upper[i:i + m_len].to_owned(),
        strand: 1,
      })
    }
  }
  // Reverse strand search
  if search_reverse {
    let rc = skbio_revcomp(seq_upper)
    for i = 0; i <= rc.length() - m_len; i = i + 1 {
      let mut is_match = true
      for j = 0; j < m_len; j = j + 1 {
        let s_char = rc[i + j].to_int()
        let m_char = motif_upper[j].to_int()
        if !skbio_iupac_match(m_char, s_char) {
          is_match = false
          break
        }
      }
      if is_match {
        results.push(SkbioMotifHit::{
          position: s_len - i - m_len,
          matched_seq: rc[i:i + m_len].to_owned(),
          strand: -1,
        })
      }
    }
  }
  results
}

///|
/// Check if two characters match under IUPAC ambiguity codes.
fn skbio_iupac_match(motif_char : Int32, seq_char : Int32) -> Bool {
  if motif_char == seq_char {
    return true
  }
  // N matches anything
  if motif_char == 78 || motif_char == 110 {
    return true
  }
  // Build IUPAC sets
  let motif_set = skbio_iupac_set(motif_char)
  let seq_set = skbio_iupac_set(seq_char)
  // Check intersection
  for i = 0; i < motif_set.length(); i = i + 1 {
    for j = 0; j < seq_set.length(); j = j + 1 {
      if motif_set[i] == seq_set[j] {
        return true
      }
    }
  }
  false
}

///|
/// Get the set of concrete bases an IUPAC code represents.
fn skbio_iupac_set(code : Int32) -> Array[Int32] {
  let c = if code >= 97 { code - 32 } else { code }
  match c {
    65 => [65] // A
    67 => [67] // C
    71 => [71] // G
    84 => [84] // T
    85 => [84] // U -> T
    82 => [65, 71] // R = A/G
    89 => [67, 84] // Y = C/T
    83 => [67, 71] // S = C/G
    87 => [65, 84] // W = A/T
    75 => [71, 84] // K = G/T
    77 => [65, 67] // M = A/C
    66 => [67, 71, 84] // B = C/G/T
    68 => [65, 71, 84] // D = A/G/T
    72 => [65, 67, 84] // H = A/C/T
    86 => [65, 67, 71] // V = A/C/G
    78 => [65, 67, 71, 84] // N = any
    _ => [c]
  }
}

///|
/// Reverse complement of a DNA sequence.
fn skbio_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i].to_int()
    let rc = match c {
      65 => 84 // A->T
      84 => 65 // T->A
      71 => 67 // G->C
      67 => 71 // C->G
      85 => 65 // U->A
      78 => 78 // N->N
      82 => 89 // R->Y
      89 => 82 // Y->R
      83 => 83 // S->S
      87 => 87 // W->W
      75 => 77 // K->M
      77 => 75 // M->K
      66 => 86 // B->V
      68 => 72 // D->H
      72 => 68 // H->D
      86 => 66 // V->B
      _ => 78
    }
    let _ = sb.write_char(rc.unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Kmer position index result.
pub struct SkbioKmerPos {
  kmer : String
  positions : Array[Int32] // all positions where this kmer occurs
}

///|
/// Build a kmer position index for a sequence.
/// Returns one entry per unique kmer with all its positions.
pub fn skbio_seq_kmer_positions(
  seq : String,
  k~ : Int32 = 3,
) -> Array[SkbioKmerPos] {
  let s = seq.to_upper()
  let n = s.length()
  if k <= 0 || k > n {
    return Array(capacity=0)
  }
  // Collect all kmers with positions
  let kmers : Array[String] = Array(capacity=0)
  let positions : Array[Int32] = Array(capacity=0)
  for i = 0; i <= n - k; i = i + 1 {
    kmers.push(s[i:i + k].to_owned())
    positions.push(i)
  }
  // Group by unique kmer
  let unique : Array[String] = Array(capacity=0)
  for i = 0; i < kmers.length(); i = i + 1 {
    let mut found = false
    for j = 0; j < unique.length(); j = j + 1 {
      if unique[j] == kmers[i] {
        found = true
      }
    }
    if !found {
      unique.push(kmers[i])
    }
  }
  // Build results
  let results : Array[SkbioKmerPos] = Array(capacity=unique.length())
  for u = 0; u < unique.length(); u = u + 1 {
    let pos : Array[Int32] = Array(capacity=0)
    for i = 0; i < kmers.length(); i = i + 1 {
      if kmers[i] == unique[u] {
        pos.push(positions[i])
      }
    }
    results.push(SkbioKmerPos::{ kmer: unique[u], positions: pos, })
  }
  results
}

///|
pub fn skbio_seqtype_dummy() -> Unit {
  ignore(SkbioRna)
  ignore(SkbioProtein)
}
