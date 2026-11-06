/// skbio_msa_subst.mbt — Advanced MSA ops, substitution matrix lookup, reverse complement
/// Items 18-20: skbio_msa_advanced, skbio_substitution_lookup, skbio_seq_reverse_complement

///|
/// MSA record: ID + aligned sequence.
pub struct SkbioMsaRecord {
  id : String
  seq : String
}

///|
/// Remove gap-only columns from an MSA.
/// Returns the filtered MSA and the list of kept column indices.
pub fn skbio_msa_remove_gap_cols(
  records : Array[SkbioMsaRecord],
) -> (Array[SkbioMsaRecord], Array[Int32]) {
  if records.length() == 0 {
    return (Array(capacity=0), Array(capacity=0))
  }
  let aln_len = records[0].seq.length()
  // Determine which columns are all gaps
  let keep_cols : Array[Int32] = Array(capacity=aln_len)
  for col = 0; col < aln_len; col = col + 1 {
    let mut all_gap = true
    for r = 0; r < records.length(); r = r + 1 {
      if col < records[r].seq.length() {
        let c = records[r].seq[col].to_int()
        if c != 45 && c != 46 { // not '-' or '.'
          all_gap = false
          break
        }
      }
    }
    if !all_gap {
      keep_cols.push(col)
    }
  }
  // Build filtered sequences
  let result : Array[SkbioMsaRecord] = Array(capacity=records.length())
  for r = 0; r < records.length(); r = r + 1 {
    let sb = StringBuilder()
    for ci = 0; ci < keep_cols.length(); ci = ci + 1 {
      let col = keep_cols[ci]
      if col < records[r].seq.length() {
        let _ = sb.write_char(records[r].seq[col].to_int().unsafe_to_char())
      }
    }
    result.push(SkbioMsaRecord::{ id: records[r].id, seq: sb.to_string(), })
  }
  (result, keep_cols)
}

///|
/// Sort MSA records by ID.
pub fn skbio_msa_sort_by_id(
  records : Array[SkbioMsaRecord],
) -> Array[SkbioMsaRecord] {
  let arr = Array(capacity=records.length())
  for i = 0; i < records.length(); i = i + 1 {
    arr.push(records[i])
  }
  arr.sort_by(fn(a, b) {
    if a.id < b.id {
      -1
    } else if a.id > b.id {
      1
    } else {
      0
    }
  })
  arr
}

///|
/// Reindex MSA: assign new sequential IDs.
pub fn skbio_msa_reindex(
  records : Array[SkbioMsaRecord],
  prefix~ : String = "seq",
) -> Array[SkbioMsaRecord] {
  let result : Array[SkbioMsaRecord] = Array(capacity=records.length())
  for i = 0; i < records.length(); i = i + 1 {
    result.push(SkbioMsaRecord::{
      id: prefix + "_" + i.to_string(),
      seq: records[i].seq,
    })
  }
  result
}

///|
/// Filter MSA to keep only conserved columns.
/// A column is conserved if the most common residue appears in >= min_fraction of sequences.
pub fn skbio_msa_filter_conserved(
  records : Array[SkbioMsaRecord],
  min_fraction~ : Double = 0.8,
) -> (Array[SkbioMsaRecord], Array[Int32]) {
  if records.length() == 0 {
    return (Array(capacity=0), Array(capacity=0))
  }
  let aln_len = records[0].seq.length()
  let n = records.length()
  let keep_cols : Array[Int32] = Array(capacity=aln_len)
  for col = 0; col < aln_len; col = col + 1 {
    // Count residues (excluding gaps)
    let counts : Array[(Int32, Int32)] = Array(capacity=0) // (char_code, count)
    let mut non_gap = 0
    for r = 0; r < n; r = r + 1 {
      if col < records[r].seq.length() {
        let c = records[r].seq[col].to_int()
        if c != 45 && c != 46 {
          non_gap = non_gap + 1
          let mut found = false
          for k = 0; k < counts.length(); k = k + 1 {
            if counts[k].0 == c {
              counts[k] = (c, counts[k].1 + 1)
              found = true
            }
          }
          if !found {
            counts.push((c, 1))
          }
        }
      }
    }
    // Find max count
    let mut max_count = 0
    for k = 0; k < counts.length(); k = k + 1 {
      if counts[k].1 > max_count {
        max_count = counts[k].1
      }
    }
    if non_gap > 0 {
      let frac = max_count.to_double() / non_gap.to_double()
      if frac >= min_fraction {
        keep_cols.push(col)
      }
    }
  }
  // Build filtered sequences
  let result : Array[SkbioMsaRecord] = Array(capacity=n)
  for r = 0; r < n; r = r + 1 {
    let sb = StringBuilder()
    for ci = 0; ci < keep_cols.length(); ci = ci + 1 {
      let col = keep_cols[ci]
      if col < records[r].seq.length() {
        let _ = sb.write_char(records[r].seq[col].to_int().unsafe_to_char())
      }
    }
    result.push(SkbioMsaRecord::{ id: records[r].id, seq: sb.to_string(), })
  }
  (result, keep_cols)
}

///|
/// Substitution matrix: lookup table for amino acid pairs.
pub struct SkbioSubstMatrix {
  name : String
  // Stored as flat array: index = aa1 * 25 + aa2 (A=0..Z=25, using ASCII-65)
  scores : Array[Int32]
}

///|
/// Get score for an amino acid pair.
pub fn SkbioSubstMatrix::get(
  self : SkbioSubstMatrix,
  aa1 : String,
  aa2 : String,
) -> Int32 {
  if aa1.length() == 0 || aa2.length() == 0 {
    return 0
  }
  let c1 = aa1[0].to_int() - 65
  let c2 = aa2[0].to_int() - 65
  if c1 < 0 || c1 >= 25 || c2 < 0 || c2 >= 25 {
    return 0
  }
  self.scores[c1 * 25 + c2]
}

///|
/// BLOSUM62 substitution matrix (partial — key amino acids).
pub fn skbio_blosum62() -> SkbioSubstMatrix {
  let scores = Array(capacity=625)
  for i = 0; i < 625; i = i + 1 {
    scores.push(0)
  }
  // A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9,
  // M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19
  // Set diagonal (self-match) values
  let diag = [4, 9, 6, 5, 6, 6, 8, 4, 5, 4, 5, 6, 7, 5, 5, 4, 5, 4, 11, 7]
  for i = 0; i < 20; i = i + 1 {
    scores[i * 25 + i] = diag[i]
  }
  // Key off-diagonal values (symmetric)
  let pairs : Array[(Int32, Int32, Int32)] = [
    (0, 5, 0), // A-G
    (0, 15, 1), // A-S
    (0, 16, 1), // A-T
    (0, 17, 0), // A-V
    (1, 11, -3), // C-N
    (1, 19, -2), // C-Y
    (2, 3, 2), // D-E
    (2, 11, 1), // D-N
    (3, 13, 2), // E-Q
    (4, 19, 3), // F-Y
    (5, 15, 0), // G-S
    (7, 17, 3), // I-V
    (8, 14, 2), // K-R
    (9, 10, 2), // L-M
    (9, 17, 1), // L-V
    (11, 15, 1), // N-S
    (13, 14, 1), // Q-R
    (15, 16, 1), // S-T
    (18, 19, 2), // W-Y
  ]
  for i = 0; i < pairs.length(); i = i + 1 {
    let (a, b, v) = pairs[i]
    scores[a * 25 + b] = v
    scores[b * 25 + a] = v
  }
  // Set mismatch defaults to -1 for known AAs
  for i = 0; i < 20; i = i + 1 {
    for j = 0; j < 20; j = j + 1 {
      if i != j && scores[i * 25 + j] == 0 {
        scores[i * 25 + j] = -1
      }
    }
  }
  SkbioSubstMatrix::{ name: "BLOSUM62", scores, }
}

///|
/// PAM250 substitution matrix (partial).
pub fn skbio_pam250() -> SkbioSubstMatrix {
  let scores = Array(capacity=625)
  for i = 0; i < 625; i = i + 1 {
    scores.push(0)
  }
  let diag = [2, 2, 2, 2, 3, 1, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 6, 3]
  for i = 0; i < 20; i = i + 1 {
    scores[i * 25 + i] = diag[i]
  }
  // Set mismatch defaults
  for i = 0; i < 20; i = i + 1 {
    for j = 0; j < 20; j = j + 1 {
      if i != j && scores[i * 25 + j] == 0 {
        scores[i * 25 + j] = -2
      }
    }
  }
  SkbioSubstMatrix::{ name: "PAM250", scores, }
}

///|
/// Look up a substitution matrix by name.
pub fn skbio_substitution_lookup(name : String) -> SkbioSubstMatrix {
  let upper = name.to_upper()
  if upper == "BLOSUM62" {
    skbio_blosum62()
  } else if upper == "PAM250" {
    skbio_pam250()
  } else {
    // Default to BLOSUM62
    skbio_blosum62()
  }
}

///|
/// Reverse complement with metadata preservation.
/// Handles IUPAC ambiguous codes and optional case preservation.
pub struct SkbioSeqWithMeta {
  id : String
  seq : String
  metadata : String
}

///|
pub fn skbio_seq_reverse_complement(
  record : SkbioSeqWithMeta,
  preserve_case~ : Bool = false,
) -> SkbioSeqWithMeta {
  let sb = StringBuilder()
  for i = record.seq.length() - 1; i >= 0; i = i - 1 {
    let c = record.seq[i].to_int()
    let is_lower = c >= 97 && c <= 122
    let upper_c = if is_lower { c - 32 } else { c }
    let rc = match upper_c {
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
      45 => 45 // - -> -
      46 => 46 // . -> .
      _ => 78
    }
    let final_c = if preserve_case && is_lower { rc + 32 } else { rc }
    let _ = sb.write_char(final_c.unsafe_to_char())
  }
  SkbioSeqWithMeta::{
    id: record.id,
    seq: sb.to_string(),
    metadata: record.metadata,
  }
}

///|
/// Reverse complement a plain DNA string (convenience wrapper).
pub fn skbio_revcomp_str(seq : String) -> String {
  let sb = StringBuilder()
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i].to_int()
    let upper_c = if c >= 97 { c - 32 } else { c }
    let rc = match upper_c {
      65 => 84
      84 => 65
      71 => 67
      67 => 71
      85 => 65
      78 => 78
      82 => 89
      89 => 82
      83 => 83
      87 => 87
      75 => 77
      77 => 75
      66 => 86
      68 => 72
      72 => 68
      86 => 66
      45 => 45
      _ => 78
    }
    let _ = sb.write_char(rc.unsafe_to_char())
  }
  sb.to_string()
}
