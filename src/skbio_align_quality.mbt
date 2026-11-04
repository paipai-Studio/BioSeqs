/// Alignment column quality and conservation scores.
///
/// scikit-bio equivalent: alignment column analysis
///
/// Per-column metrics for multiple sequence alignments:
/// conservation, Shannon entropy, gap fraction, consensus.

///|
/// Compute per-column conservation score (fraction of most common residue).
/// Returns array of scores (0.0 to 1.0) per column.
pub fn skbio_column_conservation(alignment : Array[String]) -> Array[Double] {
  if alignment.length() == 0 {
    return []
  }
  let aln_len = alignment[0].length()
  let result : Array[Double] = Array::make(aln_len, 0.0)
  for col = 0; col < aln_len; col = col + 1 {
    let counts : Map[Char, Int32] = Map([])
    for seq in alignment {
      if col < seq.length() {
        let c = seq.unsafe_get(col).unsafe_to_char()
        match counts.get(c) {
          Some(n) => counts[c] = n + 1
          None => counts[c] = 1
        }
      }
    }
    // Find max count
    let mut max_count : Int32 = 0
    let keys = counts.keys().to_array()
    for k in keys {
      if counts[k] > max_count {
        max_count = counts[k]
      }
    }
    result[col] = max_count.to_double() / alignment.length().to_double()
  }
  result
}

///|
/// Compute per-column Shannon entropy (in bits).
/// Higher entropy = more variable column.
pub fn skbio_column_entropy(alignment : Array[String]) -> Array[Double] {
  if alignment.length() == 0 {
    return []
  }
  let aln_len = alignment[0].length()
  let result : Array[Double] = Array::make(aln_len, 0.0)
  let n = alignment.length().to_double()
  for col = 0; col < aln_len; col = col + 1 {
    let counts : Map[Char, Int32] = Map([])
    for seq in alignment {
      if col < seq.length() {
        let c = seq.unsafe_get(col).unsafe_to_char()
        match counts.get(c) {
          Some(k) => counts[c] = k + 1
          None => counts[c] = 1
        }
      }
    }
    let mut entropy = 0.0
    let keys = counts.keys().to_array()
    for k in keys {
      let p = counts[k].to_double() / n
      if p > 0.0 {
        entropy = entropy - p * @math.log2(p)
      }
    }
    result[col] = entropy
  }
  result
}

///|
/// Compute per-column gap fraction.
pub fn skbio_column_gap_fraction(
  alignment : Array[String],
  gap~ : Char = '-',
) -> Array[Double] {
  if alignment.length() == 0 {
    return []
  }
  let aln_len = alignment[0].length()
  let result : Array[Double] = Array::make(aln_len, 0.0)
  let n = alignment.length().to_double()
  for col = 0; col < aln_len; col = col + 1 {
    let mut gap_count = 0.0
    for seq in alignment {
      if col < seq.length() && seq.unsafe_get(col).unsafe_to_char() == gap {
        gap_count = gap_count + 1.0
      }
    }
    result[col] = gap_count / n
  }
  result
}

///|
/// Compute per-column consensus residue (most common, ignoring gaps).
/// Returns array of (residue, frequency) pairs.
pub fn skbio_column_consensus(
  alignment : Array[String],
  gap~ : Char = '-',
) -> Array[(Char, Double)] {
  if alignment.length() == 0 {
    return []
  }
  let aln_len = alignment[0].length()
  let result : Array[(Char, Double)] = Array::make(aln_len, ('-', 0.0))
  let n = alignment.length().to_double()
  for col = 0; col < aln_len; col = col + 1 {
    let counts : Map[Char, Int32] = Map([])
    for seq in alignment {
      if col < seq.length() {
        let c = seq.unsafe_get(col).unsafe_to_char()
        if c != gap {
          match counts.get(c) {
            Some(k) => counts[c] = k + 1
            None => counts[c] = 1
          }
        }
      }
    }
    let mut best_char : Char = '-'
    let mut best_count : Int32 = 0
    let keys = counts.keys().to_array()
    for k in keys {
      if counts[k] > best_count {
        best_count = counts[k]
        best_char = k
      }
    }
    result[col] = (best_char, best_count.to_double() / n)
  }
  result
}

///|
/// Overall alignment quality: mean conservation across columns.
pub fn skbio_alignment_quality(alignment : Array[String]) -> Double {
  let conservation = skbio_column_conservation(alignment)
  if conservation.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for c in conservation {
    sum = sum + c
  }
  sum / conservation.length().to_double()
}

///|
/// Compute per-column variability (number of distinct residues, ignoring gaps).
pub fn skbio_column_variability(
  alignment : Array[String],
  gap~ : Char = '-',
) -> Array[Int32] {
  if alignment.length() == 0 {
    return []
  }
  let aln_len = alignment[0].length()
  let result : Array[Int32] = Array::make(aln_len, 0)
  for col = 0; col < aln_len; col = col + 1 {
    let seen : Set[Char] = Set([])
    for seq in alignment {
      if col < seq.length() {
        let c = seq.unsafe_get(col).unsafe_to_char()
        if c != gap {
          seen.add(c)
        }
      }
    }
    result[col] = seen.length()
  }
  result
}
