///|
/// Bio.Align.AlignInfo - Alignment statistics and consensus sequence generation.
/// Inspired by Biopython's Bio.Align.AlignInfo module.

///|
/// Column-wise statistics for a multiple sequence alignment.
pub struct AlignmentColumn {
  column_index : Int
  counts : Map[Char, Int]
  gap_count : Int
  total_count : Int
}

///|
/// Construct AlignmentColumn.
pub fn AlignmentColumn::new(
  column_index~ : Int,
  counts~ : Map[Char, Int],
  gap_count~ : Int,
  total_count~ : Int,
) -> AlignmentColumn {
  AlignmentColumn::{ column_index, counts, gap_count, total_count }
}

///|
/// Calculate the consensus residue for this column.
pub fn AlignmentColumn::consensus(self : AlignmentColumn) -> Char {
  if self.total_count == 0 {
    return '-'
  }
  let mut max_count = 0
  let mut max_char = 'N'
  for kv in self.counts {
    if kv.1 > max_count {
      max_count = kv.1
      max_char = kv.0
    }
  }
  max_char
}

///|
/// Calculate the majority residue (appears in >= 50% of sequences).
pub fn AlignmentColumn::majority(self : AlignmentColumn) -> Char? {
  if self.total_count == 0 {
    return None
  }
  let threshold = (self.total_count.to_double() / 2.0).ceil().to_int()
  for kv in self.counts {
    if kv.1 >= threshold {
      return Some(kv.0)
    }
  }
  None
}

///|
/// Calculate the fraction of the most common residue.
pub fn AlignmentColumn::consensus_fraction(self : AlignmentColumn) -> Double {
  if self.total_count == 0 {
    return 0.0
  }
  let mut max_count = 0
  for kv in self.counts {
    if kv.1 > max_count {
      max_count = kv.1
    }
  }
  max_count.to_double() / self.total_count.to_double()
}

///|
/// Calculate the Shannon entropy for this column (measure of variability).
pub fn AlignmentColumn::shannon_entropy(self : AlignmentColumn) -> Double {
  if self.total_count == 0 {
    return 0.0
  }
  let mut entropy = 0.0
  let nf = self.total_count.to_double()
  for kv in self.counts {
    let p = kv.1.to_double() / nf
    entropy = entropy - p * @math.log2(p)
  }
  entropy
}

///|
/// Check if this column is conserved (all residues are the same).
pub fn AlignmentColumn::is_conserved(
  self : AlignmentColumn,
  min_conservation? : Double = 1.0,
) -> Bool {
  self.consensus_fraction() >= min_conservation
}

///|
/// Calculate column statistics for a multiple sequence alignment.
/// alignments: Array of sequences (all must be same length).
pub fn calc_column_statistics(
  alignments : Array[String],
) -> Array[AlignmentColumn] {
  if alignments.length() == 0 {
    return []
  }
  let n_cols = alignments[0].length()
  let n_seqs = alignments.length()
  let result : Array[AlignmentColumn] = Array::new()
  let mut col_idx = 0
  while col_idx < n_cols {
    let counts : Map[Char, Int] = Map([], capacity=16)
    let mut gap_count = 0
    let mut row_idx = 0
    while row_idx < n_seqs {
      let c = alignments[row_idx].unsafe_get(col_idx).unsafe_to_char()
      if c == '-' || c == '.' || c == 'N' || c == 'n' {
        gap_count = gap_count + 1
      } else {
        let uc = c.to_ascii_uppercase()
        match counts.get(uc) {
          Some(count) => counts[uc] = count + 1
          None => counts[uc] = 1
        }
      }
      row_idx = row_idx + 1
    }
    result.push(
      AlignmentColumn::new(
        column_index=col_idx,
        counts~,
        gap_count~,
        total_count=n_seqs - gap_count,
      ),
    )
    col_idx = col_idx + 1
  }
  result
}

///|
/// Generate a consensus sequence from alignments.
pub fn alignment_consensus(alignments : Array[String]) -> String {
  let cols = calc_column_statistics(alignments)
  let mut result = ""
  for col in cols {
    result = result + col.consensus().to_string()
  }
  result
}

///|
/// Generate a majority sequence (only positions with >= 50% consensus).
/// Non-majority positions are replaced with ambiguity character.
pub fn majority_sequence(
  alignments : Array[String],
  ambiguity_char? : Char = 'N',
) -> String {
  let cols = calc_column_statistics(alignments)
  let mut result = ""
  for col in cols {
    match col.majority() {
      Some(c) => result = result + c.to_string()
      None => result = result + ambiguity_char.to_string()
    }
  }
  result
}

///|
/// Generate a strict consensus sequence (only 100% conserved positions).
/// Non-conserved positions are replaced with ambiguity character.
pub fn strict_consensus(
  alignments : Array[String],
  ambiguity_char? : Char = 'N',
) -> String {
  let cols = calc_column_statistics(alignments)
  let mut result = ""
  for col in cols {
    if col.consensus_fraction() >= 1.0 {
      result = result + col.consensus().to_string()
    } else {
      result = result + ambiguity_char.to_string()
    }
  }
  result
}

///|
/// Find conserved positions in alignment (columns with high conservation).
/// Returns array of (column_index, consensus_residue, conservation_fraction).
pub fn find_conserved_positions(
  alignments : Array[String],
  min_conservation? : Double = 0.8,
) -> Array[(Int, Char, Double)] {
  let result : Array[(Int, Char, Double)] = Array::new()
  let cols = calc_column_statistics(alignments)
  for col in cols {
    let frac = col.consensus_fraction()
    if frac >= min_conservation && col.gap_count == 0 {
      result.push((col.column_index, col.consensus(), frac))
    }
  }
  result
}

///|
/// Find variable positions in alignment (columns with low conservation).
/// Returns array of (column_index, shannon_entropy).
pub fn find_variable_positions(
  alignments : Array[String],
  max_conservation? : Double = 0.5,
) -> Array[(Int, Double)] {
  let result : Array[(Int, Double)] = Array::new()
  let cols = calc_column_statistics(alignments)
  for col in cols {
    let frac = col.consensus_fraction()
    if frac <= max_conservation && col.gap_count == 0 {
      result.push((col.column_index, col.shannon_entropy()))
    }
  }
  result
}

///|
/// Calculate the conservation profile across the alignment.
/// Returns array of conservation fractions for each column.
pub fn conservation_profile(alignments : Array[String]) -> Array[Double] {
  let cols = calc_column_statistics(alignments)
  cols.map(fn(col) { col.consensus_fraction() })
}

///|
/// Calculate the entropy profile across the alignment.
/// Returns array of Shannon entropy values for each column.
pub fn entropy_profile(alignments : Array[String]) -> Array[Double] {
  let cols = calc_column_statistics(alignments)
  cols.map(fn(col) { col.shannon_entropy() })
}

///|
/// Calculate the percentage identity between all pairs of sequences.
pub fn pairwise_identity(
  alignments : Array[String],
) -> Array[(String, String, Double)] {
  let result : Array[(String, String, Double)] = Array::new()
  let n = alignments.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let identity = sequence_identity(alignments[i], alignments[j])
      result.push(("seq" + i.to_string(), "seq" + j.to_string(), identity))
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Calculate the identity between two aligned sequences.
pub fn sequence_identity(seq1 : String, seq2 : String) -> Double {
  let min_len = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  if min_len == 0 {
    return 0.0
  }
  let mut identical = 0
  let mut total = 0
  let mut i = 0
  while i < min_len {
    let c1 = seq1.unsafe_get(i).unsafe_to_char()
    let c2 = seq2.unsafe_get(i).unsafe_to_char()
    if c1 != '-' && c1 != '.' && c2 != '-' && c2 != '.' {
      total = total + 1
      if c1.to_ascii_uppercase() == c2.to_ascii_uppercase() {
        identical = identical + 1
      }
    }
    i = i + 1
  }
  if total == 0 {
    0.0
  } else {
    identical.to_double() / total.to_double() * 100.0
  }
}

///|
/// Calculate the average pairwise identity across all sequences.
pub fn average_pairwise_identity(alignments : Array[String]) -> Double {
  let n = alignments.length()
  if n < 2 {
    return 0.0
  }
  let pairs = pairwise_identity(alignments)
  let mut total = 0.0
  for p in pairs {
    total = total + p.2
  }
  total / pairs.length().to_double()
}

///|
/// Calculate the number of parsimony-informative sites.
/// A site is parsimony-informative if it has at least 2 different residues,
/// each appearing in at least 2 sequences.
pub fn count_parsimony_informative_sites(alignments : Array[String]) -> Int {
  let cols = calc_column_statistics(alignments)
  let mut count = 0
  for col in cols {
    if col.gap_count == 0 {
      let mut num_residues = 0
      let mut has_two_or_more = false
      for kv in col.counts {
        num_residues = num_residues + 1
        if kv.1 >= 2 {
          has_two_or_more = true
        }
      }
      if num_residues >= 2 && has_two_or_more {
        count = count + 1
      }
    }
  }
  count
}

///|
/// Calculate GC content profile across alignment columns.
pub fn gc_content_profile(alignments : Array[String]) -> Array[Double] {
  let cols = calc_column_statistics(alignments)
  cols.map(fn(col) {
    let g = col.counts.get('G').unwrap_or(0)
    let c = col.counts.get('C').unwrap_or(0)
    let total = col.total_count
    if total == 0 {
      0.0
    } else {
      (g + c).to_double() / total.to_double() * 100.0
    }
  })
}

///|
/// Create a demo multiple sequence alignment.
pub fn create_demo_alignment() -> Array[String] {
  ["ATCGATCGATCG", "ATCGATCGATCG", "ATCGATCGATCG", "ATCGATCGATCG"]
}

///|
/// Create a demo alignment with some variation.
pub fn create_demo_alignment_varied() -> Array[String] {
  [
    "ATCGATCGATCG", "ATCGATCGATCG", "ATCGATGGATCG", "ATCGATCGATCA", "ATCGATCGATCG",
  ]
}
