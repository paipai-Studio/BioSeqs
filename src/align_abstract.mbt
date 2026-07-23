///|
/// Bio.Align.AlignAbstract - Abstract alignment types and validation.
/// Inspired by Biopython's Bio.Align.AlignAbstract module.
///
/// Provides abstract base classes for alignment objects, including
/// validation, comparison, and advanced alignment statistics.
/// Defines the core interface that all alignment types must implement.

///|
/// Type of alignment.
pub enum AlignAbstractType {
  /// Nucleotide alignment (DNA/RNA)
  Nucleotide
  /// Protein alignment (amino acids)
  Protein
  /// Generic/unknown alignment
  Generic
}

///|
/// Convert string to alignment type.
pub fn align_abstract_type_from_str(s : String) -> AlignAbstractType {
  match s.to_lower().trim() {
    "nucleotide" | "dna" | "rna" | "nt" => AlignAbstractType::Nucleotide
    "protein" | "amino" | "prot" | "aa" => AlignAbstractType::Protein
    _ => AlignAbstractType::Generic
  }
}

///|
/// Abstract alignment interface - defines core properties of any alignment.
pub struct AbstractAlignment {
  /// Aligned sequences
  sequences : Array[String]
  /// Sequence identifiers
  identifiers : Array[String]
  /// Type of alignment
  alignment_type : AlignAbstractType
  /// Whether alignment has been validated
  validated : Bool
  /// Number of sequences
  n_seqs : Int
  /// Alignment length
  alignment_length : Int
}

///|
/// Construct an abstract alignment from aligned sequences.
pub fn AbstractAlignment::new(
  sequences~ : Array[String],
  identifiers~ : Array[String],
  alignment_type~ : AlignAbstractType,
) -> AbstractAlignment {
  let n_seqs = sequences.length()
  let alignment_length = if n_seqs > 0 { sequences[0].length() } else { 0 }
  AbstractAlignment::{
    sequences, identifiers, alignment_type,
    validated: false, n_seqs, alignment_length
  }
}

///|
/// Validate the alignment: all sequences must have the same length.
pub fn AbstractAlignment::validate(self : AbstractAlignment) -> (Bool, String) {
  if self.n_seqs == 0 {
    return (false, "No sequences in alignment")
  }
  let first_len = self.alignment_length
  let mut i = 1
  while i < self.n_seqs {
    if self.sequences[i].length() != first_len {
      let id = self.identifiers[i]
      let len = self.sequences[i].length().to_string()
      let flen = first_len.to_string()
      let msg = "Sequence " + id + " has length " + len + ", expected " + flen
      return (false, msg)
    }
    i = i + 1
  }

  // Validate characters based on alignment type
  let valid = self.validate_characters()
  match valid {
    (true, _) => (true, "Alignment is valid")
    (false, msg) => (false, msg)
  }
}

///|
/// Validate characters in the alignment.
pub fn AbstractAlignment::validate_characters(self : AbstractAlignment) -> (Bool, String) {
  let valid_chars = get_valid_chars(self.alignment_type)

  let mut seq_idx = 0
  while seq_idx < self.n_seqs {
    let seq = self.sequences[seq_idx]
    let mut i = 0
    while i < seq.length() {
      let c = char_to_upper(seq.unsafe_get(i))
      if !char_in_array(c, valid_chars) {
        let ch = c.to_string()
        let pos = i.to_string()
        let id = self.identifiers[seq_idx]
        let msg = "Invalid character " + ch + " at position " + pos + " in sequence " + id
        return (false, msg)
      }
      i = i + 1
    }
    seq_idx = seq_idx + 1
  }
  (true, "")
}

fn get_valid_chars(align_type : AlignAbstractType) -> Array[UInt16] {
  let nucleotide_chars : Array[UInt16] = [
    65, 84, 71, 67, 85, 78, 45, 46,
    82, 89, 83, 87, 75, 77, 66, 86, 68, 72
  ]
  let protein_chars : Array[UInt16] = [
    65, 82, 78, 68, 67, 81, 69, 71, 72, 73,
    76, 75, 77, 70, 80, 83, 84, 87, 89, 86,
    66, 90, 88, 45, 46, 85, 79
  ]
  let generic_chars : Array[UInt16] = [
    97, 98, 99, 100, 101, 102, 103, 104, 105, 106,
    107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
    117, 118, 119, 120, 121, 122,
    65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
    75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
    85, 86, 87, 88, 89, 90,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    45, 46, 42
  ]

  match align_type {
    AlignAbstractType::Nucleotide => nucleotide_chars
    AlignAbstractType::Protein => protein_chars
    AlignAbstractType::Generic => generic_chars
  }
}

fn char_to_upper(c : UInt16) -> UInt16 {
  if c >= 97 && c <= 122 { c - 32 } else { c }
}

fn char_to_lower(c : UInt16) -> UInt16 {
  if c >= 65 && c <= 90 { c + 32 } else { c }
}

/// Convert a UInt16 character code to its 1-character String representation.
fn u16_to_str(c : UInt16) -> String {
  match c {
    65 => "A"
    66 => "B"
    67 => "C"
    68 => "D"
    69 => "E"
    70 => "F"
    71 => "G"
    72 => "H"
    73 => "I"
    74 => "J"
    75 => "K"
    76 => "L"
    77 => "M"
    78 => "N"
    79 => "O"
    80 => "P"
    81 => "Q"
    82 => "R"
    83 => "S"
    84 => "T"
    85 => "U"
    86 => "V"
    87 => "W"
    88 => "X"
    89 => "Y"
    90 => "Z"
    97 => "a"
    98 => "b"
    99 => "c"
    100 => "d"
    101 => "e"
    102 => "f"
    103 => "g"
    104 => "h"
    105 => "i"
    106 => "j"
    107 => "k"
    108 => "l"
    109 => "m"
    110 => "n"
    111 => "o"
    112 => "p"
    113 => "q"
    114 => "r"
    115 => "s"
    116 => "t"
    117 => "u"
    118 => "v"
    119 => "w"
    120 => "x"
    121 => "y"
    122 => "z"
    45 => "-"
    46 => "."
    42 => "*"
    43 => "+"
    47 => "/"
    48 => "0"
    49 => "1"
    50 => "2"
    51 => "3"
    52 => "4"
    53 => "5"
    54 => "6"
    55 => "7"
    56 => "8"
    57 => "9"
    _ => "?"
  }
}

fn char_in_array(c : UInt16, arr : Array[UInt16]) -> Bool {
  for ch in arr {
    if ch == c {
      return true
    }
  }
  false
}

///|
/// Get sequence count.
pub fn AbstractAlignment::abstract_n_seqs(self : AbstractAlignment) -> Int {
  self.n_seqs
}

///|
/// Get alignment length.
pub fn AbstractAlignment::abstract_alignment_length(self : AbstractAlignment) -> Int {
  self.alignment_length
}

///|
/// Get a sequence by index.
pub fn AbstractAlignment::abstract_get_seq(self : AbstractAlignment, idx : Int) -> String {
  if idx >= 0 && idx < self.n_seqs {
    self.sequences[idx]
  } else {
    ""
  }
}

///|
/// Get sequence identifier by index.
pub fn AbstractAlignment::abstract_get_id(self : AbstractAlignment, idx : Int) -> String {
  if idx >= 0 && idx < self.n_seqs {
    self.identifiers[idx]
  } else {
    ""
  }
}

///|
/// Get a column from the alignment (array of characters).
pub fn AbstractAlignment::abstract_get_column(self : AbstractAlignment, col_idx : Int) -> Array[UInt16] {
  let column : Array[UInt16] = Array::new()
  if col_idx >= 0 && col_idx < self.alignment_length {
    for seq in self.sequences {
      column.push(seq.unsafe_get(col_idx))
    }
  }
  column
}

///|
/// Get the alignment type.
pub fn AbstractAlignment::abstract_type(self : AbstractAlignment) -> AlignAbstractType {
  self.alignment_type
}

// --- Alignment Properties ---

///|
/// Column statistics for an alignment position.
pub struct AlignAbstractColumnStats {
  column_index : Int
  conservation : Double
  diversity : Double
  gap_fraction : Double
  n_unique_chars : Int
  consensus_char : UInt16
}

///|
/// Compute column statistics for a given alignment column.
pub fn abstract_column_stats(alignment : AbstractAlignment, col_idx : Int) -> AlignAbstractColumnStats {
  let column = alignment.abstract_get_column(col_idx)
  let n = column.length()
  if n == 0 {
    return AlignAbstractColumnStats::{
      column_index: col_idx, conservation: 0.0, diversity: 0.0,
      gap_fraction: 0.0, n_unique_chars: 0, consensus_char: 45
    }
  }

  // Count character frequencies
  let counts : Map[UInt16, Int] = Map([], capacity=20)
  let mut gap_count = 0
  for c in column {
    let upper = char_to_upper(c)
    if upper == 45 || upper == 46 {
      gap_count = gap_count + 1
    } else {
      match counts.get(upper) {
        Some(count) => counts[upper] = count + 1
        None => counts[upper] = 1
      }
    }
  }

  let non_gap = n - gap_count
  let gap_fraction = gap_count.to_double() / n.to_double()

  // Find consensus character
  let mut max_count = 0
  let mut consensus : UInt16 = 45
  for kv in counts {
    if kv.1 > max_count {
      max_count = kv.1
      consensus = kv.0
    }
  }

  // Conservation = most common non-gap fraction
  let conservation = if non_gap > 0 {
    max_count.to_double() / non_gap.to_double()
  } else {
    0.0
  }

  // Shannon entropy (diversity)
  let mut entropy = 0.0
  if non_gap > 0 {
    for kv in counts {
      let p = kv.1.to_double() / non_gap.to_double()
      if p > 0.0 {
        entropy = entropy - p * @math.log2(p)
      }
    }
  }

  AlignAbstractColumnStats::{
    column_index: col_idx,
    conservation,
    diversity: entropy,
    gap_fraction,
    n_unique_chars: counts.size(),
    consensus_char: consensus
  }
}

///|
/// Get an array of all column statistics for the alignment.
pub fn abstract_all_column_stats(alignment : AbstractAlignment) -> Array[AlignAbstractColumnStats] {
  let stats : Array[AlignAbstractColumnStats] = Array::new()
  let mut i = 0
  while i < alignment.alignment_length {
    stats.push(abstract_column_stats(alignment, i))
    i = i + 1
  }
  stats
}

///|
/// Get the consensus sequence of the alignment (using most common character per column).
pub fn abstract_consensus_sequence(alignment : AbstractAlignment, threshold? : Double = 0.5) -> String {
  let mut consensus = ""
  let mut i = 0
  while i < alignment.alignment_length {
    let stats = abstract_column_stats(alignment, i)
    if stats.conservation >= threshold {
      consensus = consensus + u16_to_str(stats.consensus_char)
    } else {
      consensus = consensus + "N"
    }
    i = i + 1
  }
  consensus
}

///|
/// Calculate sequence identity matrix (pairwise).
pub fn abstract_identity_matrix(alignment : AbstractAlignment) -> Array[Array[Double]] {
  let n = alignment.n_seqs
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let identity = seq_identity(alignment.sequences[i], alignment.sequences[j])
      row.push(identity)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

fn seq_identity(seq1 : String, seq2 : String) -> Double {
  if seq1.length() != seq2.length() {
    return 0.0
  }
  let n = seq1.length()
  if n == 0 {
    return 0.0
  }
  let mut matches = 0
  let mut i = 0
  while i < n {
    let c1 = char_to_upper(seq1.unsafe_get(i))
    let c2 = char_to_upper(seq2.unsafe_get(i))
    if c1 == c2 {
      matches = matches + 1
    }
    i = i + 1
  }
  matches.to_double() / n.to_double()
}

///|
/// Calculate the overall alignment identity (average of all pairwise identities).
pub fn abstract_overall_identity(alignment : AbstractAlignment) -> Double {
  let matrix = abstract_identity_matrix(alignment)
  let n = matrix.length()
  if n <= 1 {
    return 1.0
  }
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      sum = sum + matrix[i][j]
      count = count + 1
      j = j + 1
    }
    i = i + 1
  }
  if count > 0 {
    sum / count.to_double()
  } else {
    1.0
  }
}

///|
/// Calculate alignment coverage (fraction of non-gap positions).
pub fn abstract_coverage(alignment : AbstractAlignment) -> Double {
  let total = alignment.n_seqs * alignment.alignment_length
  if total == 0 {
    return 0.0
  }
  let mut gaps = 0
  for seq in alignment.sequences {
    let mut i = 0
    while i < seq.length() {
      let c = seq.unsafe_get(i)
      if c == 45 || c == 46 {
        gaps = gaps + 1
      }
      i = i + 1
    }
  }
  (total - gaps).to_double() / total.to_double()
}

///|
/// Calculate the number of parsimony-informative sites.
pub fn abstract_parsimony_sites(alignment : AbstractAlignment) -> Int {
  let mut count = 0
  let mut col = 0
  while col < alignment.alignment_length {
    if is_parsimony_informative(alignment, col) {
      count = count + 1
    }
    col = col + 1
  }
  count
}

fn is_parsimony_informative(alignment : AbstractAlignment, col : Int) -> Bool {
  let column = alignment.abstract_get_column(col)
  let counts : Map[UInt16, Int] = Map([], capacity=20)
  for c in column {
    let upper = char_to_upper(c)
    if upper != 45 && upper != 46 {
      match counts.get(upper) {
        Some(count) => counts[upper] = count + 1
        None => counts[upper] = 1
      }
    }
  }

  let mut informative_states = 0
  for kv in counts {
    if kv.1 >= 2 {
      informative_states = informative_states + 1
    }
  }
  informative_states >= 2
}

///|
/// Calculate the number of variable sites.
pub fn abstract_variable_sites(alignment : AbstractAlignment) -> Int {
  let mut count = 0
  let mut col = 0
  while col < alignment.alignment_length {
    if is_variable_site(alignment, col) {
      count = count + 1
    }
    col = col + 1
  }
  count
}

fn is_variable_site(alignment : AbstractAlignment, col : Int) -> Bool {
  let column = alignment.abstract_get_column(col)
  if column.length() == 0 {
    return false
  }
  let first = char_to_upper(column[0])
  for c in column {
    if char_to_upper(c) != first {
      return true
    }
  }
  false
}

///|
/// Calculate the number of singleton sites.
pub fn abstract_singleton_sites(alignment : AbstractAlignment) -> Int {
  let mut count = 0
  let mut col = 0
  while col < alignment.alignment_length {
    if is_singleton_site(alignment, col) {
      count = count + 1
    }
    col = col + 1
  }
  count
}

fn is_singleton_site(alignment : AbstractAlignment, col : Int) -> Bool {
  let column = alignment.abstract_get_column(col)
  let counts : Map[UInt16, Int] = Map([], capacity=20)
  for c in column {
    let upper = char_to_upper(c)
    if upper != 45 && upper != 46 {
      match counts.get(upper) {
        Some(count) => counts[upper] = count + 1
        None => counts[upper] = 1
      }
    }
  }

  if counts.size() <= 1 {
    return false
  }

  let mut states_with_one = 0
  let mut states_with_many = 0
  for kv in counts {
    if kv.1 == 1 {
      states_with_one = states_with_one + 1
    } else if kv.1 > 1 {
      states_with_many = states_with_many + 1
    }
  }
  states_with_one >= 1 && states_with_many <= 1
}

///|
/// Calculate alignment distance matrix (using simple p-distance).
pub fn abstract_distance_matrix(alignment : AbstractAlignment) -> Array[Array[Double]] {
  let n = alignment.n_seqs
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let dist = pairwise_distance(alignment.sequences[i], alignment.sequences[j])
      row.push(dist)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

fn pairwise_distance(seq1 : String, seq2 : String) -> Double {
  if seq1.length() != seq2.length() {
    return 1.0
  }
  let n = seq1.length()
  if n == 0 {
    return 0.0
  }
  let mut diffs = 0
  let mut i = 0
  while i < n {
    let c1 = char_to_upper(seq1.unsafe_get(i))
    let c2 = char_to_upper(seq2.unsafe_get(i))
    if c1 != c2 {
      diffs = diffs + 1
    }
    i = i + 1
  }
  diffs.to_double() / n.to_double()
}

///|
/// Filter alignment columns by gap fraction.
pub fn abstract_filter_gaps(alignment : AbstractAlignment, max_gap_fraction : Double) -> AbstractAlignment {
  let keep_columns : Array[Int] = Array::new()
  let mut col = 0
  while col < alignment.alignment_length {
    let stats = abstract_column_stats(alignment, col)
    if stats.gap_fraction <= max_gap_fraction {
      keep_columns.push(col)
    }
    col = col + 1
  }

  let new_seqs : Array[String] = Array::new()
  let mut seq_idx = 0
  while seq_idx < alignment.n_seqs {
    let seq = alignment.sequences[seq_idx]
    let mut new_seq = ""
    let mut k_idx = 0
    while k_idx < keep_columns.length() {
      let k = keep_columns[k_idx]
      new_seq = new_seq + u16_to_str(seq.unsafe_get(k))
      k_idx = k_idx + 1
    }
    new_seqs.push(new_seq)
    seq_idx = seq_idx + 1
  }

  AbstractAlignment::new(sequences=new_seqs, identifiers=alignment.identifiers.copy(), alignment_type=alignment.alignment_type)
}

///|
/// Filter alignment sequences by minimum coverage.
pub fn abstract_filter_coverage(alignment : AbstractAlignment, min_coverage : Double) -> AbstractAlignment {
  let keep_seqs : Array[Int] = Array::new()
  let mut idx = 0
  while idx < alignment.n_seqs {
    let seq = alignment.sequences[idx]
    let mut non_gap = 0
    let mut i = 0
    while i < seq.length() {
      let c = seq.unsafe_get(i)
      if c != 45 && c != 46 {
        non_gap = non_gap + 1
      }
      i = i + 1
    }
    let coverage = non_gap.to_double() / seq.length().to_double()
    if coverage >= min_coverage {
      keep_seqs.push(idx)
    }
    idx = idx + 1
  }

  let new_seqs : Array[String] = Array::new()
  let new_ids : Array[String] = Array::new()
  for idx2 in keep_seqs {
    new_seqs.push(alignment.sequences[idx2])
    new_ids.push(alignment.identifiers[idx2])
  }

  AbstractAlignment::new(sequences=new_seqs, identifiers=new_ids, alignment_type=alignment.alignment_type)
}

///|
/// Trim alignment to include only the region between start and end columns.
pub fn abstract_trim(alignment : AbstractAlignment, start : Int, end : Int) -> AbstractAlignment {
  let s = if start > 0 { start } else { 0 }
  let e = if end < alignment.alignment_length - 1 { end } else { alignment.alignment_length - 1 }

  let new_seqs : Array[String] = Array::new()
  for seq in alignment.sequences {
    let trimmed = seq.substring(start=s, end=e + 1)
    new_seqs.push(trimmed)
  }

  AbstractAlignment::new(sequences=new_seqs, identifiers=alignment.identifiers.copy(), alignment_type=alignment.alignment_type)
}

///|
/// Generate a summary report for the alignment.
pub fn abstract_summary(alignment : AbstractAlignment) -> String {
  let (valid, msg) = alignment.validate()
  let n_seqs = alignment.n_seqs.to_string()
  let alen = alignment.alignment_length.to_string()
  let coverage = abstract_coverage(alignment).to_string()
  let identity = abstract_overall_identity(alignment).to_string()
  let parsimony = abstract_parsimony_sites(alignment).to_string()
  let variable = abstract_variable_sites(alignment).to_string()
  let vstr = if valid { "true" } else { "false" }
  let tstr = alignment_type_str(alignment.alignment_type)

  "Alignment Summary:\n" +
  "  Sequences: " + n_seqs + "\n" +
  "  Length: " + alen + "\n" +
  "  Type: " + tstr + "\n" +
  "  Valid: " + vstr + "\n" +
  "  Coverage: " + coverage + "\n" +
  "  Overall identity: " + identity + "\n" +
  "  Variable sites: " + variable + "\n" +
  "  Parsimony-informative sites: " + parsimony + "\n"
}

pub fn alignment_type_str(t : AlignAbstractType) -> String {
  match t {
    AlignAbstractType::Nucleotide => "Nucleotide"
    AlignAbstractType::Protein => "Protein"
    AlignAbstractType::Generic => "Generic"
  }
}
