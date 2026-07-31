///|
/// Bio.motifs.Jaspar - JASPAR transcription factor binding profile parsing.
/// Inspired by Biopython's Bio.motifs.jaspar module.
///
/// JASPAR is a curated database of transcription factor binding profiles.
/// This module parses JASPAR PFM (Position Frequency Matrix) format files
/// and provides utilities for working with them.
///
/// JASPAR PFM format:
///   >MA0001.1 AP-1
///   A [1 0 0 7 0 0 3 0 2 0]
///   C [2 0 1 3 0 5 0 0 1 1]
///   G [6 0 0 5 6 0 0 0 0 1]
///   T [1 8 9 2 4 3 5 9 7 8]
///
/// Each matrix has:
///   - An ID (e.g., MA0001.1)
///   - A name (e.g., AP-1)
///   - Four rows for A, C, G, T (or A, C, G, U for RNA)
///   - Each row contains frequencies at each position

///|
/// A JASPAR position frequency matrix.
pub struct JasparMatrix {
  id : String
  name : String
  pfm : Array[Array[Int]] // 4 rows (A,C,G,T/U) x width columns
  alphabet : Array[String] // ["A", "C", "G", "T"] or ["A", "C", "G", "U"]
  width : Int
  mut total_counts : Int // Total count per position (should be constant)
  mut collection : String // e.g., "CORE", "POLII"
  species : Array[String] // Species taxonomy
  pubmed_ids : Array[String] // PubMed references
  comments : Array[String] // Additional comments
}

///|
/// Create a new JasparMatrix.
pub fn JasparMatrix::new(
  id : String,
  name : String,
  pfm : Array[Array[Int]],
  alphabet? : Array[String] = ["A", "C", "G", "T"],
) -> JasparMatrix {
  let width = if pfm.length() > 0 { pfm[0].length() } else { 0 }
  JasparMatrix::{
    id,
    name,
    pfm,
    alphabet,
    width,
    total_counts: 0,
    collection: "",
    species: Array::new(),
    pubmed_ids: Array::new(),
    comments: Array::new(),
  }
}

///|
/// Get the consensus sequence of the PFM.
pub fn JasparMatrix::consensus(self : JasparMatrix) -> String {
  let result = StringBuilder::new()
  let mut col = 0
  while col < self.width {
    let mut max_count = -1
    let mut max_base = "N"
    let mut row = 0
    while row < self.alphabet.length() {
      let count = self.pfm[row][col]
      if count > max_count {
        max_count = count
        max_base = self.alphabet[row]
      }
      row = row + 1
    }
    result.write_string(max_base)
    col = col + 1
  }
  result.to_string()
}

///|
/// Get the IUPAC consensus sequence (using ambiguity codes).
pub fn JasparMatrix::iupac_consensus(self : JasparMatrix) -> String {
  let result = StringBuilder::new()
  let iupac_map : Array[(Array[String], String)] = [
    (["A"], "A"),
    (["C"], "C"),
    (["G"], "G"),
    (["T"], "T"),
    (["A", "G"], "R"),
    (["C", "T"], "Y"),
    (["G", "C"], "S"),
    (["A", "T"], "W"),
    (["G", "T"], "K"),
    (["A", "C"], "M"),
    (["C", "G", "T"], "B"),
    (["A", "G", "T"], "D"),
    (["A", "C", "T"], "H"),
    (["A", "C", "G"], "V"),
    (["A", "C", "G", "T"], "N"),
  ]
  let mut col = 0
  while col < self.width {
    let bases_present : Array[String] = Array::new()
    let total = get_col_total(self, col)
    let mut row = 0
    while row < self.alphabet.length() {
      let count = self.pfm[row][col]
      if count > total / 5 { // More than 20% threshold
        bases_present.push(self.alphabet[row])
      }
      row = row + 1
    }
    // Find matching IUPAC code
    let mut matched = "N"
    let mut mi = 0
    while mi < iupac_map.length() {
      let entry = iupac_map[mi]
      let bases = entry.0
      let code = entry.1
      if array_set_equals(bases, bases_present) {
        matched = code
        break
      }
      mi = mi + 1
    }
    result.write_string(matched)
    col = col + 1
  }
  result.to_string()
}

///|
/// Calculate total count at a column position.
fn get_col_total(matrix : JasparMatrix, col : Int) -> Int {
  let mut total = 0
  let mut row = 0
  while row < matrix.alphabet.length() {
    total = total + matrix.pfm[row][col]
    row = row + 1
  }
  total
}

///|
/// Check if two string arrays have the same set of elements.
fn array_set_equals(a : Array[String], b : Array[String]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for item in a {
    let found = b.contains(item)
    if !found {
      return false
    }
  }
  true
}

///|
/// Find the index of a substring in a string, returns -1 if not found.
fn j_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n {
    return -1
  }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Convert PFM to a position weight matrix (PWM) with log-odds scoring.
/// Uses background frequency of 0.25 for each base.
pub fn JasparMatrix::to_pwm(self : JasparMatrix, pseudocount? : Double = 0.01) -> Array[Array[Double]] {
  let pwm = Array::new()
  let bg = 0.25
  let mut row = 0
  while row < self.alphabet.length() {
    let col_scores : Array[Double] = Array::new()
    let mut col = 0
    while col < self.width {
      let total = get_col_total(self, col).to_double()
      let freq = (self.pfm[row][col].to_double() + pseudocount) / (total + pseudocount * 4.0)
      let score = @math.ln(freq / bg) / @math.ln(2.0)
      col_scores.push(score)
      col = col + 1
    }
    pwm.push(col_scores)
    row = row + 1
  }
  pwm
}

///|
/// Calculate information content at each position (bits).
pub fn JasparMatrix::information_content(self : JasparMatrix) -> Array[Double] {
  let ic : Array[Double] = Array::new()
  let mut col = 0
  while col < self.width {
    let total = get_col_total(self, col).to_double()
    let mut entropy = 0.0
    let mut row = 0
    while row < self.alphabet.length() {
      let p = self.pfm[row][col].to_double() / total
      if p > 0.0 {
        entropy = entropy - p * @math.ln(p) / @math.ln(2.0)
      }
      row = row + 1
    }
    let max_entropy = @math.ln(self.alphabet.length().to_double()) / @math.ln(2.0)
    ic.push(max_entropy - entropy)
    col = col + 1
  }
  ic
}

///|
/// Calculate total information content of the matrix.
pub fn JasparMatrix::total_information(self : JasparMatrix) -> Double {
  let ic = self.information_content()
  let mut total = 0.0
  for v in ic {
    total = total + v
  }
  total
}

///|
/// Scan a sequence with the PWM and return the best score.
pub fn JasparMatrix::scan_sequence(self : JasparMatrix, sequence : String) -> Array[(Int, Double, String)] {
  let pwm = self.to_pwm()
  let upper_seq = sequence.to_upper()
  let seq_len = upper_seq.length()
  let results : Array[(Int, Double, String)] = Array::new()
  if seq_len < self.width {
    return results
  }
  let mut i = 0
  while i <= seq_len - self.width {
    let mut score = 0.0
    let mut valid = true
    let mut j = 0
    while j < self.width {
      let base = upper_seq.unsafe_get(i + j).unsafe_to_char().to_string()
      let base_idx = find_base_index(self, base)
      if base_idx < 0 {
        valid = false
        break
      }
      score = score + pwm[base_idx][j]
      j = j + 1
    }
    if valid {
      let subseq_view = upper_seq[i:i + self.width]
      let subseq = subseq_view.to_owned()
      let tuple : (Int, Double, String) = (i, score, subseq)
      results.push(tuple)
    }
    i = i + 1
  }
  results
}

///|
/// Find the index of a base in the alphabet.
fn find_base_index(matrix : JasparMatrix, base : String) -> Int {
  let mut i = 0
  while i < matrix.alphabet.length() {
    if matrix.alphabet[i] == base {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Sort scanning results by score (descending).
pub fn jaspar_sort_results(results : Array[(Int, Double, String)]) -> Array[(Int, Double, String)] {
  let sorted = results
  let mut i = 0
  while i < sorted.length() - 1 {
    let mut swapped = false
    let mut j = 0
    while j < sorted.length() - 1 - i {
      if sorted[j].1 < sorted[j + 1].1 {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    if !swapped {
      break
    }
    i = i + 1
  }
  sorted
}

///|
/// Parse a JASPAR PFM format string into an array of JasparMatrix.
/// Handles both single matrix and multi-matrix formats.
pub fn parse_jaspar(content : String) -> Array[JasparMatrix] {
  let matrices : Array[JasparMatrix] = Array::new()
  let lines = content.split("\n").to_array()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if line[0:1] == ">" {
      // Parse header: >ID NAME [extra info]
      let header = line[1:line.length()].trim().to_owned()
      // Split by tab or space
      let parts_views = header.split("\t").to_array()
      let parts : Array[String] = Array::new()
      for pv in parts_views {
        let trimmed = pv.trim().to_owned()
        if trimmed.length() > 0 {
          // Further split by space if needed
          let sub_parts = trimmed.split(" ").to_array()
          for sp in sub_parts {
            let sp_trimmed = sp.trim().to_owned()
            if sp_trimmed.length() > 0 {
              parts.push(sp_trimmed)
            }
          }
        }
      }
      let id = if parts.length() > 0 { parts[0] } else { "" }
      let name = if parts.length() > 1 { parts[1] } else { "" }
      // Parse the 4 rows of the matrix
      let pfm : Array[Array[Int]] = Array::new()
      let bases = ["A", "C", "G", "T"]
      let mut row_idx = 0
      let mut j = i + 1
      while j < lines.length() && row_idx < 4 {
        let row_line = lines[j].trim().to_owned()
        if row_line.length() == 0 {
          j = j + 1
          continue
        }
        // Parse the row: BASE [n1 n2 n3 ...]
        let bracket_start = j_index_of(row_line, "[")
        let bracket_end = j_index_of(row_line, "]")
        if bracket_start < 0 || bracket_end < 0 {
          break
        }
        let num_str = row_line[bracket_start + 1:bracket_end].to_owned()
        let nums = parse_int_array(num_str)
        pfm.push(nums)
        row_idx = row_idx + 1
        j = j + 1
      }
      if pfm.length() == 4 {
        let matrix = JasparMatrix::new(id, name, pfm, alphabet=bases)
        // Calculate total counts (sum of first column = total sequences)
        let mut total = 0
        let mut col = 0
        while col < pfm.length() {
          total = total + pfm[col][0]
          col = col + 1
        }
        matrix.total_counts = total
        // Try to extract collection from header
        if parts.length() > 2 {
          matrix.collection = parts[2]
        }
        matrices.push(matrix)
      }
      i = j
    } else {
      i = i + 1
    }
  }
  matrices
}

///|
/// Parse a space-separated string of integers into an array.
fn parse_int_array(s : String) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let parts = s.split(" ").to_array()
  for p in parts {
    let trimmed = p.trim().to_owned()
    if trimmed.length() > 0 {
      // Parse the integer
      let mut val = 0
      let mut neg = false
      let mut k = 0
      if trimmed[k] == '-' {
        neg = true
        k = k + 1
      }
      while k < trimmed.length() {
        let digit = trimmed[k].to_int() - '0'.to_int()
        if digit >= 0 && digit <= 9 {
          val = val * 10 + digit
        }
        k = k + 1
      }
      result.push(if neg { -val } else { val })
    }
  }
  result
}

///|
/// Get a specific matrix by ID from an array of matrices.
pub fn jaspar_get_matrix(matrices : Array[JasparMatrix], id : String) -> JasparMatrix {
  for m in matrices {
    if m.id == id {
      return m
    }
  }
  // Return a default empty matrix
  let empty_pfm : Array[Array[Int]] = Array::new()
  JasparMatrix::new("", "", empty_pfm)
}

///|
/// Filter matrices by collection.
pub fn jaspar_filter_by_collection(matrices : Array[JasparMatrix], collection : String) -> Array[JasparMatrix] {
  let result : Array[JasparMatrix] = Array::new()
  for m in matrices {
    if m.collection == collection {
      result.push(m)
    }
  }
  result
}

///|
/// Filter matrices by name pattern (case-insensitive substring match).
pub fn jaspar_filter_by_name(matrices : Array[JasparMatrix], pattern : String) -> Array[JasparMatrix] {
  let result : Array[JasparMatrix] = Array::new()
  let upper_pattern = pattern.to_upper()
  for m in matrices {
    let upper_name = m.name.to_upper()
    if j_index_of(upper_name, upper_pattern) >= 0 {
      result.push(m)
    }
  }
  result
}

///|
/// Generate a sample JASPAR PFM format string for testing.
pub fn jaspar_sample_pfm() -> String {
  ">MA0001.1 AP-1\n" +
  "A [1 0 0 7 0 0 3 0 2 0]\n" +
  "C [2 0 1 3 0 5 0 0 1 1]\n" +
  "G [6 0 0 5 6 0 0 0 0 1]\n" +
  "T [1 8 9 2 4 3 5 9 7 8]\n" +
  ">MA0002.1 GATA1\n" +
  "A [8 2 1 5 0 1 9 1 0 2]\n" +
  "C [0 1 1 2 0 0 1 0 0 1]\n" +
  "G [0 1 1 5 9 2 0 0 0 0]\n" +
  "T [2 8 7 2 10 7 0 9 10 7]"
}

///|
/// Convert a JasparMatrix to a display string.
pub fn JasparMatrix::to_string(self : JasparMatrix) -> String {
  let result = StringBuilder::new()
  result.write_string(">" + self.id + " " + self.name + "\n")
  let mut row = 0
  while row < self.alphabet.length() {
    result.write_string(self.alphabet[row] + " [")
    let mut col = 0
    while col < self.width {
      if col > 0 {
        result.write_string(" ")
      }
      result.write_string(self.pfm[row][col].to_string())
      col = col + 1
    }
    result.write_string("]\n")
    row = row + 1
  }
  result.to_string()
}
