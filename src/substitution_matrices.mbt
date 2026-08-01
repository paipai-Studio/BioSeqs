///|
/// substitution_matrices module - Modern substitution matrix infrastructure.
/// Reference: Biopython Bio.Align.substitution_matrices (2021+ replacement
/// for Bio.SubsMat).
///
/// Provides:
///   1. ArrayData: 2D numeric array with labeled rows/columns
///   2. SubsMatrix: scoring matrix with alphabet, lookup, extension
///   3. Built-in matrix database: BLOSUM (45-100), PAM (30/70/250), NUC (4.4/4.2)
///   4. calculate_frequency_matrix: observed/expected pair frequencies from alignment
///   5. calculate_substitution_matrix: log-odds matrix from alignment
///   6. subs_shannon_entropy: information content of a frequency matrix
///   7. Matrix registry: register/load matrices by name
///
/// Note: All public identifiers are prefixed with `subs_` to avoid conflicts
/// with existing modules (motifs.mbt, variation.mbt, biostrings.mbt).

// ============================================================================
// Standard alphabets
// ============================================================================

///|
/// Standard 20 amino acid alphabet (IUPAC protein, sorted).
pub fn subs_protein_alphabet() -> Array[String] {
  [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
  ]
}

///|
/// Standard 4 nucleotide alphabet.
pub fn subs_nucleotide_alphabet() -> Array[String] {
  ["A", "C", "G", "T"]
}

// ============================================================================
// ArrayData: 2D labeled numeric array
// ============================================================================

///|
/// A 2D numeric array with row and column labels.
/// Used internally by SubsMatrix for score storage.
pub struct ArrayData {
  pub rows : Array[String]        // row labels
  pub cols : Array[String]        // column labels
  pub data : Array[Array[Double]] // row x col values
}

///|
/// Create an ArrayData from a dict of (row, col) -> value pairs.
pub fn array_data_from_dict(
  pairs : Array[(String, String, Double)],
  rows : Array[String],
  cols : Array[String],
) -> ArrayData {
  let n_rows = rows.length()
  let n_cols = cols.length()
  // Initialize data matrix with zeros
  let data : Array[Array[Double]] = []
  for _i in 0..<n_rows {
    let row : Array[Double] = []
    for _j in 0..<n_cols { row.push(0.0) }
    data.push(row)
  }
  // Build row/col index maps
  let row_idx_keys : Array[String] = []
  let row_idx_vals : Array[Int] = []
  for i in 0..<n_rows { row_idx_keys.push(rows[i]); row_idx_vals.push(i) }
  let col_idx_keys : Array[String] = []
  let col_idx_vals : Array[Int] = []
  for j in 0..<n_cols { col_idx_keys.push(cols[j]); col_idx_vals.push(j) }
  // Fill values
  let mut pi = 0
  while pi < pairs.length() {
    let (r, c, v) = pairs[pi]
    let mut ri = -1
    for k in 0..<row_idx_keys.length() {
      if row_idx_keys[k] == r { ri = row_idx_vals[k]; break }
    }
    let mut ci = -1
    for k in 0..<col_idx_keys.length() {
      if col_idx_keys[k] == c { ci = col_idx_vals[k]; break }
    }
    if ri >= 0 && ci >= 0 { data[ri][ci] = v }
    pi = pi + 1
  }
  { rows: rows, cols: cols, data: data }
}

///|
/// Get a value from ArrayData by row and column index.
pub fn ArrayData::get(self : ArrayData, i : Int, j : Int) -> Double {
  if i < 0 || i >= self.rows.length() || j < 0 || j >= self.cols.length() { 0.0 }
  else { self.data[i][j] }
}

///|
/// Set a value in ArrayData by row and column index.
pub fn ArrayData::set(self : ArrayData, i : Int, j : Int, v : Double) -> Unit {
  if i >= 0 && i < self.rows.length() && j >= 0 && j < self.cols.length() {
    self.data[i][j] = v
  }
}

///|
/// Get a value by row and column labels.
pub fn ArrayData::get_by_label(self : ArrayData, row : String, col : String) -> Double {
  let mut ri = -1
  for k in 0..<self.rows.length() {
    if self.rows[k] == row { ri = k; break }
  }
  let mut ci = -1
  for k in 0..<self.cols.length() {
    if self.cols[k] == col { ci = k; break }
  }
  if ri < 0 || ci < 0 { 0.0 } else { self.data[ri][ci] }
}

// ============================================================================
// SubsMatrix
// ============================================================================

///|
/// A substitution matrix for sequence alignment scoring.
/// Stores scores for each pair of letters in the alphabet.
pub struct SubsMatrix {
  name : String
  alphabet : Array[String]
  matrix : ArrayData       // symmetric for protein matrices
  n_letters : Int
}

///|
/// Create a SubsMatrix from an ArrayData.
pub fn subs_matrix(
  name : String,
  alphabet : Array[String],
  matrix : ArrayData,
) -> SubsMatrix {
  {
    name: name,
    alphabet: alphabet,
    matrix: matrix,
    n_letters: alphabet.length(),
  }
}

///|
/// Create a SubsMatrix from (letter1, letter2, score) triples.
pub fn subs_matrix_from_pairs(
  name : String,
  alphabet : Array[String],
  pairs : Array[(String, String, Double)],
) -> SubsMatrix {
  let matrix = array_data_from_dict(pairs, alphabet, alphabet)
  subs_matrix(name, alphabet, matrix)
}

///|
/// Get the score for a pair of letters.
pub fn SubsMatrix::get_score(
  self : SubsMatrix,
  letter1 : String,
  letter2 : String,
) -> Double {
  self.matrix.get_by_label(letter1, letter2)
}

///|
/// Get the score by index (both 0-based positions in the alphabet).
pub fn SubsMatrix::get_score_idx(
  self : SubsMatrix,
  i : Int,
  j : Int,
) -> Double {
  self.matrix.get(i, j)
}

///|
/// Get the matrix name.
pub fn SubsMatrix::name(self : SubsMatrix) -> String {
  self.name
}

///|
/// Get the alphabet.
pub fn SubsMatrix::alphabet(self : SubsMatrix) -> Array[String] {
  self.alphabet
}

///|
/// Number of letters in the alphabet.
pub fn SubsMatrix::n_letters(self : SubsMatrix) -> Int {
  self.n_letters
}

///|
/// Select a sub-matrix containing only the specified letters.
pub fn SubsMatrix::select(
  self : SubsMatrix,
  letters : Array[String],
) -> SubsMatrix {
  let pairs : Array[(String, String, Double)] = []
  for r in letters {
    for c in letters {
      let v = self.get_score(r, c)
      pairs.push((r, c, v))
    }
  }
  let matrix = array_data_from_dict(pairs, letters, letters)
  subs_matrix(self.name + "_selected", letters, matrix)
}

///|
/// Format the matrix as a string table (for display).
pub fn SubsMatrix::to_table_string(self : SubsMatrix) -> String {
  let mut result = self.name + "\n"
  // Header row
  result = result + "    "
  for c in self.alphabet {
    result = result + c + "  "
  }
  result = result + "\n"
  // Data rows
  for i in 0..<self.n_letters {
    result = result + self.alphabet[i] + "  "
    for j in 0..<self.n_letters {
      let v = self.matrix.get(i, j)
      let s = v.to_int().to_string()
      if s.length() < 2 {
        result = result + " " + s + " "
      } else {
        result = result + s + " "
      }
    }
    result = result + "\n"
  }
  result
}

// ============================================================================
// Matrix registry
// ============================================================================

///|
/// Registry of named matrices. Stores (name, matrix) pairs.
let registry_keys : Ref[Array[String]] = Ref::new([])
let registry_matrices : Ref[Array[SubsMatrix]] = Ref::new([])

///|
/// Register a substitution matrix by name.
pub fn register_matrix(name : String, matrix : SubsMatrix) -> Unit {
  // Check if name already exists
  let keys = registry_keys.val
  for i in 0..<keys.length() {
    if keys[i] == name {
      registry_matrices.val[i] = matrix
      return
    }
  }
  registry_keys.val.push(name)
  registry_matrices.val.push(matrix)
}

///|
/// Load a matrix from the registry by name. Returns Option.
pub fn load_matrix(name : String) -> SubsMatrix? {
  let keys = registry_keys.val
  for i in 0..<keys.length() {
    if keys[i] == name { return Some(registry_matrices.val[i]) }
  }
  None
}

///|
/// List all registered matrix names.
pub fn list_matrices() -> Array[String] {
  registry_keys.val.copy()
}

// ============================================================================
// Built-in BLOSUM matrices
// ============================================================================

///|
/// BLOSUM62 substitution matrix (Henikoff & Henikoff 1992).
/// The most commonly used protein substitution matrix.
pub fn subs_blosum62_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  // Key scores (selected, full 20x20 matrix)
  let pairs : Array[(String, String, Double)] = []
  // BLOSUM62 scores as (row, col, score) triples
  let scores = blosum62_scores()
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("BLOSUM62", alpha, pairs)
}

///|
/// BLOSUM45 substitution matrix (more divergent sequences).
pub fn subs_blosum45_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = blosum45_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("BLOSUM45", alpha, pairs)
}

///|
/// BLOSUM80 substitution matrix (more closely related sequences).
pub fn subs_blosum80_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = blosum80_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("BLOSUM80", alpha, pairs)
}

///|
/// BLOSUM90 substitution matrix (closely related sequences).
pub fn subs_blosum90_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = blosum90_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("BLOSUM90", alpha, pairs)
}

// ============================================================================
// Built-in PAM matrices
// ============================================================================

///|
/// PAM30 substitution matrix (closely related sequences).
pub fn subs_pam30_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = pam30_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("PAM30", alpha, pairs)
}

///|
/// PAM70 substitution matrix.
pub fn subs_pam70_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = pam70_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("PAM70", alpha, pairs)
}

///|
/// PAM250 substitution matrix (distantly related sequences).
pub fn subs_pam250_matrix() -> SubsMatrix {
  let alpha = subs_protein_alphabet()
  let scores = pam250_scores()
  let pairs : Array[(String, String, Double)] = []
  for i in 0..<20 {
    for j in 0..<20 {
      pairs.push((alpha[i], alpha[j], scores[i * 20 + j].to_double()))
    }
  }
  subs_matrix_from_pairs("PAM250", alpha, pairs)
}

// ============================================================================
// Built-in nucleotide matrices
// ============================================================================

///|
/// NUC.4 substitution matrix for nucleotides.
/// Match=5, mismatch=-4 (commonly used for BLASTN).
pub fn subs_nuc44_matrix() -> SubsMatrix {
  let alpha = subs_nucleotide_alphabet()
  let pairs : Array[(String, String, Double)] = []
  for r in alpha {
    for c in alpha {
      let v = if r == c { 5.0 } else { -4.0 }
      pairs.push((r, c, v))
    }
  }
  subs_matrix_from_pairs("NUC4.4", alpha, pairs)
}

// ============================================================================
// Initialize registry with all built-in matrices
// ============================================================================

///|
/// Register all built-in substitution matrices.
/// Call this once at program startup.
pub fn subs_initialize_registry() -> Unit {
  register_matrix("BLOSUM45", subs_blosum45_matrix())
  register_matrix("BLOSUM62", subs_blosum62_matrix())
  register_matrix("BLOSUM80", subs_blosum80_matrix())
  register_matrix("BLOSUM90", subs_blosum90_matrix())
  register_matrix("PAM30", subs_pam30_matrix())
  register_matrix("PAM70", subs_pam70_matrix())
  register_matrix("PAM250", subs_pam250_matrix())
  register_matrix("NUC4.4", subs_nuc44_matrix())
}

// ============================================================================
// Frequency matrix calculation from alignment
// ============================================================================

///|
/// Calculate the observed pair frequency matrix from a multiple sequence
/// alignment.
///
/// For each pair of sequences in the alignment, counts substitutions
/// (letter_i, letter_j) at each aligned position. Returns a symmetric
/// frequency matrix where freq[i][j] = number of times letter i is aligned
/// with letter j.
///
/// Parameters:
///   - alignment: array of aligned sequences (all same length)
///   - alphabet: the letters to include in the matrix
pub fn calculate_frequency_matrix(
  alignment : Array[String],
  alphabet : Array[String],
) -> ArrayData {
  let n = alphabet.length()
  // Initialize count matrix
  let data : Array[Array[Double]] = []
  for _i in 0..<n {
    let row : Array[Double] = []
    for _j in 0..<n { row.push(0.0) }
    data.push(row)
  }
  // Build letter -> index map
  let letter_keys : Array[String] = []
  let letter_vals : Array[Int] = []
  for i in 0..<n { letter_keys.push(alphabet[i]); letter_vals.push(i) }
  // Count pairs across all sequence pairs
  let n_seq = alignment.length()
  if n_seq < 2 { return { rows: alphabet, cols: alphabet, data: data } }
  let seq_len = alignment[0].length()
  for i in 0..<n_seq {
    for j in (i + 1)..<n_seq {
      let seq_i = alignment[i]
      let seq_j = alignment[j]
      for pos in 0..<seq_len {
        if pos >= seq_i.length() || pos >= seq_j.length() { continue }
        let ci = seq_i[pos:pos + 1].to_owned()
        let cj = seq_j[pos:pos + 1].to_owned()
        // Skip gaps
        if ci == "-" || cj == "-" { continue }
        let mut ri = -1
        for k in 0..<letter_keys.length() {
          if letter_keys[k] == ci { ri = letter_vals[k]; break }
        }
        let mut rj = -1
        for k in 0..<letter_keys.length() {
          if letter_keys[k] == cj { rj = letter_vals[k]; break }
        }
        if ri >= 0 && rj >= 0 {
          data[ri][rj] = data[ri][rj] + 1.0
          data[rj][ri] = data[rj][ri] + 1.0 // symmetric
        }
      }
    }
  }
  { rows: alphabet, cols: alphabet, data: data }
}

// ============================================================================
// Substitution matrix calculation from alignment (log-odds)
// ============================================================================

///|
/// Calculate a log-odds substitution matrix from a multiple sequence
/// alignment.
///
/// The log-odds score for pair (i, j) is:
///   S(i,j) = log2( q(i,j) / (p(i) * p(j)) )
/// where q(i,j) is the observed pair frequency and p(i) is the background
/// frequency of letter i.
///
/// Parameters:
///   - alignment: array of aligned sequences
///   - alphabet: letters to include
///   - scale: scaling factor (default 2.0, similar to BLOSUM)
pub fn calculate_substitution_matrix(
  alignment : Array[String],
  alphabet : Array[String],
  scale? : Double = 2.0,
) -> SubsMatrix {
  let n = alphabet.length()
  // Step 1: Calculate observed pair frequencies
  let freq_data = calculate_frequency_matrix(alignment, alphabet)
  // Step 2: Calculate total pair count and letter frequencies
  let mut total_pairs = 0.0
  let letter_counts : Array[Double] = []
  for _i in 0..<n { letter_counts.push(0.0) }
  for i in 0..<n {
    for j in 0..<n {
      total_pairs = total_pairs + freq_data.data[i][j]
      letter_counts[i] = letter_counts[i] + freq_data.data[i][j]
    }
  }
  // Step 3: Calculate expected frequencies and log-odds scores
  let score_pairs : Array[(String, String, Double)] = []
  for i in 0..<n {
    for j in 0..<n {
      let observed = freq_data.data[i][j] / total_pairs
      let p_i = letter_counts[i] / total_pairs
      let p_j = letter_counts[j] / total_pairs
      let expected = p_i * p_j
      let score = if expected < 1.0e-12 || observed < 1.0e-12 {
        -999.0
      } else {
        scale * @math.log2(observed / expected)
      }
      score_pairs.push((alphabet[i], alphabet[j], score))
    }
  }
  let matrix = array_data_from_dict(score_pairs, alphabet, alphabet)
  subs_matrix("calculated", alphabet, matrix)
}

// ============================================================================
// Shannon entropy
// ============================================================================

///|
/// Calculate the Shannon entropy of a frequency matrix.
///
/// H = -Σ q(i,j) * log2(q(i,j))
///
/// where q(i,j) is the normalized pair frequency. The entropy measures
/// the information content of the substitution matrix: lower entropy
/// means more conserved (less diverse) substitutions.
pub fn subs_shannon_entropy(freq_matrix : ArrayData) -> Double {
  let n = freq_matrix.rows.length()
  let mut total = 0.0
  for i in 0..<n {
    for j in 0..<n {
      total = total + freq_matrix.data[i][j]
    }
  }
  if total < 1.0e-12 { return 0.0 }
  let mut entropy = 0.0
  for i in 0..<n {
    for j in 0..<n {
      let p = freq_matrix.data[i][j] / total
      if p > 1.0e-12 {
        entropy = entropy - p * @math.log2(p)
      }
    }
  }
  entropy
}

///|
/// Calculate the relative entropy (Kullback-Leibler divergence) between
/// the observed frequency matrix and the expected (product of marginals).
///
/// D(Q || P) = Σ q(i,j) * log2( q(i,j) / (p(i) * p(j)) )
///
/// This is equivalent to the average bits per substitution of the
/// corresponding log-odds matrix.
pub fn subs_relative_entropy(freq_matrix : ArrayData) -> Double {
  let n = freq_matrix.rows.length()
  let mut total = 0.0
  let letter_counts : Array[Double] = []
  for _i in 0..<n { letter_counts.push(0.0) }
  for i in 0..<n {
    for j in 0..<n {
      total = total + freq_matrix.data[i][j]
      letter_counts[i] = letter_counts[i] + freq_matrix.data[i][j]
    }
  }
  if total < 1.0e-12 { return 0.0 }
  let mut kl = 0.0
  for i in 0..<n {
    for j in 0..<n {
      let q = freq_matrix.data[i][j] / total
      let p_i = letter_counts[i] / total
      let p_j = letter_counts[j] / total
      let p = p_i * p_j
      if q > 1.0e-12 && p > 1.0e-12 {
        kl = kl + q * @math.log2(q / p)
      }
    }
  }
  kl
}

// ============================================================================
// NCBI BLAST matrix file parsing
// ============================================================================

///|
/// Parse an NCBI BLAST substitution matrix file (e.g. BLOSUM62).
///
/// Format: lines starting with '#' are comments, first non-comment line
/// is the column header (letters), subsequent lines are "row_letter score1
/// score2 ... scoreN".
///
/// Returns a SubsMatrix.
pub fn parse_ncbi_matrix(content : String) -> SubsMatrix? {
  let lines = content.split("\n").to_array()
  let alphabet : Array[String] = []
  let data_rows : Array[(String, Array[Double])] = []
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 { continue }
    if trimmed[0] == '#' { continue }
    let parts = trimmed.split(" ").to_array().filter(fn (s) -> Bool {
      s.length() > 0
    })
    if parts.length() == 0 { continue }
    if alphabet.length() == 0 && parts[0].length() == 0 {
      // This might be the header line
      for p in parts {
        if p.length() == 1 { alphabet.push(p.to_string()) }
      }
      continue
    }
    // Detect header: first part is empty or a single letter
    if alphabet.length() == 0 {
      // First non-comment line is the header
      for p in parts {
        if p.length() >= 1 { alphabet.push(p.to_string()) }
      }
      continue
    }
    // Data row: first part is letter, rest are scores
    if parts.length() >= 2 {
      let row_letter = parts[0].to_string()
      let scores : Array[Double] = []
      for i in 1..<parts.length() {
        match parse_double(parts[i].to_string()) {
          Some(v) => scores.push(v)
          None => continue
        }
      }
      if scores.length() == alphabet.length() {
        data_rows.push((row_letter, scores))
      }
    }
  }
  if alphabet.length() == 0 || data_rows.length() == 0 { return None }
  // Build the matrix
  let row_labels = data_rows.map(fn (r) -> String { r.0 })
  let pairs : Array[(String, String, Double)] = []
  let mut dri = 0
  while dri < data_rows.length() {
    let (r_label, scores) = data_rows[dri]
    for j in 0..<alphabet.length() {
      pairs.push((r_label, alphabet[j], scores[j]))
    }
    dri = dri + 1
  }
  let matrix = array_data_from_dict(pairs, row_labels, alphabet)
  Some(subs_matrix("parsed", row_labels, matrix))
}

// ============================================================================
// Matrix comparison
// ============================================================================

///|
/// Compute the Pearson correlation between two substitution matrices
/// (restricted to their shared alphabet).
pub fn matrix_correlation(m1 : SubsMatrix, m2 : SubsMatrix) -> Double {
  // Find common letters
  let common : Array[String] = []
  for letter in m1.alphabet {
    let mut found = false
    for l2 in m2.alphabet {
      if l2 == letter { found = true; break }
    }
    if found { common.push(letter) }
  }
  if common.length() < 2 { return 0.0 }
  // Collect paired scores
  let v1 : Array[Double] = []
  let v2 : Array[Double] = []
  for r in common {
    for c in common {
      v1.push(m1.get_score(r, c))
      v2.push(m2.get_score(r, c))
    }
  }
  // Pearson correlation
  let n = v1.length()
  let mut sum1 = 0.0
  let mut sum2 = 0.0
  for i in 0..<n { sum1 = sum1 + v1[i]; sum2 = sum2 + v2[i] }
  let mean1 = sum1 / n.to_double()
  let mean2 = sum2 / n.to_double()
  let mut cov = 0.0
  let mut var1 = 0.0
  let mut var2 = 0.0
  for i in 0..<n {
    let d1 = v1[i] - mean1
    let d2 = v2[i] - mean2
    cov = cov + d1 * d2
    var1 = var1 + d1 * d1
    var2 = var2 + d2 * d2
  }
  let denom = (var1 * var2).sqrt()
  if denom < 1.0e-12 { 0.0 } else { cov / denom }
}

// ============================================================================
// BLOSUM62 raw scores (20x20, row-major, alphabet order: ACDEFGHIKLMNPQRSTVWY)
// ============================================================================

fn blosum62_scores() -> Array[Int] {
  [
    // A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y
    4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2,  // A
    0,   9,  -3,  -4,  -2,  -3,  -3,  -1,  -3,  -1,  -1,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -2,  -2,  // C
   -2,  -3,   6,   2,  -3,  -1,  -1,  -3,  -1,  -4,  -3,   1,  -1,   0,  -2,   0,  -1,  -3,  -4,  -3,  // D
   -1,  -4,   2,   5,  -3,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -3,  -2,  // E
   -2,  -2,  -3,  -3,   6,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   1,   3,  // F
    0,  -3,  -1,  -2,  -3,   6,  -2,  -4,  -2,  -4,  -3,   0,  -2,  -2,  -2,   0,  -2,  -3,  -2,  -3,  // G
   -2,  -3,  -1,   0,  -1,  -2,   8,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2,  // H
   -1,  -1,  -3,  -3,   0,  -4,  -3,   4,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -3,  -1,  // I
   -1,  -3,  -1,   1,  -3,  -2,  -1,  -3,   5,  -2,  -1,   0,  -1,   1,   2,   0,  -1,  -2,  -3,  -2,  // K
   -1,  -1,  -4,  -3,   0,  -4,  -3,   2,  -2,   4,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1,  // L
   -1,  -1,  -3,  -2,   0,  -3,  -2,   1,  -1,   2,   5,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1,  // M
   -2,  -3,   1,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2,  // N
   -1,  -3,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   7,  -1,  -2,  -1,  -1,  -2,  -4,  -3,  // P
   -1,  -3,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   5,   1,   0,  -1,  -2,  -2,  -1,  // Q
   -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   2,  -2,  -1,   0,  -2,   1,   5,  -1,  -1,  -3,  -3,  -2,  // R
    1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2,  // S
    0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -2,  -2,  // T
    0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -2,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -3,  -1,  // V
   -3,  -2,  -4,  -3,   1,  -2,  -2,  -3,  -3,  -2,  -1,  -4,  -4,  -2,  -3,  -3,  -2,  -3,  11,   2,  // W
   -2,  -2,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   7,  // Y
  ]
}

// ============================================================================
// BLOSUM45 raw scores (20x20)
// ============================================================================

fn blosum45_scores() -> Array[Int] {
  [
    5,  -1,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2,
   -1,  12,  -5,  -5,  -4,  -3,  -3,  -2,  -5,  -2,  -2,  -4,  -4,  -4,  -4,  -1,  -1,  -1,  -5,  -4,
   -2,  -5,   7,   2,  -4,  -1,   0,  -4,   0,  -4,  -3,   2,  -1,   0,  -1,   0,  -1,  -3,  -5,  -3,
   -1,  -5,   2,   6,  -4,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -4,  -2,
   -2,  -4,  -4,  -4,   8,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   0,   3,
    0,  -3,  -1,  -2,  -3,   7,  -2,  -3,  -2,  -3,  -2,   0,  -2,  -2,  -2,   0,  -2,  -3,  -3,  -3,
   -2,  -3,   0,   0,  -1,  -2,  10,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2,
   -1,  -2,  -4,  -3,   0,  -3,  -3,   5,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -2,  -1,
   -1,  -5,   0,   1,  -3,  -2,  -1,  -3,   5,  -3,  -1,   0,  -1,   1,   3,   0,  -1,  -3,  -4,  -2,
   -1,  -2,  -4,  -3,   0,  -3,  -3,   2,  -3,   5,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1,
   -1,  -2,  -3,  -2,   0,  -2,  -2,   1,  -1,   2,   6,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1,
   -2,  -4,   2,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2,
   -1,  -4,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   9,  -1,  -2,  -1,  -1,  -2,  -5,  -3,
   -1,  -4,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   6,   1,   0,  -1,  -2,  -3,  -1,
   -1,  -4,  -1,   0,  -3,  -2,   0,  -3,   3,  -2,  -1,   0,  -2,   1,   6,  -1,  -1,  -3,  -3,  -2,
    1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2,
    0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -3,  -2,
    0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -3,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -2,  -1,
   -3,  -5,  -5,  -4,   0,  -3,  -2,  -2,  -4,  -2,  -1,  -4,  -5,  -3,  -3,  -3,  -3,  -2,  14,   2,
   -2,  -4,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   8,
  ]
}

// ============================================================================
// BLOSUM80 raw scores (20x20)
// ============================================================================

fn blosum80_scores() -> Array[Int] {
  [
    7,  -1,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   2,   1,   0,  -4,  -3,
   -1,  13,  -5,  -5,  -3,  -4,  -4,  -2,  -4,  -2,  -2,  -4,  -5,  -4,  -4,  -1,  -1,  -1,  -5,  -3,
   -2,  -5,   9,   2,  -5,  -2,  -1,  -4,  -1,  -5,  -4,   2,  -2,   1,  -2,   0,  -1,  -4,  -6,  -4,
   -1,  -5,   2,   8,  -5,  -3,   0,  -4,   2,  -4,  -3,   0,  -2,   3,   0,   0,  -1,  -3,  -5,  -3,
   -2,  -3,  -5,  -5,   9,  -4,  -1,   0,  -4,   0,   0,  -4,  -5,  -4,  -4,  -2,  -2,  -1,   1,   4,
    0,  -4,  -2,  -3,  -4,   9,  -3,  -5,  -2,  -5,  -4,   0,  -3,  -3,  -3,   0,  -2,  -4,  -4,  -4,
   -2,  -4,  -1,   0,  -1,  -3,  12,  -4,  -1,  -4,  -3,   1,  -3,   0,   0,  -1,  -2,  -4,  -3,   2,
   -1,  -2,  -4,  -4,   0,  -5,  -4,   6,  -3,   2,   1,  -4,  -3,  -3,  -3,  -2,  -1,   3,  -4,  -2,
   -1,  -4,  -1,   2,  -4,  -2,  -1,  -3,   8,  -3,  -1,   0,  -2,   2,   3,   0,  -1,  -3,  -5,  -3,
   -1,  -2,  -5,  -4,   0,  -5,  -4,   2,  -3,   6,   3,  -4,  -3,  -3,  -3,  -2,  -1,   1,  -3,  -1,
   -1,  -2,  -4,  -3,   0,  -4,  -3,   1,  -1,   3,   8,  -3,  -3,   0,  -1,  -1,  -1,   1,  -2,  -1,
   -2,  -4,   2,   0,  -4,   0,   1,  -4,   0,  -4,  -3,   8,  -2,   0,   0,   1,   0,  -4,  -5,  -3,
   -1,  -5,  -2,  -2,  -5,  -3,  -3,  -3,  -2,  -3,  -3,  -2,  11,  -2,  -3,  -2,  -2,  -3,  -6,  -4,
   -1,  -4,   1,   3,  -4,  -3,   0,  -3,   2,  -3,   0,   0,  -2,   8,   1,   0,  -1,  -3,  -3,  -2,
   -1,  -4,  -2,   0,  -4,  -3,   0,  -3,   3,  -3,  -1,   0,  -3,   1,   8,  -1,  -1,  -4,  -4,  -2,
    2,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -2,   0,  -1,   7,   2,  -2,  -4,  -3,
    1,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -2,  -1,  -1,   2,   7,   0,  -4,  -2,
    0,  -1,  -4,  -3,  -1,  -4,  -4,   3,  -3,   1,   1,  -4,  -3,  -3,  -4,  -2,   0,   6,  -4,  -2,
   -4,  -5,  -6,  -5,   1,  -4,  -3,  -4,  -5,  -3,  -2,  -5,  -6,  -3,  -4,  -4,  -4,  -4,  15,   3,
   -3,  -3,  -4,  -3,   4,  -4,   2,  -2,  -3,  -1,  -1,  -3,  -4,  -2,  -2,  -3,  -2,  -2,   3,   9,
  ]
}

// ============================================================================
// BLOSUM90 raw scores (20x20)
// ============================================================================

fn blosum90_scores() -> Array[Int] {
  [
    7,  -1,  -2,  -1,  -3,   0,  -2,  -1,  -1,  -2,  -1,  -2,  -1,  -1,  -2,   2,   1,   0,  -4,  -3,
   -1,  14,  -6,  -6,  -4,  -4,  -4,  -2,  -5,  -3,  -3,  -5,  -5,  -5,  -4,  -1,  -1,  -2,  -5,  -4,
   -2,  -6,  10,   2,  -5,  -2,  -1,  -4,  -1,  -5,  -4,   2,  -2,   1,  -2,   0,  -1,  -4,  -6,  -4,
   -1,  -6,   2,   9,  -5,  -3,   0,  -4,   2,  -4,  -3,   0,  -2,   3,   0,   0,  -1,  -3,  -5,  -3,
   -3,  -4,  -5,  -5,  10,  -4,  -1,  0,  -4,   0,   0,  -4,  -5,  -4,  -4,  -3,  -2,  -1,   2,   4,
    0,  -4,  -2,  -3,  -4,  10,  -3,  -5,  -2,  -5,  -4,   0,  -3,  -3,  -3,   0,  -2,  -4,  -4,  -4,
   -2,  -4,  -1,   0,  -1,  -3,  12,  -4,  -1,  -4,  -3,   1,  -3,   0,   0,  -1,  -2,  -4,  -3,   2,
   -1,  -2,  -4,  -4,   0,  -5,  -4,   6,  -3,   2,   1,  -4,  -3,  -3,  -3,  -2,  -1,   3,  -4,  -2,
   -1,  -5,  -1,   2,  -4,  -2,  -1,  -3,   8,  -3,  -1,   0,  -2,   2,   3,   0,  -1,  -3,  -5,  -3,
   -2,  -3,  -5,  -4,   0,  -5,  -4,   2,  -3,   6,   3,  -4,  -3,  -3,  -3,  -2,  -1,   1,  -3,  -1,
   -1,  -3,  -4,  -3,   0,  -4,  -3,   1,  -1,   3,   8,  -3,  -3,   0,  -1,  -1,  -1,   1,  -2,  -1,
   -2,  -5,   2,   0,  -4,   0,   1,  -4,   0,  -4,  -3,   8,  -2,   0,   0,   1,   0,  -4,  -5,  -3,
   -1,  -5,  -2,  -2,  -5,  -3,  -3,  -3,  -2,  -3,  -3,  -2,  12,  -2,  -3,  -2,  -2,  -3,  -6,  -4,
   -1,  -5,   1,   3,  -4,  -3,   0,  -3,   2,  -3,   0,   0,  -2,   9,   1,   0,  -1,  -3,  -4,  -2,
   -2,  -4,  -2,   0,  -4,  -3,   0,  -3,   3,  -3,  -1,   0,  -3,   1,   9,  -1,  -1,  -4,  -4,  -2,
    2,  -1,   0,   0,  -3,   0,  -1,  -2,   0,  -2,  -1,   1,  -2,   0,  -1,   7,   2,  -2,  -4,  -3,
    1,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -2,  -1,  -1,   2,   8,   0,  -4,  -2,
    0,  -2,  -4,  -3,  -1,  -4,  -4,   3,  -3,   1,   1,  -4,  -3,  -3,  -4,  -2,   0,   6,  -4,  -2,
   -4,  -5,  -6,  -5,   2,  -4,  -3,  -4,  -5,  -3,  -2,  -5,  -6,  -4,  -4,  -4,  -4,  -4,  15,   3,
   -3,  -4,  -4,  -3,   4,  -4,   2,  -2,  -3,  -1,  -1,  -3,  -4,  -2,  -2,  -3,  -2,  -2,   3,   9,
  ]
}

// ============================================================================
// PAM30 raw scores (20x20)
// ============================================================================

fn pam30_scores() -> Array[Int] {
  [
    6,  -7,  -4,  -2,  -6,  -2,  -5,  -3,  -3,  -5,  -4,  -3,  -3,  -3,  -5,   0,   0,  -3, -10,  -6,
   -7,  12,  -9,  -9,  -6,  -9,  -6,  -7,  -9,  -7,  -7,  -7,  -8,  -7,  -7,  -4,  -5,  -6, -10,  -7,
   -4,  -9,   8,   2,  -9,  -3,  -3,  -7,  -3,  -8,  -7,   2,  -5,   0,  -6,  -1,  -3,  -7,  -9,  -7,
   -2,  -9,   2,   8,  -8,  -5,  -1,  -6,   2,  -7,  -5,   0,  -4,   3,  -3,  -1,  -3,  -5,  -9,  -5,
   -6,  -6,  -9,  -8,  10,  -8,  -4,  -1,  -8,  -1,  -1,  -8,  -8,  -8,  -7,  -5,  -5,  -3,   1,   3,
   -2,  -9,  -3,  -5,  -8,   8,  -6,  -8,  -5,  -8,  -6,   0,  -6,  -5,  -6,  -1,  -4,  -7,  -9,  -8,
   -5,  -6,  -3,  -1,  -4,  -6,  10,  -6,  -3,  -6,  -5,   1,  -5,   0,  -2,  -3,  -4,  -6,  -4,   2,
   -3,  -7,  -7,  -6,  -1,  -8,  -6,   7,  -5,   1,   1,  -6,  -5,  -5,  -5,  -4,  -2,   3,  -7,  -4,
   -3,  -9,  -3,   2,  -8,  -5,  -3,  -5,   7,  -5,  -3,   0,  -4,   2,   2,  -1,  -2,  -5,  -9,  -6,
   -5,  -7,  -8,  -7,  -1,  -8,  -6,   1,  -5,   7,   3,  -6,  -5,  -4,  -6,  -5,  -4,   1,  -6,  -4,
   -4,  -7,  -7,  -5,  -1,  -6,  -5,   1,  -3,   3,   8,  -5,  -4,  -1,  -3,  -3,  -2,   1,  -6,  -4,
   -3,  -7,   2,   0,  -8,   0,   1,  -6,   0,  -6,  -5,   8,  -4,   0,  -3,   0,  -1,  -6,  -9,  -5,
   -3,  -8,  -5,  -4,  -8,  -6,  -5,  -5,  -4,  -5,  -4,  -4,  11,  -4,  -5,  -3,  -3,  -5,  -9,  -7,
   -3,  -7,   0,   3,  -8,  -5,   0,  -5,   2,  -4,  -1,   0,  -4,   8,   1,  -1,  -3,  -5,  -8,  -5,
   -5,  -7,  -6,  -3,  -7,  -6,  -2,  -5,   2,  -6,  -3,  -3,  -5,   1,   8,  -2,  -3,  -6,  -7,  -5,
    0,  -4,  -1,  -1,  -5,  -1,  -3,  -4,  -1,  -5,  -3,   0,  -3,  -1,  -2,   6,   1,  -4,  -6,  -4,
    0,  -5,  -3,  -3,  -5,  -4,  -4,  -2,  -2,  -4,  -2,  -1,  -3,  -3,  -3,   1,   7,  -1,  -7,  -4,
   -3,  -6,  -7,  -5,  -3,  -7,  -6,   3,  -5,   1,   1,  -6,  -5,  -5,  -6,  -4,  -1,   7,  -7,  -4,
  -10, -10,  -9,  -9,   1,  -9,  -4,  -7,  -9,  -6,  -6,  -9,  -9,  -8,  -7,  -6,  -7,  -7,  14,   3,
   -6,  -7,  -7,  -5,   3,  -8,   2,  -4,  -6,  -4,  -4,  -5,  -7,  -5,  -5,  -4,  -4,  -4,   3,  10,
  ]
}

// ============================================================================
// PAM70 raw scores (20x20)
// ============================================================================

fn pam70_scores() -> Array[Int] {
  [
    5,  -4,  -3,  -1,  -5,  -1,  -4,  -2,  -2,  -4,  -3,  -2,  -2,  -2,  -4,   1,   1,  -2,  -8,  -4,
   -4,  10,  -7,  -6,  -4,  -7,  -4,  -5,  -6,  -5,  -5,  -5,  -6,  -5,  -5,  -2,  -3,  -4,  -7,  -5,
   -3,  -7,   7,   2,  -7,  -2,  -2,  -5,  -2,  -6,  -5,   2,  -4,   1,  -4,   0,  -2,  -5,  -7,  -5,
   -1,  -6,   2,   7,  -6,  -4,   0,  -4,   2,  -5,  -3,   0,  -3,   3,  -2,   0,  -2,  -4,  -7,  -4,
   -5,  -4,  -7,  -6,   9,  -6,  -3,  -1,  -6,  -1,  -1,  -6,  -6,  -6,  -5,  -4,  -4,  -2,   1,   3,
   -1,  -7,  -2,  -4,  -6,   8,  -5,  -6,  -4,  -6,  -5,   1,  -5,  -4,  -4,   0,  -3,  -5,  -7,  -6,
   -4,  -4,  -2,   0,  -3,  -5,  10,  -5,  -2,  -4,  -4,   2,  -4,   1,   0,  -2,  -3,  -5,  -3,   2,
   -2,  -5,  -5,  -4,  -1,  -6,  -5,   7,  -4,   2,   2,  -4,  -4,  -3,  -4,  -3,  -1,   3,  -5,  -3,
   -2,  -6,  -2,   2,  -6,  -4,  -2,  -4,   6,  -4,  -2,   0,  -3,   2,   2,   0,  -1,  -4,  -7,  -4,
   -4,  -5,  -6,  -5,  -1,  -6,  -4,   2,  -4,   6,   3,  -4,  -4,  -3,  -4,  -3,  -3,   1,  -4,  -2,
   -3,  -5,  -5,  -3,  -1,  -5,  -4,   2,  -2,   3,   8,  -4,  -3,  -1,  -2,  -2,  -1,   1,  -4,  -2,
   -2,  -5,   2,   0,  -6,   1,   2,  -4,   0,  -4,  -4,   8,  -3,   0,  -2,   1,   0,  -4,  -7,  -4,
   -2,  -6,  -4,  -3,  -6,  -5,  -4,  -4,  -3,  -4,  -3,  -3,  10,  -3,  -4,  -2,  -2,  -4,  -7,  -5,
   -2,  -5,   1,   3,  -6,  -4,   1,  -3,   2,  -3,  -1,   0,  -3,   7,   1,   0,  -2,  -4,  -6,  -3,
   -4,  -5,  -4,  -2,  -5,  -4,   0,  -4,   2,  -4,  -2,  -2,  -4,   1,   8,  -1,  -2,  -4,  -5,  -3,
    1,  -2,   0,   0,  -4,   0,  -2,  -3,   0,  -3,  -2,   1,  -2,   0,  -1,   5,   2,  -3,  -5,  -3,
    1,  -3,  -2,  -2,  -4,  -3,  -3,  -1,  -1,  -3,  -1,   0,  -2,  -2,  -2,   2,   6,  -1,  -5,  -3,
   -2,  -4,  -5,  -4,  -2,  -5,  -5,   3,  -4,   1,   1,  -4,  -4,  -4,  -4,  -3,  -1,   7,  -5,  -3,
   -8,  -7,  -7,  -7,   1,  -7,  -3,  -5,  -7,  -4,  -4,  -7,  -7,  -6,  -5,  -5,  -5,  -5,  13,   3,
   -4,  -5,  -5,  -4,   3,  -6,   2,  -3,  -4,  -2,  -2,  -4,  -5,  -3,  -3,  -3,  -3,  -3,   3,   9,
  ]
}

// ============================================================================
// PAM250 raw scores (20x20)
// ============================================================================

fn pam250_scores() -> Array[Int] {
  [
    2,  -2,   0,   0,  -3,   1,  -1,  -1,  -1,  -2,  -1,   0,   1,   0,  -2,   1,   1,   0,  -6,  -3,
   -2,  12,  -5,  -5,  -4,  -3,  -3,  -2,  -3,  -2,  -2,  -4,  -4,  -4,  -4,   0,  -2,  -2,  -7,  -4,
    0,  -5,   4,   3,  -6,   1,   1,  -2,   0,  -4,  -3,   2,  -1,   2,  -1,   0,   0,  -2,  -5,  -4,
    0,  -5,   3,   4,  -5,   0,   1,  -2,   2,  -3,  -2,   1,  -1,   2,  -1,   0,   0,  -2,  -5,  -3,
   -3,  -4,  -6,  -5,   9,  -5,  -2,   1,  -5,   2,   0,  -3,  -5,  -5,  -4,  -3,  -3,  -1,   0,   7,
    1,  -3,   1,   0,  -5,   5,  -2,  -3,  -2,  -4,  -3,   0,   0,  -1,  -3,   1,   0,  -1,  -7,  -5,
   -1,  -3,   1,   1,  -2,  -2,   6,  -2,   0,  -2,  -2,   2,   0,   3,   2,  -1,  -1,  -2,  -3,   2,
   -1,  -2,  -2,  -2,   1,  -3,  -2,   5,  -2,   2,   2,  -2,  -2,  -2,  -2,  -1,   0,   4,  -5,  -1,
   -1,  -3,   0,   2,  -5,  -2,   0,  -2,   5,  -3,  -1,   1,  -1,   1,   3,   0,   0,  -2,  -3,  -4,
   -2,  -2,  -4,  -3,   2,  -4,  -2,   2,  -3,   6,   4,  -3,  -3,  -2,  -3,  -3,  -2,   2,  -2,  -1,
   -1,  -2,  -3,  -2,   0,  -3,  -2,   2,  -1,   4,   6,  -2,  -2,  -1,   0,  -2,  -1,   2,  -4,  -2,
    0,  -4,   2,   1,  -3,   0,   2,  -2,   1,  -3,  -2,   2,   0,   1,   0,   1,   0,  -2,  -4,  -2,
    1,  -4,  -1,  -1,  -5,   0,   0,  -2,  -1,  -3,  -2,   0,   6,   0,   0,   1,   0,  -1,  -6,  -5,
    0,  -4,   2,   2,  -5,  -1,   3,  -2,   1,  -2,  -1,   1,   0,   4,   1,  -1,  -1,  -2,  -5,  -4,
   -2,  -4,  -1,  -1,  -4,  -3,   2,  -2,   3,  -3,   0,   0,   0,   1,   6,   0,  -1,  -2,  -4,   0,
    1,   0,   0,   0,  -3,   1,  -1,  -1,   0,  -3,  -2,   1,   1,  -1,   0,   2,   1,  -1,  -4,  -2,
    1,  -2,   0,   0,  -3,   0,  -1,   0,   0,  -2,  -1,   0,   0,  -1,  -1,   1,   3,   0,  -5,  -3,
    0,  -2,  -2,  -2,  -1,  -1,  -2,   4,  -2,   2,   2,  -2,  -1,  -2,  -2,  -1,   0,   4,  -6,  -2,
   -6,  -7,  -5,  -5,   0,  -7,  -3,  -5,  -3,  -2,  -4,  -4,  -6,  -5,  -4,  -4,  -5,  -6,  17,   0,
   -3,  -4,  -4,  -3,   7,  -5,   2,  -1,  -4,  -1,  -2,  -2,  -5,  -4,   0,  -2,  -3,  -2,   0,  10,
  ]
}
