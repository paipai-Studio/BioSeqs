// Alignment column statistics ported from Biopython `Bio.Align`:
// MultipleSeqAlignment.substitutions (symmetric residue-pair count
// matrix with per-sequence weights), Alignment.frequencies (per-column
// letter counts/frequencies, gaps represented as a '-') and the
// relative pair frequencies q_ij / e_ij used to build log-odds
// substitution matrices from an MSA.
//
// Input rows are printed-style alignment strings of equal length, with
// '-' (or '.') denoting a gap. All arithmetic uses the sequence weights
// (Biopython `record.annotations["weight"]`, 1.0 by default).

///|
/// Error raised while computing alignment letter statistics.
pub suberror AlignmentFrequenciesError {
  AlignmentFrequenciesError(String)
}

///|
/// Per-column (optionally weighted) letter counts of an alignment,
/// equivalent to Biopython `Alignment.frequencies` for a printed
/// alignment. Gaps are stored as the letter `'-'` (the gap row is
/// present only if the alignment contains gaps). `counts[k][j]` is the
/// total weight of letter `letters[k]` in column `j`.
pub struct AlignmentLetterCounts {
  letters : Array[String]
  ncols : Int
  counts : Array[Array[Double]]
}

///|
/// Symmetric residue-substitution count matrix of a multiple sequence
/// alignment, equivalent to Biopython `MultipleSeqAlignment.substitutions`.
/// Gap columns are excluded from the alphabet; each unordered sequence
/// pair contributes `w1 * w2` to the matrix cell of its two residues,
/// split evenly across the symmetric cells (M = (M + M^T) / 2).
pub struct SubstitutionCountMatrix {
  letters : Array[String]
  matrix : Array[Array[Double]]
}

///|
fn afreq_is_gap(ch : Char) -> Bool {
  ch == '-' || ch == '.'
}

///|
fn afreq_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn afreq_contains(letters : Array[String], key : String) -> Bool {
  let mut i = 0
  while i < letters.length() {
    if letters[i] == key {
      return true
    }
    i = i + 1
  }
  false
}

///|
fn afreq_index(
  letters : Array[String],
  key : String,
) -> Int? {
  let mut i = 0
  while i < letters.length() {
    if letters[i] == key {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
fn afreq_require_index(
  letters : Array[String],
  key : String,
) -> Int raise AlignmentFrequenciesError {
  match afreq_index(letters, key) {
    Some(i) => i
    None => raise AlignmentFrequenciesError("internal error: missing letter '" + key + "'")
  }
}

///|
fn aln_prepare(
  rows : Array[String],
  weights : Array[Double],
) -> (Int, Array[Double]) raise AlignmentFrequenciesError {
  if rows.length() == 0 {
    raise AlignmentFrequenciesError("alignment must contain at least one sequence")
  }
  let ncols = rows[0].length()
  if ncols == 0 {
    raise AlignmentFrequenciesError("alignment columns must not be empty")
  }
  let mut i = 1
  while i < rows.length() {
    if rows[i].length() != ncols {
      raise AlignmentFrequenciesError(
        "alignment rows have different lengths (" +
          ncols.to_string() +
          " vs " +
          rows[i].length().to_string() +
          ")",
      )
    }
    i = i + 1
  }
  let w = if weights.length() == 0 {
    Array::make(rows.length(), 1.0)
  } else {
    if weights.length() != rows.length() {
      raise AlignmentFrequenciesError(
        "weights length (" +
          weights.length().to_string() +
          ") does not match number of sequences (" +
          rows.length().to_string() +
          ")",
      )
    }
    weights
  }
  (ncols, w)
}

///|
/// Insertion-sort helper for deterministic alphabet ordering.
fn afreq_sort_letters(letters : Array[String]) -> Unit {
  let mut i = 1
  while i < letters.length() {
    let key = letters[i]
    let mut j = i
    while j > 0 && letters[j - 1] > key {
      letters[j] = letters[j - 1]
      j = j - 1
    }
    letters[j] = key
    i = i + 1
  }
}

///|
/// Collect the sorted alphabet of an alignment. Gaps (`'-'`/`'.'`) are
/// reported as the single letter `'-'`; when `include_gaps` is false,
/// gaps are excluded entirely (Biopython substitution matrix alphabet).
fn afreq_alphabet(
  rows : Array[String],
  ncols : Int,
  include_gaps : Bool,
) -> Array[String] {
  let letters : Array[String] = []
  let mut i = 0
  while i < rows.length() {
    let row = rows[i]
    let mut j = 0
    while j < ncols {
      let ch = row.unsafe_get(j).unsafe_to_char()
      if afreq_is_gap(ch) {
        if include_gaps && not(afreq_contains(letters, "-")) {
          letters.push("-")
        }
      } else {
        let key = ch.to_string()
        if not(afreq_contains(letters, key)) {
          letters.push(key)
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  afreq_sort_letters(letters)
  letters
}

///|
/// Per-column letter counts (optionally sequence-weighted) for a
/// printed alignment, equivalent to Biopython `Alignment.frequencies`
/// (raw weighted counts; normalize with `column_frequencies`). Gaps are
/// reported under the letter `'-'`.
pub fn alignment_letter_counts(
  rows : Array[String],
  weights? : Array[Double] = [],
) -> AlignmentLetterCounts raise AlignmentFrequenciesError {
  let (ncols, w) = aln_prepare(rows, weights)
  let letters = afreq_alphabet(rows, ncols, true)
  let counts : Array[Array[Double]] = []
  let mut li = 0
  while li < letters.length() {
    counts.push(Array::make(ncols, 0.0))
    li = li + 1
  }
  let mut i = 0
  while i < rows.length() {
    let row = rows[i]
    let wi = w[i]
    let mut j = 0
    while j < ncols {
      let ch = row.unsafe_get(j).unsafe_to_char()
      let key = if afreq_is_gap(ch) { "-" } else { ch.to_string() }
      let k = afreq_require_index(letters, key)
      counts[k][j] += wi
      j = j + 1
    }
    i = i + 1
  }
  AlignmentLetterCounts::{ letters, ncols, counts }
}

///|
/// Letter index in a column-count matrix.
pub fn AlignmentLetterCounts::letter_index(
  self : AlignmentLetterCounts,
  letter : String,
) -> Int? {
  afreq_index(self.letters, letter)
}

///|
/// Weighted count of `letter` in column `col` (0 if not present).
pub fn AlignmentLetterCounts::count_of(
  self : AlignmentLetterCounts,
  letter : String,
  col : Int,
) -> Double {
  match self.letter_index(letter) {
    Some(i) => self.counts[i][col]
    None => 0.0
  }
}

///|
/// Column-normalized letter frequencies: each column sums to 1.0
/// (including the gap letter), so `result[k][j]` is the weighted fraction
/// of `letters[k]` in column `j`.
pub fn AlignmentLetterCounts::column_frequencies(
  self : AlignmentLetterCounts,
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = []
  let mut ok = 0
  while ok < self.letters.length() {
    out.push(Array::make(self.ncols, 0.0))
    ok = ok + 1
  }
  let mut j = 0
  while j < self.ncols {
    let mut total = 0.0
    let mut k = 0
    while k < self.letters.length() {
      total += self.counts[k][j]
      k = k + 1
    }
    if total > 0.0 {
      k = 0
      while k < self.letters.length() {
        out[k][j] = self.counts[k][j] / total
        k = k + 1
      }
    }
    j = j + 1
  }
  out
}

///|
/// Whole-alignment background letter frequencies (weighted), excluding
/// gaps: total count of each non-gap letter divided by the total number
/// of non-gap residues.
pub fn AlignmentLetterCounts::background_frequencies(
  self : AlignmentLetterCounts,
) -> Map[String, Double] {
  let totals : Array[Double] = Array::make(self.letters.length(), 0.0)
  let mut total = 0.0
  let mut k = 0
  while k < self.letters.length() {
    if self.letters[k] != "-" {
      let mut s = 0.0
      let mut j = 0
      while j < self.ncols {
        s += self.counts[k][j]
        j = j + 1
      }
      totals[k] = s
      total += s
    }
    k = k + 1
  }
  let out : Map[String, Double] = Map([])
  k = 0
  while k < self.letters.length() {
    if self.letters[k] != "-" {
      out.set(self.letters[k], if total > 0.0 { totals[k] / total } else { 0.0 })
    }
    k = k + 1
  }
  out
}

///|
/// Symmetric residue substitution count matrix of a multiple sequence
/// alignment (Biopython `MultipleSeqAlignment.substitutions`).
/// Gaps are excluded; sequence `i` and `j` contribute `weights[i] *
/// weights[j]` per non-gap column; the matrix is symmetrized as
/// M = (M + M^T) / 2.
pub fn alignment_substitution_matrix(
  rows : Array[String],
  weights? : Array[Double] = [],
) -> SubstitutionCountMatrix raise AlignmentFrequenciesError {
  let (ncols, w) = aln_prepare(rows, weights)
  let letters = afreq_alphabet(rows, ncols, false)
  let n = letters.length()
  let m : Array[Array[Double]] = []
  let mut mi = 0
  while mi < n {
    m.push(Array::make(n, 0.0))
    mi = mi + 1
  }
  let mut i1 = 0
  while i1 < rows.length() {
    let seq1 = rows[i1]
    let w1 = w[i1]
    let mut i2 = 0
    while i2 < i1 {
      let seq2 = rows[i2]
      let w2 = w[i2]
      let mut j = 0
      while j < ncols {
        let c1 = seq1.unsafe_get(j).unsafe_to_char()
        let c2 = seq2.unsafe_get(j).unsafe_to_char()
        if not(afreq_is_gap(c1)) && not(afreq_is_gap(c2)) {
          let k1 = afreq_require_index(letters, c1.to_string())
          let k2 = afreq_require_index(letters, c2.to_string())
          m[k1][k2] += w1 * w2
        }
        j = j + 1
      }
      i2 = i2 + 1
    }
    i1 = i1 + 1
  }
  // Symmetrize: M <- (M + M^T) / 2.
  let mut a = 0
  while a < n {
    let mut b = a
    while b < n {
      let v = (m[a][b] + m[b][a]) / 2.0
      m[a][b] = v
      m[b][a] = v
      b = b + 1
    }
    a = a + 1
  }
  SubstitutionCountMatrix::{ letters, matrix: m }
}

///|
/// Index of a residue in the substitution alphabet.
pub fn SubstitutionCountMatrix::letter_index(
  self : SubstitutionCountMatrix,
  letter : String,
) -> Int? {
  afreq_index(self.letters, letter)
}

///|
/// Symmetric substitution count between residues `a` and `b`.
pub fn SubstitutionCountMatrix::count_of(
  self : SubstitutionCountMatrix,
  a : String,
  b : String,
) -> Double raise AlignmentFrequenciesError {
  let ia = match self.letter_index(a) {
    Some(v) => v
    None => raise AlignmentFrequenciesError("unknown residue '" + a + "'")
  }
  let ib = match self.letter_index(b) {
    Some(v) => v
    None => raise AlignmentFrequenciesError("unknown residue '" + b + "'")
  }
  self.matrix[ia][ib]
}

///|
/// Total mass of the matrix (sum of all cells).
pub fn SubstitutionCountMatrix::total(self : SubstitutionCountMatrix) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < self.letters.length() {
    let mut j = 0
    while j < self.letters.length() {
      total += self.matrix[i][j]
      j = j + 1
    }
    i = i + 1
  }
  total
}

///|
/// Observed pair frequencies q_ij: the count matrix normalized so all
/// cells sum to 1.0.
pub fn SubstitutionCountMatrix::pair_frequencies(
  self : SubstitutionCountMatrix,
) -> Array[Array[Double]] {
  let total = self.total()
  let n = self.letters.length()
  let q : Array[Array[Double]] = []
  let mut qi = 0
  while qi < n {
    q.push(Array::make(n, 0.0))
    qi = qi + 1
  }
  if total > 0.0 {
    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        q[i][j] = self.matrix[i][j] / total
        j = j + 1
      }
      i = i + 1
    }
  }
  q
}

///|
/// Residue background frequencies estimated from the pair counts:
/// p_i = sum_j q_ij (row sums of the pair-frequency matrix). With
/// half-count symmetric cells this gives the standard BLOSUM-style
/// maximum-likelihood background of the alignment pairs.
pub fn SubstitutionCountMatrix::background_frequencies(
  self : SubstitutionCountMatrix,
) -> Map[String, Double] {
  let q = self.pair_frequencies()
  let out : Map[String, Double] = Map([])
  let mut i = 0
  while i < self.letters.length() {
    let mut s = 0.0
    let mut j = 0
    while j < self.letters.length() {
      s += q[i][j]
      j = j + 1
    }
    out.set(self.letters[i], s)
    i = i + 1
  }
  out
}

///|
/// Relative pair frequencies q_ij / e_ij where the expected pair
/// frequency is e_ij = p_i * p_j with `p` the maximum-likelihood
/// background estimated from the pair counts. Values of 1.0 match
/// expectation; the log-odds score used by substitution matrices is
/// log2 of this ratio.
pub fn SubstitutionCountMatrix::relative_frequencies(
  self : SubstitutionCountMatrix,
) -> Array[Array[Double]] raise AlignmentFrequenciesError {
  self.relative_frequencies_expected(self.background_frequencies())
}

///|
/// Relative pair frequencies q_ij / e_ij using a user-supplied residue
/// background `expected` (one positive probability per residue).
pub fn SubstitutionCountMatrix::relative_frequencies_expected(
  self : SubstitutionCountMatrix,
  expected : Map[String, Double],
) -> Array[Array[Double]] raise AlignmentFrequenciesError {
  let q = self.pair_frequencies()
  let n = self.letters.length()
  let p = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    match expected.get(self.letters[i]) {
      Some(v) => {
        if not(afreq_is_finite(v)) || v <= 0.0 {
          raise AlignmentFrequenciesError(
            "expected frequency of '" +
              self.letters[i] +
              "' must be a positive finite number",
          )
        }
        p[i] = v
      }
      None => raise AlignmentFrequenciesError(
        "missing expected frequency for residue '" + self.letters[i] + "'",
      )
    }
    i = i + 1
  }
  let out : Array[Array[Double]] = []
  let mut oi = 0
  while oi < n {
    out.push(Array::make(n, 0.0))
    oi = oi + 1
  }
  i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let e = p[i] * p[j]
      out[i][j] = if e > 0.0 { q[i][j] / e } else { 0.0 }
      j = j + 1
    }
    i = i + 1
  }
  out
}

///|
/// Left-pad a string with spaces to the requested width.
fn afreq_pad_left(text : String, width : Int) -> String {
  let mut out = text
  while out.length() < width {
    out = " " + out
  }
  out
}

///|
/// Fixed-point formatting helper (3 decimals) for matrix tables.
fn afreq_fmt(x : Double) -> String {
  let neg = x < 0.0
  let v = x.abs() + 0.0005
  let whole = v.floor().to_int()
  let frac = ((v - v.floor()) * 1000.0).floor().to_int()
  let frac_s = if frac >= 100 {
    frac.to_string()
  } else if frac >= 10 {
    "0" + frac.to_string()
  } else {
    "00" + frac.to_string()
  }
  (if neg { "-" } else { "" }) + whole.to_string() + "." + frac_s
}

///|
/// Render the substitution count matrix as a Biopython-style symmetric
/// table.
pub fn SubstitutionCountMatrix::to_table(self : SubstitutionCountMatrix) -> String {
  let mut line = "      "
  let mut k = 0
  while k < self.letters.length() {
    line += afreq_pad_left(self.letters[k], 6)
    k = k + 1
  }
  line += "\n"
  let mut i = 0
  while i < self.letters.length() {
    line += afreq_pad_left(self.letters[i], 5) + " "
    let mut j = 0
    while j < self.letters.length() {
      line += afreq_pad_left(afreq_fmt(self.matrix[i][j]), 6)
      j = j + 1
    }
    line += "\n"
    i = i + 1
  }
  line
}

///|
/// Render the per-column letter-count matrix as a Biopython-style table.
pub fn AlignmentLetterCounts::to_table(self : AlignmentLetterCounts) -> String {
  let mut line = "      "
  let mut j = 0
  while j < self.ncols {
    line += afreq_pad_left(j.to_string(), 7)
    j = j + 1
  }
  line += "\n"
  let mut i = 0
  while i < self.letters.length() {
    line += afreq_pad_left(self.letters[i], 5) + " "
    j = 0
    while j < self.ncols {
      line += afreq_pad_left(afreq_fmt(self.counts[i][j]), 7)
      j = j + 1
    }
    line += "\n"
    i = i + 1
  }
  line
}
