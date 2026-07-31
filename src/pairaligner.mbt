///|
/// pairaligner module - Unified pairwise sequence aligner.
/// Inspired by Biopython Bio.Align.PairwiseAligner.
///
/// Supports:
///   - Global (Needleman-Wunsch) and Local (Smith-Waterman) alignment modes
///   - Affine gap penalties (gap_open + gap_extend)
///   - DNA and PROTEIN alphabet types
///   - BLOSUM62 substitution matrix for proteins
///
/// References:
///   - Needleman SB, Wunsch CD (1970). J Mol Biol. 48(3):443-53.
///   - Smith TF, Waterman MS (1981). J Mol Biol. 147(1):195-7.
///   - Gotoh O (1982). An improved algorithm for matching biological sequences.
///     J Mol Biol. 162(3):705-8. (Affine gaps)

// ============================================================
// Enums
// ============================================================

///|
/// Alignment mode: Global (Needleman-Wunsch) or Local (Smith-Waterman).
pub enum AlignmentMode {
  Global
  Local
} derive(Eq, Debug)

///|
/// Create a Global alignment mode.
pub fn pairaligner_global() -> AlignmentMode {
  AlignmentMode::Global
}

///|
/// Create a Local alignment mode.
pub fn pairaligner_local() -> AlignmentMode {
  AlignmentMode::Local
}

///|
/// Substitution matrix choice for protein alignment.
pub enum SubstitutionMatrixChoice {
  NoSubstMatrix
  Blosum62
} derive(Eq, Debug)

///|
/// Create NoSubstMatrix choice (use match/mismatch scores).
pub fn pairaligner_no_matrix() -> SubstitutionMatrixChoice {
  SubstitutionMatrixChoice::NoSubstMatrix
}

///|
/// Create Blosum62 choice.
pub fn pairaligner_blosum62() -> SubstitutionMatrixChoice {
  SubstitutionMatrixChoice::Blosum62
}

// ============================================================
// PairwiseAlignerConfig
// ============================================================

///|
/// Configuration for the pairwise aligner.
///
/// Fields:
///   - mode: Global or Local alignment
///   - match_score: Score for identical characters (DNA/no-matrix)
///   - mismatch_score: Score for different characters (DNA/no-matrix)
///   - gap_open: Penalty for opening a gap (typically negative)
///   - gap_extend: Penalty for extending an existing gap (typically negative, >= gap_open)
///   - target_gap_open/_extend: Optional separate penalties for target gaps
///   - query_gap_open/_extend: Optional separate penalties for query gaps
///   - submatrix: Substitution matrix choice
///   - alphabet_type: "DNA" or "PROTEIN"
pub struct PairwiseAlignerConfig {
  mode : AlignmentMode
  match_score : Double
  mismatch_score : Double
  gap_open : Double
  gap_extend : Double
  target_gap_open : Double
  target_gap_extend : Double
  query_gap_open : Double
  query_gap_extend : Double
  submatrix : SubstitutionMatrixChoice
  alphabet_type : String
}

///|
/// Create a new PairwiseAlignerConfig.
/// Note: Use named argument syntax at call sites.
pub fn PairwiseAlignerConfig::new(
  mode~ : AlignmentMode,
  match_score~ : Double,
  mismatch_score~ : Double,
  gap_open~ : Double,
  gap_extend~ : Double,
  submatrix~ : SubstitutionMatrixChoice,
  alphabet_type~ : String,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode,
    match_score,
    mismatch_score,
    gap_open,
    gap_extend,
    target_gap_open: gap_open,
    target_gap_extend: gap_extend,
    query_gap_open: gap_open,
    query_gap_extend: gap_extend,
    submatrix,
    alphabet_type,
  }
}

///|
/// Default DNA config: Global mode, match=1, mismatch=-1, gap_open=-1, gap_extend=-0.5.
pub fn PairwiseAlignerConfig::default_dna() -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: AlignmentMode::Global,
    match_score: 1.0,
    mismatch_score: -1.0,
    gap_open: -1.0,
    gap_extend: -0.5,
    target_gap_open: -1.0,
    target_gap_extend: -0.5,
    query_gap_open: -1.0,
    query_gap_extend: -0.5,
    submatrix: SubstitutionMatrixChoice::NoSubstMatrix,
    alphabet_type: "DNA",
  }
}

///|
/// Default PROTEIN config: Global mode, BLOSUM62, gap_open=-10, gap_extend=-0.5.
pub fn PairwiseAlignerConfig::default_protein() -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: AlignmentMode::Global,
    match_score: 5.0,
    mismatch_score: -4.0,
    gap_open: -10.0,
    gap_extend: -0.5,
    target_gap_open: -10.0,
    target_gap_extend: -0.5,
    query_gap_open: -10.0,
    query_gap_extend: -0.5,
    submatrix: SubstitutionMatrixChoice::Blosum62,
    alphabet_type: "PROTEIN",
  }
}

// ============================================================
// Setter methods (immutable pattern, no ...self spread)
// ============================================================

///|
/// Set alignment mode, returns new config.
pub fn PairwiseAlignerConfig::set_mode(
  self : PairwiseAlignerConfig,
  mode : AlignmentMode,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode,
    match_score: self.match_score,
    mismatch_score: self.mismatch_score,
    gap_open: self.gap_open,
    gap_extend: self.gap_extend,
    target_gap_open: self.target_gap_open,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: self.query_gap_open,
    query_gap_extend: self.query_gap_extend,
    submatrix: self.submatrix,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set match score, returns new config.
pub fn PairwiseAlignerConfig::set_match_score(
  self : PairwiseAlignerConfig,
  val : Double,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: val,
    mismatch_score: self.mismatch_score,
    gap_open: self.gap_open,
    gap_extend: self.gap_extend,
    target_gap_open: self.target_gap_open,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: self.query_gap_open,
    query_gap_extend: self.query_gap_extend,
    submatrix: self.submatrix,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set mismatch score, returns new config.
pub fn PairwiseAlignerConfig::set_mismatch_score(
  self : PairwiseAlignerConfig,
  val : Double,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: self.match_score,
    mismatch_score: val,
    gap_open: self.gap_open,
    gap_extend: self.gap_extend,
    target_gap_open: self.target_gap_open,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: self.query_gap_open,
    query_gap_extend: self.query_gap_extend,
    submatrix: self.submatrix,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set gap open penalty, returns new config.
/// Also sets target_gap_open and query_gap_open to the same value.
pub fn PairwiseAlignerConfig::set_gap_open(
  self : PairwiseAlignerConfig,
  val : Double,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: self.match_score,
    mismatch_score: self.mismatch_score,
    gap_open: val,
    gap_extend: self.gap_extend,
    target_gap_open: val,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: val,
    query_gap_extend: self.query_gap_extend,
    submatrix: self.submatrix,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set gap extend penalty, returns new config.
/// Also sets target_gap_extend and query_gap_extend to the same value.
pub fn PairwiseAlignerConfig::set_gap_extend(
  self : PairwiseAlignerConfig,
  val : Double,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: self.match_score,
    mismatch_score: self.mismatch_score,
    gap_open: self.gap_open,
    gap_extend: val,
    target_gap_open: self.target_gap_open,
    target_gap_extend: val,
    query_gap_open: self.query_gap_open,
    query_gap_extend: val,
    submatrix: self.submatrix,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set substitution matrix choice, returns new config.
pub fn PairwiseAlignerConfig::set_submatrix(
  self : PairwiseAlignerConfig,
  mat : SubstitutionMatrixChoice,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: self.match_score,
    mismatch_score: self.mismatch_score,
    gap_open: self.gap_open,
    gap_extend: self.gap_extend,
    target_gap_open: self.target_gap_open,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: self.query_gap_open,
    query_gap_extend: self.query_gap_extend,
    submatrix: mat,
    alphabet_type: self.alphabet_type,
  }
}

///|
/// Set alphabet type ("DNA" or "PROTEIN"), returns new config.
pub fn PairwiseAlignerConfig::set_alphabet_type(
  self : PairwiseAlignerConfig,
  val : String,
) -> PairwiseAlignerConfig {
  PairwiseAlignerConfig::{
    mode: self.mode,
    match_score: self.match_score,
    mismatch_score: self.mismatch_score,
    gap_open: self.gap_open,
    gap_extend: self.gap_extend,
    target_gap_open: self.target_gap_open,
    target_gap_extend: self.target_gap_extend,
    query_gap_open: self.query_gap_open,
    query_gap_extend: self.query_gap_extend,
    submatrix: self.submatrix,
    alphabet_type: val,
  }
}

// ============================================================
// PairwiseAlignment
// ============================================================

///|
/// Result of a pairwise sequence alignment.
pub struct PairwiseAlignment {
  target : String
  query : String
  score : Double
  aligned_target : String
  aligned_query : String
  match_line : String
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
}

///|
/// Get aligned target sequence (aligned1).
pub fn PairwiseAlignment::aligned1(self : PairwiseAlignment) -> String {
  self.aligned_target
}

///|
/// Get aligned query sequence (aligned2).
pub fn PairwiseAlignment::aligned2(self : PairwiseAlignment) -> String {
  self.aligned_query
}

///|
/// Get alignment score.
pub fn PairwiseAlignment::score(self : PairwiseAlignment) -> Double {
  self.score
}

///|
/// Get alignment length.
pub fn PairwiseAlignment::alignment_length(self : PairwiseAlignment) -> Int {
  self.aligned_target.length()
}

///|
/// Count the number of identical positions in the alignment.
pub fn PairwiseAlignment::identities(self : PairwiseAlignment) -> Int {
  let mut count = 0
  let ml = self.match_line
  let mut i = 0
  while i < ml.length() {
    if ml.unsafe_get(i).unsafe_to_char() == '|' {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Percent identity (0.0 - 100.0).
pub fn PairwiseAlignment::identity_pct(self : PairwiseAlignment) -> Double {
  let alen = self.alignment_length()
  if alen == 0 {
    return 0.0
  }
  self.identities().to_double() / alen.to_double() * 100.0
}

///|
/// Count the total number of gap positions in both sequences.
pub fn PairwiseAlignment::gaps_count(self : PairwiseAlignment) -> Int {
  let mut count = 0
  let at = self.aligned_target
  let aq = self.aligned_query
  let mut i = 0
  while i < at.length() {
    if at.unsafe_get(i).unsafe_to_char() == '-' {
      count = count + 1
    }
    i = i + 1
  }
  let mut j = 0
  while j < aq.length() {
    if aq.unsafe_get(j).unsafe_to_char() == '-' {
      count = count + 1
    }
    j = j + 1
  }
  count
}

// ============================================================
// Internal helpers: case-insensitive comparison
// ============================================================

///|
/// Case-insensitive char comparison.
/// Convert chars to 1-char Strings and use String.to_upper().
fn pairaligner_chars_equal_ci(a : Char, b : Char) -> Bool {
  let sa = a.to_string().to_upper()
  let sb = b.to_string().to_upper()
  sa == sb
}

// ============================================================
// BLOSUM62 Internal Matrix
// ============================================================

///|
/// BLOSUM62 data for 20 standard amino acids in order:
/// A R N D C Q E G H I L K M F P S T W Y V
let pairaligner_blosum62_data : Array[Array[Int]] = [
  [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0], // A
  [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3], // R
  [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3], // N
  [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3], // D
  [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1], // C
  [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2], // Q
  [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2], // E
  [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3], // G
  [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3], // H
  [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3], // I
  [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1], // L
  [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2], // K
  [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1], // M
  [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1], // F
  [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2], // P
  [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2], // S
  [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0], // T
  [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3], // W
  [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1], // Y
  [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4], // V
]

///|
/// Internal BLOSUM62 matrix lookup using 20 standard amino acids.
/// Returns score for a pair (case-insensitive).
/// Uses mismatch_score default for unknown amino acids.
fn pairaligner_blosum62_score(a : Char, b : Char, mismatch_default : Double) -> Double {
  let ia = pairaligner_aa_index_ci(a)
  let ib = pairaligner_aa_index_ci(b)
  if ia < 0 || ib < 0 {
    return mismatch_default
  }
  pairaligner_blosum62_data[ia][ib].to_double()
}

///|
/// Convert amino acid char to index 0-19, -1 if unknown (case-insensitive).
fn pairaligner_aa_index_ci(c : Char) -> Int {
  let up = c.to_string().to_upper()
  if up.length() == 0 {
    return -1
  }
  match up.unsafe_get(0).unsafe_to_char() {
    'A' => 0
    'R' => 1
    'N' => 2
    'D' => 3
    'C' => 4
    'Q' => 5
    'E' => 6
    'G' => 7
    'H' => 8
    'I' => 9
    'L' => 10
    'K' => 11
    'M' => 12
    'F' => 13
    'P' => 14
    'S' => 15
    'T' => 16
    'W' => 17
    'Y' => 18
    'V' => 19
    _ => -1
  }
}

// ============================================================
// Substitution score helpers
// ============================================================

///|
/// Get substitution score for two characters based on config.
fn pairaligner_get_subst_score(
  a : Char,
  b : Char,
  config : PairwiseAlignerConfig,
) -> Double {
  let is_protein = config.alphabet_type == "PROTEIN"
  let use_matrix = is_protein && config.submatrix == SubstitutionMatrixChoice::Blosum62
  if use_matrix {
    pairaligner_blosum62_score(a, b, config.mismatch_score)
  } else {
    if pairaligner_chars_equal_ci(a, b) {
      config.match_score
    } else {
      config.mismatch_score
    }
  }
}

// ============================================================
// Match line builder
// ============================================================

///|
/// Build the match line: '|' match, '.' mismatch, ' ' gap.
fn pairaligner_build_match_line(s1 : String, s2 : String) -> String {
  let buf = StringBuilder::new()
  let mut i = 0
  while i < s1.length() {
    let c1 = s1.unsafe_get(i).unsafe_to_char()
    let c2 = s2.unsafe_get(i).unsafe_to_char()
    if c1 == '-' || c2 == '-' {
      buf.write_char(' ')
    } else if pairaligner_chars_equal_ci(c1, c2) {
      buf.write_char('|')
    } else {
      buf.write_char('.')
    }
    i = i + 1
  }
  buf.to_string()
}

///|
/// Reverse a string (local to pairaligner to avoid name clashes).
fn pairaligner_reverse_string(s : String) -> String {
  let n = s.length()
  let result = StringBuilder::new()
  let mut i = n - 1
  while i >= 0 {
    result.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  result.to_string()
}

// ============================================================
// Traceback direction enum (internal)
// ============================================================

priv enum TraceDir {
  Stop
  Diag
  Up
  Left
}

// ============================================================
// Main alignment function
// ============================================================

///|
/// Run pairwise sequence alignment.
///
/// Uses Gotoh-style affine gap penalties with three matrices:
///   - m_mat:  best score ending with a match/mismatch (diagonal step)
///   - ix_mat: best score ending with a gap in query (Up step, deletion)
///   - iy_mat: best score ending with a gap in target (Left step, insertion)
///
/// `target` — first sequence.
/// `query` — second sequence.
/// `config?` — optional PairwiseAlignerConfig (default: default_dna()).
pub fn pairaligner_align(
  target : String,
  query : String,
  config? : PairwiseAlignerConfig = PairwiseAlignerConfig::default_dna(),
) -> PairwiseAlignment {
  let n = target.length()
  let m = query.length()
  let is_local = config.mode == AlignmentMode::Local

  // ---- Initialize DP matrices ----
  let m_mat  : Array[Array[Double]]  = pairaligner_make_matrix(n + 1, m + 1, 0.0)
  let ix_mat : Array[Array[Double]]  = pairaligner_make_matrix(n + 1, m + 1, 0.0)
  let iy_mat : Array[Array[Double]]  = pairaligner_make_matrix(n + 1, m + 1, 0.0)
  let trace  : Array[Array[TraceDir]] = pairaligner_make_trace_matrix(n + 1, m + 1)

  let target_gap_open = config.target_gap_open
  let target_gap_extend = config.target_gap_extend
  let query_gap_open = config.query_gap_open
  let query_gap_extend = config.query_gap_extend

  let neg_inf : Double = -1000000.0

  // ---- Boundary conditions ----
  if is_local {
    // For local alignment, boundaries stay 0 and trace is Stop
    let mut i = 0
    while i <= n {
      trace[i][0] = TraceDir::Stop
      i = i + 1
    }
    let mut j = 0
    while j <= m {
      trace[0][j] = TraceDir::Stop
      j = j + 1
    }
  } else {
    // Global alignment: first row / column
    m_mat[0][0] = 0.0
    ix_mat[0][0] = neg_inf
    iy_mat[0][0] = neg_inf
    trace[0][0] = TraceDir::Stop

    let mut i = 1
    while i <= n {
      if i == 1 {
        ix_mat[1][0] = query_gap_open
      } else {
        ix_mat[i][0] = ix_mat[i - 1][0] + query_gap_extend
      }
      m_mat[i][0] = neg_inf
      iy_mat[i][0] = neg_inf
      trace[i][0] = TraceDir::Up
      i = i + 1
    }
    let mut j = 1
    while j <= m {
      if j == 1 {
        iy_mat[0][1] = target_gap_open
      } else {
        iy_mat[0][j] = iy_mat[0][j - 1] + target_gap_extend
      }
      m_mat[0][j] = neg_inf
      ix_mat[0][j] = neg_inf
      trace[0][j] = TraceDir::Left
      j = j + 1
    }
  }

  // ---- Fill DP matrices ----
  let mut max_score : Double = if is_local { 0.0 } else { neg_inf }
  let mut max_i = 0
  let mut max_j = 0

  let mut i = 1
  while i <= n {
    let mut j = 1
    while j <= m {
      let tc = target.unsafe_get(i - 1).unsafe_to_char()
      let qc = query.unsafe_get(j - 1).unsafe_to_char()
      let s = pairaligner_get_subst_score(tc, qc, config)

      // Compute m_mat[i][j]: match/mismatch
      let from_M  = m_mat[i - 1][j - 1]  + s
      let from_Ix = ix_mat[i - 1][j - 1] + s
      let from_Iy = iy_mat[i - 1][j - 1] + s
      let mut m_val = from_M
      if from_Ix > m_val { m_val = from_Ix }
      if from_Iy > m_val { m_val = from_Iy }

      // Compute ix_mat[i][j]: gap in query (Up step)
      let ix_extend       = ix_mat[i - 1][j] + query_gap_extend
      let ix_open_from_M  = m_mat[i - 1][j]  + query_gap_open
      let ix_open_from_Iy = iy_mat[i - 1][j] + query_gap_open
      let mut ix_val = ix_extend
      if ix_open_from_M  > ix_val { ix_val = ix_open_from_M  }
      if ix_open_from_Iy > ix_val { ix_val = ix_open_from_Iy }

      // Compute iy_mat[i][j]: gap in target (Left step)
      let iy_extend       = iy_mat[i][j - 1] + target_gap_extend
      let iy_open_from_M  = m_mat[i][j - 1]  + target_gap_open
      let iy_open_from_Ix = ix_mat[i][j - 1] + target_gap_open
      let mut iy_val = iy_extend
      if iy_open_from_M  > iy_val { iy_val = iy_open_from_M  }
      if iy_open_from_Ix > iy_val { iy_val = iy_open_from_Ix }

      if is_local {
        if m_val  < 0.0 { m_val  = 0.0 }
        if ix_val < 0.0 { ix_val = 0.0 }
        if iy_val < 0.0 { iy_val = 0.0 }
      }

      // Choose overall best and record trace
      let mut best_val = m_val
      let mut best_dir = TraceDir::Diag
      if ix_val > best_val {
        best_val = ix_val
        best_dir = TraceDir::Up
      }
      if iy_val > best_val {
        best_val = iy_val
        best_dir = TraceDir::Left
      }

      m_mat[i][j]  = m_val
      ix_mat[i][j] = ix_val
      iy_mat[i][j] = iy_val

      if is_local {
        if best_val == 0.0 {
          trace[i][j] = TraceDir::Stop
        } else {
          trace[i][j] = best_dir
        }
        if best_val > max_score {
          max_score = best_val
          max_i = i
          max_j = j
        }
      } else {
        trace[i][j] = best_dir
      }

      j = j + 1
    }
    i = i + 1
  }

  // ---- Determine traceback start ----
  let (start_i, start_j, final_score) : (Int, Int, Double) = if is_local {
    (max_i, max_j, max_score)
  } else {
    let m_nm  = m_mat[n][m]
    let ix_nm = ix_mat[n][m]
    let iy_nm = iy_mat[n][m]
    let mut best = m_nm
    if ix_nm > best { best = ix_nm }
    if iy_nm > best { best = iy_nm }
    (n, m, best)
  }

  // ---- Traceback ----
  let r_target = StringBuilder::new()
  let r_query = StringBuilder::new()
  let mut i_pos = start_i
  let mut j_pos = start_j

  if is_local {
    while i_pos > 0 && j_pos > 0 {
      match trace[i_pos][j_pos] {
        TraceDir::Stop => break
        TraceDir::Diag => {
          r_target.write_char(target.unsafe_get(i_pos - 1).unsafe_to_char())
          r_query.write_char(query.unsafe_get(j_pos - 1).unsafe_to_char())
          i_pos = i_pos - 1
          j_pos = j_pos - 1
        }
        TraceDir::Up => {
          r_target.write_char(target.unsafe_get(i_pos - 1).unsafe_to_char())
          r_query.write_char('-')
          i_pos = i_pos - 1
        }
        TraceDir::Left => {
          r_target.write_char('-')
          r_query.write_char(query.unsafe_get(j_pos - 1).unsafe_to_char())
          j_pos = j_pos - 1
        }
      }
    }
  } else {
    while i_pos > 0 || j_pos > 0 {
      if i_pos == 0 {
        r_target.write_char('-')
        r_query.write_char(query.unsafe_get(j_pos - 1).unsafe_to_char())
        j_pos = j_pos - 1
        continue
      }
      if j_pos == 0 {
        r_target.write_char(target.unsafe_get(i_pos - 1).unsafe_to_char())
        r_query.write_char('-')
        i_pos = i_pos - 1
        continue
      }
      match trace[i_pos][j_pos] {
        TraceDir::Stop => break
        TraceDir::Diag => {
          r_target.write_char(target.unsafe_get(i_pos - 1).unsafe_to_char())
          r_query.write_char(query.unsafe_get(j_pos - 1).unsafe_to_char())
          i_pos = i_pos - 1
          j_pos = j_pos - 1
        }
        TraceDir::Up => {
          r_target.write_char(target.unsafe_get(i_pos - 1).unsafe_to_char())
          r_query.write_char('-')
          i_pos = i_pos - 1
        }
        TraceDir::Left => {
          r_target.write_char('-')
          r_query.write_char(query.unsafe_get(j_pos - 1).unsafe_to_char())
          j_pos = j_pos - 1
        }
      }
    }
  }

  let aligned_target = pairaligner_reverse_string(r_target.to_string())
  let aligned_query = pairaligner_reverse_string(r_query.to_string())
  let match_line = pairaligner_build_match_line(aligned_target, aligned_query)

  let target_start = i_pos
  let target_end = start_i
  let query_start = j_pos
  let query_end = start_j

  PairwiseAlignment::{
    target,
    query,
    score: final_score,
    aligned_target,
    aligned_query,
    match_line,
    target_start,
    target_end,
    query_start,
    query_end,
  }
}

// ============================================================
// Matrix helpers
// ============================================================

///|
/// Create a 2D Double matrix of size rows x cols filled with val.
fn pairaligner_make_matrix(rows : Int, cols : Int, val : Double) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < rows {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < cols {
      row.push(val)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Create a 2D trace matrix filled with Stop.
fn pairaligner_make_trace_matrix(rows : Int, cols : Int) -> Array[Array[TraceDir]] {
  let result : Array[Array[TraceDir]] = Array::new()
  let mut i = 0
  while i < rows {
    let row : Array[TraceDir] = Array::new()
    let mut j = 0
    while j < cols {
      row.push(TraceDir::Stop)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

// ============================================================
// Sample data
// ============================================================

///|
/// Return a pair of sample test sequences (protein - HEAGAWGHEE vs PAWHEAE).
/// Returns (target, query).
pub fn pairaligner_sample_data() -> (String, String) {
  ("HEAGAWGHEE".to_string(), "PAWHEAE".to_string())
}
