///|
/// Needleman-Wunsch global sequence alignment via dynamic programming.
///
/// References:
/// - Needleman SB, Wunsch CD (1970). "A general method applicable to the
///   search for similarities in the amino acid sequence of two proteins."
///   J Mol Biol. 48(3):443-53.
/// - Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis.
///   Cambridge University Press, 1998. Chapter 2.

///|
/// Traceback direction codes for the NW DP matrix.
///
/// Unlike Smith-Waterman, there is no `Stop` variant because global
/// alignment always traces back from (n, m) to (0, 0).
pub enum NWTrace {
  Diag // match/mismatch: comes from (i-1, j-1)
  Up // gap in seq2 (deletion): comes from (i-1, j)
  Left // gap in seq1 (insertion): comes from (i, j-1)
  Start // only used at (0, 0) to terminate traceback
} derive(Eq, Show)

///|
/// Scoring scheme for global alignment.
///
/// Reuses the same structure concept as Smith-Waterman's `Scoring`
/// but defined independently to avoid cross-module coupling.
pub struct NWScoring {
  match_reward : Double
  mismatch_penalty : Double
  gap_penalty : Double
}

///|
/// Create a scoring scheme.
pub fn NWScoring::new(
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> NWScoring {
  NWScoring::{ match_reward, mismatch_penalty, gap_penalty }
}

///|
/// Default DNA scoring: match=2, mismatch=-1, gap=-1.
pub fn NWScoring::default_dna() -> NWScoring {
  NWScoring::{ match_reward: 2.0, mismatch_penalty: -1.0, gap_penalty: -1.0 }
}

///|
/// Default protein scoring: match=5, mismatch=-4, gap=-3.
pub fn NWScoring::default_protein() -> NWScoring {
  NWScoring::{ match_reward: 5.0, mismatch_penalty: -4.0, gap_penalty: -3.0 }
}

///|
/// Full result of a Needleman-Wunsch alignment.
pub struct NWResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  score_matrix : Array[Array[Double]]
  trace_matrix : Array[Array[NWTrace]]
}

///|
/// Get aligned sequence 1.
pub fn NWResult::aligned1(self : NWResult) -> String {
  self.aligned_seq1
}

///|
/// Get aligned sequence 2.
pub fn NWResult::aligned2(self : NWResult) -> String {
  self.aligned_seq2
}

///|
/// Get match line ('|' for match, '.' for mismatch, ' ' for gap).
pub fn NWResult::match_line(self : NWResult) -> String {
  self.match_line
}

///|
/// Get alignment score.
pub fn NWResult::score(self : NWResult) -> Double {
  self.score
}

///|
/// Get the full score matrix for inspection.
pub fn NWResult::score_matrix(self : NWResult) -> Array[Array[Double]] {
  self.score_matrix
}

///|
/// Get the full traceback matrix for inspection.
pub fn NWResult::trace_matrix(self : NWResult) -> Array[Array[NWTrace]] {
  self.trace_matrix
}

///|
/// Get alignment length (both aligned sequences have the same length).
pub fn NWResult::alignment_length(self : NWResult) -> Int {
  self.aligned_seq1.length()
}

///|
/// Count the number of matches in the alignment.
pub fn NWResult::num_matches(self : NWResult) -> Int {
  let mut count = 0
  for i = 0; i < self.match_line.length(); i = i + 1 {
    if self.match_line.unsafe_get(i).unsafe_to_char() == '|' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of mismatches (excluding gaps).
pub fn NWResult::num_mismatches(self : NWResult) -> Int {
  let mut count = 0
  for i = 0; i < self.match_line.length(); i = i + 1 {
    if self.match_line.unsafe_get(i).unsafe_to_char() == '.' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of gap positions in sequence 1.
pub fn NWResult::num_gaps_seq1(self : NWResult) -> Int {
  let mut count = 0
  for i = 0; i < self.aligned_seq1.length(); i = i + 1 {
    if self.aligned_seq1.unsafe_get(i).unsafe_to_char() == '-' {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of gap positions in sequence 2.
pub fn NWResult::num_gaps_seq2(self : NWResult) -> Int {
  let mut count = 0
  for i = 0; i < self.aligned_seq2.length(); i = i + 1 {
    if self.aligned_seq2.unsafe_get(i).unsafe_to_char() == '-' {
      count = count + 1
    }
  }
  count
}

///|
/// Calculate percent identity (0.0 – 100.0).
pub fn NWResult::percent_identity(self : NWResult) -> Double {
  let aln_len = self.alignment_length()
  if aln_len == 0 {
    return 0.0
  }
  self.num_matches().to_double() / aln_len.to_double() * 100.0
}

///|
/// Format the alignment as a readable block (3 lines).
pub fn NWResult::format(self : NWResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string(self.aligned_seq1)
  buf.write_char('\n')
  buf.write_string(self.match_line)
  buf.write_char('\n')
  buf.write_string(self.aligned_seq2)
  buf.to_string()
}

///|
/// Format the alignment with line wrapping for long sequences.
/// `width` — number of characters per line.
pub fn NWResult::format_wrapped(self : NWResult, width : Int) -> String {
  let buf = StringBuilder::new()
  let total = self.aligned_seq1.length()
  let mut start = 0
  let mut block = 1
  while start < total {
    let end = if start + width > total { total } else { start + width }
    buf.write_string("Block ")
    buf.write_string(block.to_string())
    buf.write_string(" (positions ")
    buf.write_string(start.to_string())
    buf.write_char('-')
    buf.write_string((end - 1).to_string())
    buf.write_string("):\n")
    buf.write_string(self.aligned_seq1[start:end].to_owned())
    buf.write_char('\n')
    buf.write_string(self.match_line[start:end].to_owned())
    buf.write_char('\n')
    buf.write_string(self.aligned_seq2[start:end].to_owned())
    buf.write_char('\n')
    if end < total {
      buf.write_char('\n')
    }
    start = end
    block = block + 1
  }
  buf.to_string()
}

///|
/// Run Needleman-Wunsch global alignment.
///
/// `seq1` — first sequence (string).
/// `seq2` — second sequence (string).
/// `scoring` — scoring scheme (use `NWScoring::default_dna()` etc.).
///
/// Returns an `NWResult` containing aligned sequences, score,
/// and the full DP matrices.
pub fn needleman_wunsch(
  seq1 : String,
  seq2 : String,
  scoring : NWScoring,
) -> NWResult {
  let n = seq1.length()
  let m = seq2.length()
  // Build score matrix: (n+1) x (m+1).
  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[NWTrace]] = Array::new()
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array::new()
    let trow : Array[NWTrace] = Array::new()
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(NWTrace::Start)
    }
    dp.push(drow)
    trace.push(trow)
  }
  // Initialise first row and column (gap penalties).
  // dp[0][0] = 0, trace = Start.
  for i = 1; i <= n; i = i + 1 {
    dp[i][0] = i.to_double() * scoring.gap_penalty
    trace[i][0] = NWTrace::Up
  }
  for j = 1; j <= m; j = j + 1 {
    dp[0][j] = j.to_double() * scoring.gap_penalty
    trace[0][j] = NWTrace::Left
  }
  // Fill DP matrix.
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = seq1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = seq2.unsafe_get(j - 1).unsafe_to_char()
      let sub_score = if c1 == c2 {
        scoring.match_reward
      } else {
        scoring.mismatch_penalty
      }
      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + scoring.gap_penalty
      let left = dp[i][j - 1] + scoring.gap_penalty
      // Choose the maximum (no 0-flooring — that's the NW vs SW difference).
      let (val, dir) = if diag >= up && diag >= left {
        (diag, NWTrace::Diag)
      } else if up >= left {
        (up, NWTrace::Up)
      } else {
        (left, NWTrace::Left)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }
  // Traceback from (n, m) to (0, 0).
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = n
  let mut j = m
  while i > 0 || j > 0 {
    match trace[i][j] {
      NWTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      NWTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      NWTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      NWTrace::Start => break
    }
  }
  let aligned1 = nw_reverse_string(r1.to_string())
  let aligned2 = nw_reverse_string(r2.to_string())
  let match_line = nw_build_match_line(aligned1, aligned2)
  NWResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: dp[n][m],
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Needleman-Wunsch with default DNA scoring.
pub fn needleman_wunsch_dna(seq1 : String, seq2 : String) -> NWResult {
  needleman_wunsch(seq1, seq2, NWScoring::default_dna())
}

///|
/// Convenience: run Needleman-Wunsch with default protein scoring.
pub fn needleman_wunsch_protein(seq1 : String, seq2 : String) -> NWResult {
  needleman_wunsch(seq1, seq2, NWScoring::default_protein())
}

///|
/// Run Needleman-Wunsch global alignment using a substitution matrix.
///
/// `seq1`, `seq2` — protein sequences (single-letter code strings).
/// `matrix` — substitution matrix (e.g., `blosum62()`).
/// `gap_penalty` — penalty for a gap (typically negative).
///
/// Returns an `NWResult` containing aligned sequences, score,
/// and the full DP matrices.
pub fn needleman_wunsch_with_matrix(
  seq1 : String,
  seq2 : String,
  matrix : SubstitutionMatrix,
  gap_penalty : Double,
) -> NWResult {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[NWTrace]] = Array::new()
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array::new()
    let trow : Array[NWTrace] = Array::new()
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(NWTrace::Start)
    }
    dp.push(drow)
    trace.push(trow)
  }
  for i = 1; i <= n; i = i + 1 {
    dp[i][0] = i.to_double() * gap_penalty
    trace[i][0] = NWTrace::Up
  }
  for j = 1; j <= m; j = j + 1 {
    dp[0][j] = j.to_double() * gap_penalty
    trace[0][j] = NWTrace::Left
  }
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = seq1.unsafe_get(i - 1).unsafe_to_char().to_string()
      let c2 = seq2.unsafe_get(j - 1).unsafe_to_char().to_string()
      let sub_score = matrix.get_score_case_insensitive(c1, c2).to_double()
      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + gap_penalty
      let left = dp[i][j - 1] + gap_penalty
      let (val, dir) = if diag >= up && diag >= left {
        (diag, NWTrace::Diag)
      } else if up >= left {
        (up, NWTrace::Up)
      } else {
        (left, NWTrace::Left)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = n
  let mut j = m
  while i > 0 || j > 0 {
    match trace[i][j] {
      NWTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      NWTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      NWTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
      NWTrace::Start => break
    }
  }
  let aligned1 = nw_reverse_string(r1.to_string())
  let aligned2 = nw_reverse_string(r2.to_string())
  let match_line = nw_build_match_line(aligned1, aligned2)
  NWResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: dp[n][m],
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Needleman-Wunsch with BLOSUM62 matrix.
pub fn needleman_wunsch_blosum62(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> NWResult {
  needleman_wunsch_with_matrix(seq1, seq2, blosum62(), gap_penalty)
}

///|
/// Convenience: run Needleman-Wunsch with PAM250 matrix.
pub fn needleman_wunsch_pam250(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> NWResult {
  needleman_wunsch_with_matrix(seq1, seq2, pam250(), gap_penalty)
}

///|
/// Build the match line: '|' for match, '.' for mismatch, ' ' for gap.
fn nw_build_match_line(s1 : String, s2 : String) -> String {
  let buf = StringBuilder::new()
  for i = 0; i < s1.length(); i = i + 1 {
    let c1 = s1.unsafe_get(i).unsafe_to_char()
    let c2 = s2.unsafe_get(i).unsafe_to_char()
    if c1 == '-' || c2 == '-' {
      buf.write_char(' ')
    } else if c1 == c2 {
      buf.write_char('|')
    } else {
      buf.write_char('.')
    }
  }
  buf.to_string()
}

///|
/// Reverse a string (local to avoid clashing with alignment.mbt).
fn nw_reverse_string(s : String) -> String {
  let n = s.length()
  let result = StringBuilder::new()
  let mut i = n - 1
  while i >= 0 {
    result.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  result.to_string()
}
