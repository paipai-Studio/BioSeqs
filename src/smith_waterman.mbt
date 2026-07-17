///|
/// Smith-Waterman local sequence alignment via dynamic programming.
///
/// References:
/// - Smith TF, Waterman MS (1981). "Identification of common molecular
///   subsequences." J Mol Biol. 147(1):195-7.
/// - Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis.
///   Cambridge University Press, 1998. Chapter 2.

///|
/// Traceback direction codes for the DP matrix.
pub enum TraceDirection {
  Stop // score is 0, alignment ends here
  Diag // match/mismatch: comes from (i-1, j-1)
  Up // gap in seq2 (deletion): comes from (i-1, j)
  Left // gap in seq1 (insertion): comes from (i, j-1)
} derive(Eq, Show)

///|
/// Scoring scheme for sequence alignment.
///
/// `match_reward` — score for identical characters.
/// `mismatch_penalty` — score (typically negative) for different characters.
/// `gap_penalty` — score (typically negative) for a gap.
pub struct Scoring {
  match_reward : Double
  mismatch_penalty : Double
  gap_penalty : Double
}

///|
/// Create a simple scoring scheme.
pub fn Scoring::new(
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> Scoring {
  Scoring::{ match_reward, mismatch_penalty, gap_penalty }
}

///|
/// Default DNA scoring: match=2, mismatch=-1, gap=-1.
pub fn Scoring::default_dna() -> Scoring {
  Scoring::{ match_reward: 2.0, mismatch_penalty: -1.0, gap_penalty: -1.0 }
}

///|
/// Default protein scoring: match=5, mismatch=-4, gap=-3
/// (simplified BLOSUM-like scheme without a full substitution matrix).
pub fn Scoring::default_protein() -> Scoring {
  Scoring::{ match_reward: 5.0, mismatch_penalty: -4.0, gap_penalty: -3.0 }
}

///|
/// Full result of a Smith-Waterman alignment.
pub struct SWResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  start_pos1 : (Int, Int) // (start, end) in seq1, 0-based, end exclusive
  start_pos2 : (Int, Int) // (start, end) in seq2
  score_matrix : Array[Array[Double]]
  trace_matrix : Array[Array[TraceDirection]]
}

///|
/// Get aligned sequence 1.
pub fn SWResult::aligned1(self : SWResult) -> String {
  self.aligned_seq1
}

///|
/// Get aligned sequence 2.
pub fn SWResult::aligned2(self : SWResult) -> String {
  self.aligned_seq2
}

///|
/// Get match line ('|' for match, '.' for mismatch, ' ' for gap).
pub fn SWResult::match_line(self : SWResult) -> String {
  self.match_line
}

///|
/// Get alignment score.
pub fn SWResult::score(self : SWResult) -> Double {
  self.score
}

///|
/// Get start/end positions in sequence 1 (0-based, end exclusive).
pub fn SWResult::pos1(self : SWResult) -> (Int, Int) {
  self.start_pos1
}

///|
/// Get start/end positions in sequence 2 (0-based, end exclusive).
pub fn SWResult::pos2(self : SWResult) -> (Int, Int) {
  self.start_pos2
}

///|
/// Get the full score matrix for inspection.
pub fn SWResult::score_matrix(self : SWResult) -> Array[Array[Double]] {
  self.score_matrix
}

///|
/// Get the full traceback matrix for inspection.
pub fn SWResult::trace_matrix(self : SWResult) -> Array[Array[TraceDirection]] {
  self.trace_matrix
}

///|
/// Get alignment length.
pub fn SWResult::alignment_length(self : SWResult) -> Int {
  self.aligned_seq1.length()
}

///|
/// Count the number of matches in the alignment.
pub fn SWResult::num_matches(self : SWResult) -> Int {
  let mut count = 0
  for i = 0; i < self.match_line.length(); i = i + 1 {
    if self.match_line.unsafe_get(i).unsafe_to_char() == '|' {
      count = count + 1
    }
  }
  count
}

///|
/// Calculate percent identity (0.0 – 100.0).
pub fn SWResult::percent_identity(self : SWResult) -> Double {
  let aln_len = self.alignment_length()
  if aln_len == 0 {
    return 0.0
  }
  self.num_matches().to_double() / aln_len.to_double() * 100.0
}

///|
/// Format the alignment as a readable block (3 lines).
pub fn SWResult::format(self : SWResult) -> String {
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
pub fn SWResult::format_wrapped(self : SWResult, width : Int) -> String {
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
/// Run Smith-Waterman local alignment.
///
/// `seq1` — first sequence (string).
/// `seq2` — second sequence (string).
/// `scoring` — scoring scheme (use `Scoring::default_dna()` etc.).
///
/// Returns a `SWResult` containing aligned sequences, score, positions,
/// and the full DP matrices.
pub fn smith_waterman(
  seq1 : String,
  seq2 : String,
  scoring : Scoring,
) -> SWResult {
  let n = seq1.length()
  let m = seq2.length()
  // Build score matrix: (n+1) x (m+1), initialised to 0.
  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[TraceDirection]] = Array::new()
  for i = 0; i <= n; i = i + 1 {
    let row = FixedArray::make(m + 1, 0.0)
    let trow = FixedArray::make(m + 1, TraceDirection::Stop)
    let drow : Array[Double] = Array::new()
    let trrow : Array[TraceDirection] = Array::new()
    for j = 0; j <= m; j = j + 1 {
      drow.push(row[j])
      trrow.push(trow[j])
    }
    dp.push(drow)
    trace.push(trrow)
  }
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0
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
      // Choose the best of diag, up, left, or 0 (Smith-Waterman key step).
      let (val, dir) = if diag >= up && diag >= left && diag > 0.0 {
        (diag, TraceDirection::Diag)
      } else if up >= left && up > 0.0 {
        (up, TraceDirection::Up)
      } else if left > 0.0 {
        (left, TraceDirection::Left)
      } else {
        (0.0, TraceDirection::Stop)
      }
      dp[i][j] = val
      trace[i][j] = dir
      if val > max_score {
        max_score = val
        max_i = i
        max_j = j
      }
    }
  }
  // Traceback from (max_i, max_j) until a Stop cell.
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = max_i
  let mut j = max_j
  while i > 0 && j > 0 {
    match trace[i][j] {
      TraceDirection::Stop => break
      TraceDirection::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      TraceDirection::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      TraceDirection::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  SWResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: max_score,
    start_pos1: (i, max_i),
    start_pos2: (j, max_j),
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Smith-Waterman with default DNA scoring.
pub fn smith_waterman_dna(seq1 : String, seq2 : String) -> SWResult {
  smith_waterman(seq1, seq2, Scoring::default_dna())
}

///|
/// Convenience: run Smith-Waterman with default protein scoring.
pub fn smith_waterman_protein(seq1 : String, seq2 : String) -> SWResult {
  smith_waterman(seq1, seq2, Scoring::default_protein())
}

///|
/// Run Smith-Waterman local alignment using a substitution matrix.
///
/// `seq1`, `seq2` — protein sequences (single-letter code strings).
/// `matrix` — substitution matrix (e.g., `blosum62()`).
/// `gap_penalty` — penalty for a gap (typically negative).
///
/// Returns a `SWResult` containing aligned sequences, score, positions,
/// and the full DP matrices.
pub fn smith_waterman_with_matrix(
  seq1 : String,
  seq2 : String,
  matrix : SubstitutionMatrix,
  gap_penalty : Double,
) -> SWResult {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array::new()
  let trace : Array[Array[TraceDirection]] = Array::new()
  for i = 0; i <= n; i = i + 1 {
    let row = FixedArray::make(m + 1, 0.0)
    let trow = FixedArray::make(m + 1, TraceDirection::Stop)
    let drow : Array[Double] = Array::new()
    let trrow : Array[TraceDirection] = Array::new()
    for j = 0; j <= m; j = j + 1 {
      drow.push(row[j])
      trrow.push(trow[j])
    }
    dp.push(drow)
    trace.push(trrow)
  }
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = seq1.unsafe_get(i - 1).unsafe_to_char().to_string()
      let c2 = seq2.unsafe_get(j - 1).unsafe_to_char().to_string()
      let sub_score = matrix.get_score_case_insensitive(c1, c2).to_double()
      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + gap_penalty
      let left = dp[i][j - 1] + gap_penalty
      let (val, dir) = if diag >= up && diag >= left && diag > 0.0 {
        (diag, TraceDirection::Diag)
      } else if up >= left && up > 0.0 {
        (up, TraceDirection::Up)
      } else if left > 0.0 {
        (left, TraceDirection::Left)
      } else {
        (0.0, TraceDirection::Stop)
      }
      dp[i][j] = val
      trace[i][j] = dir
      if val > max_score {
        max_score = val
        max_i = i
        max_j = j
      }
    }
  }
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()
  let mut i = max_i
  let mut j = max_j
  while i > 0 && j > 0 {
    match trace[i][j] {
      TraceDirection::Stop => break
      TraceDirection::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      TraceDirection::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      TraceDirection::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  SWResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: max_score,
    start_pos1: (i, max_i),
    start_pos2: (j, max_j),
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: run Smith-Waterman with BLOSUM62 matrix.
pub fn smith_waterman_blosum62(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> SWResult {
  smith_waterman_with_matrix(seq1, seq2, blosum62(), gap_penalty)
}

///|
/// Convenience: run Smith-Waterman with PAM250 matrix.
pub fn smith_waterman_pam250(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> SWResult {
  smith_waterman_with_matrix(seq1, seq2, pam250(), gap_penalty)
}

///|
/// Build the match line: '|' for match, '.' for mismatch, ' ' for gap.
fn build_match_line(s1 : String, s2 : String) -> String {
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
/// `reverse_string` is defined in `alignment.mbt` and reused here.
