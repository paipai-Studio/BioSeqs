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
#derive(Eq, Debug)
pub enum TraceDirection {
  Stop // score is 0, alignment ends here
  Diag // match/mismatch: comes from (i-1, j-1)
  Up // gap in seq2 (deletion): comes from (i-1, j)
  Left // gap in seq1 (insertion): comes from (i, j-1)
}

///|
pub extend TraceDirection with Eq::{equal, not_equal}

///|
pub extend TraceDirection with @moonbitlang/core/debug.Debug::{to_repr}

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
  Scoring::{ match_reward, mismatch_penalty, gap_penalty, }
}

///|
/// Default DNA scoring: match=2, mismatch=-1, gap=-1.
pub fn Scoring::default_dna() -> Scoring {
  Scoring::{ match_reward: 2.0, mismatch_penalty: -1.0, gap_penalty: -1.0, }
}

///|
/// Default protein scoring: match=5, mismatch=-4, gap=-3
/// (simplified BLOSUM-like scheme without a full substitution matrix).
pub fn Scoring::default_protein() -> Scoring {
  Scoring::{ match_reward: 5.0, mismatch_penalty: -4.0, gap_penalty: -3.0, }
}

///|
/// Full result of a Smith-Waterman alignment.
pub struct SWResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  start_pos1 : (Int32, Int32) // (start, end) in seq1, 0-based, end exclusive
  start_pos2 : (Int32, Int32) // (start, end) in seq2
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
pub fn SWResult::pos1(self : SWResult) -> (Int32, Int32) {
  self.start_pos1
}

///|
/// Get start/end positions in sequence 2 (0-based, end exclusive).
pub fn SWResult::pos2(self : SWResult) -> (Int32, Int32) {
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
pub fn SWResult::alignment_length(self : SWResult) -> Int32 {
  self.aligned_seq1.length()
}

///|
/// Count the number of matches in the alignment.
pub fn SWResult::num_matches(self : SWResult) -> Int32 {
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
  let buf = StringBuilder()
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
pub fn SWResult::format_wrapped(self : SWResult, width : Int32) -> String {
  let buf = StringBuilder()
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

// ──────────────────────────────────────────────────────────────────────
// Shared alignment engine (used by both Smith-Waterman and Needleman-Wunsch)
// ──────────────────────────────────────────────────────────────────────

///|
/// Alignment mode: local (Smith-Waterman) or global (Needleman-Wunsch).
#derive(Eq)
enum AlignMode {
  Local
  Global
}

///|
pub extend AlignMode with Eq::{equal, not_equal}

///|
/// Unified traceback direction codes used internally by the shared engine.
priv enum AlignTrace {
  Stop // local: score floored to 0; global: only at (0,0)
  Diag // match/mismatch: from (i-1, j-1)
  Up // gap in seq2: from (i-1, j)
  Left // gap in seq1: from (i, j-1)
}

///|
/// Internal alignment output produced by the shared DP engine.
priv struct AlignOutput {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  start1 : Int32 // start index in seq1 (local mode)
  end1 : Int32 // end index in seq1 (local mode)
  start2 : Int32
  end2 : Int32
  score_matrix : Array[Array[Double]]
  trace_matrix : Array[Array[AlignTrace]]
}

///|
/// Convert an AlignTrace code to a Smith-Waterman TraceDirection.
fn align_trace_to_sw(t : AlignTrace) -> TraceDirection {
  match t {
    AlignTrace::Stop => TraceDirection::Stop
    AlignTrace::Diag => TraceDirection::Diag
    AlignTrace::Up => TraceDirection::Up
    AlignTrace::Left => TraceDirection::Left
  }
}

///|
/// Shared DP engine for both simple (match/mismatch) scoring modes.
fn align_core_simple(
  seq1 : String,
  seq2 : String,
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
  mode : AlignMode,
) -> AlignOutput {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array(capacity=0)
  let trace : Array[Array[AlignTrace]] = Array(capacity=0)
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array(capacity=0)
    let trow : Array[AlignTrace] = Array(capacity=0)
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(AlignTrace::Stop)
    }
    dp.push(drow)
    trace.push(trow)
  }
  // Global mode: initialise first row / column with gap penalties.
  if mode == AlignMode::Global {
    for i = 1; i <= n; i = i + 1 {
      dp[i][0] = i.to_double() * gap_penalty
      trace[i][0] = AlignTrace::Up
    }
    for j = 1; j <= m; j = j + 1 {
      dp[0][j] = j.to_double() * gap_penalty
      trace[0][j] = AlignTrace::Left
    }
  }
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0
  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = seq1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = seq2.unsafe_get(j - 1).unsafe_to_char()
      let sub_score = if c1 == c2 { match_reward } else { mismatch_penalty }
      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + gap_penalty
      let left = dp[i][j - 1] + gap_penalty
      let (val, dir) = match mode {
        AlignMode::Local =>
          if diag >= up && diag >= left && diag > 0.0 {
            (diag, AlignTrace::Diag)
          } else if up >= left && up > 0.0 {
            (up, AlignTrace::Up)
          } else if left > 0.0 {
            (left, AlignTrace::Left)
          } else {
            (0.0, AlignTrace::Stop)
          }
        AlignMode::Global =>
          if diag >= up && diag >= left {
            (diag, AlignTrace::Diag)
          } else if up >= left {
            (up, AlignTrace::Up)
          } else {
            (left, AlignTrace::Left)
          }
      }
      dp[i][j] = val
      trace[i][j] = dir
      if mode == AlignMode::Local && val > max_score {
        max_score = val
        max_i = i
        max_j = j
      }
    }
  }
  // Traceback start depends on mode.
  let (start_i, start_j, score) = match mode {
    AlignMode::Local => (max_i, max_j, max_score)
    AlignMode::Global => (n, m, dp[n][m])
  }
  let mut i = start_i
  let mut j = start_j
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  while i > 0 || j > 0 {
    match trace[i][j] {
      AlignTrace::Stop => break
      AlignTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      AlignTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      AlignTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  AlignOutput::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score,
    start1: i,
    end1: match mode {
      AlignMode::Local => max_i
      AlignMode::Global => n
    },
    start2: j,
    end2: match mode {
      AlignMode::Local => max_j
      AlignMode::Global => m
    },
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Shared DP engine for substitution-matrix scoring modes.
fn align_core_matrix(
  seq1 : String,
  seq2 : String,
  matrix : @data.SubstitutionMatrix,
  gap_penalty : Double,
  mode : AlignMode,
) -> AlignOutput {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array(capacity=0)
  let trace : Array[Array[AlignTrace]] = Array(capacity=0)
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array(capacity=0)
    let trow : Array[AlignTrace] = Array(capacity=0)
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(AlignTrace::Stop)
    }
    dp.push(drow)
    trace.push(trow)
  }
  if mode == AlignMode::Global {
    for i = 1; i <= n; i = i + 1 {
      dp[i][0] = i.to_double() * gap_penalty
      trace[i][0] = AlignTrace::Up
    }
    for j = 1; j <= m; j = j + 1 {
      dp[0][j] = j.to_double() * gap_penalty
      trace[0][j] = AlignTrace::Left
    }
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
      let (val, dir) = match mode {
        AlignMode::Local =>
          if diag >= up && diag >= left && diag > 0.0 {
            (diag, AlignTrace::Diag)
          } else if up >= left && up > 0.0 {
            (up, AlignTrace::Up)
          } else if left > 0.0 {
            (left, AlignTrace::Left)
          } else {
            (0.0, AlignTrace::Stop)
          }
        AlignMode::Global =>
          if diag >= up && diag >= left {
            (diag, AlignTrace::Diag)
          } else if up >= left {
            (up, AlignTrace::Up)
          } else {
            (left, AlignTrace::Left)
          }
      }
      dp[i][j] = val
      trace[i][j] = dir
      if mode == AlignMode::Local && val > max_score {
        max_score = val
        max_i = i
        max_j = j
      }
    }
  }
  let (start_i, start_j, score) = match mode {
    AlignMode::Local => (max_i, max_j, max_score)
    AlignMode::Global => (n, m, dp[n][m])
  }
  let mut i = start_i
  let mut j = start_j
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  while i > 0 || j > 0 {
    match trace[i][j] {
      AlignTrace::Stop => break
      AlignTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      AlignTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      AlignTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = reverse_string(r1.to_string())
  let aligned2 = reverse_string(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  AlignOutput::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score,
    start1: i,
    end1: match mode {
      AlignMode::Local => max_i
      AlignMode::Global => n
    },
    start2: j,
    end2: match mode {
      AlignMode::Local => max_j
      AlignMode::Global => m
    },
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Build the match line: '|' for match, '.' for mismatch, ' ' for gap.
fn build_match_line(s1 : String, s2 : String) -> String {
  let buf = StringBuilder()
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
  let out = align_core_simple(
    seq1,
    seq2,
    scoring.match_reward,
    scoring.mismatch_penalty,
    scoring.gap_penalty,
    AlignMode::Local,
  )
  let trace : Array[Array[TraceDirection]] = out.trace_matrix.map(fn(row) {
    row.map(align_trace_to_sw)
  })
  SWResult::{
    aligned_seq1: out.aligned_seq1,
    aligned_seq2: out.aligned_seq2,
    match_line: out.match_line,
    score: out.score,
    start_pos1: (out.start1, out.end1),
    start_pos2: (out.start2, out.end2),
    score_matrix: out.score_matrix,
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
/// `matrix` — substitution matrix (e.g., `@data.blosum62()`).
/// `gap_penalty` — penalty for a gap (typically negative).
///
/// Returns a `SWResult` containing aligned sequences, score, positions,
/// and the full DP matrices.
pub fn smith_waterman_with_matrix(
  seq1 : String,
  seq2 : String,
  matrix : @data.SubstitutionMatrix,
  gap_penalty : Double,
) -> SWResult {
  let out = align_core_matrix(seq1, seq2, matrix, gap_penalty, AlignMode::Local)
  let trace : Array[Array[TraceDirection]] = out.trace_matrix.map(fn(row) {
    row.map(align_trace_to_sw)
  })
  SWResult::{
    aligned_seq1: out.aligned_seq1,
    aligned_seq2: out.aligned_seq2,
    match_line: out.match_line,
    score: out.score,
    start_pos1: (out.start1, out.end1),
    start_pos2: (out.start2, out.end2),
    score_matrix: out.score_matrix,
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
  smith_waterman_with_matrix(seq1, seq2, @data.blosum62(), gap_penalty)
}

///|
/// Convenience: run Smith-Waterman with PAM250 matrix.
pub fn smith_waterman_pam250(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> SWResult {
  smith_waterman_with_matrix(seq1, seq2, @data.pam250(), gap_penalty)
}

///|
/// `reverse_string` is defined in `alignment.mbt` and reused here.
