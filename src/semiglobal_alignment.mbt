///|
/// Semiglobal (glocal) sequence alignment via dynamic programming.
///
/// Semiglobal alignment is "global in the small, local in the large":
/// the entire of both sequences must be consumed, but gaps at the
/// beginning and end of either sequence are free. It is the natural
/// mode for aligning a read to a reference genome (leading/trailing
/// clipping) and for overlap-layout-consensus assembly steps.
///
/// Unlike the banded approximation in `align_banded`
/// (BandedAlignMode::Semiglobal, which was previously approximated by
/// global alignment), this module implements the exact DP:
///   - initialisation: dp[0][j] = dp[i][0] = 0 (free leading gaps)
///   - recurrence: standard match/mismatch/gap (linear gap penalty)
///   - termination: best score along the last row or last column
///     (free trailing gaps)
///
/// References:
/// - Durbin R, Eddy SR, Krogh A, Mitchison G. Biological Sequence Analysis.
///   Cambridge University Press, 1998. Chapter 2.6 exercises.
/// - Gibbons FD. "Semiglobal alignment of sequences." UMD technical note.

///|
/// Traceback direction codes for the semiglobal DP matrix.
#derive(Eq, Debug)
pub enum SGTrace {
  Diag // match/mismatch: comes from (i-1, j-1)
  Up // gap in seq2 (deletion): comes from (i-1, j)
  Left // gap in seq1 (insertion): comes from (i, j-1)
  Start // terminates traceback (reached row 0 or column 0)
}

///|
pub extend SGTrace with Eq::{equal, not_equal}

///|
pub extend SGTrace with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Scoring scheme for semiglobal alignment.
///
/// Same shape as `NWScoring`: positive match reward, negative
/// mismatch penalty, negative linear gap penalty.
pub struct SGScoring {
  match_reward : Double
  mismatch_penalty : Double
  gap_penalty : Double
}

///|
/// Create a scoring scheme.
pub fn SGScoring::new(
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> SGScoring {
  SGScoring::{ match_reward, mismatch_penalty, gap_penalty, }
}

///|
/// Default DNA scoring: match=2, mismatch=-1, gap=-1.
pub fn SGScoring::default_dna() -> SGScoring {
  SGScoring::{ match_reward: 2.0, mismatch_penalty: -1.0, gap_penalty: -1.0, }
}

///|
/// Default protein scoring: match=5, mismatch=-4, gap=-3.
pub fn SGScoring::default_protein() -> SGScoring {
  SGScoring::{ match_reward: 5.0, mismatch_penalty: -4.0, gap_penalty: -3.0, }
}

///|
/// Full result of a semiglobal alignment.
pub struct SGResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  start1 : Int32 // index in seq1 where the aligned (non-free-gap) region starts
  start2 : Int32 // index in seq2 where the aligned region starts
  end1 : Int32 // index in seq1 where the aligned region ends (exclusive)
  end2 : Int32 // index in seq2 where the aligned region ends (exclusive)
  score_matrix : Array[Array[Double]]
  trace_matrix : Array[Array[SGTrace]]
}

///|
/// Get aligned sequence 1.
pub fn SGResult::aligned1(self : SGResult) -> String {
  self.aligned_seq1
}

///|
/// Get aligned sequence 2.
pub fn SGResult::aligned2(self : SGResult) -> String {
  self.aligned_seq2
}

///|
/// Get match line ('|' for match, '.' for mismatch, ' ' for gap).
pub fn SGResult::match_line(self : SGResult) -> String {
  self.match_line
}

///|
/// Get alignment score.
pub fn SGResult::score(self : SGResult) -> Double {
  self.score
}

///|
/// Get the start index of the aligned region in seq1.
pub fn SGResult::start1(self : SGResult) -> Int32 {
  self.start1
}

///|
/// Get the start index of the aligned region in seq2.
pub fn SGResult::start2(self : SGResult) -> Int32 {
  self.start2
}

///|
/// Get the end index (exclusive) of the aligned region in seq1.
pub fn SGResult::end1(self : SGResult) -> Int32 {
  self.end1
}

///|
/// Get the end index (exclusive) of the aligned region in seq2.
pub fn SGResult::end2(self : SGResult) -> Int32 {
  self.end2
}

///|
/// Get the full score matrix for inspection.
pub fn SGResult::score_matrix(self : SGResult) -> Array[Array[Double]] {
  self.score_matrix
}

///|
/// Get the full traceback matrix for inspection.
pub fn SGResult::trace_matrix(self : SGResult) -> Array[Array[SGTrace]] {
  self.trace_matrix
}

///|
/// Get alignment length (both aligned sequences have the same length).
pub fn SGResult::alignment_length(self : SGResult) -> Int32 {
  self.aligned_seq1.length()
}

///|
/// Count the number of matches in the alignment.
pub fn SGResult::num_matches(self : SGResult) -> Int32 {
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
pub fn SGResult::num_mismatches(self : SGResult) -> Int32 {
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
pub fn SGResult::num_gaps_seq1(self : SGResult) -> Int32 {
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
pub fn SGResult::num_gaps_seq2(self : SGResult) -> Int32 {
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
pub fn SGResult::percent_identity(self : SGResult) -> Double {
  let aln_len = self.alignment_length()
  if aln_len == 0 {
    return 0.0
  }
  self.num_matches().to_double() / aln_len.to_double() * 100.0
}

///|
/// Format the alignment as a readable block (3 lines).
pub fn SGResult::format(self : SGResult) -> String {
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
pub fn SGResult::format_wrapped(self : SGResult, width : Int32) -> String {
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

///|
/// Run exact semiglobal alignment with linear gap penalty.
///
/// Leading and trailing gaps are free:
///   - dp[0][j] = dp[i][0] = 0 for all i, j
///   - the alignment terminates at the best cell on the last row or
///     the last column of the DP matrix
///
/// `seq1` — first sequence (string).
/// `seq2` — second sequence (string).
/// `scoring` — scoring scheme (use `SGScoring::default_dna()` etc.).
///
/// Returns an `SGResult` containing aligned sequences, score, the
/// aligned region coordinates (start1/start2/end1/end2, with the
/// free-gap flanks excluded), and the full DP matrices.
pub fn semiglobal_align(
  seq1 : String,
  seq2 : String,
  scoring : SGScoring,
) -> SGResult {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array(capacity=0)
  let trace : Array[Array[SGTrace]] = Array(capacity=0)
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array(capacity=0)
    let trow : Array[SGTrace] = Array(capacity=0)
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(SGTrace::Start)
    }
    dp.push(drow)
    trace.push(trow)
  }
  // 初始化：首行首列全 0 —— 起始空位不计罚（semiglobal 的关键）
  for i = 1; i <= n; i = i + 1 {
    trace[i][0] = SGTrace::Start
  }
  for j = 1; j <= m; j = j + 1 {
    trace[0][j] = SGTrace::Start
  }
  for i = 1; i <= n; i = i + 1 {
    let c1 = seq1.unsafe_get(i - 1).unsafe_to_char()
    for j = 1; j <= m; j = j + 1 {
      let c2 = seq2.unsafe_get(j - 1).unsafe_to_char()
      let sub_score = if c1 == c2 {
        scoring.match_reward
      } else {
        scoring.mismatch_penalty
      }
      let diag = dp[i - 1][j - 1] + sub_score
      let up = dp[i - 1][j] + scoring.gap_penalty
      let left = dp[i][j - 1] + scoring.gap_penalty
      let (val, dir) = if diag >= up && diag >= left {
        (diag, SGTrace::Diag)
      } else if up >= left {
        (up, SGTrace::Up)
      } else {
        (left, SGTrace::Left)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }
  // 终止：末行或末列上的最大分值 —— 末端空位不计罚
  let mut best_score = dp[n][0]
  let mut best_i = n
  let mut best_j = 0
  for i = 0; i <= n; i = i + 1 {
    // 末列 j = m
    if dp[i][m] > best_score || (i == 0 && m == 0) {
      best_score = dp[i][m]
      best_i = i
      best_j = m
    }
  }
  for j = 0; j <= m; j = j + 1 {
    // 末行 i = n
    if dp[n][j] > best_score {
      best_score = dp[n][j]
      best_i = n
      best_j = j
    }
  }
  // 回溯：从 (best_i, best_j) 走到某行/列的 0 边界
  let mut i = best_i
  let mut j = best_j
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  while i > 0 && j > 0 {
    match trace[i][j] {
      SGTrace::Start => break
      SGTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      SGTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      SGTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = sg_reverse(r1.to_string())
  let aligned2 = sg_reverse(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  // 对齐区间：回溯结束的位置即前导自由空位的分界
  let (start1, start2) = (i, j)
  let (end1, end2) = (best_i, best_j)
  SGResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: best_score,
    start1,
    start2,
    end1,
    end2,
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: semiglobal alignment with default DNA scoring.
pub fn semiglobal_align_dna(seq1 : String, seq2 : String) -> SGResult {
  semiglobal_align(seq1, seq2, SGScoring::default_dna())
}

///|
/// Convenience: semiglobal alignment with default protein scoring.
pub fn semiglobal_align_protein(seq1 : String, seq2 : String) -> SGResult {
  semiglobal_align(seq1, seq2, SGScoring::default_protein())
}

///|
/// Semiglobal alignment using a substitution matrix.
///
/// `seq1`, `seq2` — protein sequences (single-letter code strings).
/// `matrix` — substitution matrix (e.g., `@data.blosum62()`).
/// `gap_penalty` — penalty for a gap (typically negative).
pub fn semiglobal_align_with_matrix(
  seq1 : String,
  seq2 : String,
  matrix : @data.SubstitutionMatrix,
  gap_penalty : Double,
) -> SGResult {
  let n = seq1.length()
  let m = seq2.length()
  let dp : Array[Array[Double]] = Array(capacity=0)
  let trace : Array[Array[SGTrace]] = Array(capacity=0)
  for i = 0; i <= n; i = i + 1 {
    let drow : Array[Double] = Array(capacity=0)
    let trow : Array[SGTrace] = Array(capacity=0)
    for j = 0; j <= m; j = j + 1 {
      drow.push(0.0)
      trow.push(SGTrace::Start)
    }
    dp.push(drow)
    trace.push(trow)
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
        (diag, SGTrace::Diag)
      } else if up >= left {
        (up, SGTrace::Up)
      } else {
        (left, SGTrace::Left)
      }
      dp[i][j] = val
      trace[i][j] = dir
    }
  }
  let mut best_score = dp[0][m]
  let mut best_i = 0
  let mut best_j = m
  for i = 0; i <= n; i = i + 1 {
    if dp[i][m] > best_score {
      best_score = dp[i][m]
      best_i = i
      best_j = m
    }
  }
  for j = 0; j <= m; j = j + 1 {
    if dp[n][j] > best_score {
      best_score = dp[n][j]
      best_i = n
      best_j = j
    }
  }
  let mut i = best_i
  let mut j = best_j
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  while i > 0 && j > 0 {
    match trace[i][j] {
      SGTrace::Start => break
      SGTrace::Diag => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        i = i - 1
        j = j - 1
      }
      SGTrace::Up => {
        r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
        r2.write_char('-')
        i = i - 1
      }
      SGTrace::Left => {
        r1.write_char('-')
        r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
        j = j - 1
      }
    }
  }
  let aligned1 = sg_reverse(r1.to_string())
  let aligned2 = sg_reverse(r2.to_string())
  let match_line = build_match_line(aligned1, aligned2)
  SGResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    match_line,
    score: best_score,
    start1: i,
    start2: j,
    end1: best_i,
    end2: best_j,
    score_matrix: dp,
    trace_matrix: trace,
  }
}

///|
/// Convenience: semiglobal alignment with BLOSUM62 matrix.
pub fn semiglobal_align_blosum62(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> SGResult {
  semiglobal_align_with_matrix(seq1, seq2, @data.blosum62(), gap_penalty)
}

///|
/// Convenience: semiglobal alignment with PAM250 matrix.
pub fn semiglobal_align_pam250(
  seq1 : String,
  seq2 : String,
  gap_penalty : Double,
) -> SGResult {
  semiglobal_align_with_matrix(seq1, seq2, @data.pam250(), gap_penalty)
}

///|
/// Semiglobal alignment of a short sequence (read) onto a long
/// reference (e.g., a genome), returning the zero-based reference
/// interval covered by the alignment.
///
/// This is the typical read-mapping query: leading and trailing
/// reference flanks get free gaps, so the read can sit anywhere.
pub fn semiglobal_map(
  read : String,
  reference : String,
  scoring~ : SGScoring = SGScoring::default_dna(),
) -> SGResult {
  // 参考序列作为 seq1（行维度），read 作为 seq2（列维度），
  // 使自由空位出现在参考侧的首尾
  semiglobal_align(reference, read, scoring)
}

///|
fn sg_reverse(s : String) -> String {
  let chars = s.to_array()
  let n = chars.length()
  let result = StringBuilder()
  for i = n - 1; i >= 0; i = i - 1 {
    result.write_char(chars[i])
  }
  result.to_string()
}
