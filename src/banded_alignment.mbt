///|
/// Banded global alignment — restricts DP to a diagonal band of width 2*k+1.
///
/// When two sequences are known to be similar (e.g., a protein vs its
/// close homolog, or a read vs its reference), the optimal alignment
/// path stays close to the main diagonal. By only computing cells where
/// |i - j| <= k, we reduce time from O(m*n) to O(m*k) and space similarly.
///
/// If the band is too narrow to contain the optimal path, the result is
/// still a valid alignment but may be suboptimal — the caller should
/// increase k if needed.
///
/// Reference:
/// - Chao KM, Pearson GL, Miller W. "Aligning two sequences within a
///   specified diagonal band." CABIOS 8(6):481-487, 1992.

///|
/// Result of a banded alignment — same shape as NWResult but without
/// the full score/trace matrices (they are band-only and not useful
/// to expose).
pub struct BandedAlignResult {
  aligned_seq1 : String
  aligned_seq2 : String
  match_line : String
  score : Double
  band_used : Int // actual k used (may be widened if sequences differ in length)
}

///|
/// Banded Needleman-Wunsch alignment.
///
/// `k` is the half-bandwidth: only cells where |i - j| <= k are computed.
/// If |n - m| > k, the band is automatically widened to |n - m| + k so
/// the alignment can still reach (n, m).
pub fn banded_needleman_wunsch(
  seq1 : String,
  seq2 : String,
  k : Int,
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> BandedAlignResult {
  let n = seq1.length()
  let m = seq2.length()
  // Widen band if sequences differ in length by more than k
  let len_diff = if n > m { n - m } else { m - n }
  let band = if k >= len_diff { k } else { len_diff + k }
  // For very short sequences or wide bands, fall back to full DP
  if band >= n || band >= m {
    let scoring = NWScoring::{
      match_reward,
      mismatch_penalty,
      gap_penalty,
    }
    let full = needleman_wunsch(seq1, seq2, scoring)
    return BandedAlignResult::{
      aligned_seq1: full.aligned_seq1,
      aligned_seq2: full.aligned_seq2,
      match_line: full.match_line,
      score: full.score,
      band_used: band,
    }
  }
  // DP arrays: two rows (current and previous), each of width 2*band+1
  // Column j in the band for row i ranges from max(0, i-band) to min(m, i+band)
  // We store dp as a full (n+1) x (m+1) matrix but only fill band cells.
  // For simplicity and correctness, use a sparse approach: store only
  // the band cells in a 1D array per row.
  //
  // Row i has columns from lo_i = max(0, i-band) to hi_i = min(m, i+band)
  // Width = hi_i - lo_i + 1
  // We use two arrays: prev_dp (row i-1) and cur_dp (row i)
  // Traceback: store direction per cell
  let trace : Array[Array[Int]] = Array(capacity=n + 1)
  for i = 0; i <= n; i = i + 1 {
    let lo = if i > band { i - band } else { 0 }
    let hi = if i + band < m { i + band } else { m }
    let width = hi - lo + 1
    let row : Array[Int] = Array(capacity=width)
    for dummy = 0; dummy < width; dummy = dummy + 1 {
      row.push(0) // 0=Start, 1=Diag, 2=Up, 3=Left
    }
    trace.push(row)
  }
  // Initialize row 0: dp[0][j] = j * gap_penalty for j in [0, min(m, band)]
  let mut prev_dp : Array[Double] = Array(capacity=2 * band + 1)
  let hi0 = if band < m { band } else { m }
  for j = 0; j <= hi0; j = j + 1 {
    prev_dp.push(j.to_double() * gap_penalty)
  }
  // Fill DP row by row
  for i = 1; i <= n; i = i + 1 {
    let lo_i = if i > band { i - band } else { 0 }
    let hi_i = if i + band < m { i + band } else { m }
    let width_i = hi_i - lo_i + 1
    let cur_dp : Array[Double] = Array(capacity=width_i)
    // Previous row's range
    let lo_prev = if i - 1 > band { i - 1 - band } else { 0 }
    let hi_prev = if i - 1 + band < m { i - 1 + band } else { m }
    for j_idx = 0; j_idx < width_i; j_idx = j_idx + 1 {
      let j = lo_i + j_idx
      // Option 1: gap in seq2 (from row i-1, same column j)
      let up = if j >= lo_prev && j <= hi_prev {
        let prev_jdx = j - lo_prev
        prev_dp[prev_jdx] + gap_penalty
      } else {
        // Out of band — use -inf
        -1.0e18
      }
      // Option 2: gap in seq1 (from row i, column j-1)
      let left = if j > lo_i && j_idx > 0 {
        cur_dp[j_idx - 1] + gap_penalty
      } else if j == 0 {
        // First column: dp[i][0] = i * gap_penalty
        i.to_double() * gap_penalty
      } else {
        -1.0e18
      }
      // Option 3: diagonal (from row i-1, column j-1)
      let diag = if j > 0 && j - 1 >= lo_prev && j - 1 <= hi_prev {
        let prev_jdx = j - 1 - lo_prev
        let c1 = seq1.unsafe_get(i - 1).unsafe_to_char()
        let c2 = seq2.unsafe_get(j - 1).unsafe_to_char()
        let sub = if c1 == c2 { match_reward } else { mismatch_penalty }
        prev_dp[prev_jdx] + sub
      } else {
        -1.0e18
      }
      // Choose best
      let (val, dir) = if diag >= up && diag >= left {
        (diag, 1)
      } else if up >= left {
        (up, 2)
      } else {
        (left, 3)
      }
      cur_dp.push(val)
      trace[i][j_idx] = dir
    }
    prev_dp = cur_dp
  }
  // Traceback from (n, m) to (0, 0)
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  let mut i = n
  let mut j = m
  while i > 0 || j > 0 {
    let lo_i = if i > band { i - band } else { 0 }
    let j_idx = j - lo_i
    let dir = trace[i][j_idx]
    if dir == 1 {
      // Diag
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      i = i - 1
      j = j - 1
    } else if dir == 2 {
      // Up: gap in seq2
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if dir == 3 {
      // Left: gap in seq1
      r1.write_char('-')
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      j = j - 1
    } else {
      // Start (dir == 0): only at (0, 0)
      break
    }
  }
  let aligned1 = banded_reverse(r1.to_string())
  let aligned2 = banded_reverse(r2.to_string())
  let match_line = banded_build_match_line(aligned1, aligned2)
  // Final score is dp[n][m]
  let lo_n = if n > band { n - band } else { 0 }
  let final_score = prev_dp[m - lo_n]
  BandedAlignResult::{ aligned_seq1: aligned1, aligned_seq2: aligned2, match_line, score: final_score, band_used: band }
}

///|
/// Banded Smith-Waterman local alignment.
///
/// Same band restriction as banded NW, but with 0-flooring (local alignment).
pub fn banded_smith_waterman(
  seq1 : String,
  seq2 : String,
  k : Int,
  match_reward : Double,
  mismatch_penalty : Double,
  gap_penalty : Double,
) -> BandedAlignResult {
  let n = seq1.length()
  let m = seq2.length()
  let len_diff = if n > m { n - m } else { m - n }
  let band = if k >= len_diff { k } else { len_diff + k }
  if band >= n || band >= m {
    let scoring = Scoring::{
      match_reward,
      mismatch_penalty,
      gap_penalty,
    }
    let full = smith_waterman(seq1, seq2, scoring)
    return BandedAlignResult::{
      aligned_seq1: full.aligned_seq1,
      aligned_seq2: full.aligned_seq2,
      match_line: full.match_line,
      score: full.score,
      band_used: band,
    }
  }
  // DP with 0-flooring
  let trace : Array[Array[Int]] = Array(capacity=n + 1)
  for i = 0; i <= n; i = i + 1 {
    let lo = if i > band { i - band } else { 0 }
    let hi = if i + band < m { i + band } else { m }
    let width = hi - lo + 1
    let row : Array[Int] = Array(capacity=width)
    for dummy = 0; dummy < width; dummy = dummy + 1 {
      row.push(0)
    }
    trace.push(row)
  }
  let mut prev_dp : Array[Double] = Array(capacity=2 * band + 1)
  let hi0 = if band < m { band } else { m }
  for j = 0; j <= hi0; j = j + 1 {
    prev_dp.push(0.0) // local: dp[0][j] = 0
  }
  let mut best_score = 0.0
  let mut best_i = 0
  let mut best_j = 0
  for i = 1; i <= n; i = i + 1 {
    let lo_i = if i > band { i - band } else { 0 }
    let hi_i = if i + band < m { i + band } else { m }
    let width_i = hi_i - lo_i + 1
    let cur_dp : Array[Double] = Array(capacity=width_i)
    let lo_prev = if i - 1 > band { i - 1 - band } else { 0 }
    let hi_prev = if i - 1 + band < m { i - 1 + band } else { m }
    for j_idx = 0; j_idx < width_i; j_idx = j_idx + 1 {
      let j = lo_i + j_idx
      let up = if j >= lo_prev && j <= hi_prev {
        let prev_jdx = j - lo_prev
        prev_dp[prev_jdx] + gap_penalty
      } else {
        -1.0e18
      }
      let left = if j > lo_i && j_idx > 0 {
        cur_dp[j_idx - 1] + gap_penalty
      } else {
        -1.0e18
      }
      let diag = if j > 0 && j - 1 >= lo_prev && j - 1 <= hi_prev {
        let prev_jdx = j - 1 - lo_prev
        let c1 = seq1.unsafe_get(i - 1).unsafe_to_char()
        let c2 = seq2.unsafe_get(j - 1).unsafe_to_char()
        let sub = if c1 == c2 { match_reward } else { mismatch_penalty }
        prev_dp[prev_jdx] + sub
      } else {
        -1.0e18
      }
      // 0-flooring for local alignment
      let (val, dir) = if diag >= up && diag >= left && diag > 0.0 {
        (diag, 1)
      } else if up >= left && up > 0.0 {
        (up, 2)
      } else if left > 0.0 {
        (left, 3)
      } else {
        (0.0, 0)
      }
      cur_dp.push(val)
      trace[i][j_idx] = dir
      if val > best_score {
        best_score = val
        best_i = i
        best_j = j
      }
    }
    prev_dp = cur_dp
  }
  // Traceback from best cell
  let r1 = StringBuilder()
  let r2 = StringBuilder()
  let mut i = best_i
  let mut j = best_j
  while i > 0 || j > 0 {
    let lo_i = if i > band { i - band } else { 0 }
    let j_idx = j - lo_i
    if j_idx < 0 || j_idx >= trace[i].length() {
      break
    }
    let dir = trace[i][j_idx]
    if dir == 1 {
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      i = i - 1
      j = j - 1
    } else if dir == 2 {
      r1.write_char(seq1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if dir == 3 {
      r1.write_char('-')
      r2.write_char(seq2.unsafe_get(j - 1).unsafe_to_char())
      j = j - 1
    } else {
      break
    }
  }
  let aligned1 = banded_reverse(r1.to_string())
  let aligned2 = banded_reverse(r2.to_string())
  let match_line = banded_build_match_line(aligned1, aligned2)
  BandedAlignResult::{ aligned_seq1: aligned1, aligned_seq2: aligned2, match_line, score: best_score, band_used: band }
}

///|
fn banded_reverse(s : String) -> String {
  let chars = s.to_array()
  let n = chars.length()
  let result = StringBuilder()
  for i = n - 1; i >= 0; i = i - 1 {
    result.write_char(chars[i])
  }
  result.to_string()
}

///|
fn banded_build_match_line(s1 : String, s2 : String) -> String {
  let n = s1.length()
  let result = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let c1 = s1.unsafe_get(i).unsafe_to_char()
    let c2 = s2.unsafe_get(i).unsafe_to_char()
    if c1 == c2 && c1 != '-' {
      result.write_char('|')
    } else if c1 == '-' || c2 == '-' {
      result.write_char(' ')
    } else {
      result.write_char('.')
    }
  }
  result.to_string()
}
