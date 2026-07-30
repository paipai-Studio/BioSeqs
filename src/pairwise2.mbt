///|
/// Bio.Pairwise2 module - Flexible pairwise sequence alignment.
/// Inspired by Biopython's Bio.Pairwise2 package.
/// Provides customizable global/local alignment with arbitrary scoring matrices.

///|
/// Global vs local alignment mode.
pub enum PairwiseMode {
  Global
  Local
  GlobalMS
  LocalMS
}

///|
/// Alignment result from Pairwise2.
pub struct PairwiseAlignResult {
  aligned_seq1 : String
  aligned_seq2 : String
  score : Double
  mode : String
  seq1_start : Int
  seq1_end : Int
  seq2_start : Int
  seq2_end : Int
}

///|
/// Get the aligned first sequence.
pub fn PairwiseAlignResult::get_aligned_seq1(self : PairwiseAlignResult) -> String {
  self.aligned_seq1
}

///|
/// Get the aligned second sequence.
pub fn PairwiseAlignResult::get_aligned_seq2(self : PairwiseAlignResult) -> String {
  self.aligned_seq2
}

///|
/// Get the alignment score.
pub fn PairwiseAlignResult::get_score(self : PairwiseAlignResult) -> Double {
  self.score
}

///|
/// Get the alignment mode.
pub fn PairwiseAlignResult::get_mode(self : PairwiseAlignResult) -> String {
  self.mode
}

///|
/// Get the alignment as a formatted string.
pub fn PairwiseAlignResult::format_alignment(self : PairwiseAlignResult, line_width : Int) -> String {
  let w = if line_width <= 0 { 60 } else { line_width }
  let mut result = ""
  let len = self.aligned_seq1.length()

  let mut pos = 0
  while pos < len {
    let end = if pos + w < len { pos + w } else { len }
    let s1 = self.aligned_seq1.substring(start=pos, end=end)
    let s2 = self.aligned_seq2.substring(start=pos, end=end)
    let pos1 = self.seq1_start + pos
    let pos2 = self.seq2_start + pos

    let mut mid = ""
    let mut i = pos
    while i < end {
      let c1 = self.aligned_seq1.substring(start=i, end=i + 1)
      let c2 = self.aligned_seq2.substring(start=i, end=i + 1)
      if c1 == "-" || c2 == "-" {
        mid = mid + " "
      } else if c1 == c2 {
        mid = mid + "|"
      } else {
        mid = mid + "."
      }
      i = i + 1
    }

    result = result + "seq1: \{pos1 + 1}\{s1}\n"
    result = result + "      \{mid}\n"
    result = result + "seq2: \{pos2 + 1}\{s2}\n\n"
    pos = end
  }

  result = result + "Score: \{self.score.to_string()}\n"
  result
}

///|
/// Create a new PairwiseAlignResult.
pub fn pairwise_result_new(
  aligned_seq1 : String,
  aligned_seq2 : String,
  score : Double,
  mode : String,
  seq1_start : Int,
  seq1_end : Int,
  seq2_start : Int,
  seq2_end : Int,
) -> PairwiseAlignResult {
  PairwiseAlignResult::{
    aligned_seq1,
    aligned_seq2,
    score,
    mode,
    seq1_start,
    seq1_end,
    seq2_start,
    seq2_end,
  }
}

///|
/// Create a simple match/mismatch scoring function.
/// Returns a function that scores based on match/mismatch values.
pub fn simple_score(match_sc : Double, mismatch_sc : Double) -> (String, String) -> Double {
  fn(c1 : String, c2 : String) -> Double {
    if c1.to_upper() == c2.to_upper() { match_sc } else { mismatch_sc }
  }
}

///|
/// Create an identity-based scoring function (match = 1, mismatch = 0).
pub fn identity_score() -> (String, String) -> Double {
  simple_score(1.0, 0.0)
}

///|
/// Create a scoring function using a substitution matrix (Map-based).
pub fn matrix_score(matrix : Map[String, Double], default : Double) -> (String, String) -> Double {
  fn(c1 : String, c2 : String) -> Double {
    let k = if c1 <= c2 { c1.to_upper() + c2.to_upper() } else { c2.to_upper() + c1.to_upper() }
    match matrix.get(k) {
      Some(v) => v
      None => default
    }
  }
}

///|
/// Create a standard DNA scoring matrix (A/T/C/G pairwise).
pub fn dna_matrix(match_sc : Double, mismatch_sc : Double) -> Map[String, Double] {
  let m : Map[String, Double] = Map([], capacity=10)
  m["AA"] = match_sc
  m["AC"] = mismatch_sc
  m["AG"] = mismatch_sc
  m["AT"] = mismatch_sc
  m["CC"] = match_sc
  m["CG"] = mismatch_sc
  m["CT"] = mismatch_sc
  m["GG"] = match_sc
  m["GT"] = mismatch_sc
  m["TT"] = match_sc
  m
}

/// Helper: get element from flat Double matrix
fn get_flat_double(arr : Array[Double], i : Int, j : Int, m_plus_1 : Int) -> Double {
  arr[i * m_plus_1 + j]
}

/// Helper: get element from flat Int matrix
fn get_flat_int(arr : Array[Int], i : Int, j : Int, m_plus_1 : Int) -> Int {
  arr[i * m_plus_1 + j]
}

/// Helper: initialize flat Double matrix of size (n+1) x (m+1) with fill_val
fn init_flat_double(n : Int, m : Int, fill_val : Double) -> Array[Double] {
  let size = (n + 1) * (m + 1)
  let arr : Array[Double] = Array::new()
  let mut i = 0
  while i < size {
    arr.push(fill_val)
    i = i + 1
  }
  arr
}

/// Helper: initialize flat Int matrix of size (n+1) x (m+1) with fill_val
fn init_flat_int(n : Int, m : Int, fill_val : Int) -> Array[Int] {
  let size = (n + 1) * (m + 1)
  let arr : Array[Int] = Array::new()
  let mut i = 0
  while i < size {
    arr.push(fill_val)
    i = i + 1
  }
  arr
}

///|
/// Global alignment with custom scoring function and gap penalties.
/// Uses Needleman-Wunsch algorithm with configurable match/mismatch/gap.
pub fn pairwise_globalxx(
  seq1 : String,
  seq2 : String,
  match_sc : Double,
  mismatch_sc : Double,
  gap_open : Double,
  gap_extend : Double,
) -> PairwiseAlignResult {
  let n = seq1.length()
  let m = seq2.length()

  if n == 0 || m == 0 {
    return PairwiseAlignResult::{
      aligned_seq1: seq1,
      aligned_seq2: seq2,
      score: 0.0,
      mode: "globalxx",
      seq1_start: 0,
      seq1_end: n,
      seq2_start: 0,
      seq2_end: m,
    }
  }

  let m_plus_1 = m + 1
  let sc = init_flat_double(n, m, 0.0)
  let tb = init_flat_int(n, m, 0)

  // Set up first row and column
  sc.set(0, 0.0)
  tb.set(0, 0)

  let mut i = 1
  while i <= n {
    sc.set(i * m_plus_1 + 0, gap_open + (i.to_double() - 1.0) * gap_extend)
    tb.set(i * m_plus_1 + 0, 1)
    i = i + 1
  }

  let mut j = 1
  while j <= m {
    sc.set(0 * m_plus_1 + j, gap_open + (j.to_double() - 1.0) * gap_extend)
    tb.set(0 * m_plus_1 + j, 2)
    j = j + 1
  }

  // Fill DP matrix
  i = 1
  while i <= n {
    let c1 = seq1.substring(start=i - 1, end=i).to_upper()
    j = 1
    while j <= m {
      let c2 = seq2.substring(start=j - 1, end=j).to_upper()
      let s = if c1 == c2 { match_sc } else { mismatch_sc }

      let diag = get_flat_double(sc, i - 1, j - 1, m_plus_1) + s
      let up = get_flat_double(sc, i - 1, j, m_plus_1) + gap_extend
      let left = get_flat_double(sc, i, j - 1, m_plus_1) + gap_extend

      let mut best = diag
      let mut best_dir = 0

      if up > best {
        best = up
        best_dir = 1
      }
      if left > best {
        best = left
        best_dir = 2
      }

      sc.set(i * m_plus_1 + j, best)
      tb.set(i * m_plus_1 + j, best_dir)

      j = j + 1
    }
    i = i + 1
  }

  // Traceback
  let mut aligned1 = ""
  let mut aligned2 = ""
  let mut p = n
  let mut q = m

  while p > 0 && q > 0 {
    let dir = get_flat_int(tb, p, q, m_plus_1)
    if dir == 0 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      p = p - 1
      q = q - 1
    } else if dir == 1 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = "-" + aligned2
      p = p - 1
    } else {
      aligned1 = "-" + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      q = q - 1
    }
  }

  while p > 0 {
    aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
    aligned2 = "-" + aligned2
    p = p - 1
  }

  while q > 0 {
    aligned1 = "-" + aligned1
    aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
    q = q - 1
  }

  let final_sc = get_flat_double(sc, n, m, m_plus_1)

  PairwiseAlignResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    score: final_sc,
    mode: "globalxx",
    seq1_start: 0,
    seq1_end: n,
    seq2_start: 0,
    seq2_end: m,
  }
}

///|
/// Local alignment with custom scoring function and gap penalties.
/// Uses Smith-Waterman algorithm with configurable match/mismatch/gap.
pub fn pairwise_localxx(
  seq1 : String,
  seq2 : String,
  match_sc : Double,
  mismatch_sc : Double,
  gap_open : Double,
  gap_extend : Double,
) -> PairwiseAlignResult {
  let n = seq1.length()
  let m = seq2.length()

  if n == 0 || m == 0 {
    return PairwiseAlignResult::{
      aligned_seq1: "",
      aligned_seq2: "",
      score: 0.0,
      mode: "localxx",
      seq1_start: 0,
      seq1_end: 0,
      seq2_start: 0,
      seq2_end: 0,
    }
  }

  let m_plus_1 = m + 1
  let sc = init_flat_double(n, m, 0.0)
  let tb = init_flat_int(n, m, -1)

  let mut best_score = 0.0
  let mut best_i = 0
  let mut best_j = 0

  let mut i = 1
  while i <= n {
    let c1 = seq1.substring(start=i - 1, end=i).to_upper()
    let mut j = 1
    while j <= m {
      let c2 = seq2.substring(start=j - 1, end=j).to_upper()
      let s = if c1 == c2 { match_sc } else { mismatch_sc }

      let diag = get_flat_double(sc, i - 1, j - 1, m_plus_1) + s
      let up = get_flat_double(sc, i - 1, j, m_plus_1) + gap_extend
      let left = get_flat_double(sc, i, j - 1, m_plus_1) + gap_extend

      let mut best = 0.0
      let mut best_dir = -1

      if diag > best {
        best = diag
        best_dir = 0
      }
      if up > best {
        best = up
        best_dir = 1
      }
      if left > best {
        best = left
        best_dir = 2
      }

      sc.set(i * m_plus_1 + j, best)
      tb.set(i * m_plus_1 + j, best_dir)

      if best > best_score {
        best_score = best
        best_i = i
        best_j = j
      }

      j = j + 1
    }
    i = i + 1
  }

  if best_score == 0.0 {
    return PairwiseAlignResult::{
      aligned_seq1: "",
      aligned_seq2: "",
      score: 0.0,
      mode: "localxx",
      seq1_start: 0,
      seq1_end: 0,
      seq2_start: 0,
      seq2_end: 0,
    }
  }

  // Traceback
  let mut aligned1 = ""
  let mut aligned2 = ""
  let mut p = best_i
  let mut q = best_j
  let mut start1 = p
  let mut start2 = q

  while p > 0 && q > 0 && get_flat_int(tb, p, q, m_plus_1) != -1 && get_flat_double(sc, p, q, m_plus_1) > 0.0 {
    let dir = get_flat_int(tb, p, q, m_plus_1)
    if dir == 0 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      start1 = p - 1
      start2 = q - 1
      p = p - 1
      q = q - 1
    } else if dir == 1 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = "-" + aligned2
      start1 = p - 1
      p = p - 1
    } else {
      aligned1 = "-" + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      start2 = q - 1
      q = q - 1
    }
  }

  PairwiseAlignResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    score: best_score,
    mode: "localxx",
    seq1_start: start1,
    seq1_end: best_i,
    seq2_start: start2,
    seq2_end: best_j,
  }
}

///|
/// Global alignment using a substitution matrix (GlobalMS).
pub fn pairwise_globalms(
  seq1 : String,
  seq2 : String,
  matrix : Map[String, Double],
  gap_open : Double,
  gap_extend : Double,
) -> PairwiseAlignResult {
  let scorer = matrix_score(matrix, 0.0)
  let n = seq1.length()
  let m = seq2.length()

  if n == 0 || m == 0 {
    return PairwiseAlignResult::{
      aligned_seq1: seq1,
      aligned_seq2: seq2,
      score: 0.0,
      mode: "globalms",
      seq1_start: 0,
      seq1_end: n,
      seq2_start: 0,
      seq2_end: m,
    }
  }

  let m_plus_1 = m + 1
  let sc = init_flat_double(n, m, 0.0)
  let tb = init_flat_int(n, m, 0)

  // Set up first row and column
  sc.set(0, 0.0)
  tb.set(0, 0)

  let mut i = 1
  while i <= n {
    sc.set(i * m_plus_1 + 0, gap_open + (i.to_double() - 1.0) * gap_extend)
    tb.set(i * m_plus_1 + 0, 1)
    i = i + 1
  }

  let mut j = 1
  while j <= m {
    sc.set(0 * m_plus_1 + j, gap_open + (j.to_double() - 1.0) * gap_extend)
    tb.set(0 * m_plus_1 + j, 2)
    j = j + 1
  }

  // Fill DP matrix
  i = 1
  while i <= n {
    let c1 = seq1.substring(start=i - 1, end=i)
    j = 1
    while j <= m {
      let c2 = seq2.substring(start=j - 1, end=j)
      let s = scorer(c1, c2)

      let diag = get_flat_double(sc, i - 1, j - 1, m_plus_1) + s
      let up = get_flat_double(sc, i - 1, j, m_plus_1) + gap_extend
      let left = get_flat_double(sc, i, j - 1, m_plus_1) + gap_extend

      let mut best = diag
      let mut best_dir = 0

      if up > best {
        best = up
        best_dir = 1
      }
      if left > best {
        best = left
        best_dir = 2
      }

      sc.set(i * m_plus_1 + j, best)
      tb.set(i * m_plus_1 + j, best_dir)

      j = j + 1
    }
    i = i + 1
  }

  // Traceback
  let mut aligned1 = ""
  let mut aligned2 = ""
  let mut p = n
  let mut q = m

  while p > 0 && q > 0 {
    let dir = get_flat_int(tb, p, q, m_plus_1)
    if dir == 0 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      p = p - 1
      q = q - 1
    } else if dir == 1 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = "-" + aligned2
      p = p - 1
    } else {
      aligned1 = "-" + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      q = q - 1
    }
  }

  while p > 0 {
    aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
    aligned2 = "-" + aligned2
    p = p - 1
  }

  while q > 0 {
    aligned1 = "-" + aligned1
    aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
    q = q - 1
  }

  let final_sc = get_flat_double(sc, n, m, m_plus_1)

  PairwiseAlignResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    score: final_sc,
    mode: "globalms",
    seq1_start: 0,
    seq1_end: n,
    seq2_start: 0,
    seq2_end: m,
  }
}

///|
/// Local alignment using a substitution matrix (LocalMS).
pub fn pairwise_localms(
  seq1 : String,
  seq2 : String,
  matrix : Map[String, Double],
  gap_open : Double,
  gap_extend : Double,
) -> PairwiseAlignResult {
  let scorer = matrix_score(matrix, 0.0)
  let n = seq1.length()
  let m = seq2.length()

  if n == 0 || m == 0 {
    return PairwiseAlignResult::{
      aligned_seq1: "",
      aligned_seq2: "",
      score: 0.0,
      mode: "localms",
      seq1_start: 0,
      seq1_end: 0,
      seq2_start: 0,
      seq2_end: 0,
    }
  }

  let m_plus_1 = m + 1
  let sc = init_flat_double(n, m, 0.0)
  let tb = init_flat_int(n, m, -1)

  let mut best_score = 0.0
  let mut best_i = 0
  let mut best_j = 0

  let mut i = 1
  while i <= n {
    let c1 = seq1.substring(start=i - 1, end=i)
    let mut j = 1
    while j <= m {
      let c2 = seq2.substring(start=j - 1, end=j)
      let s = scorer(c1, c2)

      let diag = get_flat_double(sc, i - 1, j - 1, m_plus_1) + s
      let up = get_flat_double(sc, i - 1, j, m_plus_1) + gap_extend
      let left = get_flat_double(sc, i, j - 1, m_plus_1) + gap_extend

      let mut best = 0.0
      let mut best_dir = -1

      if diag > best {
        best = diag
        best_dir = 0
      }
      if up > best {
        best = up
        best_dir = 1
      }
      if left > best {
        best = left
        best_dir = 2
      }

      sc.set(i * m_plus_1 + j, best)
      tb.set(i * m_plus_1 + j, best_dir)

      if best > best_score {
        best_score = best
        best_i = i
        best_j = j
      }

      j = j + 1
    }
    i = i + 1
  }

  if best_score == 0.0 {
    return PairwiseAlignResult::{
      aligned_seq1: "",
      aligned_seq2: "",
      score: 0.0,
      mode: "localms",
      seq1_start: 0,
      seq1_end: 0,
      seq2_start: 0,
      seq2_end: 0,
    }
  }

  // Traceback
  let mut aligned1 = ""
  let mut aligned2 = ""
  let mut p = best_i
  let mut q = best_j
  let mut start1 = p
  let mut start2 = q

  while p > 0 && q > 0 && get_flat_int(tb, p, q, m_plus_1) != -1 && get_flat_double(sc, p, q, m_plus_1) > 0.0 {
    let dir = get_flat_int(tb, p, q, m_plus_1)
    if dir == 0 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      start1 = p - 1
      start2 = q - 1
      p = p - 1
      q = q - 1
    } else if dir == 1 {
      aligned1 = seq1.substring(start=p - 1, end=p) + aligned1
      aligned2 = "-" + aligned2
      start1 = p - 1
      p = p - 1
    } else {
      aligned1 = "-" + aligned1
      aligned2 = seq2.substring(start=q - 1, end=q) + aligned2
      start2 = q - 1
      q = q - 1
    }
  }

  PairwiseAlignResult::{
    aligned_seq1: aligned1,
    aligned_seq2: aligned2,
    score: best_score,
    mode: "localms",
    seq1_start: start1,
    seq1_end: best_i,
    seq2_start: start2,
    seq2_end: best_j,
  }
}

///|
/// Convenience: global alignment with match=2, mismatch=-1, gap_open=-5, gap_extend=-0.5
pub fn pairwise_global(seq1 : String, seq2 : String) -> PairwiseAlignResult {
  pairwise_globalxx(seq1, seq2, 2.0, -1.0, -5.0, -0.5)
}

///|
/// Convenience: local alignment with match=2, mismatch=-1, gap_open=-5, gap_extend=-0.5
pub fn pairwise_local(seq1 : String, seq2 : String) -> PairwiseAlignResult {
  pairwise_localxx(seq1, seq2, 2.0, -1.0, -5.0, -0.5)
}

///|
/// Generate a summary of alignment statistics.
pub fn alignment_summary(result : PairwiseAlignResult) -> String {
  let s1 = result.aligned_seq1
  let s2 = result.aligned_seq2
  let len = s1.length()
  let mut matches = 0
  let mut mismatches = 0
  let mut gaps = 0

  let mut i = 0
  while i < len {
    let c1 = s1.substring(start=i, end=i + 1)
    let c2 = s2.substring(start=i, end=i + 1)
    if c1 == "-" || c2 == "-" {
      gaps = gaps + 1
    } else if c1 == c2 {
      matches = matches + 1
    } else {
      mismatches = mismatches + 1
    }
    i = i + 1
  }

  let identity = if len > 0 { (matches.to_double() / len.to_double()) * 100.0 } else { 0.0 }
  let similarity = if len > 0 { ((matches + mismatches).to_double() / len.to_double()) * 100.0 } else { 0.0 }

  let mut summary = "=== Pairwise Alignment Summary ===\n".to_string()
  summary = summary + "Mode: \{result.mode}\n"
  summary = summary + "Alignment length: \{len}\n"
  summary = summary + "Matches: \{matches}\n"
  summary = summary + "Mismatches: \{mismatches}\n"
  summary = summary + "Gaps: \{gaps}\n"
  summary = summary + "Score: \{result.score.to_string()}\n"
  summary = summary + "Identity: \{identity.to_string()}%\n"
  summary = summary + "Similarity: \{similarity.to_string()}%\n"
  summary
}