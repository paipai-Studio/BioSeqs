///|
/// bp3_pairwise2.mbt — Biopython round 3: pairwise2 with custom gap functions
/// #20 pairwise2 custom gap functions — position-specific gap penalties

///|
pub struct Bp3Pairwise2Result {
  score : Double
  aligned_seq1 : String
  aligned_seq2 : String
  start1 : Int32
  start2 : Int32
  end1 : Int32
  end2 : Int32
}

///|
/// Gap function type: position -> gap penalty
pub fn bp3_pairwise2_align_global(
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_open : Double,
  gap_extend : Double,
  gap_fn1~ : ((Int32) -> Double)? = None,
  gap_fn2~ : ((Int32) -> Double)? = None,
) -> Bp3Pairwise2Result {
  let m = seq1.length()
  let n = seq2.length()
  let neg_inf = -100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.0
  let dp : Array[Array[Double]] = Array(capacity=m + 1)
  let _m1 = m + 1
  for _i in 0..<_m1 {
    let row = Array(capacity=n + 1)
    let _n1 = n + 1
    for _j in 0..<_n1 {
      row.push(neg_inf)
    }
    dp.push(row)
  }
  dp[0][0] = 0.0
  let _m1b = m + 1
  for i in 1..<_m1b {
    let gap_pen = match gap_fn1 {
      Some(f) => f(i)
      None => gap_open + gap_extend * i.to_double()
    }
    dp[i][0] = -gap_pen
  }
  let _n1b = n + 1
  for j in 1..<_n1b {
    let gap_pen = match gap_fn2 {
      Some(f) => f(j)
      None => gap_open + gap_extend * j.to_double()
    }
    dp[0][j] = -gap_pen
  }
  let _m1b = m + 1
  for i in 1..<_m1b {
    let _n1b = n + 1
    for j in 1..<_n1b {
      let c1 = seq1[i - 1:i].to_owned()
      let c2 = seq2[j - 1:j].to_owned()
      let s = if c1 == c2 { match_score } else { mismatch_score }
      let diag = dp[i - 1][j - 1] + s
      let gap_up = match gap_fn1 {
        Some(f) => dp[i - 1][j] - f(i)
        None => dp[i - 1][j] - gap_open - gap_extend
      }
      let gap_left = match gap_fn2 {
        Some(f) => dp[i][j - 1] - f(j)
        None => dp[i][j - 1] - gap_open - gap_extend
      }
      let best = if diag >= gap_up && diag >= gap_left {
        diag
      } else if gap_up >= gap_left {
        gap_up
      } else {
        gap_left
      }
      dp[i][j] = best
    }
  }
  let mut i = m
  let mut j = n
  let aln1 = Array(capacity=0)
  let aln2 = Array(capacity=0)
  while i > 0 || j > 0 {
    if i > 0 && j > 0 {
      let c1 = seq1[i - 1:i].to_owned()
      let c2 = seq2[j - 1:j].to_owned()
      let s = if c1 == c2 { match_score } else { mismatch_score }
      if dp[i][j] - dp[i - 1][j - 1] - s < 0.001 &&
        dp[i][j] - dp[i - 1][j - 1] - s > -0.001 {
        aln1.push(c1)
        aln2.push(c2)
        i = i - 1
        j = j - 1
        continue
      }
    }
    if i > 0 {
      aln1.push(seq1[i - 1:i].to_owned())
      aln2.push("-")
      i = i - 1
      continue
    }
    if j > 0 {
      aln1.push("-")
      aln2.push(seq2[j - 1:j].to_owned())
      j = j - 1
      continue
    }
  }
  let sb1 = StringBuilder()
  let sb2 = StringBuilder()
  for k in 0..<aln1.length() {
    let _ = sb1.write(aln1[aln1.length() - 1 - k])
    let _ = sb2.write(aln2[aln2.length() - 1 - k])
  }
  {
    score: dp[m][n],
    aligned_seq1: sb1.to_string(),
    aligned_seq2: sb2.to_string(),
    start1: 0,
    start2: 0,
    end1: m,
    end2: n,
  }
}

///|
/// Position-specific gap penalty: higher penalty in conserved regions
pub fn bp3_position_specific_gap_penalty(
  conserved_positions : Array[Int32],
  base_penalty : Double,
  conserved_penalty : Double,
) -> (Int32) -> Double {
  fn(pos : Int32) -> Double {
    let mut is_conserved = false
    for i in 0..<conserved_positions.length() {
      if conserved_positions[i] == pos {
        is_conserved = true
        break
      }
    }
    if is_conserved {
      conserved_penalty
    } else {
      base_penalty
    }
  }
}

///|
/// Affine gap function with position-dependent opening
pub fn bp3_affine_gap_fn(open : Double, extend_ : Double) -> (Int32) -> Double {
  fn(length_ : Int32) -> Double {
    if length_ <= 0 {
      0.0
    } else {
      open + extend_ * (length_.to_double() - 1.0)
    }
  }
}

///|
/// Local alignment (Smith-Waterman) with custom gap functions
pub fn bp3_pairwise2_align_local(
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_open : Double,
  gap_extend : Double,
) -> Bp3Pairwise2Result {
  let m = seq1.length()
  let n = seq2.length()
  let dp : Array[Array[Double]] = Array(capacity=m + 1)
  let _m1 = m + 1
  for _i in 0..<_m1 {
    let row = Array(capacity=n + 1)
    let _n1 = n + 1
    for _j in 0..<_n1 {
      row.push(0.0)
    }
    dp.push(row)
  }
  let mut best_score = 0.0
  let mut best_i = 0
  let mut best_j = 0
  let _m1b = m + 1
  for i in 1..<_m1b {
    let _n1b = n + 1
    for j in 1..<_n1b {
      let c1 = seq1[i - 1:i].to_owned()
      let c2 = seq2[j - 1:j].to_owned()
      let s = if c1 == c2 { match_score } else { mismatch_score }
      let diag = dp[i - 1][j - 1] + s
      let gap_up = dp[i - 1][j] - gap_open - gap_extend
      let gap_left = dp[i][j - 1] - gap_open - gap_extend
      let best = if diag >= 0.0 && diag >= gap_up && diag >= gap_left {
        diag
      } else if gap_up >= 0.0 && gap_up >= gap_left {
        gap_up
      } else if gap_left >= 0.0 {
        gap_left
      } else {
        0.0
      }
      dp[i][j] = best
      if best > best_score {
        best_score = best
        best_i = i
        best_j = j
      }
    }
  }
  let mut i = best_i
  let mut j = best_j
  let aln1 = Array(capacity=0)
  let aln2 = Array(capacity=0)
  while i > 0 && j > 0 && dp[i][j] > 0.0 {
    let c1 = seq1[i - 1:i].to_owned()
    let c2 = seq2[j - 1:j].to_owned()
    let s = if c1 == c2 { match_score } else { mismatch_score }
    if dp[i][j] - dp[i - 1][j - 1] - s < 0.001 &&
      dp[i][j] - dp[i - 1][j - 1] - s > -0.001 {
      aln1.push(c1)
      aln2.push(c2)
      i = i - 1
      j = j - 1
    } else if dp[i][j] - dp[i - 1][j] + gap_open + gap_extend < 0.001 &&
      dp[i][j] - dp[i - 1][j] + gap_open + gap_extend > -0.001 {
      aln1.push(c1)
      aln2.push("-")
      i = i - 1
    } else {
      aln1.push("-")
      aln2.push(c2)
      j = j - 1
    }
  }
  let sb1 = StringBuilder()
  let sb2 = StringBuilder()
  for k in 0..<aln1.length() {
    let _ = sb1.write(aln1[aln1.length() - 1 - k])
    let _ = sb2.write(aln2[aln2.length() - 1 - k])
  }
  {
    score: best_score,
    aligned_seq1: sb1.to_string(),
    aligned_seq2: sb2.to_string(),
    start1: i,
    start2: j,
    end1: best_i,
    end2: best_j,
  }
}

///|
pub fn _bp3_pairwise2_dummy() -> Unit {
  let _ = Bp3Pairwise2Result::{
    score: 0.0,
    aligned_seq1: "",
    aligned_seq2: "",
    start1: 0,
    start2: 0,
    end1: 0,
    end2: 0,
  }
}
