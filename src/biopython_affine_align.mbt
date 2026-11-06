/// biopython_affine_align.mbt — Affine gap penalty alignment (Gotoh algorithm)
/// Inspired by Biopython Bio.Align.PairwiseAligner with affine gap

///|
/// Affine gap alignment result.
pub struct BiopythonAffineAlignResult {
  score : Double
  aligned1 : String
  aligned2 : String
  start1 : Int32
  end1 : Int32
  start2 : Int32
  end2 : Int32
}

///|
/// Global alignment with affine gap penalties (Gotoh algorithm).
/// Gap cost = gap_open + gap_extend * (gap_length - 1)
pub fn biopython_affine_align_global(
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_open : Double,
  gap_extend : Double,
) -> BiopythonAffineAlignResult {
  let m = seq1.length()
  let n = seq2.length()
  // Three matrices: M (match/mismatch), Ix (gap in seq2), Iy (gap in seq1)
  let neg_inf = -1000000000000000.0
  let m_mat : Array[Array[Double]] = Array::make(
    m + 1,
    Array::make(n + 1, neg_inf),
  )
  let ix : Array[Array[Double]] = Array::make(
    m + 1,
    Array::make(n + 1, neg_inf),
  )
  let iy : Array[Array[Double]] = Array::make(
    m + 1,
    Array::make(n + 1, neg_inf),
  )
  m_mat[0][0] = 0.0
  for i = 1; i <= m; i = i + 1 {
    ix[i][0] = neg_inf
    iy[i][0] = gap_open + gap_extend * i.to_double()
  }
  for j = 1; j <= n; j = j + 1 {
    ix[0][j] = gap_open + gap_extend * j.to_double()
    iy[0][j] = neg_inf
  }
  for i = 1; i <= m; i = i + 1 {
    for j = 1; j <= n; j = j + 1 {
      let sub = if seq1[i - 1] == seq2[j - 1] {
        match_score
      } else {
        mismatch_score
      }
      m_mat[i][j] = sub +
        biopython_max3(m_mat[i - 1][j - 1], ix[i - 1][j - 1], iy[i - 1][j - 1])
      ix[i][j] = biopython_max2(
        m_mat[i - 1][j] + gap_open + gap_extend,
        ix[i - 1][j] + gap_extend,
      )
      iy[i][j] = biopython_max2(
        m_mat[i][j - 1] + gap_open + gap_extend,
        iy[i][j - 1] + gap_extend,
      )
    }
  }
  let score = biopython_max3(m_mat[m][n], ix[m][n], iy[m][n])
  // Traceback
  let (a1, a2) = biopython_affine_traceback(
    m_mat, ix, iy, seq1, seq2, match_score, mismatch_score, gap_open, gap_extend,
  )
  { score, aligned1: a1, aligned2: a2, start1: 0, end1: m, start2: 0, end2: n, }
}

///|
/// Local alignment with affine gap penalties.
pub fn biopython_affine_align_local(
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_open : Double,
  gap_extend : Double,
) -> BiopythonAffineAlignResult {
  let m = seq1.length()
  let n = seq2.length()
  let neg_inf = -1000000000000000.0
  let m_mat : Array[Array[Double]] = Array::make(m + 1, Array::make(n + 1, 0.0))
  let ix : Array[Array[Double]] = Array::make(
    m + 1,
    Array::make(n + 1, neg_inf),
  )
  let iy : Array[Array[Double]] = Array::make(
    m + 1,
    Array::make(n + 1, neg_inf),
  )
  let mut best_score = 0.0
  let mut best_i = 0
  let mut best_j = 0
  for i = 1; i <= m; i = i + 1 {
    for j = 1; j <= n; j = j + 1 {
      let sub = if seq1[i - 1] == seq2[j - 1] {
        match_score
      } else {
        mismatch_score
      }
      m_mat[i][j] = biopython_max4(
        0.0,
        sub + m_mat[i - 1][j - 1],
        sub + ix[i - 1][j - 1],
        sub + iy[i - 1][j - 1],
      )
      ix[i][j] = biopython_max2(
        m_mat[i - 1][j] + gap_open + gap_extend,
        ix[i - 1][j] + gap_extend,
      )
      iy[i][j] = biopython_max2(
        m_mat[i][j - 1] + gap_open + gap_extend,
        iy[i][j - 1] + gap_extend,
      )
      if m_mat[i][j] > best_score {
        best_score = m_mat[i][j]
        best_i = i
        best_j = j
      }
    }
  }
  {
    score: best_score,
    aligned1: seq1[0:best_i].to_owned(),
    aligned2: seq2[0:best_j].to_owned(),
    start1: 0,
    end1: best_i,
    start2: 0,
    end2: best_j,
  }
}

///|
fn biopython_affine_traceback(
  m_mat : Array[Array[Double]],
  _ix : Array[Array[Double]],
  _iy : Array[Array[Double]],
  seq1 : String,
  seq2 : String,
  match_score : Double,
  mismatch_score : Double,
  _gap_open : Double,
  _gap_extend : Double,
) -> (String, String) {
  // Simplified traceback: reconstruct from M matrix diagonal
  let m = seq1.length()
  let n = seq2.length()
  let a1 : Array[Int32] = Array(capacity=0)
  let a2 : Array[Int32] = Array(capacity=0)
  let mut i = m
  let mut j = n
  while i > 0 && j > 0 {
    let sub = if seq1[i - 1] == seq2[j - 1] {
      match_score
    } else {
      mismatch_score
    }
    if m_mat[i][j] == m_mat[i - 1][j - 1] + sub {
      a1.push(seq1[i - 1].to_int())
      a2.push(seq2[j - 1].to_int())
      i = i - 1
      j = j - 1
    } else if i > 0 &&
      m_mat[i][j] == m_mat[i - 1][j] + match_score + mismatch_score {
      a1.push(seq1[i - 1].to_int())
      a2.push(45) // '-'
      i = i - 1
    } else {
      a1.push(45) // '-'
      a2.push(seq2[j - 1].to_int())
      j = j - 1
    }
  }
  while i > 0 {
    a1.push(seq1[i - 1].to_int())
    a2.push(45)
    i = i - 1
  }
  while j > 0 {
    a1.push(45)
    a2.push(seq2[j - 1].to_int())
    j = j - 1
  }
  let s1 = StringBuilder()
  let s2 = StringBuilder()
  for k = a1.length() - 1; k >= 0; k = k - 1 {
    let _ = s1.write_char(a1[k].unsafe_to_char())
    let _ = s2.write_char(a2[k].unsafe_to_char())
  }
  (s1.to_string(), s2.to_string())
}

///|
fn biopython_max2(a : Double, b : Double) -> Double {
  if a >= b {
    a
  } else {
    b
  }
}

///|
fn biopython_max3(a : Double, b : Double, c : Double) -> Double {
  biopython_max2(biopython_max2(a, b), c)
}

///|
fn biopython_max4(a : Double, b : Double, c : Double, d : Double) -> Double {
  biopython_max2(biopython_max3(a, b, c), d)
}
