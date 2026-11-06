/// biotite_evalue.mbt — E-value estimation for alignment significance
/// Inspired by biotite.sequence.align.EValueEstimator (Karlin-Altschul)

///|
/// E-value estimator using Karlin-Altschul statistics.
/// E = K * m * n * exp(-lambda * S)
/// where K and lambda are parameters estimated from random alignments,
/// m and n are sequence lengths, and S is the alignment score.
pub struct BiotiteEValueEstimator {
  lambda : Double // Karlin-Altschul lambda parameter
  k : Double // Karlin-Altschul K parameter
  h : Double // relative entropy
}

///|
/// Create an E-value estimator with known parameters.
pub fn BiotiteEValueEstimator::new(
  lambda~ : Double = 0.318,
  k~ : Double = 0.134,
  h~ : Double = 0.483,
) -> BiotiteEValueEstimator {
  { lambda, k, h, }
}

///|
/// Estimate Karlin-Altschul parameters from a substitution matrix
/// and background frequencies.
/// This uses the method of iterating to find lambda such that
/// sum(p_i * p_j * exp(lambda * s_ij)) = 1
pub fn BiotiteEValueEstimator::estimate(
  sub_matrix : Array[Array[Double]],
  freqs : Array[Double],
) -> BiotiteEValueEstimator {
  let n = freqs.length()
  // Binary search for lambda
  let mut lo = 0.001
  let mut hi = 1.0
  for _iter = 0; _iter < 100; _iter = _iter + 1 {
    let mid = (lo + hi) / 2.0
    let mut total = 0.0
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        if i < sub_matrix.length() && j < sub_matrix[i].length() {
          total = total +
            freqs[i] * freqs[j] * @math.exp(mid * sub_matrix[i][j])
        }
      }
    }
    if total > 1.0 {
      hi = mid
    } else {
      lo = mid
    }
  }
  let lambda = (lo + hi) / 2.0
  // Estimate H (relative entropy)
  let mut h = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i < sub_matrix.length() && j < sub_matrix[i].length() {
        let p = freqs[i] * freqs[j] * @math.exp(lambda * sub_matrix[i][j])
        if p > 0.0 && sub_matrix[i][j] != 0.0 {
          h = h + p * sub_matrix[i][j]
        }
      }
    }
  }
  // Estimate K (simplified: K ≈ lambda * H)
  let k = lambda * h
  { lambda, k, h, }
}

///|
/// Calculate E-value for an alignment score.
/// E = K * m * n * exp(-lambda * S)
pub fn BiotiteEValueEstimator::evalue(
  self : BiotiteEValueEstimator,
  score : Double,
  m : Int32,
  n : Int32,
) -> Double {
  self.k * m.to_double() * n.to_double() * @math.exp(-self.lambda * score)
}

///|
/// Calculate bit score from raw score.
/// S' = (lambda * S - ln(K)) / ln(2)
pub fn BiotiteEValueEstimator::bit_score(
  self : BiotiteEValueEstimator,
  score : Double,
) -> Double {
  (self.lambda * score - @math.ln(self.k)) / @math.ln(2.0)
}

///|
/// Calculate p-value from E-value.
/// p = 1 - exp(-E)
pub fn BiotiteEValueEstimator::pvalue(
  self : BiotiteEValueEstimator,
  score : Double,
  m : Int32,
  n : Int32,
) -> Double {
  let e = self.evalue(score, m, n)
  1.0 - @math.exp(-e)
}

///|
/// Calculate score threshold for a given E-value.
/// S = ln(K * m * n / E) / lambda
pub fn BiotiteEValueEstimator::score_from_evalue(
  self : BiotiteEValueEstimator,
  evalue : Double,
  m : Int32,
  n : Int32,
) -> Double {
  @math.ln(self.k * m.to_double() * n.to_double() / evalue) / self.lambda
}
