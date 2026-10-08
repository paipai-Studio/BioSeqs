///|
/// Alignment statistics — E-value and p-value computation.
/// Ported from biotite.sequence.align.statistics.
///
/// Karlin-Altschul theory: E = K * m * n * exp(-lambda * S)
/// where K and lambda are parameters estimated from the scoring system.

///|
/// Alignment statistics parameters.
pub struct AlignStats {
  lambda : Double    // Karlin-Altschul lambda
  k : Double         // Karlin-Altschul K
  h : Double         // relative entropy
  score_mean : Double  // expected score per aligned position
}

///|
/// Compute E-value from alignment score using Karlin-Altschul formula.
///
/// Parameters:
///   `score` — alignment score
///   `m` — effective query length
///   `n` — effective database length
///   `stats` — pre-computed lambda, K, H
pub fn evalue(
  score : Double,
  m : Int,
  n : Int,
  stats : AlignStats,
) -> Double {
  // E = K * m * n * exp(-lambda * S)
  let m_d = m.to_double()
  let n_d = n.to_double()
  stats.k * m_d * n_d * @math.exp(-stats.lambda * score)
}

///|
/// Compute p-value from E-value.
/// p = 1 - exp(-E)
pub fn pvalue(e_value : Double) -> Double {
  1.0 - @math.exp(-e_value)
}

///|
/// Compute bit-score from raw score.
/// S' = (lambda * S - ln(K)) / ln(2)
pub fn bit_score(
  score : Double,
  stats : AlignStats,
) -> Double {
  (stats.lambda * score - @math.ln(stats.k)) / 0.6931471805599453
}

///|
/// Compute E-value from bit-score.
/// E = m * n / 2^S'
pub fn evalue_from_bitscore(
  bitscore : Double,
  m : Int,
  n : Int,
) -> Double {
  let m_d = m.to_double()
  let n_d = n.to_double()
  m_d * n_d / pow2(bitscore)
}

///|
/// Estimate Karlin-Altschul parameters from a scoring matrix.
///
/// Uses the method of Altschul & Erickson (1986) for a simple
/// scoring system. For protein matrices (BLOSUM, PAM), pre-computed
/// values are typically used.
///
/// Parameters:
///   `scores` — scoring matrix values (flat array)
///   `freqs` — background frequencies (must sum to 1.0)
pub fn estimate_karlin_altschul(
  scores : Array[Double],
  freqs : Array[Double],
) -> AlignStats? {
  // This is a simplified estimation.
  // Full estimation requires solving for lambda in:
  //   sum_i sum_j freq_i * freq_j * exp(lambda * s_ij) = 1
  //
  // We use binary search to find lambda.
  let n = scores.length()
  if n == 0 || freqs.length() != n {
    return None
  }
  // Check that expected score is negative (required for local alignment)
  let mut expected_score = 0.0
  for i = 0; i < n; i = i + 1 {
    expected_score = expected_score + freqs[i] * scores[i]
  }
  if expected_score >= 0.0 {
    return None
  }
  // Binary search for lambda
  let mut lo = 0.0
  let mut hi = 1.0
  // First, find an upper bound where the sum > 1
  let mut found_upper = false
  for _i = 0; _i < 100; _i = _i + 1 {
    let mut sum_val = 0.0
    for j = 0; j < n; j = j + 1 {
      sum_val = sum_val + freqs[j] * @math.exp(hi * scores[j])
    }
    if sum_val > 1.0 {
      found_upper = true
    } else {
      break
    }
    hi = hi * 2.0
  }
  if not(found_upper) {
    return None
  }
  // Binary search
  for _i = 0; _i < 200; _i = _i + 1 {
    let mid = (lo + hi) / 2.0
    let mut sum_val = 0.0
    for j = 0; j < n; j = j + 1 {
      sum_val = sum_val + freqs[j] * @math.exp(mid * scores[j])
    }
    if sum_val > 1.0 {
      hi = mid
    } else {
      lo = mid
    }
  }
  let lambda = (lo + hi) / 2.0
  // Compute H (relative entropy)
  let mut h = 0.0
  for j = 0; j < n; j = j + 1 {
    let p = freqs[j] * @math.exp(lambda * scores[j])
    if p > 0.0 {
      h = h + p * (lambda * scores[j])
    }
  }
  // K estimation (simplified): K ≈ 1.0 (conservative)
  // Full K estimation is complex; this is a rough approximation
  let k = 1.0
  Some(AlignStats::{
    lambda,
    k,
    h,
    score_mean: expected_score,
  })
}

///|
/// Pre-computed BLASTN parameters (NCBI defaults).
/// lambda = 1.37, K = 0.711, H = 1.31
pub fn blastn_stats() -> AlignStats {
  AlignStats::{ lambda: 1.37, k: 0.711, h: 1.31, score_mean: -1.0 }
}

///|
/// Pre-computed BLOSUM62 parameters (NCBI defaults).
/// lambda = 0.318, K = 0.134, H = 0.481
pub fn blosum62_stats() -> AlignStats {
  AlignStats::{ lambda: 0.318, k: 0.134, h: 0.481, score_mean: -0.52 }
}

///|
/// 2^x for Double.
fn pow2(x : Double) -> Double {
  @math.exp(x * 0.6931471805599453)
}
