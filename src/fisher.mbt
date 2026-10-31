///|
/// Fisher's exact test (two-tailed, hypergeometric) and the `bedtools fisher`
/// interval-overlap enrichment wrapper. The statistical core mirrors
/// bedtools2 `fisher/kfunc.cpp` (`kf_lgamma` / `lbinom` / `hypergeo` /
/// `kt_fisher_exact`).

///|
pub struct FisherResult {
  left : Double
  right : Double
  two_tail : Double
  ratio : Double
}

///|
/// log binomial coefficient  ln C(n, k).
fn fisher_lbinom(n : Int32, k : Int32) -> Double {
  if k == 0 || n == k {
    return 0.0
  }
  lgamma(n.to_double() + 1.0) -
  lgamma(k.to_double() + 1.0) -
  lgamma((n - k).to_double() + 1.0)
}

///|
/// Hypergeometric probability P(X = n11) given margins n1_, n_1 and total n.
fn fisher_hypergeo(n11 : Int32, n1_ : Int32, n_1 : Int32, n : Int32) -> Double {
  @math.exp(
    fisher_lbinom(n1_, n11) +
    fisher_lbinom(n - n1_, n_1 - n11) -
    fisher_lbinom(n, n_1),
  )
}

///|
fn fisher_ratio(n11 : Int32, n12 : Int32, n21 : Int32, n22 : Int32) -> Double {
  if n12 == 0 || n21 == 0 || n22 == 0 {
    return 0.0 / 0.0
  }
  n11.to_double() / n12.to_double() / (n21.to_double() / n22.to_double())
}

///|
fn fi_abs(v : Int32) -> Int32 {
  if v < 0 {
    0 - v
  } else {
    v
  }
}

///|
/// Two-tailed Fisher exact test on a 2×2 contingency table
/// `[[n11, n12], [n21, n22]]`, returning left-tail, right-tail and two-tailed
/// p-values plus the odds-ratio-style `(n11/n12)/(n21/n22)` statistic (NaN on
/// division by zero).
pub fn fisher_exact(
  n11 : Int32,
  n12 : Int32,
  n21 : Int32,
  n22 : Int32,
) -> FisherResult {
  let n1_ = n11 + n12
  let n_1 = n11 + n21
  let n = n11 + n12 + n21 + n22
  let max = if n_1 < n1_ { n_1 } else { n1_ }
  let mut min = n1_ + n_1 - n
  if min < 0 {
    min = 0
  }
  let ratio = fisher_ratio(n11, n12, n21, n22)
  if min == max {
    return FisherResult::{ left: 1.0, right: 1.0, two_tail: 1.0, ratio, }
  }
  let q = fisher_hypergeo(n11, n1_, n_1, n)
  if q == 0.0 {
    // the table's probability underflowed: pick the side based on the mode
    if n11 * (n + 2) < (n_1 + 1) * (n1_ + 1) {
      return FisherResult::{ left: 0.0, right: 1.0, two_tail: 0.0, ratio, }
    }
    return FisherResult::{ left: 1.0, right: 0.0, two_tail: 0.0, ratio, }
  }
  // left tail
  let mut p = fisher_hypergeo(min, n1_, n_1, n)
  let mut left = 0.0
  let mut i = min + 1
  while p < 0.99999999 * q && i <= max {
    left = left + p
    p = fisher_hypergeo(i, n1_, n_1, n)
    i = i + 1
  }
  i = i - 1
  if p < 1.00000001 * q {
    left = left + p
  } else {
    i = i - 1
  }
  // right tail
  p = fisher_hypergeo(max, n1_, n_1, n)
  let mut right = 0.0
  let mut j = max - 1
  while p < 0.99999999 * q && j >= 0 {
    right = right + p
    p = fisher_hypergeo(j, n1_, n_1, n)
    j = j - 1
  }
  j = j + 1
  if p < 1.00000001 * q {
    right = right + p
  } else {
    j = j + 1
  }
  let mut two = left + right
  if two > 1.0 {
    two = 1.0
  }
  if fi_abs(i - n11) < fi_abs(j - n11) {
    right = 1.0 - left + q
  } else {
    left = 1.0 - right + q
  }
  FisherResult::{ left, right, two_tail: two, ratio, }
}

///|
/// Fisher exact test on the interval-overlap enrichment between A and B
/// (`bedtools fisher`): builds a 2×2 table whose cells count intervals
/// (overlap, A-only, B-only) with the "neither" cell estimated from the total
/// `genome_size` via bedtools' bMean heuristic, then runs `fisher_exact`.
pub fn bedtools_fisher(
  a : Array[BedInterval],
  b : Array[BedInterval],
  genome_size : Int32,
) -> FisherResult {
  let a_count = a.length()
  let b_count = b.length()
  let mut a_union = 0
  for ai in a {
    a_union = a_union + ai.end - ai.start
  }
  let mut b_union = 0
  for bi in b {
    b_union = b_union + bi.end - bi.start
  }
  let mut overlap_count = 0
  for ai in a {
    for bi in b {
      if ai.chrom == bi.chrom && ai.start < bi.end && bi.start < ai.end {
        overlap_count = overlap_count + 1
      }
    }
  }
  let q_mean = if a_count > 0 {
    1.0 + a_union.to_double() / a_count.to_double()
  } else {
    1.0
  }
  let d_mean = if b_count > 0 {
    1.0 + b_union.to_double() / b_count.to_double()
  } else {
    1.0
  }
  let b_mean = q_mean + d_mean
  let n11 = overlap_count
  let n12 = if a_count - overlap_count > 0 {
    a_count - overlap_count
  } else {
    0
  }
  let n21 = if b_count - overlap_count > 0 {
    b_count - overlap_count
  } else {
    0
  }
  let n22_full_base = n21 + n12 + n11
  let n22_full_est = (genome_size.to_double() / b_mean).to_int()
  let n22_full = if n22_full_base > n22_full_est {
    n22_full_base
  } else {
    n22_full_est
  }
  let n22 = if n22_full - n12 - n21 - n11 > 0 {
    n22_full - n12 - n21 - n11
  } else {
    0
  }
  fisher_exact(n11, n12, n21, n22)
}
