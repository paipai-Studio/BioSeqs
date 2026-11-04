/// biopython_multiple_testing.mbt — Multiple testing correction
/// Inspired by Biopython Bio.Statistics: Bonferroni, BH, Holm, BY

///|
/// Bonferroni correction: multiply each p-value by the number of tests.
pub fn biopython_bonferroni_correction(
  pvalues : Array[Double],
) -> Array[Double] {
  let m = pvalues.length().to_double()
  let result : Array[Double] = Array::make(pvalues.length(), 0.0)
  for i = 0; i < pvalues.length(); i = i + 1 {
    result[i] = biopython_min(pvalues[i] * m, 1.0)
  }
  result
}

///|
/// Benjamini-Hochberg FDR correction.
/// Sorts p-values, then adjusts: p_adj = p * m / rank
pub fn biopython_benjamini_hochberg(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length()
  if m == 0 {
    return []
  }
  // Create (p-value, original index) pairs
  let indexed : Array[(Double, Int32)] = Array::make(m, (0.0, 0))
  for i = 0; i < m; i = i + 1 {
    indexed[i] = (pvalues[i], i)
  }
  // Sort by p-value ascending
  indexed.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  // Adjust: p_adj = p * m / rank, then enforce monotonicity from bottom
  let adjusted : Array[Double] = Array::make(m, 0.0)
  let m_d = m.to_double()
  adjusted[m - 1] = biopython_min(indexed[m - 1].0, 1.0)
  for i = m - 2; i >= 0; i = i - 1 {
    let rank = (i + 1).to_double()
    let raw = indexed[i].0 * m_d / rank
    adjusted[i] = biopython_min(raw, adjusted[i + 1])
  }
  // Map back to original order
  let result : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    result[indexed[i].1] = adjusted[i]
  }
  result
}

///|
/// Holm-Bonferroni correction (step-down).
pub fn biopython_holm_correction(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length()
  if m == 0 {
    return []
  }
  let indexed : Array[(Double, Int32)] = Array::make(m, (0.0, 0))
  for i = 0; i < m; i = i + 1 {
    indexed[i] = (pvalues[i], i)
  }
  indexed.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  let result : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    let remaining = (m - i).to_double()
    result[indexed[i].1] = biopython_min(indexed[i].0 * remaining, 1.0)
  }
  // Enforce monotonicity
  let sorted_result : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    sorted_result[i] = result[indexed[i].1]
  }
  for i = 1; i < m; i = i + 1 {
    if sorted_result[i] < sorted_result[i - 1] {
      sorted_result[i] = sorted_result[i - 1]
    }
  }
  for i = 0; i < m; i = i + 1 {
    result[indexed[i].1] = sorted_result[i]
  }
  result
}

///|
/// Benjamini-Yekutieli FDR correction (more conservative than BH).
pub fn biopython_benjamini_yekutieli(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length()
  if m == 0 {
    return []
  }
  // Calculate harmonic number H(m) = sum(1/i for i=1..m)
  let mut harmonic = 0.0
  for i = 1; i <= m; i = i + 1 {
    harmonic = harmonic + 1.0 / i.to_double()
  }
  // BH correction with additional factor c(m) = H(m)
  let bh = biopython_benjamini_hochberg(pvalues)
  let result : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    result[i] = biopython_min(bh[i] * harmonic, 1.0)
  }
  result
}

///|
/// Helper: min of two doubles
fn biopython_min(a : Double, b : Double) -> Double {
  if a < b {
    a
  } else {
    b
  }
}
