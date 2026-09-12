///|
/// Biopython ``Bio.motifs.thresholds`` port: approximate score distribution
/// of a DNA motif PSSM via dynamic programming.
///
/// Given a log-odds position-specific scoring matrix (rows A, C, G, T), the
/// distribution of scores of all length-L sequences under (i) the motif model
/// and (ii) the background model is computed on a discrete grid with a fixed
/// precision, following Norbert Dojer's ScoreDistribution implementation in
/// Biopython. Threshold selectors target false positive / false negative
/// rates, a balanced error rate, and the patser criterion (Hertz & Stormo
/// 1999). All numeric semantics mirror Biopython 1.88 bit-for-bit (including
/// the positive-index initialization quirk of the original Python code).
pub suberror MotifThresholdError {
  MotifThresholdError(String)
}

///|
/// Approximate discrete score distribution of a motif PSSM.
///
/// Grid bin ``i`` represents scores around ``min_score + i * step``.
/// ``motif_density[i]`` is the probability mass under the motif model,
/// ``background_density[i]`` under the background model.
pub struct MotifScoreDistribution {
  min_score : Double
  step : Double
  n_points : Int
  information_content : Double
  motif_density : Array[Double]
  background_density : Array[Double]
}

///|
pub fn MotifScoreDistribution::from_pssm(
  pssm : Array[Array[Double]],
  background~ : Array[Double] = [0.25, 0.25, 0.25, 0.25],
  precision~ : Int = 1000
) -> MotifScoreDistribution raise MotifThresholdError {
  if pssm.length() != 4 {
    raise MotifThresholdError(
      "PSSM must contain exactly 4 rows ordered A, C, G, T"
    )
  }
  let length = pssm[0].length()
  if length == 0 {
    raise MotifThresholdError("PSSM must contain at least one column")
  }
  for row in pssm {
    if row.length() != length {
      raise MotifThresholdError("PSSM rows must all have the same length")
    }
    for value in row {
      if !motif_threshold_finite(value) {
        raise MotifThresholdError(
          "PSSM contains a non-finite value (NaN or infinity); add pseudocounts before building the PSSM"
        )
      }
    }
  }
  if background.length() != 4 {
    raise MotifThresholdError("Background must contain 4 probabilities (A, C, G, T)")
  }
  let bg_total = {
    let mut total = 0.0
    for value in background {
      if !motif_threshold_finite(value) || value < 0.0 {
        raise MotifThresholdError("Background probabilities must be finite and non-negative")
      }
      total += value
    }
    total
  }
  if bg_total <= 0.0 {
    raise MotifThresholdError("Background probabilities must sum to a positive value")
  }
  if precision < 1 {
    raise MotifThresholdError("Precision must be at least 1")
  }
  let n_points = precision * length
  if n_points < 2 {
    raise MotifThresholdError("Need at least two grid points (increase precision)")
  }
  let bg = background.map(fn(value : Double) -> Double { value / bg_total })
  // pssm.min / pssm.max in Biopython are the minimum and maximum *total*
  // scores obtainable against the matrix: sums of the per-column extrema.
  let mut min_total = 0.0
  let mut max_total = 0.0
  for position in 0..<length {
    let mut column_min = pssm[0][position]
    let mut column_max = pssm[0][position]
    for letter in 1..<4 {
      let value = pssm[letter][position]
      if value < column_min {
        column_min = value
      }
      if value > column_max {
        column_max = value
      }
    }
    min_total += column_min
    max_total += column_max
  }
  let min_score = if min_total < 0.0 { min_total } else { 0.0 }
  let max_score = max_total
  let interval = (if max_score > 0.0 { max_score } else { 0.0 }) - min_score
  let step = interval / (n_points - 1).to_double()
  let index_diff =
    fn(x : Double) -> Int {
      ((x + 0.5 * step) / step).floor().to_int()
    }
  let mut motif_density = Array::make(n_points, 0.0)
  let mut background_density = Array::make(n_points, 0.0)
  // Mirror Biopython: densities[-index_diff(min_score)] = 1.0, i.e. the
  // negated (positive) discretization index of min_score.
  let initial_index = 0 - index_diff(min_score)
  if initial_index < 0 || initial_index >= n_points {
    raise MotifThresholdError("Internal error: initial grid index out of range")
  }
  motif_density[initial_index] = 1.0
  background_density[initial_index] = 1.0
  let add_index = fn(i : Int, j : Int) -> Int {
    let k = i + j
    if k < 0 {
      0
    } else if k > n_points - 1 {
      n_points - 1
    } else {
      k
    }
  }
  let mut information_content = 0.0
  for position in 0..<length {
    let motif_new = Array::make(n_points, 0.0)
    let background_new = Array::make(n_points, 0.0)
    for letter in 0..<4 {
      let score = pssm[letter][position]
      let b = bg[letter]
      let motif_prob = @math.pow(2.0, score) * b
      information_content += motif_prob * score
      let shift = index_diff(score)
      for i in 0..<n_points {
        let target = add_index(i, shift)
        motif_new[target] += motif_density[i] * motif_prob
        background_new[target] += background_density[i] * b
      }
    }
    motif_density = motif_new
    background_density = background_new
  }
  {
    min_score,
    step,
    n_points,
    information_content,
    motif_density,
    background_density
  }
}

///|
/// Build a score distribution directly from nucleotide counts.
///
/// Rows are ordered A, C, G, T. A scalar ``pseudocount`` (added to every
/// count, mirroring ``Bio.motifs`` scalar pseudocounts) avoids -infinite
/// log-odds for unobserved letters; the position weight is
/// ``(count + pseudocount) / (column_total + 4 * pseudocount)`` and the
/// PSSM holds base-2 log-odds against ``background``.
pub fn MotifScoreDistribution::from_counts(
  counts : Array[Array[Int]],
  background~ : Array[Double] = [0.25, 0.25, 0.25, 0.25],
  pseudocount~ : Double = 0.0,
  precision~ : Int = 1000
) -> MotifScoreDistribution raise MotifThresholdError {
  if counts.length() != 4 {
    raise MotifThresholdError(
      "Counts must contain exactly 4 rows ordered A, C, G, T"
    )
  }
  let length = counts[0].length()
  if length == 0 {
    raise MotifThresholdError("Counts must contain at least one column")
  }
  for row in counts {
    if row.length() != length {
      raise MotifThresholdError("Count rows must all have the same length")
    }
    for value in row {
      if value < 0 {
        raise MotifThresholdError("Counts must be non-negative")
      }
    }
  }
  if !motif_threshold_finite(pseudocount) || pseudocount < 0.0 {
    raise MotifThresholdError("Pseudocount must be a finite non-negative number")
  }
  let bg_total = {
    let mut total = 0.0
    for value in background {
      total += value
    }
    total
  }
  if background.length() != 4 || bg_total <= 0.0 {
    raise MotifThresholdError("Background must contain 4 positive probabilities")
  }
  let pssm : Array[Array[Double]] = []
  for letter in 0..<4 {
    let b = background[letter] / bg_total
    let row : Array[Double] = []
    for position in 0..<length {
      let mut column_total = 0
      for r in 0..<4 {
        column_total += counts[r][position]
      }
      let probability =
        (counts[letter][position].to_double() + pseudocount) /
        (column_total.to_double() + 4.0 * pseudocount)
      if probability <= 0.0 {
        raise MotifThresholdError(
          "Zero probability cell yields -infinite log-odds; pass a positive pseudocount"
        )
      }
      row.push(@math.log2(probability / b))
    }
    pssm.push(row)
  }
  MotifScoreDistribution::from_pssm(
    pssm,
    background=background,
    precision=precision
  )
}

///|
/// Lowest grid score (``min(0, PSSM minimum)``).
pub fn MotifScoreDistribution::min_score(self : MotifScoreDistribution) -> Double {
  self.min_score
}

///|
/// Width of one discretization bin.
pub fn MotifScoreDistribution::step(self : MotifScoreDistribution) -> Double {
  self.step
}

///|
/// Number of grid bins (``precision * motif_length``).
pub fn MotifScoreDistribution::n_points(self : MotifScoreDistribution) -> Int {
  self.n_points
}

///|
/// Total information content of the PSSM under the given background
/// (expected log-odds score of a true motif occurrence).
pub fn MotifScoreDistribution::information_content(
  self : MotifScoreDistribution
) -> Double {
  self.information_content
}

///|
/// Motif-model probability mass on the discrete score grid (copy).
pub fn MotifScoreDistribution::motif_density(
  self : MotifScoreDistribution
) -> Array[Double] {
  self.motif_density.copy()
}

///|
/// Background-model probability mass on the discrete score grid (copy).
pub fn MotifScoreDistribution::background_density(
  self : MotifScoreDistribution
) -> Array[Double] {
  self.background_density.copy()
}

///|
/// Score represented by grid bin ``index`` (``min_score + index * step``).
pub fn MotifScoreDistribution::grid_score(
  self : MotifScoreDistribution,
  index : Int
) -> Double {
  self.min_score + index.to_double() * self.step
}

///|
/// Log-odds threshold achieving (approximately) the requested false positive
/// rate: the fraction of background sequences scoring at least the threshold
/// is at most ``fpr``.
pub fn MotifScoreDistribution::threshold_fpr(
  self : MotifScoreDistribution,
  fpr : Double
) -> Double raise MotifThresholdError {
  motif_threshold_check_rate(fpr, "false positive rate")
  let mut i = self.n_points
  let mut prob = 0.0
  while prob < fpr {
    i -= 1
    if i < 0 {
      raise MotifThresholdError(
        "False positive rate larger than the discretized background mass"
      )
    }
    prob += self.background_density[i]
  }
  self.min_score + i.to_double() * self.step
}

///|
/// Log-odds threshold achieving (approximately) the requested false negative
/// rate: the fraction of motif occurrences scoring below the threshold is at
/// most ``fnr``.
pub fn MotifScoreDistribution::threshold_fnr(
  self : MotifScoreDistribution,
  fnr : Double
) -> Double raise MotifThresholdError {
  motif_threshold_check_rate(fnr, "false negative rate")
  let mut i = -1
  let mut prob = 0.0
  while prob < fnr {
    i += 1
    if i >= self.n_points {
      raise MotifThresholdError(
        "False negative rate larger than the discretized motif mass"
      )
    }
    prob += self.motif_density[i]
  }
  self.min_score + i.to_double() * self.step
}

///|
/// Log-odds threshold balancing the two error rates so that
/// ``FNR == FPR * rate_proportion`` (default equal rates).
pub fn MotifScoreDistribution::threshold_balanced(
  self : MotifScoreDistribution,
  rate_proportion~ : Double = 1.0
) -> Double raise MotifThresholdError {
  let (threshold, _) = self.threshold_balanced_with_rate(
    rate_proportion=rate_proportion
  )
  threshold
}

///|
/// Same as ``threshold_balanced`` but also returns the actual false positive
/// rate at the selected threshold.
pub fn MotifScoreDistribution::threshold_balanced_with_rate(
  self : MotifScoreDistribution,
  rate_proportion~ : Double = 1.0
) -> (Double, Double) raise MotifThresholdError {
  if !motif_threshold_finite(rate_proportion) || rate_proportion < 0.0 {
    raise MotifThresholdError("Rate proportion must be a finite non-negative number")
  }
  let mut i = self.n_points
  let mut fpr = 0.0
  let mut fnr = 1.0
  while fpr * rate_proportion < fnr {
    i -= 1
    if i < 0 {
      raise MotifThresholdError(
        "Cannot balance error rates at the requested rate proportion"
      )
    }
    fpr += self.background_density[i]
    fnr -= self.motif_density[i]
  }
  (self.min_score + i.to_double() * self.step, fpr)
}

///|
/// patser threshold (Hertz & Stormo 1999): FPR = 2 ** (-IC(M)). The original
/// patser uses natural logarithms, so values are not directly comparable.
pub fn MotifScoreDistribution::threshold_patser(
  self : MotifScoreDistribution
) -> Double raise MotifThresholdError {
  self.threshold_fpr(@math.pow(2.0, -self.information_content))
}

fn motif_threshold_check_rate(
  rate : Double,
  label : String
) -> Unit raise MotifThresholdError {
  if !motif_threshold_finite(rate) || rate < 0.0 || rate > 1.0 {
    raise MotifThresholdError(
      "Requested " + label + " must be a finite number within [0, 1]"
    )
  }
}

fn motif_threshold_finite(value : Double) -> Bool {
  !value.is_nan() && value.abs() <= 1.7976931348623157e308
}
