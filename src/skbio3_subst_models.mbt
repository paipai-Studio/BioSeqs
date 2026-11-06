// Scikit-bio third-round: evolutionary substitution models (items 1-5)
// JC69, K80, HKY85, GTR, TN93 distance corrections

///|
pub fn skbio3_jc69_distance(
  p : Double,
) -> Double {
  // Jukes-Cantor 1969: d = -3/4 * ln(1 - 4p/3)
  // p = observed proportion of differences
  let val = 1.0 - 4.0 * p / 3.0
  if val <= 0.0 {
    return 0.0
  }
  -0.75 * @math.ln(val)
}

///|
pub fn skbio3_jc69_from_sequences(
  seq1 : String,
  seq2 : String,
) -> Double {
  // Calculate JC69 corrected distance from two aligned sequences
  let n = seq1.length()
  if n == 0 {
    return 0.0
  }
  let mut diffs = 0
  for i in 0..<n {
    if seq1[i] != seq2[i] {
      diffs = diffs + 1
    }
  }
  let p = diffs.to_double() / n.to_double()
  skbio3_jc69_distance(p)
}

///|
pub fn skbio3_k80_distance(
  transitions : Double,
  transversions : Double,
  n : Int,
) -> Double {
  // Kimura 2-parameter: d = -1/2 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
  // P = transition fraction, Q = transversion fraction
  if n == 0 {
    return 0.0
  }
  let p = transitions / n.to_double()
  let q = transversions / n.to_double()
  let term1 = 1.0 - 2.0 * p - q
  let term2 = 1.0 - 2.0 * q
  if term1 <= 0.0 || term2 <= 0.0 {
    return 0.0
  }
  -0.5 * @math.ln(term1 * term2.sqrt())
}

///|
pub fn skbio3_k80_from_sequences(
  seq1 : String,
  seq2 : String,
) -> Double {
  // Count transitions (A<->G, C<->T) and transversions separately
  let n = seq1.length()
  if n == 0 {
    return 0.0
  }
  let mut transitions = 0.0
  let mut transversions = 0.0
  for i in 0..<n {
    let c1 = seq1[i]
    let c2 = seq2[i]
    if c1 != c2 {
      let is_purine1 = c1 == 'A' || c1 == 'G'
      let is_purine2 = c2 == 'A' || c2 == 'G'
      if is_purine1 == is_purine2 {
        transitions = transitions + 1.0
      } else {
        transversions = transversions + 1.0
      }
    }
  }
  skbio3_k80_distance(transitions, transversions, n)
}

///|
pub fn skbio3_hky85_distance(
  transitions : Double,
  transversions : Double,
  n : Int,
  freq_a : Double,
  freq_g : Double,
  freq_c : Double,
  freq_t : Double,
) -> Double {
  // HKY85: considers base frequencies and ts/tv ratio
  // Simplified: use Kimura-like correction with frequency adjustment
  if n == 0 {
    return 0.0
  }
  let p = transitions / n.to_double()
  let q = transversions / n.to_double()
  let pi_r = freq_a + freq_g  // purine frequency
  let pi_y = freq_c + freq_t  // pyrimidine frequency
  if pi_r <= 0.0 || pi_y <= 0.0 {
    return skbio3_k80_distance(transitions, transversions, n)
  }
  let term1 = 1.0 - p / (2.0 * pi_r) - q / (2.0 * pi_y)
  let term2 = 1.0 - q / (2.0 * pi_r * pi_y)
  if term1 <= 0.0 || term2 <= 0.0 {
    return 0.0
  }
  -pi_r * @math.ln(term1) - pi_y * @math.ln(term2)
}

///|
pub fn skbio3_gtr_distance(
  p_ag : Double,
  p_at : Double,
  p_ac : Double,
  p_gc : Double,
  p_gt : Double,
  p_ct : Double,
  n : Int,
  freq_a : Double,
  freq_g : Double,
  freq_c : Double,
  freq_t : Double,
) -> Double {
  // GTR: general time reversible — 6 substitution rates + 4 frequencies
  // Simplified: weighted sum of substitutions normalized by frequencies
  if n == 0 {
    return 0.0
  }
  let total_sub = p_ag + p_at + p_ac + p_gc + p_gt + p_ct
  let p = total_sub / n.to_double()
  // Expected substitution rate under equilibrium
  let expected_rate = 2.0 * (
    freq_a * freq_g + freq_a * freq_t + freq_a * freq_c +
    freq_g * freq_c + freq_g * freq_t + freq_c * freq_t
  )
  if expected_rate <= 0.0 {
    return 0.0
  }
  let val = 1.0 - p / expected_rate
  if val <= 0.0 {
    return 0.0
  }
  -expected_rate * @math.ln(val)
}

///|
pub fn skbio3_tn93_distance(
  transitions1 : Double,
  transitions2 : Double,
  transversions : Double,
  n : Int,
  freq_a : Double,
  freq_g : Double,
  freq_c : Double,
  freq_t : Double,
) -> Double {
  // Tamura-Nei 1993: two transition classes (A<->G and C<->T) with different rates
  if n == 0 {
    return 0.0
  }
  let p1 = transitions1 / n.to_double()  // A<->G transitions
  let p2 = transitions2 / n.to_double()  // C<->T transitions
  let q = transversions / n.to_double()
  let pi_r = freq_a + freq_g
  let pi_y = freq_c + freq_t
  if pi_r <= 0.0 || pi_y <= 0.0 || freq_a <= 0.0 || freq_g <= 0.0 || freq_c <= 0.0 || freq_t <= 0.0 {
    return skbio3_k80_distance(transitions1 + transitions2, transversions, n)
  }
  let term1 = 1.0 - p1 / (2.0 * freq_a * freq_g / pi_r) - q / (2.0 * pi_y)
  let term2 = 1.0 - p2 / (2.0 * freq_c * freq_t / pi_y) - q / (2.0 * pi_r)
  let term3 = 1.0 - q / (2.0 * pi_r * pi_y)
  if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 {
    return 0.0
  }
  -2.0 * freq_a * freq_g * @math.ln(term1) / pi_r -
    2.0 * freq_c * freq_t * @math.ln(term2) / pi_y -
    2.0 * pi_r * pi_y * @math.ln(term3)
}
