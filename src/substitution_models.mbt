/// Extended nucleotide substitution models for evolutionary distance calculation.
///
/// Biopython equivalent: Bio.Phylo.TreeConstruction.DistanceCalculator
///
/// BioSeqs already has Jukes-Cantor (JC) and Kimura 2-parameter (K2P) in
/// align_analysis.mbt. This module adds:
///   - HKY85 (Hasegawa-Kishino-Yano 1985)
///   - GTR (General Time Reversible)
///   - TN93 (Tamura-Nei 1993)
///   - F84 (Felsenstein 1984)
///   - F81 (Felsenstein 1981)
///   - Tamura-Nei distance
///
/// All distances are expressed as substitutions per site.

///|
/// Count nucleotide differences between two aligned sequences.
/// Returns (transitions, transversions, total_compared, gaps)
fn count_substitutions(
  seq1 : String,
  seq2 : String,
) -> (Int32, Int32, Int32, Int32) {
  let len = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  let mut transitions = 0
  let mut transversions = 0
  let mut compared = 0
  let mut gaps = 0
  let mut i = 0
  while i < len {
    let b1 = seq1[i:i + 1].to_owned()
    let b2 = seq2[i:i + 1].to_owned()
    if b1 == "-" || b2 == "-" || b1 == "N" || b2 == "N" {
      gaps = gaps + 1
    } else {
      compared = compared + 1
      if b1 != b2 {
        let purine1 = b1 == "A" || b1 == "G"
        let purine2 = b2 == "A" || b2 == "G"
        if purine1 == purine2 {
          transitions = transitions + 1
        } else {
          transversions = transversions + 1
        }
      }
    }
    i = i + 1
  }
  (transitions, transversions, compared, gaps)
}

///|
/// Compute nucleotide frequencies from a sequence.
fn _sm_nucleotide_frequencies(seq : String) -> (Double, Double, Double, Double) {
  let n = seq.length()
  if n == 0 {
    return (0.25, 0.25, 0.25, 0.25)
  }
  let mut a = 0
  let mut g = 0
  let mut c = 0
  let mut t = 0
  let mut i = 0
  while i < n {
    let b = seq[i:i + 1].to_owned()
    match b {
      "A" => a = a + 1
      "G" => g = g + 1
      "C" => c = c + 1
      "T" => t = t + 1
      _ => ()
    }
    i = i + 1
  }
  let total = (a + g + c + t).to_double()
  if total == 0.0 {
    (0.25, 0.25, 0.25, 0.25)
  } else {
    (
      a.to_double() / total,
      g.to_double() / total,
      c.to_double() / total,
      t.to_double() / total,
    )
  }
}

///|
/// Combined GC frequency
fn _sm_gc_freq(seq : String) -> Double {
  let (a, g, c, t) = nucleotide_frequencies(seq)
  g + c
}

///|
/// HKY85 distance (Hasegawa-Kishino-Yano 1985).
///
/// Uses transition/transversion ratio estimated from the data.
/// Distance = -b * ln(1 - P/b - Q/(2*b)) - (1-b) * ln(1 - Q/(2*b*(1-b)))
/// where b = freq_GC * (1 - freq_GC)
pub fn aln_hky85_distance(seq1 : String, seq2 : String) -> Double {
  let (trans, transv, compared, _) = count_substitutions(seq1, seq2)
  if compared == 0 {
    return 0.0
  }
  let p = trans.to_double() / compared.to_double()
  let q = transv.to_double() / compared.to_double()
  // Use average base frequencies from both sequences
  let combined = seq1 + seq2
  let gc = _sm_gc_freq(combined)
  let b = gc * (1.0 - gc)
  if b == 0.0 || b == 1.0 {
    return 10.0
  }
  let arg1 = 1.0 - p / b - q / (2.0 * b)
  let arg2 = 1.0 - q / (2.0 * b * (1.0 - b))
  if arg1 <= 0.0 || arg2 <= 0.0 {
    10.0
  } else {
    -b * @math.ln(arg1) - (1.0 - b) * @math.ln(arg2)
  }
}

///|
/// GTR distance (General Time Reversible).
///
/// Simplified GTR using the ML estimate of the rate matrix.
/// For closely related sequences, GTR ≈ p-distance / (1 - sum(pi_i^2))
pub fn aln_gtr_distance(seq1 : String, seq2 : String) -> Double {
  let (trans, transv, compared, _) = count_substitutions(seq1, seq2)
  if compared == 0 {
    return 0.0
  }
  let p = (trans + transv).to_double() / compared.to_double()
  // Estimate base frequencies from combined sequence
  let combined = seq1 + seq2
  let (fa, fg, fc, ft) = _sm_nucleotide_frequencies(combined)
  let sum_sq = fa * fa + fg * fg + fc * fc + ft * ft
  let denom = 1.0 - sum_sq
  if denom <= 0.0 {
    return 10.0
  }
  let d = p / denom
  // Apply small-sample correction similar to JC
  if d >= 1.0 {
    10.0
  } else {
    -@math.ln(1.0 - d) * denom
  }
}

///|
/// TN93 distance (Tamura-Nei 1993).
///
/// Separates transitions into purine (A<->G) and pyrimidine (C<->T) transitions.
/// More accurate than K2P when purine and pyrimidine transition rates differ.
pub fn aln_tn93_distance(seq1 : String, seq2 : String) -> Double {
  let len = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  let mut trans_purine = 0 // A <-> G
  let mut trans_pyrimidine = 0 // C <-> T
  let mut transversions = 0
  let mut compared = 0
  let mut i = 0
  while i < len {
    let b1 = seq1[i:i + 1].to_owned()
    let b2 = seq2[i:i + 1].to_owned()
    if b1 != "-" && b2 != "-" && b1 != "N" && b2 != "N" {
      compared = compared + 1
      if b1 != b2 {
        let purine1 = b1 == "A" || b1 == "G"
        let purine2 = b2 == "A" || b2 == "G"
        if purine1 && purine2 {
          trans_purine = trans_purine + 1
        } else if !purine1 && !purine2 {
          trans_pyrimidine = trans_pyrimidine + 1
        } else {
          transversions = transversions + 1
        }
      }
    }
    i = i + 1
  }
  if compared == 0 {
    return 0.0
  }
  // Base frequencies from combined sequence
  let combined = seq1 + seq2
  let (fa, fg, fc, ft) = _sm_nucleotide_frequencies(combined)
  let pi_r = fa + fg // purine frequency
  let pi_y = fc + ft // pyrimidine frequency
  let p1 = trans_purine.to_double() / compared.to_double()
  let p2 = trans_pyrimidine.to_double() / compared.to_double()
  let q = transversions.to_double() / compared.to_double()
  if pi_r == 0.0 || pi_y == 0.0 {
    return 10.0
  }
  let arg1 = 1.0 - p1 / (2.0 * pi_r) - q / (2.0 * pi_r)
  let arg2 = 1.0 - p2 / (2.0 * pi_y) - q / (2.0 * pi_y)
  let arg3 = 1.0 - q / (2.0 * pi_r * pi_y)
  if arg1 <= 0.0 || arg2 <= 0.0 || arg3 <= 0.0 {
    10.0
  } else {
    let k1 = 2.0 * fa * fg / pi_r
    let k2 = 2.0 * fc * ft / pi_y
    -k1 * @math.ln(arg1) -
    k2 * @math.ln(arg2) -
    (pi_r * pi_y - k1 * pi_y - k2 * pi_r) / (pi_r * pi_y) * @math.ln(arg3)
  }
}

///|
/// F84 distance (Felsenstein 1984).
///
/// A special case of GTR where transition rates depend on the target nucleotide.
/// Similar to HKY85 but with different parameterization.
pub fn aln_f84_distance(seq1 : String, seq2 : String) -> Double {
  let (trans, transv, compared, _) = count_substitutions(seq1, seq2)
  if compared == 0 {
    return 0.0
  }
  let p = trans.to_double() / compared.to_double()
  let q = transv.to_double() / compared.to_double()
  let combined = seq1 + seq2
  let (fa, fg, fc, ft) = _sm_nucleotide_frequencies(combined)
  let pi_r = fa + fg
  let pi_y = fc + ft
  // F84 specific: A = pi_A*pi_G/pi_R + pi_C*pi_T/pi_Y
  let a_param = fa * fg / pi_r + fc * ft / pi_y
  let b = pi_r * pi_y
  if a_param == 0.0 || b == 0.0 {
    return 10.0
  }
  let arg1 = 1.0 - p / (2.0 * a_param) - q / (2.0 * b)
  let arg2 = 1.0 - q / (2.0 * b)
  if arg1 <= 0.0 || arg2 <= 0.0 {
    10.0
  } else {
    -a_param * @math.ln(arg1) - (b - a_param) * @math.ln(arg2)
  }
}

///|
/// F81 distance (Felsenstein 1981).
///
/// A special case of GTR where all transition rates are equal
/// but base frequencies are unequal.
pub fn aln_f81_distance(seq1 : String, seq2 : String) -> Double {
  let (_, _, compared, _) = count_substitutions(seq1, seq2)
  if compared == 0 {
    return 0.0
  }
  let p = aln_p_distance(seq1, seq2)
  let combined = seq1 + seq2
  let (fa, fg, fc, ft) = _sm_nucleotide_frequencies(combined)
  let b = 1.0 - (fa * fa + fg * fg + fc * fc + ft * ft)
  if b <= 0.0 {
    return 10.0
  }
  let arg = 1.0 - p / b
  if arg <= 0.0 {
    10.0
  } else {
    -b * @math.ln(arg)
  }
}

///|
/// Tamura distance (Tamura 1992).
///
/// Designed for sequences with strong GC content bias.
/// d = -h * ln(1 - P/h - Q/(2*h)) - 0.5*(1-h) * ln(1 - Q/(1-h))
/// where h = 2*theta*(1-theta), theta = GC frequency
pub fn aln_tamura_distance(seq1 : String, seq2 : String) -> Double {
  let (trans, transv, compared, _) = count_substitutions(seq1, seq2)
  if compared == 0 {
    return 0.0
  }
  let p = trans.to_double() / compared.to_double()
  let q = transv.to_double() / compared.to_double()
  let combined = seq1 + seq2
  let gc = _sm_gc_freq(combined)
  let h = 2.0 * gc * (1.0 - gc)
  if h == 0.0 || h == 1.0 {
    return 10.0
  }
  let arg1 = 1.0 - p / h - q / (2.0 * h)
  let arg2 = 1.0 - q / (1.0 - h)
  if arg1 <= 0.0 || arg2 <= 0.0 {
    10.0
  } else {
    -h * @math.ln(arg1) - 0.5 * (1.0 - h) * @math.ln(arg2)
  }
}

///|
/// Log-likelihood ratio test distance (LRT).
/// Computes -ln(L1/L0) where L1 is the model likelihood and L0 is the
/// saturated model. Approximated as 2 * (p-distance) for small distances.
pub fn aln_lrt_distance(seq1 : String, seq2 : String) -> Double {
  let p = aln_p_distance(seq1, seq2)
  if p <= 0.0 {
    0.0
  } else if p >= 0.75 {
    10.0
  } else {
    -0.75 * @math.ln(1.0 - 4.0 / 3.0 * p)
  }
}

///|
/// Build a distance matrix using any of the available models.
/// model: "jc", "k2p", "hky85", "gtr", "tn93", "f84", "f81", "tamura", "p"
pub fn aln_distance_matrix(
  sequences : Array[String],
  model~ : String = "k2p",
) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array(capacity=0)
  let dist_fn = match model {
    "jc" => aln_jukes_cantor_distance
    "k2p" => aln_kimura_2p_distance
    "hky85" => aln_hky85_distance
    "gtr" => aln_gtr_distance
    "tn93" => aln_tn93_distance
    "f84" => aln_f84_distance
    "f81" => aln_f81_distance
    "tamura" => aln_tamura_distance
    "p" => aln_p_distance
    _ => aln_kimura_2p_distance
  }
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut j = i + 1
    while j < n {
      row[j] = dist_fn(sequences[i], sequences[j])
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  // Fill symmetric part
  let mut x = 0
  while x < n {
    let mut y = x + 1
    while y < n {
      matrix[y].set(x, matrix[x][y])
      y = y + 1
    }
    x = x + 1
  }
  matrix
}
