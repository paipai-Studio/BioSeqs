/// Relative Codon Deoptimization (RCD) analysis.
///
/// Used in attenuated vaccine design to measure how much a coding sequence
/// has been deoptimized relative to a reference (e.g., wild-type) codon usage.
///
/// Key concepts:
///   - Codon Adaptation Index (CAI) measures optimization (higher = more optimized)
///   - RCD measures deoptimization (higher = more deoptimized)
///   - RCD = -ln(CAI) or equivalently the sum of negative codon weights
///
/// References:
///   - Coleman et al. (2008) Science 320:1784
///   - Mueller et al. (2006) J Virol 80:3319
///   - Burns et al. (2006) J Virol 80:5769

///|
/// RCD result
pub struct RcdResult {
  sequence : String
  rcd : Double
  relative_deoptimization : Double
  codon_weights : Map[String, Double]
  num_codons : Int32
}

///|
pub fn RcdResult::sequence(self : RcdResult) -> String {
  self.sequence
}

///|
pub fn RcdResult::rcd(self : RcdResult) -> Double {
  self.rcd
}

///|
pub fn RcdResult::relative_deoptimization(self : RcdResult) -> Double {
  self.relative_deoptimization
}

///|
pub fn RcdResult::codon_weights(self : RcdResult) -> Map[String, Double] {
  self.codon_weights
}

///|
pub fn RcdResult::num_codons(self : RcdResult) -> Int32 {
  self.num_codons
}

///|
/// Compute codon deoptimization weights from a reference codon usage table.
///
/// For each codon, weight = ln(frequency / max_frequency_for_same_AA)
/// - Most optimal codon: weight = 0
/// - Less optimal codons: weight < 0
///
/// These are the negative of CAI weights.
pub fn compute_deoptimization_weights(
  reference_table : Map[String, Double],
) -> Map[String, Double] {
  let reverse_map = codon_table_reverse(table_id=1)
  let weights : Map[String, Double] = Map([])
  let aa_keys = reverse_map.keys().to_array()
  let mut k = 0
  while k < aa_keys.length() {
    let aa = aa_keys[k]
    match reverse_map.get(aa) {
      Some(codons) => {
        // Find max frequency for this amino acid
        let mut max_freq = 0.0
        let mut j = 0
        while j < codons.length() {
          let freq = match reference_table.get(codons[j]) {
            Some(f) => f
            None => 0.0
          }
          if freq > max_freq {
            max_freq = freq
          }
          j = j + 1
        }
        // Compute weight for each codon
        let mut m = 0
        while m < codons.length() {
          let freq = match reference_table.get(codons[m]) {
            Some(f) => f
            None => 0.0
          }
          if max_freq > 0.0 && freq > 0.0 {
            weights[codons[m]] = @math.ln(freq / max_freq)
          } else {
            weights[codons[m]] = -10.0
            // Highly deoptimized if frequency is zero
          }
          m = m + 1
        }
      }
      None => ()
    }
    k = k + 1
  }
  weights
}

///|
/// Compute Relative Codon Deoptimization (RCD) for a coding sequence.
///
/// RCD = -mean(codon_weights) = -sum(weights) / num_codons
/// Higher RCD = more deoptimized relative to reference.
pub fn relative_codon_deoptimization(
  seq : String,
  reference_table : Map[String, Double],
) -> Double {
  let weights = compute_deoptimization_weights(reference_table)
  let n = seq.length()
  let mut sum_weights = 0.0
  let mut num_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    match weights.get(codon) {
      Some(w) => {
        sum_weights = sum_weights + w
        num_codons = num_codons + 1
      }
      None => ()
    }
    i = i + 3
  }
  if num_codons == 0 {
    0.0
  } else {
    -sum_weights / num_codons.to_double()
  }
}

///|
/// Compute RCD using pre-computed deoptimization weights.
pub fn rcd_with_weights(seq : String, weights : Map[String, Double]) -> Double {
  let n = seq.length()
  let mut sum_weights = 0.0
  let mut num_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    match weights.get(codon) {
      Some(w) => {
        sum_weights = sum_weights + w
        num_codons = num_codons + 1
      }
      None => ()
    }
    i = i + 3
  }
  if num_codons == 0 {
    0.0
  } else {
    -sum_weights / num_codons.to_double()
  }
}

///|
/// Full RCD analysis with detailed result.
///
/// relative_deoptimization compares the sequence RCD to the reference RCD:
///   - 0.0 = same as reference
///   - >0 = more deoptimized than reference
///   - <0 = more optimized than reference
pub fn rcd_analysis(
  seq : String,
  reference_table : Map[String, Double],
) -> RcdResult {
  let weights = compute_deoptimization_weights(reference_table)
  let n = seq.length()
  let mut sum_weights = 0.0
  let mut num_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    match weights.get(codon) {
      Some(w) => {
        sum_weights = sum_weights + w
        num_codons = num_codons + 1
      }
      None => ()
    }
    i = i + 3
  }
  let rcd = if num_codons == 0 {
    0.0
  } else {
    -sum_weights / num_codons.to_double()
  }
  // Compute reference RCD (expected RCD if using the reference frequencies)
  let mut ref_sum = 0.0
  let mut ref_count = 0
  let weight_keys = weights.keys().to_array()
  let mut k = 0
  while k < weight_keys.length() {
    let codon = weight_keys[k]
    let w = weights[codon]
    let freq = match reference_table.get(codon) {
      Some(f) => f
      None => 0.0
    }
    ref_sum = ref_sum + w * freq
    ref_count = ref_count + 1
    k = k + 1
  }
  let ref_rcd = if ref_count == 0 { 0.0 } else { -ref_sum }
  let relative = rcd - ref_rcd
  {
    sequence: seq,
    rcd,
    relative_deoptimization: relative,
    codon_weights: weights,
    num_codons,
  }
}

///|
/// Deoptimize a coding sequence by selecting the least optimal codon
/// for each amino acid (maximizes RCD).
///
/// This is used in synthetic attenuated virus engineering (SAVE).
pub fn deoptimize_sequence(
  seq : String,
  reference_table : Map[String, Double],
  table_id~ : Int32 = 1,
) -> String {
  let weights = compute_deoptimization_weights(reference_table)
  let reverse_map = codon_table_reverse(table_id~)
  let ct = codon_table(table_id)
  let n = seq.length()
  let buf = StringBuilder()
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    let aa = match ct.get(codon) {
      Some(a) => a
      None => "X"
    }
    if aa == "X" || aa == "*" {
      buf.write_string(codon)
    } else {
      match reverse_map.get(aa) {
        Some(syn_codons) => {
          // Pick the codon with the most negative weight (most deoptimized)
          let mut worst_codon = codon
          let mut worst_weight = 1.0
          let mut j = 0
          while j < syn_codons.length() {
            let w = match weights.get(syn_codons[j]) {
              Some(v) => v
              None => 0.0
            }
            if w < worst_weight {
              worst_weight = w
              worst_codon = syn_codons[j]
            }
            j = j + 1
          }
          buf.write_string(worst_codon)
        }
        None => buf.write_string(codon)
      }
    }
    i = i + 3
  }
  buf.to_string()
}

///|
/// Compute the codon deoptimization index (CDI) as defined by
/// Mueller et al. (2006).
///
/// CDI = exp(sum(f_i * ln(f_i / f_ref_i)))
/// where f_i is the frequency of codon i in the query sequence
/// and f_ref_i is the frequency in the reference.
/// CDI > 1 means deoptimized, CDI < 1 means optimized.
pub fn codon_deoptimization_index(
  seq : String,
  reference_table : Map[String, Double],
) -> Double {
  let codon_counts : Map[String, Int32] = Map([])
  let n = seq.length()
  let mut total_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned().to_upper()
    match codon_counts.get(codon) {
      Some(c) => codon_counts[codon] = c + 1
      None => codon_counts[codon] = 1
    }
    total_codons = total_codons + 1
    i = i + 3
  }
  if total_codons == 0 {
    return 1.0
  }
  let mut sum = 0.0
  let codon_keys = codon_counts.keys().to_array()
  let mut k = 0
  while k < codon_keys.length() {
    let codon = codon_keys[k]
    let f_query = match codon_counts.get(codon) {
      Some(c) => c.to_double() / total_codons.to_double()
      None => 0.0
    }
    let f_ref = match reference_table.get(codon) {
      Some(f) => f
      None => 0.0001
      // Small pseudo-frequency
    }
    if f_query > 0.0 && f_ref > 0.0 {
      sum = sum + f_query * @math.ln(f_query / f_ref)
    }
    k = k + 1
  }
  @math.exp(sum)
}
