/// Codon Pair Bias (CPB) analysis.
///
/// Biopython equivalent: Not directly in Biopython, but a well-known metric
/// from Coleman et al. (2008) Science, used in mRNA vaccine design.
///
/// CPB measures the non-random usage of adjacent codon pairs. A negative CPB
/// indicates that underrepresented codon pairs are used more than expected,
/// which is associated with reduced translation efficiency — a key mechanism
/// in attenuated vaccine design (e.g., polio, influenza).
///
/// Key formulas:
///   - CPS(codon_i, codon_j) = ln(observed_pair / expected_pair)
///   - expected_pair = (freq_codon_i * freq_codon_j) * total_codon_pairs
///   - CPB = mean of CPS over all consecutive codon pairs in a sequence

///|
/// Codon pair bias result
pub struct CodonPairBiasResult {
  sequence : String
  cpb : Double
  cps_values : Array[Double]
  num_pairs : Int32
}

///|
pub fn CodonPairBiasResult::sequence(self : CodonPairBiasResult) -> String {
  self.sequence
}

///|
pub fn CodonPairBiasResult::cpb(self : CodonPairBiasResult) -> Double {
  self.cpb
}

///|
pub fn CodonPairBiasResult::cps_values(
  self : CodonPairBiasResult,
) -> Array[Double] {
  self.cps_values
}

///|
pub fn CodonPairBiasResult::num_pairs(self : CodonPairBiasResult) -> Int32 {
  self.num_pairs
}

///|
/// Count all codon pairs in a coding sequence.
/// Returns a map from "codon1codon2" -> count.
fn count_codon_pairs(seq : String) -> Map[String, Int32] {
  let n = seq.length()
  let pairs : Map[String, Int32] = Map([])
  let mut i = 0
  while i + 6 <= n {
    let codon1 = seq[i:i + 3].to_owned()
    let codon2 = seq[i + 3:i + 6].to_owned()
    let pair = codon1 + codon2
    match pairs.get(pair) {
      Some(c) => pairs[pair] = c + 1
      None => pairs[pair] = 1
    }
    i = i + 3
  }
  pairs
}

///|
/// Count individual codon frequencies in a sequence.
fn count_codons_in_seq(seq : String) -> Map[String, Int32] {
  let n = seq.length()
  let counts : Map[String, Int32] = Map([])
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned()
    match counts.get(codon) {
      Some(c) => counts[codon] = c + 1
      None => counts[codon] = 1
    }
    i = i + 3
  }
  counts
}

///|
/// Compute codon pair score (CPS) table from a reference coding sequence.
///
/// CPS(codon_i, codon_j) = ln(observed(i,j) / expected(i,j))
/// where expected(i,j) = count(i) * count(j) / total_pairs
///
/// Returns a map from "codon1codon2" -> CPS value.
pub fn compute_cps_table(reference_seq : String) -> Map[String, Double] {
  let pair_counts = count_codon_pairs(reference_seq)
  let codon_counts = count_codons_in_seq(reference_seq)
  let total_codons = codon_counts
    .values()
    .to_array()
    .fold(init=0, fn(acc : Int32, v : Int32) -> Int32 { acc + v })
  let total_pairs = pair_counts
    .values()
    .to_array()
    .fold(init=0, fn(acc : Int32, v : Int32) -> Int32 { acc + v })
  let cps_table : Map[String, Double] = Map([])
  if total_pairs == 0 || total_codons == 0 {
    return cps_table
  }
  let pair_keys = pair_counts.keys().to_array()
  let mut k = 0
  while k < pair_keys.length() {
    let pair = pair_keys[k]
    let observed = match pair_counts.get(pair) {
      Some(c) => c.to_double()
      None => 0.0
    }
    // pair = codon1 + codon2 (each 3 chars)
    let codon1 = pair[0:3].to_owned()
    let codon2 = pair[3:6].to_owned()
    let count1 = match codon_counts.get(codon1) {
      Some(c) => c.to_double()
      None => 1.0
    }
    let count2 = match codon_counts.get(codon2) {
      Some(c) => c.to_double()
      None => 1.0
    }
    let expected = count1 * count2 / total_codons.to_double()
    if expected > 0.0 && observed > 0.0 {
      cps_table[pair] = @math.ln(observed / expected)
    } else {
      cps_table[pair] = 0.0
    }
    k = k + 1
  }
  cps_table
}

///|
/// Compute Codon Pair Bias (CPB) for a coding sequence.
///
/// CPB = average CPS over all consecutive codon pairs.
/// A CPB of 0 means random codon pair usage.
/// Negative CPB indicates underrepresented pairs are overused.
pub fn codon_pair_bias(seq : String) -> Double {
  let n = seq.length()
  if n < 6 {
    return 0.0
  }
  // Build CPS table from the sequence itself
  let cps_table = compute_cps_table(seq)
  let mut sum_cps = 0.0
  let mut num_pairs = 0
  let mut i = 0
  while i + 6 <= n {
    let codon1 = seq[i:i + 3].to_owned()
    let codon2 = seq[i + 3:i + 6].to_owned()
    let pair = codon1 + codon2
    match cps_table.get(pair) {
      Some(cps) => {
        sum_cps = sum_cps + cps
        num_pairs = num_pairs + 1
      }
      None => ()
    }
    i = i + 3
  }
  if num_pairs == 0 {
    0.0
  } else {
    sum_cps / num_pairs.to_double()
  }
}

///|
/// Compute CPB using a pre-computed CPS table from a reference.
pub fn codon_pair_bias_with_reference(
  seq : String,
  cps_table : Map[String, Double],
) -> Double {
  let n = seq.length()
  if n < 6 {
    return 0.0
  }
  let mut sum_cps = 0.0
  let mut num_pairs = 0
  let mut i = 0
  while i + 6 <= n {
    let codon1 = seq[i:i + 3].to_owned()
    let codon2 = seq[i + 3:i + 6].to_owned()
    let pair = codon1 + codon2
    match cps_table.get(pair) {
      Some(cps) => {
        sum_cps = sum_cps + cps
        num_pairs = num_pairs + 1
      }
      None => ()
    }
    i = i + 3
  }
  if num_pairs == 0 {
    0.0
  } else {
    sum_cps / num_pairs.to_double()
  }
}

///|
/// Full CPB analysis with detailed result.
pub fn codon_pair_bias_analysis(seq : String) -> CodonPairBiasResult {
  let n = seq.length()
  let cps_table = compute_cps_table(seq)
  let cps_values : Array[Double] = []
  let mut sum_cps = 0.0
  let mut num_pairs = 0
  let mut i = 0
  while i + 6 <= n {
    let codon1 = seq[i:i + 3].to_owned()
    let codon2 = seq[i + 3:i + 6].to_owned()
    let pair = codon1 + codon2
    match cps_table.get(pair) {
      Some(cps) => {
        cps_values.push(cps)
        sum_cps = sum_cps + cps
        num_pairs = num_pairs + 1
      }
      None => ()
    }
    i = i + 3
  }
  let cpb = if num_pairs == 0 { 0.0 } else { sum_cps / num_pairs.to_double() }
  { sequence: seq, cpb, cps_values, num_pairs, }
}

///|
/// Minimize CPB by recoding a coding sequence with synonymous codons
/// that produce the most underrepresented codon pairs.
///
/// This is the principle behind CPB-deoptimized vaccines (Coleman et al. 2008).
/// Returns the recoded sequence and its CPB.
pub fn deoptimize_cpb(
  seq : String,
  cps_table : Map[String, Double],
  table_id~ : Int32 = 1,
) -> String {
  let reverse_map = codon_table_reverse(table_id~)
  let n = seq.length()
  let buf = StringBuilder()
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned()
    // Get the amino acid for this codon
    let ct = codon_table(table_id)
    let aa = match ct.get(codon) {
      Some(a) => a
      None => "X"
    }
    if aa == "X" || aa == "*" {
      buf.write_string(codon)
    } else {
      // Find synonymous codons
      match reverse_map.get(aa) {
        Some(syn_codons) => {
          // Pick the codon that, combined with the previous codon,
          // gives the lowest CPS
          let prev_codon = if i >= 3 {
            Some(buf.to_string()[-3:].to_owned())
          } else {
            None
          }
          let mut worst_codon = codon
          let mut worst_score = 999.0
          let mut j = 0
          while j < syn_codons.length() {
            let candidate = syn_codons[j]
            let score = match prev_codon {
              Some(prev) => {
                let pair = prev + candidate
                match cps_table.get(pair) {
                  Some(s) => s
                  None => 0.0
                }
              }
              None => 0.0
            }
            if score < worst_score {
              worst_score = score
              worst_codon = candidate
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
