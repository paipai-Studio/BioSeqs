///|
/// motif_matrix - Unified motif matrix abstraction (PFM/PPM/PWM/PSSM).
/// Inspired by Biopython Bio.motifs.matrix module.
/// Provides Position-Frequency/Probability/Weight/Scoring matrices with conversions.

///|
/// Position Frequency Matrix: raw counts of each base at each position.
pub struct MotifPFM {
  alphabet : Array[String]
  counts : Array[Array[Double]]
  n_seqs : Double
}

///|
/// Position Probability Matrix: normalized frequencies (0..1).
pub struct MotifPPM {
  alphabet : Array[String]
  values : Array[Array[Double]]
  pseudocount : Double
}

///|
/// Position Weight Matrix: log-odds weights relative to background.
pub struct MotifPWM {
  alphabet : Array[String]
  values : Array[Array[Double]]
  background : Array[Double]
  pseudocount : Double
}

///|
/// Position Specific Scoring Matrix: PWM rounded to integers for fast scanning.
pub struct MotifPSSM {
  alphabet : Array[String]
  values : Array[Array[Int]]
  background : Array[Double]
  pseudocount : Double
  scale : Double
}

///|
/// Construct a PFM from raw count matrix (rows=positions, cols=alphabet order).
pub fn MotifPFM::new(
  alphabet : Array[String],
  counts : Array[Array[Double]]
) -> MotifPFM {
  let n_seqs = if counts.length() > 0 { counts[0][0] + counts[0][1] + counts[0][2] + counts[0][3] } else { 0.0 }
  MotifPFM::{ alphabet, counts, n_seqs }
}

///|
/// Create a 2D Double array with independent rows.
fn motif_matrix_make2d(rows : Int, cols : Int) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < rows {
    result.push(Array::make(cols, 0.0))
    i = i + 1
  }
  result
}

///|
/// Create a 2D Int array with independent rows.
fn motif_matrix_make2d_int(rows : Int, cols : Int) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = []
  let mut i = 0
  while i < rows {
    result.push(Array::make(cols, 0))
    i = i + 1
  }
  result
}

///|
/// Construct a PFM from a list of binding sites (strings over alphabet).
pub fn motif_pfm_from_sites(alphabet : Array[String], sites : Array[String]) -> MotifPFM {
  let n_alpha = alphabet.length()
  let n_pos = if sites.length() > 0 { sites[0].length() } else { 0 }
  let counts : Array[Array[Double]] = []
  let mut pi = 0
  while pi < n_pos {
    counts.push(Array::make(n_alpha, 0.0))
    pi = pi + 1
  }
  let mut s = 0
  while s < sites.length() {
    let site = sites[s]
    let mut p = 0
    while p < site.length() {
      let ch = site[p].unsafe_to_char().to_string()
      let mut a = 0
      while a < n_alpha {
        if ch == alphabet[a] {
          let row = counts[p]
          row[a] = row[a] + 1.0
          counts[p] = row
        }
        a = a + 1
      }
      p = p + 1
    }
    s = s + 1
  }
  MotifPFM::new(alphabet, counts)
}

///|
/// Convert PFM to PPM by normalizing counts to probabilities.
pub fn motif_pfm_to_ppm(pfm : MotifPFM, pseudocount~ : Double = 0.0) -> MotifPPM {
  let n_pos = pfm.counts.length()
  let n_alpha = pfm.alphabet.length()
  let values = motif_matrix_make2d(n_pos, n_alpha)
  let total = pfm.n_seqs + pseudocount * n_alpha.to_double()
  let mut p = 0
  while p < n_pos {
    let mut a = 0
    while a < n_alpha {
      let row = values[p]
      row[a] = (pfm.counts[p][a] + pseudocount) / total
      values[p] = row
      a = a + 1
    }
    p = p + 1
  }
  MotifPPM::{ alphabet: pfm.alphabet, values, pseudocount }
}

///|
/// Convert PPM to PWM using log-odds with background frequencies.
pub fn motif_ppm_to_pwm(
  ppm : MotifPPM,
  background~ : Array[Double] = []
) -> MotifPWM {
  let n_pos = ppm.values.length()
  let n_alpha = ppm.alphabet.length()
  let bg = if background.length() == n_alpha {
    background
  } else {
    Array::make(n_alpha, 1.0 / n_alpha.to_double())
  }
  let values = motif_matrix_make2d(n_pos, n_alpha)
  let mut p = 0
  while p < n_pos {
    let mut a = 0
    while a < n_alpha {
      let row = values[p]
      row[a] = @math.ln(ppm.values[p][a] / bg[a])
      values[p] = row
      a = a + 1
    }
    p = p + 1
  }
  MotifPWM::{ alphabet: ppm.alphabet, values, background: bg, pseudocount: ppm.pseudocount }
}

///|
/// Convert PFM directly to PWM (convenience: PFM→PPM→PWM).
pub fn motif_pfm_to_pwm(
  pfm : MotifPFM,
  background~ : Array[Double] = [],
  pseudocount~ : Double = 1.0
) -> MotifPWM {
  let ppm = motif_pfm_to_ppm(pfm, pseudocount~)
  motif_ppm_to_pwm(ppm, background~)
}

///|
/// Convert PWM to PSSM by scaling and rounding to integers.
pub fn motif_pwm_to_pssm(pwm : MotifPWM, scale~ : Double = 100.0) -> MotifPSSM {
  let n_pos = pwm.values.length()
  let n_alpha = pwm.alphabet.length()
  let values = motif_matrix_make2d_int(n_pos, n_alpha)
  let mut p = 0
  while p < n_pos {
    let mut a = 0
    while a < n_alpha {
      let row = values[p]
      row[a] = (pwm.values[p][a] * scale).round().to_int()
      values[p] = row
      a = a + 1
    }
    p = p + 1
  }
  MotifPSSM::{ alphabet: pwm.alphabet, values, background: pwm.background, pseudocount: pwm.pseudocount, scale }
}

///|
/// Full pipeline: sites → PFM → PPM → PWM → PSSM.
pub fn motif_build_pssm(
  alphabet : Array[String],
  sites : Array[String],
  background~ : Array[Double] = [],
  pseudocount~ : Double = 1.0,
  scale~ : Double = 100.0
) -> MotifPSSM {
  let pfm = motif_pfm_from_sites(alphabet, sites)
  let pwm = motif_pfm_to_pwm(pfm, background~, pseudocount~)
  motif_pwm_to_pssm(pwm, scale~)
}

///|
/// Information content per position (in bits) from PPM.
pub fn motif_ppm_ic_per_position(ppm : MotifPPM) -> Array[Double] {
  let n_pos = ppm.values.length()
  let n_alpha = ppm.alphabet.length()
  let log2 = @math.ln(2.0)
  let bg = 1.0 / n_alpha.to_double()
  let result : Array[Double] = Array::make(n_pos, 0.0)
  let mut p = 0
  while p < n_pos {
    let mut ic = 0.0
    let mut a = 0
    while a < n_alpha {
      let prob = ppm.values[p][a]
      if prob > 0.0 {
        ic = ic + prob * @math.ln(prob / bg) / log2
      }
      a = a + 1
    }
    result[p] = ic
    p = p + 1
  }
  result
}

///|
/// Total information content (sum over positions, in bits).
pub fn motif_ppm_total_ic(ppm : MotifPPM) -> Double {
  let ic_arr = motif_ppm_ic_per_position(ppm)
  let mut total = 0.0
  let mut i = 0
  while i < ic_arr.length() {
    total = total + ic_arr[i]
    i = i + 1
  }
  total
}

///|
/// Score a sequence against a PWM (log-odds score).
pub fn motif_pwm_score(pwm : MotifPWM, seq : String) -> Double {
  let n_pos = pwm.values.length()
  if seq.length() < n_pos {
    return -1.0e30
  }
  let mut score = 0.0
  let mut p = 0
  while p < n_pos {
    let ch = seq[p].unsafe_to_char().to_string()
    let mut a = 0
    let mut found = false
    while a < pwm.alphabet.length() {
      if ch == pwm.alphabet[a] {
        score = score + pwm.values[p][a]
        found = true
      }
      a = a + 1
    }
    if not(found) {
      return -1.0e30
    }
    p = p + 1
  }
  score
}

///|
/// Score a sequence against a PSSM (integer score).
pub fn motif_pssm_score(pssm : MotifPSSM, seq : String) -> Int {
  let n_pos = pssm.values.length()
  if seq.length() < n_pos {
    return -2147483648
  }
  let mut score = 0
  let mut p = 0
  while p < n_pos {
    let ch = seq[p].unsafe_to_char().to_string()
    let mut a = 0
    let mut found = false
    while a < pssm.alphabet.length() {
      if ch == pssm.alphabet[a] {
        score = score + pssm.values[p][a]
        found = true
      }
      a = a + 1
    }
    if not(found) {
      return -2147483648
    }
    p = p + 1
  }
  score
}

///|
/// Find all positions in a sequence where PWM score exceeds threshold.
pub fn motif_pwm_scan(
  pwm : MotifPWM,
  seq : String,
  threshold : Double
) -> Array[(Int, Double)] {
  let n_pos = pwm.values.length()
  let results : Array[(Int, Double)] = []
  let mut i = 0
  while i <= seq.length() - n_pos {
    let sub = seq[i:i + n_pos].to_string()
    let score = motif_pwm_score(pwm, sub)
    if score >= threshold {
      results.push((i, score))
    }
    i = i + 1
  }
  results
}

///|
/// Consensus sequence from PPM (most likely base at each position).
pub fn motif_ppm_consensus(ppm : MotifPPM) -> String {
  let n_pos = ppm.values.length()
  let sb = StringBuilder()
  let mut p = 0
  while p < n_pos {
    let mut best_a = 0
    let mut best_v = -1.0
    let mut a = 0
    while a < ppm.alphabet.length() {
      if ppm.values[p][a] > best_v {
        best_v = ppm.values[p][a]
        best_a = a
      }
      a = a + 1
    }
    sb.write_string(ppm.alphabet[best_a])
    p = p + 1
  }
  sb.to_string()
}

///|
/// Reverse complement of a DNA PWM (for scanning both strands).
pub fn motif_pwm_reverse_complement(pwm : MotifPWM) -> MotifPWM {
  let comp_map : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=4)
    m["A"] = "T"
    m["T"] = "A"
    m["C"] = "G"
    m["G"] = "C"
    m
  }
  let n_pos = pwm.values.length()
  let n_alpha = pwm.alphabet.length()
  let rev_values = motif_matrix_make2d(n_pos, n_alpha)
  let mut p = 0
  while p < n_pos {
    let orig_p = n_pos - 1 - p
    let mut a = 0
    while a < n_alpha {
      let base = pwm.alphabet[a]
      match comp_map.get(base) {
        Some(comp_base) => {
          let mut c = 0
          while c < n_alpha {
            if pwm.alphabet[c] == comp_base {
              let row = rev_values[p]
              row[a] = pwm.values[orig_p][c]
              rev_values[p] = row
            }
            c = c + 1
          }
        }
        None => {
          let row = rev_values[p]
          row[a] = pwm.values[orig_p][a]
          rev_values[p] = row
        }
      }
      a = a + 1
    }
    p = p + 1
  }
  MotifPWM::{ alphabet: pwm.alphabet, values: rev_values, background: pwm.background, pseudocount: pwm.pseudocount }
}

///|
/// Format a PFM as a human-readable table.
pub fn motif_pfm_to_string(pfm : MotifPFM) -> String {
  let sb = StringBuilder()
  sb.write_string("PFM (")
  sb.write_string(pfm.n_seqs.to_string())
  sb.write_string(" sequences):\n")
  let mut a = 0
  while a < pfm.alphabet.length() {
    sb.write_string(pfm.alphabet[a])
    sb.write_string(": ")
    let mut p = 0
    while p < pfm.counts.length() {
      sb.write_string(pfm.counts[p][a].to_string())
      sb.write_string(" ")
      p = p + 1
    }
    sb.write_string("\n")
    a = a + 1
  }
  sb.to_string()
}

///|
/// Format a PWM as a human-readable table.
pub fn motif_pwm_to_string(pwm : MotifPWM) -> String {
  let sb = StringBuilder()
  sb.write_string("PWM (log-odds):\n")
  let mut a = 0
  while a < pwm.alphabet.length() {
    sb.write_string(pwm.alphabet[a])
    sb.write_string(": ")
    let mut p = 0
    while p < pwm.values.length() {
      sb.write_string(pwm.values[p][a].to_string())
      sb.write_string(" ")
      p = p + 1
    }
    sb.write_string("\n")
    a = a + 1
  }
  sb.to_string()
}

///|
/// Example: build a PWM from sample binding sites.
pub fn motif_matrix_example() -> MotifPWM {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGTT", "ACGTT", "ACGTA", "TCGTT", "ACGCT"]
  let pfm = motif_pfm_from_sites(alphabet, sites)
  motif_pfm_to_pwm(pfm)
}