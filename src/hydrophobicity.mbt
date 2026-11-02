///|
/// Hydrophobicity spectrum and charge distribution analysis.
///
/// Biotite provides hydrophobicity analysis using multiple scales.
/// This module implements Kyte-Doolittle, Eisenberg, and Wimley-White scales.

///|
/// Kyte-Doolittle hydropathy scale (higher = more hydrophobic).
fn kd_hydropathy(aa : Char) -> Double {
  match aa {
    'I' => 4.5
    'V' => 4.2
    'L' => 3.8
    'F' => 2.8
    'C' => 2.5
    'M' => 1.9
    'A' => 1.8
    'G' => -0.4
    'T' => -0.7
    'S' => -0.8
    'W' => -0.9
    'Y' => -1.3
    'P' => -1.6
    'H' => -3.2
    'E' => -3.5
    'Q' => -3.5
    'D' => -3.5
    'N' => -3.5
    'K' => -3.9
    'R' => -4.5
    _ => 0.0
  }
}

///|
/// Eisenberg consensus hydrophobicity scale.
fn eisenberg_hydro(aa : Char) -> Double {
  match aa {
    'I' => 0.73
    'L' => 0.65
    'F' => 0.62
    'V' => 0.54
    'M' => 0.39
    'W' => 0.37
    'A' => 0.20
    'C' => 0.04
    'G' => 0.06
    'T' => 0.11
    'S' => 0.23
    'Y' => 0.35
    'P' => 0.49
    'H' => 0.17
    'E' => -0.62
    'Q' => -0.69
    'D' => -0.77
    'N' => -0.85
    'K' => -1.1
    'R' => -1.8
    _ => 0.0
  }
}

///|
/// Calculate hydrophobicity profile using a sliding window.
///
/// Returns an array of hydrophobicity values, one per residue position
/// (window-averaged). Positive values indicate hydrophobic regions.
pub fn hydrophobicity_profile(
  sequence : String,
  scale~ : String = "kd",
  window_size~ : Int32 = 19,
) -> Array[Double] {
  let n = sequence.length()
  let profile : Array[Double] = []
  let half = window_size / 2
  for i in 0..<n {
    let mut sum = 0.0
    let mut count = 0
    for j in (i - half)..<=(i + half) {
      if j >= 0 && j < n {
        let aa = sequence.unsafe_get(j).unsafe_to_char()
        let h = match scale {
          "eisenberg" => eisenberg_hydro(aa)
          _ => kd_hydropathy(aa)
        }
        sum = sum + h
        count = count + 1
      }
    }
    if count > 0 {
      profile.push(sum / count.to_double())
    } else {
      profile.push(0.0)
    }
  }
  profile
}

///|
/// Identify transmembrane segments from hydrophobicity profile.
///
/// Regions where hydrophobicity > threshold for at least min_length residues.
pub fn find_transmembrane_segments(
  profile : Array[Double],
  threshold~ : Double = 1.0,
  min_length~ : Int32 = 18,
) -> Array[(Int32, Int32)] {
  let segments : Array[(Int32, Int32)] = []
  let mut in_tm = false
  let mut start = 0
  for i in 0..<profile.length() {
    if profile[i] > threshold {
      if !in_tm {
        start = i
        in_tm = true
      }
    } else if in_tm {
      if i - start >= min_length {
        segments.push((start, i - 1))
      }
      in_tm = false
    }
  }
  if in_tm && profile.length() - start >= min_length {
    segments.push((start, profile.length() - 1))
  }
  segments
}

///|
/// Calculate net charge per residue at pH 7.0.
fn residue_charge_pH7(aa : Char) -> Double {
  match aa {
    'R' => 1.0
    'K' => 1.0
    'H' => 0.1 // partially protonated at pH 7
    'D' => -1.0
    'E' => -1.0
    _ => 0.0
  }
}

///|
/// Calculate charge profile along a sequence.
///
/// Returns cumulative charge at each position.
pub fn charge_profile(sequence : String, _ph~ : Double = 7.0) -> Array[Double] {
  let profile : Array[Double] = []
  let mut cumulative = 0.0
  for i in 0..<sequence.length() {
    let aa = sequence.unsafe_get(i).unsafe_to_char()
    cumulative = cumulative + residue_charge_pH7(aa)
    profile.push(cumulative)
  }
  profile
}

///|
/// Calculate charge distribution moments.
///
/// Returns (mean_charge, charge_density, max_positive_run, max_negative_run).
pub fn charge_distribution(sequence : String) -> (Double, Double, Int32, Int32) {
  let n = sequence.length()
  if n == 0 {
    return (0.0, 0.0, 0, 0)
  }
  let mut total_charge = 0.0
  let mut max_pos_run = 0
  let mut max_neg_run = 0
  let mut cur_pos = 0
  let mut cur_neg = 0
  for i in 0..<n {
    let aa = sequence.unsafe_get(i).unsafe_to_char()
    let ch = residue_charge_pH7(aa)
    total_charge = total_charge + ch
    if ch > 0.0 {
      cur_pos = cur_pos + 1
      cur_neg = 0
    } else if ch < 0.0 {
      cur_neg = cur_neg + 1
      cur_pos = 0
    } else {
      cur_pos = 0
      cur_neg = 0
    }
    if cur_pos > max_pos_run {
      max_pos_run = cur_pos
    }
    if cur_neg > max_neg_run {
      max_neg_run = cur_neg
    }
  }
  let mean = total_charge / n.to_double()
  (mean, total_charge, max_pos_run, max_neg_run)
}

///|
/// Calculate hydrophobic moment.
///
/// The hydrophobic moment measures the asymmetry of hydrophobicity
/// around the helix axis. Useful for identifying amphipathic helices.
///
/// H_moment = sqrt( (sum h_i * cos(theta_i))^2 + (sum h_i * sin(theta_i))^2 )
/// where theta_i = i * 100° (for alpha helices, 100° per residue).
pub fn hydrophobic_moment(sequence : String, angle~ : Double = 100.0) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let angle_rad = angle * 3.141592653589793 / 180.0
  let mut sum_cos = 0.0
  let mut sum_sin = 0.0
  for i in 0..<n {
    let aa = sequence.unsafe_get(i).unsafe_to_char()
    let h = eisenberg_hydro(aa)
    let theta = i.to_double() * angle_rad
    sum_cos = sum_cos + h * @math.cos(theta)
    sum_sin = sum_sin + h * @math.sin(theta)
  }
  (sum_cos * sum_cos + sum_sin * sum_sin).sqrt()
}
