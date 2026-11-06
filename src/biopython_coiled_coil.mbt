/// biopython_coiled_coil.mbt — Coiled coil prediction
/// Inspired by Biopython COILS-like heuristics

///|
/// Coiled coil prediction result.
pub struct BiopythonCoiledCoil {
  start : Int32
  end : Int32
  length : Int32
  score : Double
  heptad_register : String // heptad repeat positions (a-g)
}

///|
/// Heptad repeat position-specific scoring.
/// Positions a and d are hydrophobic core positions.
/// Positions e and g are charged/polar positions.
fn biopython_cc_position_score(aa : String, pos : Int32) -> Double {
  // pos is 0-6 (a=0, b=1, c=2, d=3, e=4, f=5, g=6)
  let kd = biopython_cc_kd(aa)
  match pos {
    0 | 3 =>
      // a, d positions: prefer hydrophobic
      if kd > 0.0 {
        kd * 1.5
      } else {
        kd * 0.5
      }
    4 | 6 =>
      // e, g positions: prefer charged
      if aa == "K" || aa == "R" || aa == "E" || aa == "D" {
        1.0
      } else {
        kd * 0.3
      }
    _ =>
      // b, c, f positions: neutral
      kd * 0.5
  }
}

///|
/// Kyte-Doolittle hydrophobicity.
fn biopython_cc_kd(aa : String) -> Double {
  match aa {
    "A" => 1.8
    "R" => -4.5
    "N" => -3.5
    "D" => -3.5
    "C" => 2.5
    "Q" => -3.5
    "E" => -3.5
    "G" => -0.4
    "H" => -3.2
    "I" => 4.5
    "L" => 3.8
    "K" => -3.9
    "M" => 1.9
    "F" => 2.8
    "P" => -1.6
    "S" => -0.8
    "T" => -0.7
    "W" => -0.9
    "Y" => -1.3
    "V" => 4.2
    _ => 0.0
  }
}

///|
/// Predict coiled coil regions in a protein sequence.
/// Uses a sliding window and heptad repeat scoring.
pub fn biopython_predict_coiled_coil(
  sequence : String,
  window~ : Int32 = 28, // typical coiled coil window (4 heptads)
  threshold~ : Double = 0.5,
) -> Array[BiopythonCoiledCoil] {
  let n = sequence.length()
  if n < 7 {
    return []
  }
  // Score each position for each heptad register offset (0-6)
  let scores : Array[Double] = Array::make(n, 0.0)
  let registers : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let mut best_score = -1000000000000000.0
    let mut best_offset = 0
    for offset = 0; offset < 7; offset = offset + 1 {
      let half = window / 2
      let start = if i > half { i - half } else { 0 }
      let end = if i + half + 1 < n { i + half + 1 } else { n }
      let mut sum = 0.0
      let count = end - start
      for j = start; j < end; j = j + 1 {
        let pos = (j - i + offset + 7) % 7
        let aa = sequence[j:j + 1].to_owned()
        sum = sum + biopython_cc_position_score(aa, pos)
      }
      let avg = if count > 0 { sum / count.to_double() } else { 0.0 }
      if avg > best_score {
        best_score = avg
        best_offset = offset
      }
    }
    scores[i] = best_score
    registers[i] = best_offset
  }
  // Find contiguous regions above threshold
  let results : Array[BiopythonCoiledCoil] = Array(capacity=0)
  let mut i = 0
  while i < n {
    if scores[i] > threshold {
      let start = i
      let mut end = i
      while end < n && scores[end] > threshold {
        end = end + 1
      }
      // Calculate average score
      let mut avg = 0.0
      for j = start; j < end; j = j + 1 {
        avg = avg + scores[j]
      }
      avg = avg / (end - start).to_double()
      // Build heptad register string
      let register_sb = StringBuilder()
      for j = start; j < end; j = j + 1 {
        let pos = (j - start + registers[start]) % 7
        let ch = match pos {
          0 => 'a'
          1 => 'b'
          2 => 'c'
          3 => 'd'
          4 => 'e'
          5 => 'f'
          _ => 'g'
        }
        let _ = register_sb.write_char(ch)
      }
      results.push({
        start,
        end,
        length: end - start,
        score: avg,
        heptad_register: register_sb.to_string(),
      })
      i = end
    } else {
      i = i + 1
    }
  }
  results
}

///|
/// Check if a sequence is likely a coiled coil.
pub fn biopython_is_coiled_coil(
  sequence : String,
  threshold~ : Double = 0.5,
) -> Bool {
  biopython_predict_coiled_coil(sequence, threshold~).length() > 0
}

///|
/// Predict oligomerization state from coiled coil.
/// Counts residues at a and d positions to estimate 2-mer vs 3-mer.
pub fn biopython_coiled_coil_oligomer(
  sequence : String,
  heptad_register : String,
) -> Int32 {
  // Count hydrophobic residues at a/d positions
  let mut hydrophobic_ad = 0
  let mut total_ad = 0
  for i = 0; i < heptad_register.length() && i < sequence.length(); i = i + 1 {
    let pos = heptad_register[i]
    if pos == 97 || pos == 100 { // 'a' or 'd'
      total_ad = total_ad + 1
      let aa = sequence[i:i + 1].to_owned()
      if aa == "I" || aa == "L" || aa == "V" || aa == "M" || aa == "F" {
        hydrophobic_ad = hydrophobic_ad + 1
      }
    }
  }
  if total_ad == 0 {
    return 2 // default dimer
  }
  let ratio = hydrophobic_ad.to_double() / total_ad.to_double()
  // High hydrophobicity at a/d -> dimer (2)
  // Moderate -> trimer (3)
  // Lower -> tetramer (4)
  if ratio > 0.7 {
    2
  } else if ratio > 0.5 {
    3
  } else {
    4
  }
}
