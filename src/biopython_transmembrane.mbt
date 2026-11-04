/// biopython_transmembrane.mbt — Transmembrane helix prediction
/// Inspired by Biopython TMHMM-like heuristics

///|
/// Transmembrane segment prediction.
pub struct BiopythonTMSegment {
  start : Int32 // start position (0-based)
  end : Int32 // end position (exclusive)
  length : Int32 // segment length
  score : Double // hydrophobicity score
  is_tm : Bool // is transmembrane
}

///|
/// Kyte-Doolittle hydrophobicity.
fn biopython_tm_kd(aa : String) -> Double {
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
/// Calculate hydrophobicity profile using sliding window.
pub fn biopython_hydrophobicity_profile(
  sequence : String,
  window~ : Int32 = 19, // typical TM helix window
) -> Array[Double] {
  let n = sequence.length()
  if n == 0 {
    return []
  }
  let profile : Array[Double] = Array::make(n, 0.0)
  let half = window / 2
  for i = 0; i < n; i = i + 1 {
    let start = if i > half { i - half } else { 0 }
    let end = if i + half + 1 < n { i + half + 1 } else { n }
    let mut sum = 0.0
    let count = end - start
    for j = start; j < end; j = j + 1 {
      sum = sum + biopython_tm_kd(sequence[j:j + 1].to_owned())
    }
    profile[i] = if count > 0 { sum / count.to_double() } else { 0.0 }
  }
  profile
}

///|
/// Predict transmembrane segments.
/// Uses hydrophobicity threshold and minimum segment length.
pub fn biopython_predict_transmembrane(
  sequence : String,
  threshold~ : Double = 0.0, // hydrophobicity threshold
  min_length~ : Int32 = 18, // minimum TM helix length
  max_length~ : Int32 = 33, // maximum TM helix length
) -> Array[BiopythonTMSegment] {
  let profile = biopython_hydrophobicity_profile(sequence)
  let n = profile.length()
  let segments : Array[BiopythonTMSegment] = Array(capacity=0)
  let mut i = 0
  while i < n {
    if profile[i] > threshold {
      // Start of hydrophobic region
      let start = i
      let mut end = i
      while end < n && profile[end] > threshold {
        end = end + 1
      }
      let length = end - start
      // Calculate average score
      let mut avg = 0.0
      for j = start; j < end; j = j + 1 {
        avg = avg + profile[j]
      }
      avg = avg / length.to_double()
      // Check length constraints
      if length >= min_length && length <= max_length {
        segments.push({ start, end, length, score: avg, is_tm: true, })
      } else if length > max_length {
        // Split long segments
        let mut s = start
        while s + min_length <= end {
          let e = if s + max_length < end { s + max_length } else { end }
          let mut seg_avg = 0.0
          for j = s; j < e; j = j + 1 {
            seg_avg = seg_avg + profile[j]
          }
          seg_avg = seg_avg / (e - s).to_double()
          segments.push({
            start: s,
            end: e,
            length: e - s,
            score: seg_avg,
            is_tm: true,
          })
          s = e
        }
      }
      i = end
    } else {
      i = i + 1
    }
  }
  segments
}

///|
/// Count number of transmembrane segments.
pub fn biopython_tm_count(
  sequence : String,
  threshold~ : Double = 0.0,
) -> Int32 {
  biopython_predict_transmembrane(sequence, threshold~).length()
}

///|
/// Classify protein as membrane/non-membrane.
pub fn biopython_tm_classify(
  sequence : String,
  threshold~ : Double = 0.0,
) -> String {
  let count = biopython_tm_count(sequence, threshold~)
  if count == 0 {
    "non-membrane"
  } else if count == 1 {
    "single-pass"
  } else {
    "multi-pass"
  }
}

///|
/// Get inside/outside topology prediction.
/// Alternates inside/outside for each TM segment.
pub fn biopython_tm_topology(
  sequence : String,
  threshold~ : Double = 0.0,
) -> Array[(Int32, Int32, String)] {
  let segments = biopython_predict_transmembrane(sequence, threshold~)
  let topology : Array[(Int32, Int32, String)] = Array(capacity=0)
  let mut inside = true // N-terminus starts inside
  for i = 0; i < segments.length(); i = i + 1 {
    let seg = segments[i]
    // Non-TM region before this segment
    let region = if inside { "inside" } else { "outside" }
    if i == 0 {
      if seg.start > 0 {
        topology.push((0, seg.start, region))
      }
    } else {
      let prev_end = segments[i - 1].end
      if seg.start > prev_end {
        topology.push((prev_end, seg.start, region))
      }
    }
    // TM segment
    topology.push((seg.start, seg.end, "TM"))
    inside = !inside
  }
  // Last non-TM region
  if segments.length() > 0 {
    let last_end = segments[segments.length() - 1].end
    if last_end < sequence.length() {
      let region = if inside { "inside" } else { "outside" }
      topology.push((last_end, sequence.length(), region))
    }
  }
  topology
}
