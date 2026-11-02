///|
/// Sequence complexity and low-complexity region detection.
///
/// Complements existing `seq_complexity.mbt` and `biostrings.mbt` with
/// SEG masking, compositional complexity, and sequence logo information.

///|
/// Log base 2 helper.
fn sc_log2(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  @math.ln(x) / 0.6931471805599453
}

///|
/// Extract a substring of length k starting at position i.
fn sc_substring(s : String, start : Int32, length : Int32) -> String {
  let mut result = ""
  for i in 0..<length {
    if start + i < s.length() {
      result = result + s.unsafe_get(start + i).unsafe_to_char().to_string()
    }
  }
  result
}

///|
/// Detect low-complexity regions using a sliding window.
///
/// A region is low-complexity if its Shannon entropy is below threshold.
/// Returns (start, end) pairs for low-complexity segments.
pub fn find_low_complexity_regions(
  sequence : String,
  window_size~ : Int32 = 12,
  threshold~ : Double = 2.0,
) -> Array[(Int32, Int32)] {
  let n = sequence.length()
  let regions : Array[(Int32, Int32)] = []
  let mut in_region = false
  let mut region_start = 0
  for i in 0..<=(n - window_size) {
    let window = sc_substring(sequence, i, window_size)
    let entropy = sc_window_entropy(window)
    if entropy < threshold {
      if !in_region {
        region_start = i
        in_region = true
      }
    } else if in_region {
      regions.push((region_start, i + window_size - 2))
      in_region = false
    }
  }
  if in_region {
    regions.push((region_start, n - 1))
  }
  regions
}

///|
/// Calculate Shannon entropy of a window.
fn sc_window_entropy(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let counts : Map[UInt16, Int32] = Map([], capacity=32)
  for i in 0..<n {
    let c = sequence.unsafe_get(i)
    counts[c] = counts.get(c).unwrap_or(0) + 1
  }
  let mut entropy = 0.0
  let keys = counts.keys().collect()
  for k in keys {
    let p = counts[k].to_double() / n.to_double()
    entropy = entropy - p * sc_log2(p)
  }
  entropy
}

///|
/// Calculate compositional complexity (Wootton-Federhen complexity score).
///
/// This is the basis of the SEG algorithm.
pub fn compositional_complexity(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let counts : Map[UInt16, Int32] = Map([], capacity=32)
  for i in 0..<n {
    let c = sequence.unsafe_get(i)
    counts[c] = counts.get(c).unwrap_or(0) + 1
  }
  let mut score = 0.0
  let keys = counts.keys().collect()
  for k in keys {
    let ni = counts[k].to_double()
    let ratio = n.to_double() / ni
    score = score + ni * sc_log2(ratio)
  }
  score / n.to_double()
}

///|
/// SEG algorithm for masking low-complexity regions.
///
/// Returns a sequence with low-complexity regions replaced by 'X'.
pub fn seg_mask(
  sequence : String,
  window_size~ : Int32 = 12,
  threshold~ : Double = 2.2,
) -> String {
  let regions = find_low_complexity_regions(sequence, window_size~, threshold~)
  let chars : Array[Char] = []
  for i in 0..<sequence.length() {
    chars.push(sequence.unsafe_get(i).unsafe_to_char())
  }
  for region in regions {
    for i in region.0..<=region.1 {
      if i < chars.length() {
        chars[i] = 'X'
      }
    }
  }
  let mut result = ""
  for c in chars {
    result = result + c.to_string()
  }
  result
}

///|
/// Calculate sequence logo information content.
///
/// For each position, computes R_i = log2(alphabet_size) - H_i.
pub fn sequence_logo_information(
  alignment : Array[String],
  alphabet_size~ : Int32 = 20,
) -> Array[Double] {
  if alignment.length() == 0 {
    return []
  }
  let n = alignment[0].length()
  let n_seqs = alignment.length()
  let result : Array[Double] = []
  for pos in 0..<n {
    let counts : Map[UInt16, Int32] = Map([], capacity=32)
    for seq in alignment {
      if pos < seq.length() {
        let c = seq.unsafe_get(pos)
        counts[c] = counts.get(c).unwrap_or(0) + 1
      }
    }
    let mut entropy = 0.0
    let keys = counts.keys().collect()
    for k in keys {
      let p = counts[k].to_double() / n_seqs.to_double()
      entropy = entropy - p * sc_log2(p)
    }
    let info = sc_log2(alphabet_size.to_double()) - entropy
    result.push(info)
  }
  result
}
