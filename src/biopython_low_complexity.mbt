/// biopython_low_complexity.mbt — Low complexity region detection
/// Inspired by Biopython SEG algorithm

///|
/// Low complexity segment.
pub struct BiopythonLowComplexity {
  start : Int32
  end : Int32
  length : Int32
  complexity : Double // Shannon entropy
  composition : Array[(String, Double)] // (residue, frequency)
}

///|
/// Calculate Shannon entropy of a sequence region.
fn biopython_lc_entropy(region : String) -> Double {
  let n = region.length()
  if n == 0 {
    return 0.0
  }
  // Count residue frequencies
  let chars : Array[String] = Array(capacity=0)
  let counts : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let ch = region[i:i + 1].to_owned()
    let mut found = false
    for j = 0; j < chars.length(); j = j + 1 {
      if chars[j] == ch {
        counts[j] = counts[j] + 1
        found = true
        break
      }
    }
    if !found {
      chars.push(ch)
      counts.push(1)
    }
  }
  // Shannon entropy: H = -sum(p * log2(p))
  let mut entropy = 0.0
  for j = 0; j < chars.length(); j = j + 1 {
    let p = counts[j].to_double() / n.to_double()
    if p > 0.0 {
      entropy = entropy - p * @math.log2(p)
    }
  }
  entropy
}

///|
/// Get composition of a region.
fn biopython_lc_composition(region : String) -> Array[(String, Double)] {
  let n = region.length()
  let chars : Array[String] = Array(capacity=0)
  let counts : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let ch = region[i:i + 1].to_owned()
    let mut found = false
    for j = 0; j < chars.length(); j = j + 1 {
      if chars[j] == ch {
        counts[j] = counts[j] + 1
        found = true
        break
      }
    }
    if !found {
      chars.push(ch)
      counts.push(1)
    }
  }
  let result : Array[(String, Double)] = Array(capacity=0)
  for j = 0; j < chars.length(); j = j + 1 {
    result.push((chars[j], counts[j].to_double() / n.to_double()))
  }
  // Sort by frequency descending
  result.sort_by(fn(a, b) {
    if b.1 < a.1 {
      -1
    } else if b.1 > a.1 {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Find low complexity regions using a sliding window approach.
/// Regions with entropy below threshold are flagged.
pub fn biopython_find_low_complexity(
  sequence : String,
  window~ : Int32 = 12,
  threshold~ : Double = 2.0, // max entropy for "low complexity"
  min_length~ : Int32 = 6,
) -> Array[BiopythonLowComplexity] {
  let n = sequence.length()
  if n < window {
    return []
  }
  // Calculate entropy at each position
  let is_low : Array[Bool] = Array::make(n, false)
  let half = window / 2
  for i = 0; i < n; i = i + 1 {
    let start = if i > half { i - half } else { 0 }
    let end = if i + half + 1 < n { i + half + 1 } else { n }
    let region = sequence[start:end].to_owned()
    let entropy = biopython_lc_entropy(region)
    if entropy < threshold {
      is_low[i] = true
    }
  }
  // Merge contiguous low complexity positions
  let results : Array[BiopythonLowComplexity] = Array(capacity=0)
  let mut i = 0
  while i < n {
    if is_low[i] {
      let start = i
      let mut end = i
      while end < n && is_low[end] {
        end = end + 1
      }
      let length = end - start
      if length >= min_length {
        let region = sequence[start:end].to_owned()
        results.push({
          start,
          end,
          length,
          complexity: biopython_lc_entropy(region),
          composition: biopython_lc_composition(region),
        })
      }
      i = end
    } else {
      i = i + 1
    }
  }
  results
}

///|
/// Mask low complexity regions (replace with 'X').
pub fn biopython_mask_low_complexity(
  sequence : String,
  window~ : Int32 = 12,
  threshold~ : Double = 2.0,
) -> String {
  let regions = biopython_find_low_complexity(sequence, window~, threshold~)
  let chars : Array[Int32] = Array::make(sequence.length(), 0)
  for i = 0; i < sequence.length(); i = i + 1 {
    chars[i] = sequence[i].to_int()
  }
  for r = 0; r < regions.length(); r = r + 1 {
    for i = regions[r].start; i < regions[r].end; i = i + 1 {
      chars[i] = 88 // 'X'
    }
  }
  let sb = StringBuilder()
  for i = 0; i < chars.length(); i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Overall sequence complexity score.
pub fn biopython_sequence_complexity(sequence : String) -> Double {
  biopython_lc_entropy(sequence)
}

///|
/// Check if a sequence is a tandem repeat.
/// Simple check: does the sequence repeat itself?
pub fn biopython_is_tandem_repeat(sequence : String) -> (Bool, Int32) {
  let n = sequence.length()
  for unit_len = 1; unit_len <= n / 2; unit_len = unit_len + 1 {
    if n % unit_len == 0 {
      let unit = sequence[0:unit_len].to_owned()
      let mut is_repeat = true
      for i = unit_len; i < n; i = i + unit_len {
        if sequence[i:i + unit_len].to_owned() != unit {
          is_repeat = false
          break
        }
      }
      if is_repeat {
        return (true, n / unit_len)
      }
    }
  }
  (false, 1)
}
