///|
/// biotite_application_tantan.mbt — low-complexity masking in the spirit
/// of TANTAN, implemented with a sliding-window Shannon entropy.

///|
/// Shannon entropy (bits) of a window, ignoring gaps.
fn biotite_tantan_window_entropy(
  sequence : String,
  start : Int32,
  stop : Int32,
) -> Double {
  let counts : Array[Int32] = Array::make(128, 0)
  let mut total = 0
  for i in start..<stop {
    if i >= 0 && i < sequence.length() {
      let code = sequence[i].to_int()
      if code >= 0 && code < 128 {
        counts[code] = counts[code] + 1
        total = total + 1
      }
    }
  }
  if total == 0 {
    return 0.0
  }
  let mut h = 0.0
  for i in 0..<128 {
    if counts[i] > 0 {
      let p = counts[i].to_double() / total.to_double()
      h = h - p * @math.ln(p) / @math.ln(2.0)
    }
  }
  h
}

///|
/// Per-position window entropy (bits).
pub fn biotite_tantan_entropy_profile(
  sequence : String,
  window~ : Int32 = 20,
) -> Array[Double] {
  let n = sequence.length()
  let result : Array[Double] = Array(capacity=n)
  let w = if window < 1 { 1 } else { window }
  let half = w / 2
  for i in 0..<n {
    let start = if i - half < 0 { 0 } else { i - half }
    let stop = if start + w > n { n } else { start + w }
    result.push(biotite_tantan_window_entropy(sequence, start, stop))
  }
  result
}

///|
/// Mask low-complexity positions with 'N'.
pub fn biotite_tantan_mask(
  sequence : String,
  window~ : Int32 = 20,
  threshold~ : Double = 1.0,
) -> String {
  let profile = biotite_tantan_entropy_profile(sequence, window~)
  let sb = StringBuilder()
  for i in 0..<sequence.length() {
    if profile[i] < threshold {
      let _ = sb.write_char('N')
    } else {
      let _ = sb.write(sequence[i:i + 1])
    }
  }
  sb.to_string()
}

///|
/// Command-line construction for the external TANTAN binary.
pub fn biotite_tantan_command(
  input_file : String,
  output_file : String,
  window~ : Int32 = 20,
) -> String {
  "tantan -in " +
  input_file +
  " -out " +
  output_file +
  " -window " +
  window.to_string()
}
