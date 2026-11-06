/// Sequence masking operations.
///
/// scikit-bio equivalent: `soft_mask()`, `hard_mask()`
///
/// Soft mask: convert specified regions to lowercase
/// Hard mask: replace specified regions with N (or X for protein)

///|
/// Helper: convert char to lowercase
fn skbio_char_to_lower(c : Char) -> Char {
  if c >= 'A' && c <= 'Z' {
    (c.to_int() + 32).unsafe_to_char()
  } else {
    c
  }
}

///|
/// Helper: convert char to uppercase
fn skbio_char_to_upper(c : Char) -> Char {
  if c >= 'a' && c <= 'z' {
    (c.to_int() - 32).unsafe_to_char()
  } else {
    c
  }
}

///|
/// Soft mask a sequence: convert regions to lowercase.
/// regions: array of (start, end) pairs (0-based, end exclusive)
pub fn skbio_soft_mask(seq : String, regions : Array[(Int32, Int32)]) -> String {
  let mask_arr : Array[Bool] = Array::make(seq.length(), false)
  for r in regions {
    let s = if r.0 < 0 { 0 } else { r.0 }
    let e = if r.1 > seq.length() { seq.length() } else { r.1 }
    for i = s; i < e; i = i + 1 {
      mask_arr[i] = true
    }
  }
  let sb = StringBuilder()
  let mut i = 0
  for c in seq {
    if mask_arr[i] {
      sb.write_char(skbio_char_to_lower(c))
    } else {
      sb.write_char(skbio_char_to_upper(c))
    }
    i = i + 1
  }
  sb.to_string()
}

///|
/// Hard mask a sequence: replace masked regions with a specified character.
/// Default mask_char is 'N' for DNA, 'X' for protein.
pub fn skbio_hard_mask(
  seq : String,
  regions : Array[(Int32, Int32)],
  mask_char~ : Char = 'N',
) -> String {
  let mask_arr : Array[Bool] = Array::make(seq.length(), false)
  for r in regions {
    let s = if r.0 < 0 { 0 } else { r.0 }
    let e = if r.1 > seq.length() { seq.length() } else { r.1 }
    for i = s; i < e; i = i + 1 {
      mask_arr[i] = true
    }
  }
  let sb = StringBuilder()
  let mut i = 0
  for c in seq {
    if mask_arr[i] {
      sb.write_char(mask_char)
    } else {
      sb.write_char(c)
    }
    i = i + 1
  }
  sb.to_string()
}

///|
/// Unmask a sequence: convert all lowercase to uppercase.
pub fn skbio_unmask(seq : String) -> String {
  let sb = StringBuilder()
  for c in seq {
    sb.write_char(skbio_char_to_upper(c))
  }
  sb.to_string()
}

///|
/// Get mask regions from a sequence (positions where char is lowercase).
pub fn skbio_get_mask_regions(seq : String) -> Array[(Int32, Int32)] {
  let regions : Array[(Int32, Int32)] = Array(capacity=0)
  let mut in_mask = false
  let mut start = 0
  let mut i = 0
  for c in seq {
    let is_lower = c >= 'a' && c <= 'z'
    if is_lower && !in_mask {
      start = i
      in_mask = true
    } else if !is_lower && in_mask {
      regions.push((start, i))
      in_mask = false
    }
    i = i + 1
  }
  if in_mask {
    regions.push((start, i))
  }
  regions
}
