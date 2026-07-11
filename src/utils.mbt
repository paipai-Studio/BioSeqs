///|
/// Utility functions for string parsing and manipulation.

///|
/// Split string by character code.
pub fn split_by_char(s : String, delim : Int) -> Array[String] {
  let n = s.length()
  let mut count = 1
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == delim {
      count = count + 1
    }
    i = i + 1
  }

  let result = Array::make(count, "")
  let mut idx = 0
  let mut start = 0
  i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == delim {
      result[idx] = s[start:i].to_owned()
      idx = idx + 1
      start = i + 1
    }
    i = i + 1
  }
  if start < n {
    result[idx] = s[start:n].to_owned()
  }
  result
}

///|
/// Parse integer from String.
pub fn parse_int(s : String) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      break
    }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Parse double from String.
pub fn parse_double(s : String) -> Double? {
  if s == "." || s.length() == 0 {
    return None
  }
  let n = s.length()
  let mut result = 0.0
  let mut decimal = 0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()

  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }

  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      decimal = 1
      i = i + 1
      continue
    }
    if c < zero || c > zero + 9 {
      break
    }
    let digit = (c - zero).to_double()
    if decimal == 0 {
      result = result * 10.0 + digit
    } else {
      result = result + digit / @math.pow(10.0, decimal.to_double())
      decimal = decimal + 1
    }
    i = i + 1
  }

  if i > 0 || (sign < 0.0 && i > 1) {
    Some(sign * result)
  } else {
    None
  }
}

///|
/// Split string by tab character.
pub fn split_by_tab(s : String) -> Array[String] {
  split_by_char(s, 9)
}

///|
/// Parse comma-separated string into array.
pub fn parse_comma_list(s : String) -> Array[String] {
  if s == "." || s.length() == 0 {
    return []
  }
  split_by_char(s, ','.to_int())
}

///|
/// Parse semicolon-separated string into array.
pub fn parse_semicolon_list(s : String) -> Array[String] {
  if s == "." || s.length() == 0 {
    return []
  }
  split_by_char(s, ';'.to_int())
}
