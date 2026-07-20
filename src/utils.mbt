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
      result[idx] = s[start:i].to_string()
      idx = idx + 1
      start = i + 1
    }
    i = i + 1
  }
  if start < n {
    result[idx] = s[start:n].to_string()
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
/// Parse double from String, supporting scientific notation (e.g., 1e-50, 1e-30).
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
    if c == 'e'.to_int() || c == 'E'.to_int() {
      break
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

  let mut exponent = 0
  let mut exp_sign = 1
  if i < n &&
    (
      s.unsafe_get(i).to_int() == 'e'.to_int() ||
      s.unsafe_get(i).to_int() == 'E'.to_int()
    ) {
    i = i + 1
    if i < n && s.unsafe_get(i).to_int() == '-'.to_int() {
      exp_sign = -1
      i = i + 1
    } else if i < n && s.unsafe_get(i).to_int() == '+'.to_int() {
      i = i + 1
    }
    while i < n {
      let c = s.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 {
        break
      }
      exponent = exponent * 10 + (c - zero)
      i = i + 1
    }
    exponent = exponent * exp_sign
  }

  if i > 0 || (sign < 0.0 && i > 1) {
    Some(sign * result * @math.pow(10.0, exponent.to_double()))
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

///|
/// Get character at position.
pub fn char_at(s : String, idx : Int) -> String {
  let c = s.unsafe_get(idx).to_int()
  if c >= 'A'.to_int() && c <= 'Z'.to_int() {
    let mut result = ""
    for i = 0; i < 26; i = i + 1 {
      if c == 'A'.to_int() + i {
        result = result +
          [
            "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
            "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
          ][i]
        break
      }
    }
    result
  } else if c >= 'a'.to_int() && c <= 'z'.to_int() {
    let mut result = ""
    for i = 0; i < 26; i = i + 1 {
      if c == 'a'.to_int() + i {
        result = result +
          [
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
            "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
          ][i]
        break
      }
    }
    result
  } else if c >= '0'.to_int() && c <= '9'.to_int() {
    let mut result = ""
    for i = 0; i < 10; i = i + 1 {
      if c == '0'.to_int() + i {
        result = result + ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"][i]
        break
      }
    }
    result
  } else if c == ' '.to_int() {
    " "
  } else if c == '\t'.to_int() {
    "\t"
  } else if c == '\n'.to_int() {
    "\n"
  } else if c == '\r'.to_int() {
    "\r"
  } else if c == '>'.to_int() {
    ">"
  } else if c == '#'.to_int() {
    "#"
  } else if c == '/'.to_int() {
    "/"
  } else if c == ':'.to_int() {
    ":"
  } else if c == ';'.to_int() {
    ";"
  } else if c == ','.to_int() {
    ","
  } else if c == '('.to_int() {
    "("
  } else if c == ')'.to_int() {
    ")"
  } else if c == '.'.to_int() {
    "."
  } else if c == '-'.to_int() {
    "-"
  } else if c == '_'.to_int() {
    "_"
  } else if c == '['.to_int() {
    "["
  } else if c == ']'.to_int() {
    "]"
  } else if c == '@'.to_int() {
    "@"
  } else if c == '&'.to_int() {
    "&"
  } else {
    c.to_string()
  }
}

///|
/// Substring function.
pub fn substring(s : String, start : Int, len : Int) -> String {
  let mut result = ""
  for i = 0; i < len && start + i < s.length(); i = i + 1 {
    result = result + char_at(s, start + i)
  }
  result
}

///|
/// Check if string starts with prefix at position.
pub fn starts_with(s : String, pos : Int, prefix : String) -> Bool {
  if pos + prefix.length() > s.length() {
    return false
  }

  for i = 0; i < prefix.length(); i = i + 1 {
    if s.unsafe_get(pos + i) != prefix.unsafe_get(i) {
      return false
    }
  }

  true
}

///|
/// Convert string to lowercase.
pub fn to_lowercase(s : String) -> String {
  let mut result = ""
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c >= 'A'.to_int() && c <= 'Z'.to_int() {
      result = result + (c + 32).to_string()
    } else {
      result = result + c.to_string()
    }
    i = i + 1
  }
  result
}

///|
/// Check if string contains substring.
pub fn contains(s : String, substr : String) -> Bool {
  if substr.length() == 0 {
    return true
  }
  if s.length() < substr.length() {
    return false
  }
  let mut i = 0
  while i <= s.length() - substr.length() {
    let mut match_found = true
    let mut j = 0
    while j < substr.length() {
      if s.unsafe_get(i + j) != substr.unsafe_get(j) {
        match_found = false
        break
      }
      j = j + 1
    }
    if match_found {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// Trim whitespace from both ends of string.
pub fn trim(s : String) -> String {
  let mut start = 0
  while start < s.length() && s.unsafe_get(start).to_int() == ' '.to_int() {
    start = start + 1
  }
  let mut end = s.length()
  while end > start && s.unsafe_get(end - 1).to_int() == ' '.to_int() {
    end = end - 1
  }
  s[start:end].to_string()
}
