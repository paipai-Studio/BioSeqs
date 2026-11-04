/// AnnData type coercion for obs/var columns.
///
/// Convert string columns to numeric, boolean, or categorical.

///|
/// Try to parse a string as a double.
fn anndata_try_parse_double(s : String) -> Double? {
  let trimmed = s.trim()
  if trimmed.length() == 0 {
    return None
  }
  let chars = trimmed.to_array()
  let n = chars.length()
  if n == 0 {
    return None
  }
  let mut i = 0
  let mut sign = 1.0
  if chars[0] == '-' {
    sign = -1.0
    i = 1
  } else if chars[0] == '+' {
    i = 1
  }
  if i >= n {
    return None
  }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut has_digit = false
  let mut has_dot = false
  let mut in_frac = false
  let mut frac_div = 1.0
  while i < n {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      has_digit = true
      if in_frac {
        frac_div = frac_div * 10.0
        frac_part = frac_part + (c.to_int() - 48).to_double() / frac_div
      } else {
        int_part = int_part * 10.0 + (c.to_int() - 48).to_double()
      }
    } else if c == '.' && !has_dot {
      has_dot = true
      in_frac = true
    } else {
      return None
    }
    i = i + 1
  }
  if !has_digit {
    return None
  }
  Some(sign * (int_part + frac_part))
}

///|
/// Try to parse a string as boolean.
fn anndata_try_parse_bool(s : String) -> Bool? {
  let trimmed = s.trim().to_lower()
  if trimmed == "true" || trimmed == "1" || trimmed == "yes" {
    return Some(true)
  }
  if trimmed == "false" || trimmed == "0" || trimmed == "no" {
    return Some(false)
  }
  None
}

///|
/// Coerce an obs column to numeric (Double array).
/// Returns None if any value cannot be parsed.
pub fn anndata_obs_to_numeric(adata : AnnData, key : String) -> Array[Double]? {
  match adata.obs_column(key) {
    None => None
    Some(col) => {
      let out : Array[Double] = Array(capacity=col.length())
      for s in col {
        match anndata_try_parse_double(s) {
          Some(v) => out.push(v)
          None => return None
        }
      }
      Some(out)
    }
  }
}

///|
/// Coerce a var column to numeric.
pub fn anndata_var_to_numeric(adata : AnnData, key : String) -> Array[Double]? {
  match adata.var_column(key) {
    None => None
    Some(col) => {
      let out : Array[Double] = Array(capacity=col.length())
      for s in col {
        match anndata_try_parse_double(s) {
          Some(v) => out.push(v)
          None => return None
        }
      }
      Some(out)
    }
  }
}

///|
/// Coerce an obs column to boolean.
pub fn anndata_obs_to_bool(adata : AnnData, key : String) -> Array[Bool]? {
  match adata.obs_column(key) {
    None => None
    Some(col) => {
      let out : Array[Bool] = Array(capacity=col.length())
      for s in col {
        match anndata_try_parse_bool(s) {
          Some(v) => out.push(v)
          None => return None
        }
      }
      Some(out)
    }
  }
}

///|
/// Coerce a var column to boolean.
pub fn anndata_var_to_bool(adata : AnnData, key : String) -> Array[Bool]? {
  match adata.var_column(key) {
    None => None
    Some(col) => {
      let out : Array[Bool] = Array(capacity=col.length())
      for s in col {
        match anndata_try_parse_bool(s) {
          Some(v) => out.push(v)
          None => return None
        }
      }
      Some(out)
    }
  }
}

///|
/// Coerce a numeric array to string array (for storage in obs/var).
pub fn anndata_numeric_to_strings(values : Array[Double]) -> Array[String] {
  let out : Array[String] = Array(capacity=values.length())
  for v in values {
    out.push(v.to_string())
  }
  out
}

///|
/// Coerce a boolean array to string array.
pub fn anndata_bool_to_strings(values : Array[Bool]) -> Array[String] {
  let out : Array[String] = Array(capacity=values.length())
  for v in values {
    out.push(if v { "true" } else { "false" })
  }
  out
}

///|
/// Infer the type of a string column.
/// Returns: 0=numeric, 1=boolean, 2=categorical, 3=string
pub fn anndata_infer_type(values : Array[String]) -> Int32 {
  if values.length() == 0 {
    return 3
  }
  // Try numeric
  let mut all_numeric = true
  for s in values {
    match anndata_try_parse_double(s) {
      Some(_) => ()
      None => {
        all_numeric = false
        break
      }
    }
  }
  if all_numeric {
    return 0
  }
  // Try boolean
  let mut all_bool = true
  for s in values {
    match anndata_try_parse_bool(s) {
      Some(_) => ()
      None => {
        all_bool = false
        break
      }
    }
  }
  if all_bool {
    return 1
  }
  // Check if categorical (few unique values relative to length)
  let cat = anndata_encode_categorical(values)
  let n_unique = cat.n_categories()
  if n_unique <= 50 && n_unique < values.length() / 2 {
    return 2
  }
  3
}

///|
/// Get type name string.
pub fn anndata_type_name(type_id : Int32) -> String {
  match type_id {
    0 => "numeric"
    1 => "boolean"
    2 => "categorical"
    _ => "string"
  }
}
