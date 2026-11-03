/// AnnData CSV I/O.
///
/// `anndata.read_csv()` / `adata.write_csvs()` equivalent.
///
/// Reads expression matrix from CSV, writes obs/var as CSV.

///|
/// Parse a CSV string into a 2D array of strings.
fn anndata_parse_csv(content : String) -> Array[Array[String]] {
  let lines = content.split("\n").to_array()
  let out : Array[Array[String]] = Array(capacity=0)
  for line in lines {
    let trimmed = line.to_string().trim()
    if trimmed.length() == 0 {
      continue
    }
    let fields = trimmed.split(",").to_array()
    let row : Array[String] = Array(capacity=0)
    for f in fields {
      row.push(f.to_string().trim().to_string())
    }
    out.push(row)
  }
  out
}

///|
/// Read AnnData from a CSV expression matrix.
/// First row: header with var_names (gene names).
/// First column: obs_names (cell names).
/// Remaining cells: expression values.
pub fn anndata_read_csv(content : String) -> AnnData? {
  let rows = anndata_parse_csv(content)
  if rows.length() < 2 {
    return None
  }
  if rows[0].length() < 2 {
    return None
  }
  // Header: var_names (skip first column)
  let var_names : Array[String] = Array(capacity=0)
  for j = 1; j < rows[0].length(); j = j + 1 {
    var_names.push(rows[0][j])
  }
  let n_vars = var_names.length()
  let obs_names : Array[String] = Array(capacity=0)
  let x : Array[Array[Double]] = Array(capacity=0)
  for i = 1; i < rows.length(); i = i + 1 {
    if rows[i].length() < 1 {
      continue
    }
    obs_names.push(rows[i][0])
    let row : Array[Double] = Array::make(n_vars, 0.0)
    for j = 0; j < n_vars && j + 1 < rows[i].length(); j = j + 1 {
      row[j] = anndata_parse_double_or_zero(rows[i][j + 1])
    }
    x.push(row)
  }
  Some(AnnData::new(x, obs_names, var_names))
}

///|
/// Write obs annotation as CSV.
pub fn anndata_write_obs_csv(adata : AnnData) -> String {
  let obs_keys = adata.obs.keys().collect()
  let mut header = "obs_name"
  for k in obs_keys {
    header = header + "," + k
  }
  header = header + "\n"
  let mut out = header
  let obs_names = adata.get_obs_names()
  for i = 0; i < obs_names.length(); i = i + 1 {
    out = out + obs_names[i]
    for k in obs_keys {
      match adata.obs_column(k) {
        Some(col) =>
          if i < col.length() {
            out = out + "," + col[i]
          } else {
            out = out + ","
          }
        None => out = out + ","
      }
    }
    out = out + "\n"
  }
  out
}

///|
/// Write var annotation as CSV.
pub fn anndata_write_var_csv(adata : AnnData) -> String {
  let var_keys = adata.variables.keys().collect()
  let mut header = "var_name"
  for k in var_keys {
    header = header + "," + k
  }
  header = header + "\n"
  let mut out = header
  let var_names = adata.get_var_names()
  for i = 0; i < var_names.length(); i = i + 1 {
    out = out + var_names[i]
    for k in var_keys {
      match adata.var_column(k) {
        Some(col) =>
          if i < col.length() {
            out = out + "," + col[i]
          } else {
            out = out + ","
          }
        None => out = out + ","
      }
    }
    out = out + "\n"
  }
  out
}

///|
/// Write X matrix as CSV.
pub fn anndata_write_x_csv(adata : AnnData) -> String {
  let var_names = adata.get_var_names()
  let mut header = ""
  for j = 0; j < var_names.length(); j = j + 1 {
    if j > 0 {
      header = header + ","
    }
    header = header + var_names[j]
  }
  header = header + "\n"
  let mut out = header
  for i = 0; i < adata.n_obs(); i = i + 1 {
    for j = 0; j < adata.n_vars(); j = j + 1 {
      if j > 0 {
        out = out + ","
      }
      out = out + adata.get_x(i, j).to_string()
    }
    out = out + "\n"
  }
  out
}

///|
/// Write all components as separate CSV strings.
/// Returns (obs_csv, var_csv, x_csv).
pub fn anndata_write_csvs(adata : AnnData) -> (String, String, String) {
  (
    anndata_write_obs_csv(adata),
    anndata_write_var_csv(adata),
    anndata_write_x_csv(adata),
  )
}

///|
/// Parse a double from string, return 0.0 on failure.
fn anndata_parse_double_or_zero(s : String) -> Double {
  let trimmed = s.trim()
  if trimmed.length() == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let chars = trimmed.to_array()
  let n = chars.length()
  if n == 0 {
    return 0.0
  }
  if chars[0] == '-' {
    sign = -1.0
    i = 1
  } else if chars[0] == '+' {
    i = 1
  }
  let mut has_digit = false
  let mut has_dot = false
  while i < n {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      result = result * 10.0 + (c.to_int() - 48).to_double()
      has_digit = true
    } else if c == '.' && !has_dot {
      has_dot = true
    } else {
      return 0.0
    }
    i = i + 1
  }
  if !has_digit {
    return 0.0
  }
  // Handle decimal part
  if has_dot {
    let mut int_part = 0.0
    let mut frac_part = 0.0
    let mut in_frac = false
    let mut frac_div = 1.0
    let mut j = 0
    if chars[0] == '-' || chars[0] == '+' {
      j = 1
    }
    while j < n {
      let c = chars[j]
      if c == '.' {
        in_frac = true
      } else if c >= '0' && c <= '9' {
        if in_frac {
          frac_div = frac_div * 10.0
          frac_part = frac_part + (c.to_int() - 48).to_double() / frac_div
        } else {
          int_part = int_part * 10.0 + (c.to_int() - 48).to_double()
        }
      }
      j = j + 1
    }
    result = int_part + frac_part
  }
  sign * result
}
