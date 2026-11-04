/// Scanpy-inspired: read_csv — read CSV as AnnData
/// Enhanced version of anndata_read_csv with obs/var parsing.

///|
/// Read a CSV file as AnnData.
/// First column = row names (obs_names), remaining columns = var_names.
/// Values are parsed as doubles.
///
/// Parameters:
/// - `content`: CSV file content
/// - `delimiter`: column delimiter (default ",")
/// - `first_col_is_index`: if true, first column is obs names (default true)
pub fn scanpy_read_csv(
  content : String,
  delimiter~ : String = ",",
  first_col_is_index~ : Bool = true,
) -> AnnData {
  let lines = content.split("\n").to_array()
  if lines.length() == 0 {
    return AnnData::new([], [], [])
  }
  let header_parts = lines[0]
    .to_owned()
    .trim()
    .to_owned()
    .split(delimiter)
    .to_array()
  let var_names : Array[String] = []
  let mut data_start = 0
  if header_parts.length() > 0 {
    let first = header_parts[0].to_owned().trim().to_owned()
    let is_header = !scanpy_csv_is_number(first)
    if is_header {
      data_start = 1
      let start = if first_col_is_index { 1 } else { 0 }
      for i in start..<header_parts.length() {
        var_names.push(header_parts[i].to_owned().trim().to_owned())
      }
    }
  }
  let x : Array[Array[Double]] = []
  let obs_names : Array[String] = []
  for line_idx in data_start..<lines.length() {
    let trimmed = lines[line_idx].to_owned().trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    let parts = trimmed.split(delimiter).to_array()
    let mut col_start = 0
    if first_col_is_index && parts.length() > 0 {
      obs_names.push(parts[0].to_owned().trim().to_owned())
      col_start = 1
    } else {
      obs_names.push("row_" + line_idx.to_string())
    }
    let row : Array[Double] = []
    for j in col_start..<parts.length() {
      let val = match parse_double(parts[j].to_owned().trim().to_owned()) {
        Some(v) => v
        None => 0.0
      }
      row.push(val)
    }
    x.push(row)
  }
  if var_names.length() == 0 && x.length() > 0 {
    for j in 0..<x[0].length() {
      var_names.push("col_" + j.to_string())
    }
  }
  AnnData::new(x, obs_names, var_names)
}

///|
fn scanpy_csv_is_number(s : String) -> Bool {
  match parse_double(s) {
    Some(_) => true
    None => false
  }
}

///|
/// Read CSV with obs annotations (metadata) separately.
/// Parameters:
/// - `x_content`: CSV with expression matrix
/// - `obs_content`: CSV with cell metadata
/// - `var_content`: CSV with gene metadata
pub fn scanpy_read_csv_with_annotations(
  x_content : String,
  obs_content : String,
  var_content : String,
) -> AnnData {
  let adata = scanpy_read_csv(x_content)
  let obs_lines = obs_content.split("\n").to_array()
  if obs_lines.length() > 1 {
    let header = obs_lines[0].to_owned().trim().to_owned().split(",").to_array()
    for col_idx in 0..<header.length() {
      let col_name = header[col_idx].to_owned().trim().to_owned()
      let vals : Array[String] = []
      for line_idx in 1..<obs_lines.length() {
        let parts = obs_lines[line_idx]
          .to_owned()
          .trim()
          .to_owned()
          .split(",")
          .to_array()
        if col_idx < parts.length() {
          vals.push(parts[col_idx].to_owned().trim().to_owned())
        } else {
          vals.push("")
        }
      }
      ignore(adata.set_obs(col_name, vals))
    }
  }
  let var_lines = var_content.split("\n").to_array()
  if var_lines.length() > 1 {
    let header = var_lines[0].to_owned().trim().to_owned().split(",").to_array()
    for col_idx in 0..<header.length() {
      let col_name = header[col_idx].to_owned().trim().to_owned()
      let vals : Array[String] = []
      for line_idx in 1..<var_lines.length() {
        let parts = var_lines[line_idx]
          .to_owned()
          .trim()
          .to_owned()
          .split(",")
          .to_array()
        if col_idx < parts.length() {
          vals.push(parts[col_idx].to_owned().trim().to_owned())
        } else {
          vals.push("")
        }
      }
      ignore(adata.set_var(col_name, vals))
    }
  }
  adata
}
