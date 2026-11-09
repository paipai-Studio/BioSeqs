/// anndata_io_excel.mbt — Read AnnData from Excel-like tabular data
/// Item 3: read_excel

///|
/// Read AnnData from Excel-like tabular data (TSV proxy for xlsx).
/// Format: first sheet = X matrix, second sheet = obs, third sheet = var.
/// Each sheet is a tab-separated block separated by "---SHEET---".
pub fn anndata_read_excel(
  content : String,
  sheet_x~ : Int32 = 0,
  sheet_obs~ : Int32 = 1,
  sheet_var~ : Int32 = 2,
) -> AnnData? {
  let sheets = content.split("---SHEET---").to_array()
  if sheets.length() == 0 {
    return None
  }
  // Parse X sheet
  let x_lines = if sheet_x < sheets.length() {
    sheets[sheet_x].trim().split("\n").to_array()
  } else {
    Array(capacity=0)
  }
  let x : Array[Array[Double]] = Array(capacity=x_lines.length())
  let mut n_vars = 0
  for line in x_lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    let fields = trimmed.split("\t").to_array()
    let row : Array[Double] = Array(capacity=fields.length())
    for f in fields {
      row.push(parse_double(f.to_owned()).unwrap_or(0.0))
    }
    if n_vars == 0 {
      n_vars = row.length()
    }
    x.push(row)
  }
  let n_obs = x.length()
  // Build names
  let obs_names : Array[String] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    obs_names.push("cell_" + i.to_string())
  }
  let var_names : Array[String] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    var_names.push("gene_" + j.to_string())
  }
  let adata = AnnData::new(x, obs_names, var_names)
  // Parse obs sheet
  if sheet_obs < sheets.length() {
    let obs_lines = sheets[sheet_obs].trim().split("\n").to_array()
    if obs_lines.length() > 1 {
      let header = obs_lines[0].split("\t").to_array()
      for c = 0; c < header.length(); c = c + 1 {
        let col_name = header[c].to_owned()
        let vals : Array[String] = Array(capacity=obs_lines.length() - 1)
        for r = 1; r < obs_lines.length(); r = r + 1 {
          let fields = obs_lines[r].split("\t").to_array()
          if c < fields.length() {
            vals.push(fields[c].to_owned())
          } else {
            vals.push("")
          }
        }
        adata.obs[col_name] = vals
      }
    }
  }
  // Parse var sheet
  if sheet_var < sheets.length() {
    let var_lines = sheets[sheet_var].trim().split("\n").to_array()
    if var_lines.length() > 1 {
      let header = var_lines[0].split("\t").to_array()
      for c = 0; c < header.length(); c = c + 1 {
        let col_name = header[c].to_owned()
        let vals : Array[String] = Array(capacity=var_lines.length() - 1)
        for r = 1; r < var_lines.length(); r = r + 1 {
          let fields = var_lines[r].split("\t").to_array()
          if c < fields.length() {
            vals.push(fields[c].to_owned())
          } else {
            vals.push("")
          }
        }
        adata.variables[col_name] = vals
      }
    }
  }
  Some(adata)
}

///|
/// Write AnnData to Excel-like multi-sheet TSV.
pub fn anndata_write_excel(adata : AnnData) -> String {
  let sb = StringBuilder()
  // Sheet 0: X matrix
  let x = adata.x_matrix()
  for i = 0; i < x.length(); i = i + 1 {
    for j = 0; j < x[i].length(); j = j + 1 {
      if j > 0 {
        let _ = sb.write("\t")
      }
      let _ = sb.write(x[i][j].to_string())
    }
    let _ = sb.write("\n")
  }
  let _ = sb.write("---SHEET---\n")
  // Sheet 1: obs
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
  }
  let _ = sb.write("\n")
  for i = 0; i < adata.n_obs(); i = i + 1 {
    for k in obs_keys {
      let vals = adata.obs[k]
      if i < vals.length() {
        let _ = sb.write(vals[i])
      }
      let _ = sb.write("\t")
    }
    let _ = sb.write("\n")
  }
  let _ = sb.write("---SHEET---\n")
  // Sheet 2: var
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
  }
  let _ = sb.write("\n")
  for j = 0; j < adata.n_vars(); j = j + 1 {
    for k in var_keys {
      let vals = adata.variables[k]
      if j < vals.length() {
        let _ = sb.write(vals[j])
      }
      let _ = sb.write("\t")
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}
