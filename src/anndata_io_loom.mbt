/// anndata_io_loom.mbt — Loom file format read/write
/// Item 1: read_loom / write_loom

///|
/// Loom file format: HDF5-based with /matrix, /row_attrs, /col_attrs groups.
/// We implement a text-based proxy (TSV-like) since we don't have HDF5.

///|
pub struct LoomFile {
  matrix : Array[Array[Double]]
  row_attrs : Map[String, Array[String]]
  col_attrs : Map[String, Array[String]]
  row_names : Array[String]
  col_names : Array[String]
}

///|
/// Write AnnData to Loom-format text representation.
/// Loom convention: rows = genes (var), columns = cells (obs).
pub fn anndata_write_loom(adata : AnnData) -> String {
  let sb = StringBuilder()
  let _ = sb.write("loom\t")
  let _ = sb.write(adata.n_obs().to_string())
  let _ = sb.write("\t")
  let _ = sb.write(adata.n_vars().to_string())
  let _ = sb.write("\n")
  // Row attributes (var)
  let var_keys = adata.variables.keys().collect()
  let _ = sb.write("row_attrs\t")
  let _ = sb.write(var_keys.length().to_string())
  let _ = sb.write("\n")
  for k in var_keys {
    let _ = sb.write(k)
    let vals = adata.variables[k]
    for v in vals {
      let _ = sb.write("\t")
      let _ = sb.write(v)
    }
    let _ = sb.write("\n")
  }
  // Col attributes (obs)
  let obs_keys = adata.obs.keys().collect()
  let _ = sb.write("col_attrs\t")
  let _ = sb.write(obs_keys.length().to_string())
  let _ = sb.write("\n")
  for k in obs_keys {
    let _ = sb.write(k)
    let vals = adata.obs[k]
    for v in vals {
      let _ = sb.write("\t")
      let _ = sb.write(v)
    }
    let _ = sb.write("\n")
  }
  // Matrix (transposed: var × obs in loom)
  let _ = sb.write("matrix\n")
  let x = adata.x_matrix()
  for j = 0; j < adata.n_vars(); j = j + 1 {
    for i = 0; i < adata.n_obs(); i = i + 1 {
      if i > 0 {
        let _ = sb.write("\t")
      }
      let _ = sb.write(x[i][j].to_string())
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
/// Read Loom-format text into AnnData.
/// Loom: rows = genes (var), columns = cells (obs). We transpose to AnnData convention.
pub fn anndata_read_loom(content : String) -> AnnData? {
  let lines = content.split("\n").to_array()
  if lines.length() < 1 {
    return None
  }
  let header = lines[0].split("\t").to_array()
  if header.length() < 3 || header[0].to_owned() != "loom" {
    return None
  }
  let n_obs = parse_int(header[1].to_owned())
  let n_vars = parse_int(header[2].to_owned())
  let mut idx = 1
  // Parse row_attrs (var attributes)
  let var_attrs : Map[String, Array[String]] = Map([], capacity=4)
  if idx < lines.length() {
    let parts = lines[idx].split("\t").to_array()
    if parts.length() >= 2 && parts[0].to_owned() == "row_attrs" {
      let n_row_attrs = parse_int(parts[1].to_owned())
      idx = idx + 1
      for _r = 0; _r < n_row_attrs; _r = _r + 1 {
        if idx >= lines.length() {
          break
        }
        let fields = lines[idx].split("\t").to_array()
        if fields.length() >= 1 {
          let key = fields[0].to_owned()
          let vals : Array[String] = Array(capacity=fields.length() - 1)
          for f = 1; f < fields.length(); f = f + 1 {
            vals.push(fields[f].to_owned())
          }
          var_attrs[key] = vals
        }
        idx = idx + 1
      }
    }
  }
  // Parse col_attrs (obs attributes)
  let obs_attrs : Map[String, Array[String]] = Map([], capacity=4)
  if idx < lines.length() {
    let parts = lines[idx].split("\t").to_array()
    if parts.length() >= 2 && parts[0].to_owned() == "col_attrs" {
      let n_col_attrs = parse_int(parts[1].to_owned())
      idx = idx + 1
      for _c = 0; _c < n_col_attrs; _c = _c + 1 {
        if idx >= lines.length() {
          break
        }
        let fields = lines[idx].split("\t").to_array()
        if fields.length() >= 1 {
          let key = fields[0].to_owned()
          let vals : Array[String] = Array(capacity=fields.length() - 1)
          for f = 1; f < fields.length(); f = f + 1 {
            vals.push(fields[f].to_owned())
          }
          obs_attrs[key] = vals
        }
        idx = idx + 1
      }
    }
  }
  // Skip "matrix" line
  if idx < lines.length() && lines[idx].to_owned() == "matrix" {
    idx = idx + 1
  }
  // Parse matrix (var × obs in loom → transpose to obs × var)
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for _i = 0; _i < n_obs; _i = _i + 1 {
    let row : Array[Double] = Array(capacity=n_vars)
    for _j = 0; _j < n_vars; _j = _j + 1 {
      row.push(0.0)
    }
    x.push(row)
  }
  for j = 0; j < n_vars && idx + j < lines.length(); j = j + 1 {
    let fields = lines[idx + j].split("\t").to_array()
    for i = 0; i < n_obs && i < fields.length(); i = i + 1 {
      x[i][j] = parse_double(fields[i].to_owned()).unwrap_or(0.0)
    }
  }
  // Build obs_names and var_names
  let obs_names : Array[String] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    obs_names.push("cell_" + i.to_string())
  }
  let var_names : Array[String] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    var_names.push("gene_" + j.to_string())
  }
  let adata = AnnData::new(x, obs_names, var_names)
  // Copy attributes
  for k in obs_attrs.keys().collect() {
    adata.obs[k] = obs_attrs[k]
  }
  for k in var_attrs.keys().collect() {
    adata.variables[k] = var_attrs[k]
  }
  Some(adata)
}

///|
pub fn anndata_io_loom_dummy() -> Unit {
  ignore(LoomFile::{
    matrix: Array(capacity=0),
    row_attrs: Map([], capacity=0),
    col_attrs: Map([], capacity=0),
    row_names: Array(capacity=0),
    col_names: Array(capacity=0),
  })
}
