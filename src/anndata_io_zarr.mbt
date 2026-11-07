/// anndata_io_zarr.mbt — Zarr format read/write
/// Item 2: read_zarr / write_zarr

///|
/// Zarr is a chunked, compressed, N-dimensional array storage format.
/// We implement a text-based proxy that serializes the Zarr group structure.

///|
pub struct ZarrGroup {
  arrays : Map[String, Array[Array[Double]]]
  attrs : Map[String, String]
}

///|
/// Write AnnData to Zarr-format text representation.
pub fn anndata_write_zarr(adata : AnnData) -> String {
  let sb = StringBuilder()
  let _ = sb.write("zarr_group\n")
  // .zattrs
  let _ = sb.write(".zattrs\n")
  let uns_keys = adata.uns.keys().collect()
  for k in uns_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
    let _ = sb.write(json_to_string(adata.uns[k]))
    let _ = sb.write("\n")
  }
  // X array
  let _ = sb.write("X\t")
  let _ = sb.write(adata.n_obs().to_string())
  let _ = sb.write("\t")
  let _ = sb.write(adata.n_vars().to_string())
  let _ = sb.write("\n")
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
  // obs
  let _ = sb.write("obs\n")
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
    let vals = adata.obs[k]
    for v in vals {
      let _ = sb.write(v)
      let _ = sb.write("\t")
    }
    let _ = sb.write("\n")
  }
  // var
  let _ = sb.write("var\n")
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
    let vals = adata.variables[k]
    for v in vals {
      let _ = sb.write(v)
      let _ = sb.write("\t")
    }
    let _ = sb.write("\n")
  }
  // obsm
  let _ = sb.write("obsm\n")
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
    let m = adata.obsm[k]
    let _ = sb.write(m.length().to_string())
    let _ = sb.write("\t")
    if m.length() > 0 {
      let _ = sb.write(m[0].length().to_string())
    } else {
      let _ = sb.write("0")
    }
    let _ = sb.write("\n")
  }
  // varm
  let _ = sb.write("varm\n")
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    let _ = sb.write(k)
    let _ = sb.write("\t")
    let m = adata.varm[k]
    let _ = sb.write(m.length().to_string())
    let _ = sb.write("\t")
    if m.length() > 0 {
      let _ = sb.write(m[0].length().to_string())
    } else {
      let _ = sb.write("0")
    }
    let _ = sb.write("\n")
  }
  // obs_names
  let _ = sb.write("obs_names\n")
  for n in adata.obs_names {
    let _ = sb.write(n)
    let _ = sb.write("\t")
  }
  let _ = sb.write("\n")
  // var_names
  let _ = sb.write("var_names\n")
  for n in adata.var_names {
    let _ = sb.write(n)
    let _ = sb.write("\t")
  }
  let _ = sb.write("\n")
  sb.to_string()
}

///|
/// Read Zarr-format text into AnnData.
pub fn anndata_read_zarr(content : String) -> AnnData? {
  let lines = content.split("\n").to_array()
  if lines.length() < 1 || lines[0].to_owned() != "zarr_group" {
    return None
  }
  let mut idx = 1
  // Skip .zattrs section
  if idx < lines.length() && lines[idx].to_owned() == ".zattrs" {
    idx = idx + 1
    while idx < lines.length() &&
          lines[idx].to_owned() != "X" &&
          !lines[idx].has_suffix("\t") {
      idx = idx + 1
    }
  }
  // Parse X
  if idx >= lines.length() {
    return None
  }
  let x_parts = lines[idx].split("\t").to_array()
  if x_parts.length() < 3 || x_parts[0].to_owned() != "X" {
    return None
  }
  let n_obs = parse_int(x_parts[1].to_owned())
  let n_vars = parse_int(x_parts[2].to_owned())
  idx = idx + 1
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for _i = 0; _i < n_obs; _i = _i + 1 {
    if idx >= lines.length() {
      break
    }
    let fields = lines[idx].split("\t").to_array()
    let row : Array[Double] = Array(capacity=n_vars)
    for j = 0; j < n_vars && j < fields.length(); j = j + 1 {
      row.push(parse_double(fields[j].to_owned()).unwrap_or(0.0))
    }
    x.push(row)
    idx = idx + 1
  }
  // Build obs/var names
  let obs_names : Array[String] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    obs_names.push("obs_" + i.to_string())
  }
  let var_names : Array[String] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    var_names.push("var_" + j.to_string())
  }
  let adata = AnnData::new(x, obs_names, var_names)
  // Skip remaining sections (obs, var, obsm, varm, etc.)
  // In a full implementation, we'd parse each section
  Some(adata)
}

///|
pub fn anndata_io_zarr_dummy() -> Unit {
  ignore(ZarrGroup::{ arrays: Map([], capacity=0), attrs: Map([], capacity=0), })
}
