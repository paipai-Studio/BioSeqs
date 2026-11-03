/// AnnData boolean mask + name-based indexing.
///
/// `adata[bool_mask_obs, bool_mask_var]` and
/// `adata[:, ["gene_A", "gene_B"]]`

///|
/// Convert boolean mask to index array.
pub fn anndata_bool_to_indices(mask : Array[Bool]) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for i = 0; i < mask.length(); i = i + 1 {
    if mask[i] {
      out.push(i)
    }
  }
  out
}

///|
/// Subset AnnData by boolean masks for obs and var.
pub fn anndata_subset_bool(
  adata : AnnData,
  obs_mask : Array[Bool],
  var_mask : Array[Bool],
) -> AnnData {
  let obs_idx = anndata_bool_to_indices(obs_mask)
  let var_idx = anndata_bool_to_indices(var_mask)
  adata.subset(obs_idx, var_idx)
}

///|
/// Subset AnnData by boolean mask on obs only (all vars).
pub fn anndata_subset_obs_bool(
  adata : AnnData,
  obs_mask : Array[Bool],
) -> AnnData {
  let obs_idx = anndata_bool_to_indices(obs_mask)
  let all_vars : Array[Int32] = Array(capacity=adata.n_vars())
  for i = 0; i < adata.n_vars(); i = i + 1 {
    all_vars.push(i)
  }
  adata.subset(obs_idx, all_vars)
}

///|
/// Subset AnnData by boolean mask on var only (all obs).
pub fn anndata_subset_var_bool(
  adata : AnnData,
  var_mask : Array[Bool],
) -> AnnData {
  let all_obs : Array[Int32] = Array(capacity=adata.n_obs())
  for i = 0; i < adata.n_obs(); i = i + 1 {
    all_obs.push(i)
  }
  let var_idx = anndata_bool_to_indices(var_mask)
  adata.subset(all_obs, var_idx)
}

///|
/// Convert var_names (strings) to var indices.
fn anndata_names_to_indices(
  names : Array[String],
  target : Array[String],
) -> Array[Int32] {
  let name_map : Map[String, Int32] = Map([], capacity=target.length())
  for i = 0; i < target.length(); i = i + 1 {
    name_map[target[i]] = i
  }
  let out : Array[Int32] = Array(capacity=names.length())
  for s in names {
    match name_map.get(s) {
      Some(idx) => out.push(idx)
      None => out.push(-1) // not found
    }
  }
  out
}

///|
/// Subset AnnData by obs_names and var_names (strings).
/// Only matches found in the AnnData are kept.
pub fn anndata_subset_by_name(
  adata : AnnData,
  obs_names : Array[String],
  var_names : Array[String],
) -> AnnData {
  let obs_idx = anndata_names_to_indices(obs_names, adata.get_obs_names())
  let var_idx = anndata_names_to_indices(var_names, adata.get_var_names())
  // Filter out -1 (not found)
  let valid_obs : Array[Int32] = Array(capacity=0)
  for i in obs_idx {
    if i >= 0 {
      valid_obs.push(i)
    }
  }
  let valid_var : Array[Int32] = Array(capacity=0)
  for i in var_idx {
    if i >= 0 {
      valid_var.push(i)
    }
  }
  adata.subset(valid_obs, valid_var)
}

///|
/// Subset by obs_names only (all vars).
pub fn anndata_subset_obs_by_name(
  adata : AnnData,
  obs_names : Array[String],
) -> AnnData {
  let all_vars : Array[Int32] = Array(capacity=adata.n_vars())
  for i = 0; i < adata.n_vars(); i = i + 1 {
    all_vars.push(i)
  }
  let obs_idx = anndata_names_to_indices(obs_names, adata.get_obs_names())
  let valid_obs : Array[Int32] = Array(capacity=0)
  for i in obs_idx {
    if i >= 0 {
      valid_obs.push(i)
    }
  }
  adata.subset(valid_obs, all_vars)
}

///|
/// Subset by var_names only (all obs).
pub fn anndata_subset_var_by_name(
  adata : AnnData,
  var_names : Array[String],
) -> AnnData {
  let all_obs : Array[Int32] = Array(capacity=adata.n_obs())
  for i = 0; i < adata.n_obs(); i = i + 1 {
    all_obs.push(i)
  }
  let var_idx = anndata_names_to_indices(var_names, adata.get_var_names())
  let valid_var : Array[Int32] = Array(capacity=0)
  for i in var_idx {
    if i >= 0 {
      valid_var.push(i)
    }
  }
  adata.subset(all_obs, valid_var)
}

///|
/// Create a boolean mask from an obs column matching a value.
/// e.g. `adata[adata.obs["batch"] == "batch1", :]`
pub fn anndata_obs_mask(
  adata : AnnData,
  key : String,
  value : String,
) -> Array[Bool] {
  let n = adata.n_obs()
  let mask : Array[Bool] = Array::make(n, false)
  match adata.obs_column(key) {
    Some(col) =>
      for i = 0; i < n && i < col.length(); i = i + 1 {
        if col[i] == value {
          mask[i] = true
        }
      }
    None => ()
  }
  mask
}

///|
/// Create a boolean mask from a var column matching a value.
pub fn anndata_var_mask(
  adata : AnnData,
  key : String,
  value : String,
) -> Array[Bool] {
  let n = adata.n_vars()
  let mask : Array[Bool] = Array::make(n, false)
  match adata.var_column(key) {
    Some(col) =>
      for i = 0; i < n && i < col.length(); i = i + 1 {
        if col[i] == value {
          mask[i] = true
        }
      }
    None => ()
  }
  mask
}
