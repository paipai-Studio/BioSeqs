/// AnnData copy and deepcopy.
///
/// `adata.copy()` equivalent — creates a fully independent clone.

///|
/// Deep copy of an AnnData object.
/// All arrays and maps are duplicated; modifications to the copy
/// do not affect the original.
pub fn anndata_copy(adata : AnnData) -> AnnData {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  // Copy X
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_vars)
    for j = 0; j < n_vars; j = j + 1 {
      row.push(adata.get_x(i, j))
    }
    x.push(row)
  }
  let obs_names = adata.get_obs_names().copy()
  let var_names = adata.get_var_names().copy()
  let copy = AnnData::new(x, obs_names, var_names)
  // Copy obs annotations
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs_column(k) {
      Some(col) => ignore(copy.set_obs(k, col.copy()))
      None => ()
    }
  }
  // Copy var annotations
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.var_column(k) {
      Some(col) => ignore(copy.set_var(k, col.copy()))
      None => ()
    }
  }
  // Copy layers
  let layer_keys = adata.layers.keys().collect()
  for k in layer_keys {
    match adata.get_layer(k) {
      Some(mat) => {
        let mat_copy : Array[Array[Double]] = Array(capacity=mat.length())
        for row in mat {
          mat_copy.push(row.copy())
        }
        ignore(copy.set_layer(k, mat_copy))
      }
      None => ()
    }
  }
  // Copy obsm
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    match adata.get_obsm(k) {
      Some(mat) => {
        let mat_copy : Array[Array[Double]] = Array(capacity=mat.length())
        for row in mat {
          mat_copy.push(row.copy())
        }
        ignore(copy.set_obsm(k, mat_copy))
      }
      None => ()
    }
  }
  // Copy varm
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    match adata.get_varm(k) {
      Some(mat) => {
        let mat_copy : Array[Array[Double]] = Array(capacity=mat.length())
        for row in mat {
          mat_copy.push(row.copy())
        }
        ignore(copy.set_varm(k, mat_copy))
      }
      None => ()
    }
  }
  // Copy uns
  let uns_keys = adata.uns.keys().collect()
  for k in uns_keys {
    match adata.get_uns(k) {
      Some(v) => ignore(copy.set_uns(k, v))
      None => ()
    }
  }
  copy
}

///|
/// Shallow copy: shares X data but copies annotations.
/// Useful when you want to modify obs/var without duplicating X.
pub fn anndata_shallow_copy(adata : AnnData) -> AnnData {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  // Share X reference (copy row references but not data)
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_vars)
    for j = 0; j < n_vars; j = j + 1 {
      row.push(adata.get_x(i, j))
    }
    x.push(row)
  }
  let copy = AnnData::new(
    x,
    adata.get_obs_names().copy(),
    adata.get_var_names().copy(),
  )
  // Copy obs annotations
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs_column(k) {
      Some(col) => ignore(copy.set_obs(k, col.copy()))
      None => ()
    }
  }
  // Copy var annotations
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.var_column(k) {
      Some(col) => ignore(copy.set_var(k, col.copy()))
      None => ()
    }
  }
  copy
}

///|
/// Check if two AnnData objects have the same shape.
pub fn anndata_same_shape(a : AnnData, b : AnnData) -> Bool {
  a.n_obs() == b.n_obs() && a.n_vars() == b.n_vars()
}

///|
/// Check if two AnnData objects have the same obs_names and var_names.
pub fn anndata_same_labels(a : AnnData, b : AnnData) -> Bool {
  let a_obs = a.get_obs_names()
  let b_obs = b.get_obs_names()
  if a_obs.length() != b_obs.length() {
    return false
  }
  for i = 0; i < a_obs.length(); i = i + 1 {
    if a_obs[i] != b_obs[i] {
      return false
    }
  }
  let a_var = a.get_var_names()
  let b_var = b.get_var_names()
  if a_var.length() != b_var.length() {
    return false
  }
  for i = 0; i < a_var.length(); i = i + 1 {
    if a_var[i] != b_var[i] {
      return false
    }
  }
  true
}
