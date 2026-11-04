/// AnnData transpose — `adata.T` equivalent.
///
/// Swaps obs ↔ var, including all annotations:
/// obs ↔ var, obsm ↔ varm, obsp ↔ varp, obs_names ↔ var_names.
/// X matrix is transposed.

///|
/// Transpose an AnnData: swap observations and variables.
pub fn anndata_transpose(adata : AnnData) -> AnnData {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  // Transpose X
  let x_t : Array[Array[Double]] = Array::make(n_vars, Array::make(n_obs, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      x_t[j][i] = adata.get_x(i, j)
    }
  }
  let result = AnnData::new(
    x_t,
    adata.get_var_names().copy(),
    adata.get_obs_names().copy(),
  )
  // Swap obs ↔ var (string annotations)
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs_column(k) {
      Some(v) => ignore(result.set_var(k, v))
      None => ()
    }
  }
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.var_column(k) {
      Some(v) => ignore(result.set_obs(k, v))
      None => ()
    }
  }
  // Swap obsm ↔ varm
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    match adata.get_obsm(k) {
      Some(m) => ignore(result.set_varm(k, m))
      None => ()
    }
  }
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    match adata.get_varm(k) {
      Some(m) => ignore(result.set_obsm(k, m))
      None => ()
    }
  }
  // Swap obsp ↔ varp
  let obsp_keys = adata.obsp.keys().collect()
  for k in obsp_keys {
    match adata.obsp.get(k) {
      Some(m) => ignore(result.set_varp(k, m))
      None => ()
    }
  }
  let varp_keys = adata.varp.keys().collect()
  for k in varp_keys {
    match adata.varp.get(k) {
      Some(m) => ignore(result.set_obsp(k, m))
      None => ()
    }
  }
  // Transpose layers
  let layer_keys = adata.layers.keys().collect()
  for k in layer_keys {
    match adata.get_layer(k) {
      Some(mat) => {
        let lt : Array[Array[Double]] = Array::make(
          n_vars,
          Array::make(n_obs, 0.0),
        )
        for i = 0; i < mat.length(); i = i + 1 {
          for j = 0; j < mat[i].length(); j = j + 1 {
            lt[j][i] = mat[i][j]
          }
        }
        ignore(result.set_layer(k, lt))
      }
      None => ()
    }
  }
  // Carry over uns
  let uns_keys = adata.uns.keys().collect()
  for k in uns_keys {
    match adata.get_uns(k) {
      Some(v) => ignore(result.set_uns(k, v))
      None => ()
    }
  }
  result
}
