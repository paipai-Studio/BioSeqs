/// anndata_equals.mbt — Deep equality comparison for AnnData
/// Item 12: equals / __eq__

///|
/// Deep equality check: compare X, obs, var, obs_names, var_names.
/// Optionally compare layers, obsm, varm, obsp, varp, uns.
pub fn anndata_equals(
  a : AnnData,
  b : AnnData,
  check_layers~ : Bool = true,
  check_obsm~ : Bool = true,
  check_varm~ : Bool = true,
  check_obsp~ : Bool = true,
  check_varp~ : Bool = true,
  check_uns~ : Bool = true,
  tol~ : Double = 0.0000000001,
) -> Bool {
  // Shape check
  if a.n_obs() != b.n_obs() || a.n_vars() != b.n_vars() {
    return false
  }
  // obs_names
  if a.obs_names.length() != b.obs_names.length() {
    return false
  }
  for i = 0; i < a.obs_names.length(); i = i + 1 {
    if a.obs_names[i] != b.obs_names[i] {
      return false
    }
  }
  // var_names
  if a.var_names.length() != b.var_names.length() {
    return false
  }
  for i = 0; i < a.var_names.length(); i = i + 1 {
    if a.var_names[i] != b.var_names[i] {
      return false
    }
  }
  // X matrix
  let ax = a.x_matrix()
  let bx = b.x_matrix()
  for i = 0; i < ax.length(); i = i + 1 {
    for j = 0; j < ax[i].length(); j = j + 1 {
      if ad_eq_abs(ax[i][j], bx[i][j], tol) == false {
        return false
      }
    }
  }
  // obs
  let a_obs_keys = a.obs.keys().collect()
  let b_obs_keys = b.obs.keys().collect()
  if a_obs_keys.length() != b_obs_keys.length() {
    return false
  }
  for k in a_obs_keys {
    if not(b.obs.contains(k)) {
      return false
    }
    let av = a.obs[k]
    let bv = b.obs[k]
    if av.length() != bv.length() {
      return false
    }
    for i = 0; i < av.length(); i = i + 1 {
      if av[i] != bv[i] {
        return false
      }
    }
  }
  // var
  let a_var_keys = a.variables.keys().collect()
  let b_var_keys = b.variables.keys().collect()
  if a_var_keys.length() != b_var_keys.length() {
    return false
  }
  for k in a_var_keys {
    if not(b.variables.contains(k)) {
      return false
    }
    let av = a.variables[k]
    let bv = b.variables[k]
    if av.length() != bv.length() {
      return false
    }
    for i = 0; i < av.length(); i = i + 1 {
      if av[i] != bv[i] {
        return false
      }
    }
  }
  // Optional checks
  if check_layers {
    if not(ad_eq_maps_2d(a.layers, b.layers, tol)) {
      return false
    }
  }
  if check_obsm {
    if not(ad_eq_maps_2d(a.obsm, b.obsm, tol)) {
      return false
    }
  }
  if check_varm {
    if not(ad_eq_maps_2d(a.varm, b.varm, tol)) {
      return false
    }
  }
  if check_obsp {
    if not(ad_eq_maps_2d(a.obsp, b.obsp, tol)) {
      return false
    }
  }
  if check_varp {
    if not(ad_eq_maps_2d(a.varp, b.varp, tol)) {
      return false
    }
  }
  if check_uns {
    let a_uns_keys = a.uns.keys().collect()
    let b_uns_keys = b.uns.keys().collect()
    if a_uns_keys.length() != b_uns_keys.length() {
      return false
    }
    for k in a_uns_keys {
      if not(b.uns.contains(k)) {
        return false
      }
    }
  }
  true
}

///|
/// Shallow equality: only shape and names.
pub fn anndata_equals_shallow(a : AnnData, b : AnnData) -> Bool {
  a.n_obs() == b.n_obs() && a.n_vars() == b.n_vars()
}

///|
fn ad_eq_abs(x : Double, y : Double, tol : Double) -> Bool {
  let d = x - y
  let abs_d = if d < 0.0 { -d } else { d }
  abs_d < tol
}

///|
fn ad_eq_maps_2d(
  a : Map[String, Array[Array[Double]]],
  b : Map[String, Array[Array[Double]]],
  tol : Double,
) -> Bool {
  let a_keys = a.keys().collect()
  let b_keys = b.keys().collect()
  if a_keys.length() != b_keys.length() {
    return false
  }
  for k in a_keys {
    if not(b.contains(k)) {
      return false
    }
    let am = a[k]
    let bm = b[k]
    if am.length() != bm.length() {
      return false
    }
    for i = 0; i < am.length(); i = i + 1 {
      if am[i].length() != bm[i].length() {
        return false
      }
      for j = 0; j < am[i].length(); j = j + 1 {
        if ad_eq_abs(am[i][j], bm[i][j], tol) == false {
          return false
        }
      }
    }
  }
  true
}
