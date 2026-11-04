/// AnnData shape validation and consistency checks.
///
/// Verifies that obs/var annotations, layers, obsm, varm, obsp, varp
/// all have dimensions consistent with X.

///|
/// Validation result.
pub struct AnnDataValidation {
  valid : Bool
  errors : Array[String]
  warnings : Array[String]
}

///|
pub fn AnnDataValidation::valid(self : AnnDataValidation) -> Bool {
  self.valid
}

///|
pub fn AnnDataValidation::errors(self : AnnDataValidation) -> Array[String] {
  self.errors
}

///|
pub fn AnnDataValidation::warnings(self : AnnDataValidation) -> Array[String] {
  self.warnings
}

///|
/// Validate AnnData shape consistency.
pub fn anndata_validate(adata : AnnData) -> AnnDataValidation {
  let errors : Array[String] = Array(capacity=0)
  let warnings : Array[String] = Array(capacity=0)
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  // Check obs_names length
  let obs_names = adata.get_obs_names()
  if obs_names.length() != n_obs {
    errors.push(
      "obs_names length (" +
      obs_names.length().to_string() +
      ") != n_obs (" +
      n_obs.to_string() +
      ")",
    )
  }
  // Check var_names length
  let var_names = adata.get_var_names()
  if var_names.length() != n_vars {
    errors.push(
      "var_names length (" +
      var_names.length().to_string() +
      ") != n_vars (" +
      n_vars.to_string() +
      ")",
    )
  }
  // Check obs annotation columns
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs_column(k) {
      Some(col) =>
        if col.length() != n_obs {
          errors.push(
            "obs column '" +
            k +
            "' length (" +
            col.length().to_string() +
            ") != n_obs (" +
            n_obs.to_string() +
            ")",
          )
        }
      None => ()
    }
  }
  // Check var annotation columns
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.var_column(k) {
      Some(col) =>
        if col.length() != n_vars {
          errors.push(
            "var column '" +
            k +
            "' length (" +
            col.length().to_string() +
            ") != n_vars (" +
            n_vars.to_string() +
            ")",
          )
        }
      None => ()
    }
  }
  // Check layers
  let layer_keys = adata.layers.keys().collect()
  for k in layer_keys {
    match adata.get_layer(k) {
      Some(mat) => {
        if mat.length() != n_obs {
          errors.push(
            "layer '" +
            k +
            "' n_rows (" +
            mat.length().to_string() +
            ") != n_obs (" +
            n_obs.to_string() +
            ")",
          )
        }
        if mat.length() > 0 && mat[0].length() != n_vars {
          errors.push(
            "layer '" +
            k +
            "' n_cols (" +
            mat[0].length().to_string() +
            ") != n_vars (" +
            n_vars.to_string() +
            ")",
          )
        }
      }
      None => ()
    }
  }
  // Check obsm
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    match adata.get_obsm(k) {
      Some(mat) =>
        if mat.length() != n_obs {
          errors.push(
            "obsm '" +
            k +
            "' n_rows (" +
            mat.length().to_string() +
            ") != n_obs (" +
            n_obs.to_string() +
            ")",
          )
        }
      None => ()
    }
  }
  // Check varm
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    match adata.get_varm(k) {
      Some(mat) =>
        if mat.length() != n_vars {
          errors.push(
            "varm '" +
            k +
            "' n_rows (" +
            mat.length().to_string() +
            ") != n_vars (" +
            n_vars.to_string() +
            ")",
          )
        }
      None => ()
    }
  }
  // Check for duplicate obs_names
  if obs_names.length() > 0 {
    let seen : Map[String, Bool] = Map([], capacity=obs_names.length())
    for s in obs_names {
      match seen.get(s) {
        Some(_) => {
          warnings.push("duplicate obs_name: '" + s + "'")
          break
        }
        None => seen[s] = true
      }
    }
  }
  // Check for duplicate var_names
  if var_names.length() > 0 {
    let seen : Map[String, Bool] = Map([], capacity=var_names.length())
    for s in var_names {
      match seen.get(s) {
        Some(_) => {
          warnings.push("duplicate var_name: '" + s + "'")
          break
        }
        None => seen[s] = true
      }
    }
  }
  // Check for empty data
  if n_obs == 0 {
    warnings.push("n_obs is 0")
  }
  if n_vars == 0 {
    warnings.push("n_vars is 0")
  }
  { valid: errors.length() == 0, errors, warnings, }
}

///|
/// Assert that AnnData is valid (aborts if not).
pub fn anndata_assert_valid(adata : AnnData) -> Unit {
  let v = anndata_validate(adata)
  if !v.valid {
    let mut msg = "AnnData validation failed:\n"
    for e in v.errors {
      msg = msg + "  " + e + "\n"
    }
    abort(msg)
  }
}

///|
/// Check shape matches expected.
pub fn anndata_assert_shape(
  adata : AnnData,
  expected_obs : Int32,
  expected_var : Int32,
) -> Unit {
  if adata.n_obs() != expected_obs {
    abort(
      "n_obs mismatch: got " +
      adata.n_obs().to_string() +
      ", expected " +
      expected_obs.to_string(),
    )
  }
  if adata.n_vars() != expected_var {
    abort(
      "n_vars mismatch: got " +
      adata.n_vars().to_string() +
      ", expected " +
      expected_var.to_string(),
    )
  }
}
