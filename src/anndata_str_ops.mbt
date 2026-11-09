/// anndata_str_ops.mbt — String operations on obs/var
/// Items 18-20: obs_str/var_str, obs_to_str/var_to_str, obs_filter/var_filter

///|
/// Get a string column from obs.
pub fn anndata_obs_str(adata : AnnData, key : String) -> Array[String]? {
  if adata.obs.contains(key) {
    Some(adata.obs[key])
  } else {
    None
  }
}

///|
/// Get a string column from var.
pub fn anndata_var_str(adata : AnnData, key : String) -> Array[String]? {
  if adata.variables.contains(key) {
    Some(adata.variables[key])
  } else {
    None
  }
}

///|
/// Convert a numeric obs column to string representation.
pub fn anndata_obs_to_str(adata : AnnData, key : String) -> Array[String]? {
  // Check if it's already a string column
  if adata.obs.contains(key) {
    Some(adata.obs[key])
  } else {
    // Try to get from numeric annotation
    let num_key = key + "_num"
    if adata.obs.contains(num_key) {
      let vals = adata.obs[num_key]
      let result : Array[String] = Array(capacity=vals.length())
      for v in vals {
        result.push(v)
      }
      Some(result)
    } else {
      None
    }
  }
}

///|
/// Convert a numeric var column to string representation.
pub fn anndata_var_to_str(adata : AnnData, key : String) -> Array[String]? {
  if adata.variables.contains(key) {
    Some(adata.variables[key])
  } else {
    let num_key = key + "_num"
    if adata.variables.contains(num_key) {
      let vals = adata.variables[num_key]
      let result : Array[String] = Array(capacity=vals.length())
      for v in vals {
        result.push(v)
      }
      Some(result)
    } else {
      None
    }
  }
}

///|
/// Filter observations by a string condition on a column.
/// Returns indices of rows where col matches the value.
pub fn anndata_obs_filter(
  adata : AnnData,
  key : String,
  value : String,
) -> Array[Int32] {
  if not(adata.obs.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.obs[key]
  let result : Array[Int32] = Array(capacity=vals.length())
  for i = 0; i < vals.length(); i = i + 1 {
    if vals[i] == value {
      result.push(i)
    }
  }
  result
}

///|
/// Filter variables by a string condition on a column.
pub fn anndata_var_filter(
  adata : AnnData,
  key : String,
  value : String,
) -> Array[Int32] {
  if not(adata.variables.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.variables[key]
  let result : Array[Int32] = Array(capacity=vals.length())
  for i = 0; i < vals.length(); i = i + 1 {
    if vals[i] == value {
      result.push(i)
    }
  }
  result
}

///|
/// Filter obs by a predicate function (encoded as fn_id).
/// fn_id 0: equals, 1: contains, 2: starts_with, 3: ends_with
pub fn anndata_obs_filter_pred(
  adata : AnnData,
  key : String,
  value : String,
  fn_id~ : Int32 = 0,
) -> Array[Int32] {
  if not(adata.obs.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.obs[key]
  let result : Array[Int32] = Array(capacity=vals.length())
  for i = 0; i < vals.length(); i = i + 1 {
    let matched = match fn_id {
      0 => vals[i] == value
      1 => vals[i].contains(value)
      2 => vals[i].has_prefix(value)
      3 => vals[i].has_suffix(value)
      _ => false
    }
    if matched {
      result.push(i)
    }
  }
  result
}

///|
/// Filter var by a predicate function.
pub fn anndata_var_filter_pred(
  adata : AnnData,
  key : String,
  value : String,
  fn_id~ : Int32 = 0,
) -> Array[Int32] {
  if not(adata.variables.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.variables[key]
  let result : Array[Int32] = Array(capacity=vals.length())
  for i = 0; i < vals.length(); i = i + 1 {
    let matched = match fn_id {
      0 => vals[i] == value
      1 => vals[i].contains(value)
      2 => vals[i].has_prefix(value)
      3 => vals[i].has_suffix(value)
      _ => false
    }
    if matched {
      result.push(i)
    }
  }
  result
}

///|
/// Get unique values from an obs column.
pub fn anndata_obs_unique(adata : AnnData, key : String) -> Array[String] {
  if not(adata.obs.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.obs[key]
  let seen : Set[String] = Set([])
  let result : Array[String] = Array(capacity=vals.length())
  for v in vals {
    if not(seen.contains(v)) {
      seen.add(v)
      result.push(v)
    }
  }
  result
}

///|
/// Get unique values from a var column.
pub fn anndata_var_unique(adata : AnnData, key : String) -> Array[String] {
  if not(adata.variables.contains(key)) {
    return Array(capacity=0)
  }
  let vals = adata.variables[key]
  let seen : Set[String] = Set([])
  let result : Array[String] = Array(capacity=vals.length())
  for v in vals {
    if not(seen.contains(v)) {
      seen.add(v)
      result.push(v)
    }
  }
  result
}
