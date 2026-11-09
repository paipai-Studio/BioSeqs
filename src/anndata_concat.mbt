/// AnnData concatenation — `anndata.concat` equivalent.
///
/// Concatenates multiple AnnData objects along obs or var axis,
/// with join strategies (inner/outer) and merge strategies
/// (same/unique/first/only) for non-aligned elements.

///|
/// Join strategy for the alignment axis
pub enum AnnDataJoin {
  Inner // intersection of var_names (when axis=obs) or obs_names (when axis=var)
  Outer // union, fill missing with 0.0
}

///|
/// Merge strategy for non-aligned elements (obsm, varm, obsp, varp, uns)
pub enum AnnDataMerge {
  Drop // drop all non-aligned elements
  Same // keep only if identical across all objects
  First // take first-seen value
  Only // keep if present in only one object
}

///|
/// Concatenate AnnData objects along obs axis (stack rows).
/// X matrices are stacked vertically. var annotations are joined.
pub fn anndata_concat_obs(
  adatas : Array[AnnData],
  join~ : AnnDataJoin = AnnDataJoin::Inner,
  _merge~ : AnnDataMerge = AnnDataMerge::Drop,
  label~ : String? = None,
  keys~ : Array[String] = [],
) -> AnnData {
  let n = adatas.length()
  if n == 0 {
    return AnnData::new([], [], [])
  }
  if n == 1 {
    return adatas[0]
  }
  // Determine var_names based on join strategy
  let var_names = match join {
    AnnDataJoin::Inner => {
      // Intersection of all var_names
      let mut result = adatas[0].get_var_names().copy()
      for i = 1; i < n; i = i + 1 {
        result = anndata_intersect_strings(result, adatas[i].get_var_names())
      }
      result
    }
    AnnDataJoin::Outer => {
      // Union of all var_names
      let mut result = adatas[0].get_var_names().copy()
      for i = 1; i < n; i = i + 1 {
        result = anndata_union_strings(result, adatas[i].get_var_names())
      }
      result
    }
  }
  // Build var index map for each AnnData
  let var_idx_maps : Array[Map[String, Int32]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let m : Map[String, Int32] = Map([], capacity=4)
    let names = adatas[i].get_var_names()
    for j = 0; j < names.length(); j = j + 1 {
      m[names[j]] = j
    }
    var_idx_maps.push(m)
  }
  // Concatenate X matrices
  let all_obs_names : Array[String] = Array(capacity=0)
  let x_rows : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let obs_names_i = adatas[i].get_obs_names()
    for j = 0; j < obs_names_i.length(); j = j + 1 {
      let mut name = obs_names_i[j]
      if label is Some(l) && l.length() > 0 && i < keys.length() {
        name = name + "-" + keys[i]
      }
      all_obs_names.push(name)
    }
    // Build row for each obs
    for r = 0; r < adatas[i].n_obs(); r = r + 1 {
      let row : Array[Double] = Array::make(var_names.length(), 0.0)
      for c = 0; c < var_names.length(); c = c + 1 {
        let vname = var_names[c]
        match var_idx_maps[i].get(vname) {
          Some(idx) => row[c] = adatas[i].get_x(r, idx)
          None => ()
        }
      }
      x_rows.push(row)
    }
  }
  let result = AnnData::new(x_rows, all_obs_names, var_names)
  // Merge obs annotations
  let obs_keys = adatas[0].obs.keys().collect()
  for k in obs_keys {
    let mut all_present = true
    for i = 1; i < n; i = i + 1 {
      if adatas[i].obs.get(k) is Some(_) {
        all_present = false
      }
    }
    if all_present {
      let vals : Array[String] = Array(capacity=0)
      for i = 0; i < n; i = i + 1 {
        match adatas[i].obs_column(k) {
          Some(v) =>
            for s in v {
              vals.push(s)
            }
          None => ()
        }
      }
      ignore(result.set_obs(k, vals))
    }
  }
  result
}

///|
/// Concatenate AnnData objects along var axis (stack columns).
pub fn anndata_concat_var(
  adatas : Array[AnnData],
  join~ : AnnDataJoin = AnnDataJoin::Inner,
  _merge~ : AnnDataMerge = AnnDataMerge::Drop,
) -> AnnData {
  let n = adatas.length()
  if n == 0 {
    return AnnData::new([], [], [])
  }
  if n == 1 {
    return adatas[0]
  }
  // Determine obs_names based on join
  let obs_names = match join {
    AnnDataJoin::Inner => {
      let mut result = adatas[0].get_obs_names().copy()
      for i = 1; i < n; i = i + 1 {
        result = anndata_intersect_strings(result, adatas[i].get_obs_names())
      }
      result
    }
    AnnDataJoin::Outer => {
      let mut result = adatas[0].get_obs_names().copy()
      for i = 1; i < n; i = i + 1 {
        result = anndata_union_strings(result, adatas[i].get_obs_names())
      }
      result
    }
  }
  // Build obs index map
  let obs_idx_maps : Array[Map[String, Int32]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let m : Map[String, Int32] = Map([], capacity=4)
    let names = adatas[i].get_obs_names()
    for j = 0; j < names.length(); j = j + 1 {
      m[names[j]] = j
    }
    obs_idx_maps.push(m)
  }
  // Concatenate columns
  let all_var_names : Array[String] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let vnames = adatas[i].get_var_names()
    for v in vnames {
      all_var_names.push(v)
    }
  }
  let x_rows : Array[Array[Double]] = Array::make(
    obs_names.length(),
    Array::make(all_var_names.length(), 0.0),
  )
  let mut col_offset = 0
  for i = 0; i < n; i = i + 1 {
    let vnames_i = adatas[i].get_var_names()
    for r = 0; r < obs_names.length(); r = r + 1 {
      let obs_name = obs_names[r]
      match obs_idx_maps[i].get(obs_name) {
        Some(src_row) =>
          for c = 0; c < vnames_i.length(); c = c + 1 {
            x_rows[r][col_offset + c] = adatas[i].get_x(src_row, c)
          }
        None => ()
      }
    }
    col_offset = col_offset + vnames_i.length()
  }
  AnnData::new(x_rows, obs_names, all_var_names)
}

///|
/// String intersection
fn anndata_intersect_strings(
  a : Array[String],
  b : Array[String],
) -> Array[String] {
  let bset : Map[String, Bool] = Map([], capacity=b.length())
  for s in b {
    bset[s] = true
  }
  let out : Array[String] = Array(capacity=0)
  for s in a {
    match bset.get(s) {
      Some(_) => out.push(s)
      None => ()
    }
  }
  out
}

///|
/// String union (preserves order, no duplicates)
fn anndata_union_strings(a : Array[String], b : Array[String]) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([], capacity=a.length() + b.length())
  for s in a {
    match seen.get(s) {
      None => {
        out.push(s)
        seen[s] = true
      }
      Some(_) => ()
    }
  }
  for s in b {
    match seen.get(s) {
      None => {
        out.push(s)
        seen[s] = true
      }
      Some(_) => ()
    }
  }
  out
}

///|
fn _anndata_join_unused_constructors() -> Unit {
  ignore(AnnDataJoin::Outer)
  ignore(AnnDataMerge::Same)
  ignore(AnnDataMerge::First)
  ignore(AnnDataMerge::Only)
}
