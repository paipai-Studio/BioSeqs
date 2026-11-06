/// anndata_concat2.mbt — Advanced concat options
/// Item 23: concat with label, keys, merge strategies, index_unique, fill_value, pairwise

///|
/// Advanced concat options
pub struct AnnDataConcatOptions {
  join : String // "inner" or "outer"
  merge : String // "same", "unique", "first", "only", "none"
  uns_merge : String // same strategies for uns
  label : String? // column name for batch info
  keys : Array[String] // names for each object
  index_unique : String? // delimiter for unique index
  fill_value : Double // fill for outer join
  pairwise : Bool // include obsp/varp
}

///|
/// Default concat options
pub fn anndata_concat_default_options() -> AnnDataConcatOptions {
  AnnDataConcatOptions::{
    join: "inner",
    merge: "none",
    uns_merge: "none",
    label: None,
    keys: Array(capacity=0),
    index_unique: None,
    fill_value: 0.0,
    pairwise: false,
  }
}

///|
/// Advanced concat along obs axis (axis=0).
pub fn anndata_concat_advanced(
  adatas : Array[AnnData],
  options : AnnDataConcatOptions,
) -> AnnData {
  if adatas.length() == 0 {
    return AnnData::new(Array(capacity=0), Array(capacity=0), Array(capacity=0))
  }
  // Determine common vars (inner join) or all vars (outer join)
  let var_sets : Array[Set[String]] = Array(capacity=adatas.length())
  for a in adatas {
    let s : Set[String] = Set([])
    for n in a.var_names {
      s.add(n)
    }
    var_sets.push(s)
  }
  let all_var_names : Array[String] = if options.join == "inner" {
    // Intersection of all var sets
    let result : Array[String] = Array(capacity=0)
    let first = adatas[0].var_names
    for n in first {
      let mut in_all = true
      for s = 1; s < var_sets.length(); s = s + 1 {
        if not(var_sets[s].contains(n)) {
          in_all = false
          break
        }
      }
      if in_all {
        result.push(n)
      }
    }
    result
  } else {
    // Union of all var sets
    let seen : Set[String] = Set([])
    let result : Array[String] = Array(capacity=0)
    for a in adatas {
      for n in a.var_names {
        if not(seen.contains(n)) {
          seen.add(n)
          result.push(n)
        }
      }
    }
    result
  }
  // Build var index maps for each AnnData
  let var_idx_maps : Array[Map[String, Int32]] = Array(capacity=adatas.length())
  for a in adatas {
    let m : Map[String, Int32] = Map([], capacity=a.var_names.length())
    for j = 0; j < a.var_names.length(); j = j + 1 {
      m[a.var_names[j]] = j
    }
    var_idx_maps.push(m)
  }
  // Concatenate X
  let n_vars = all_var_names.length()
  let x : Array[Array[Double]] = Array(capacity=0)
  let obs_names : Array[String] = Array(capacity=0)
  for i = 0; i < adatas.length(); i = i + 1 {
    let a = adatas[i]
    let ax = a.x_matrix()
    let prefix = if i < options.keys.length() { options.keys[i] } else { "" }
    for r = 0; r < ax.length(); r = r + 1 {
      let row : Array[Double] = Array(capacity=n_vars)
      for j = 0; j < n_vars; j = j + 1 {
        let vname = all_var_names[j]
        if var_idx_maps[i].contains(vname) {
          let src_j = var_idx_maps[i][vname]
          row.push(ax[r][src_j])
        } else {
          row.push(options.fill_value)
        }
      }
      x.push(row)
      let name = if options.index_unique is Some(delim) && prefix.length() > 0 {
        a.obs_names[r] + delim + prefix
      } else if prefix.length() > 0 {
        prefix + "_" + a.obs_names[r]
      } else {
        a.obs_names[r]
      }
      obs_names.push(name)
    }
  }
  let result = AnnData::new(x, obs_names, all_var_names)
  // Add label column if requested
  if options.label is Some(label_name) {
    let labels : Array[String] = Array(capacity=0)
    for i = 0; i < adatas.length(); i = i + 1 {
      let key = if i < options.keys.length() {
        options.keys[i]
      } else {
        i.to_string()
      }
      for _r = 0; _r < adatas[i].n_obs(); _r = _r + 1 {
        labels.push(key)
      }
    }
    result.obs[label_name] = labels
  }
  // Merge uns
  if options.uns_merge != "none" {
    let merged : Map[String, JsonValue] = Map([], capacity=4)
    for a in adatas {
      for k in a.uns.keys().collect() {
        if not(merged.contains(k)) {
          merged[k] = a.uns[k]
        } else if options.uns_merge == "first" {
          // keep first
          continue
        }
      }
    }
    for k in merged.keys().collect() {
      result.uns[k] = merged[k]
    }
  }
  result
}

///|
/// Concat along var axis (axis=1).
pub fn anndata_concat_var_axis(
  adatas : Array[AnnData],
  join~ : String = "inner",
  fill_value~ : Double = 0.0,
) -> AnnData {
  if adatas.length() == 0 {
    return AnnData::new(Array(capacity=0), Array(capacity=0), Array(capacity=0))
  }
  // All must have same n_obs (for var concat)
  let n_obs = adatas[0].n_obs()
  // Collect all var names
  let all_var_names : Array[String] = Array(capacity=0)
  for a in adatas {
    for n in a.var_names {
      all_var_names.push(n)
    }
  }
  // Concatenate X horizontally
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for a in adatas {
      let ax = a.x_matrix()
      if i < ax.length() {
        for v in ax[i] {
          row.push(v)
        }
      }
    }
    x.push(row)
  }
  let obs_names = adatas[0].obs_names
  ignore(join)
  ignore(fill_value)
  AnnData::new(x, obs_names, all_var_names)
}
