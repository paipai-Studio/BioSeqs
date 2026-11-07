/// anndata_mapping.mbt — Mapping/iteration interface for AnnData annotations
/// Items 6-10: obsm_keys, varm_keys, obsp_keys, varp_keys, layer_keys,
///   has_layer/has_obsm/has_varm/has_obsp/has_varp,
///   obsm_names/varm_names,
///   del_layer/del_obsm/del_varm/del_obsp/del_varp

///|
// === Key listing ===

///|
/// List all obsm keys.
pub fn anndata_obsm_keys(adata : AnnData) -> Array[String] {
  adata.obsm.keys().collect()
}

///|
/// List all varm keys.
pub fn anndata_varm_keys(adata : AnnData) -> Array[String] {
  adata.varm.keys().collect()
}

///|
/// List all obsp keys.
pub fn anndata_obsp_keys(adata : AnnData) -> Array[String] {
  adata.obsp.keys().collect()
}

///|
/// List all varp keys.
pub fn anndata_varp_keys(adata : AnnData) -> Array[String] {
  adata.varp.keys().collect()
}

///|
/// List all layer keys.
pub fn anndata_layer_keys(adata : AnnData) -> Array[String] {
  adata.layers.keys().collect()
}

///|
// === Has checks ===

///|
pub fn anndata_has_layer(adata : AnnData, key : String) -> Bool {
  adata.layers.contains(key)
}

///|
pub fn anndata_has_obsm(adata : AnnData, key : String) -> Bool {
  adata.obsm.contains(key)
}

///|
pub fn anndata_has_varm(adata : AnnData, key : String) -> Bool {
  adata.varm.contains(key)
}

///|
pub fn anndata_has_obsp(adata : AnnData, key : String) -> Bool {
  adata.obsp.contains(key)
}

///|
pub fn anndata_has_varp(adata : AnnData, key : String) -> Bool {
  adata.varp.contains(key)
}

///|
// === Dimension names ===

///|
/// Get dimension names for an obsm array.
pub fn anndata_obsm_names(adata : AnnData, key : String) -> Array[String]? {
  if adata.obsm.contains(key) {
    let m = adata.obsm[key]
    let names : Array[String] = Array(capacity=m[0].length())
    for j = 0; j < m[0].length(); j = j + 1 {
      names.push("obsm_" + key + "_" + j.to_string())
    }
    Some(names)
  } else {
    None
  }
}

///|
/// Get dimension names for a varm array.
pub fn anndata_varm_names(adata : AnnData, key : String) -> Array[String]? {
  if adata.varm.contains(key) {
    let m = adata.varm[key]
    let names : Array[String] = Array(capacity=m[0].length())
    for j = 0; j < m[0].length(); j = j + 1 {
      names.push("varm_" + key + "_" + j.to_string())
    }
    Some(names)
  } else {
    None
  }
}

///|
// === Delete operations ===

///|
/// Remove a layer by key. Returns true if key existed.
pub fn anndata_del_layer(adata : AnnData, key : String) -> Bool {
  if adata.layers.contains(key) {
    let _ = adata.layers.remove(key)
    true
  } else {
    false
  }
}

///|
/// Remove an obsm entry by key.
pub fn anndata_del_obsm(adata : AnnData, key : String) -> Bool {
  if adata.obsm.contains(key) {
    let _ = adata.obsm.remove(key)
    true
  } else {
    false
  }
}

///|
/// Remove a varm entry by key.
pub fn anndata_del_varm(adata : AnnData, key : String) -> Bool {
  if adata.varm.contains(key) {
    let _ = adata.varm.remove(key)
    true
  } else {
    false
  }
}

///|
/// Remove an obsp entry by key.
pub fn anndata_del_obsp(adata : AnnData, key : String) -> Bool {
  if adata.obsp.contains(key) {
    let _ = adata.obsp.remove(key)
    true
  } else {
    false
  }
}

///|
/// Remove a varp entry by key.
pub fn anndata_del_varp(adata : AnnData, key : String) -> Bool {
  if adata.varp.contains(key) {
    let _ = adata.varp.remove(key)
    true
  } else {
    false
  }
}
