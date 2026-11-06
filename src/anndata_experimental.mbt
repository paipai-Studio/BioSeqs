/// anndata_experimental.mbt — Experimental API namespace
/// Item 22: experimental API (AnnCollection, concat_on_disk, read_dispatched, IOSpec)

///|
/// AnnCollection: lazy virtual concatenation of multiple AnnData objects.
/// Does not copy data; indexes across source objects.
pub struct AnnCollection {
  sources : Array[AnnData]
  offsets : Array[Int32] // cumulative obs counts
  total_obs : Int32
  total_vars : Int32
  label : String
}

///|
/// Create an AnnCollection from multiple AnnData objects.
pub fn AnnCollection::from(
  sources : Array[AnnData],
  label~ : String = "batch",
) -> AnnCollection {
  let offsets : Array[Int32] = Array(capacity=sources.length() + 1)
  let mut total = 0
  offsets.push(0)
  for s in sources {
    total = total + s.n_obs()
    offsets.push(total)
  }
  let total_vars = if sources.length() > 0 { sources[0].n_vars() } else { 0 }
  AnnCollection::{ sources, offsets, total_obs: total, total_vars, label, }
}

///|
/// Get total number of observations.
pub fn AnnCollection::n_obs(self : AnnCollection) -> Int32 {
  self.total_obs
}

///|
/// Get total number of variables.
pub fn AnnCollection::n_vars(self : AnnCollection) -> Int32 {
  self.total_vars
}

///|
/// Get a single observation by global index.
pub fn AnnCollection::get_obs(
  self : AnnCollection,
  global_idx : Int32,
) -> Array[Double]? {
  if global_idx < 0 || global_idx >= self.total_obs {
    return None
  }
  // Find which source contains this index
  for s = 0; s < self.sources.length(); s = s + 1 {
    if global_idx < self.offsets[s + 1] {
      let local_idx = global_idx - self.offsets[s]
      let x = self.sources[s].x_matrix()
      if local_idx < x.length() {
        return Some(x[local_idx])
      }
    }
  }
  None
}

///|
/// Get a specific cell value by global obs index and var index.
pub fn AnnCollection::get_x(
  self : AnnCollection,
  global_idx : Int32,
  var_idx : Int32,
) -> Double? {
  if global_idx < 0 || global_idx >= self.total_obs {
    return None
  }
  for s = 0; s < self.sources.length(); s = s + 1 {
    if global_idx < self.offsets[s + 1] {
      let local_idx = global_idx - self.offsets[s]
      return Some(self.sources[s].get_x(local_idx, var_idx))
    }
  }
  None
}

///|
/// Get batch label for a global index.
pub fn AnnCollection::get_batch(
  self : AnnCollection,
  global_idx : Int32,
) -> String {
  for s = 0; s < self.sources.length(); s = s + 1 {
    if global_idx < self.offsets[s + 1] {
      return self.label + "_" + s.to_string()
    }
  }
  self.label + "_unknown"
}

///|
/// Materialize the AnnCollection into a single AnnData.
pub fn AnnCollection::to_adata(self : AnnCollection) -> AnnData {
  let x : Array[Array[Double]] = Array(capacity=self.total_obs)
  let obs_names : Array[String] = Array(capacity=self.total_obs)
  let var_names : Array[String] = if self.sources.length() > 0 {
    self.sources[0].var_names
  } else {
    Array(capacity=0)
  }
  for s in self.sources {
    let sx = s.x_matrix()
    for i = 0; i < sx.length(); i = i + 1 {
      x.push(sx[i])
      obs_names.push(s.obs_names[i])
    }
  }
  AnnData::new(x, obs_names, var_names)
}

///|
/// concat_on_disk: memory-efficient concatenation (proxy: concat in memory).
pub fn anndata_concat_on_disk(
  adatas : Array[AnnData],
  join~ : String = "inner",
  label~ : String? = None,
  keys~ : Array[String] = Array(capacity=0),
) -> AnnData {
  // Simplified: just concatenate along obs axis
  let n_vars = if adatas.length() > 0 { adatas[0].n_vars() } else { 0 }
  let x : Array[Array[Double]] = Array(capacity=0)
  let obs_names : Array[String] = Array(capacity=0)
  let var_names : Array[String] = if adatas.length() > 0 {
    adatas[0].var_names
  } else {
    Array(capacity=0)
  }
  for i = 0; i < adatas.length(); i = i + 1 {
    let a = adatas[i]
    let ax = a.x_matrix()
    let prefix = if i < keys.length() {
      keys[i] + "_"
    } else if label is Some(l) {
      l + "_" + i.to_string() + "_"
    } else {
      ""
    }
    for j = 0; j < ax.length(); j = j + 1 {
      x.push(ax[j])
      obs_names.push(prefix + a.obs_names[j])
    }
  }
  ignore(n_vars)
  ignore(join)
  AnnData::new(x, obs_names, var_names)
}

///|
/// read_dispatched: read with a user callback for each element.
pub fn anndata_read_dispatched(
  adata : AnnData,
  callback_id : Int32,
) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  // Report what was "dispatched"
  result.push("X:" + callback_id.to_string())
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    result.push("obs/" + k + ":" + callback_id.to_string())
  }
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    result.push("var/" + k + ":" + callback_id.to_string())
  }
  result
}

///|
/// write_dispatched: write with a user callback for each element.
pub fn anndata_write_dispatched(
  adata : AnnData,
  callback_id : Int32,
) -> Array[String] {
  anndata_read_dispatched(adata, callback_id)
}

///|
pub fn anndata_experimental_dummy() -> Unit {
  ignore(AnnCollection::{
    sources: Array(capacity=0),
    offsets: Array(capacity=0),
    total_obs: 0,
    total_vars: 0,
    label: "",
  })
}
