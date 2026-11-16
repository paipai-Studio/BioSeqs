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

// ============================================================================
// AnnCollectionView — AnnCollection 的轻量 view 包装层
// ============================================================================

///|
/// AnnCollection 的轻量 view：持有 AnnCollection + 可选 obs/var 索引过滤器。
/// 不复制数据，仅通过索引映射访问底层 AnnCollection。
pub struct AnnCollectionView {
  coll : AnnCollection
  obs_idx : Array[Int32] // 若非空 → 只暴露这些 obs；空 → 全部
  var_idx : Array[Int32] // 若非空 → 只暴露这些 var；空 → 全部
}

///|
/// 从 AnnCollection 创建完整 view（无过滤）。
pub fn AnnCollectionView::from(coll : AnnCollection) -> AnnCollectionView {
  AnnCollectionView::{
    coll,
    obs_idx: Array(capacity=0),
    var_idx: Array(capacity=0),
  }
}

///|
/// 创建 view 并指定 obs 过滤索引。
pub fn AnnCollectionView::with_obs(
  coll : AnnCollection,
  obs_idx : Array[Int32],
) -> AnnCollectionView {
  AnnCollectionView::{ coll, obs_idx, var_idx: Array(capacity=0), }
}

///|
/// 创建 view 并指定 var 过滤索引。
pub fn AnnCollectionView::with_vars(
  coll : AnnCollection,
  var_idx : Array[Int32],
) -> AnnCollectionView {
  AnnCollectionView::{ coll, obs_idx: Array(capacity=0), var_idx, }
}

///|
/// 创建 view 并同时指定 obs + var 过滤索引。
pub fn AnnCollectionView::with_both(
  coll : AnnCollection,
  obs_idx : Array[Int32],
  var_idx : Array[Int32],
) -> AnnCollectionView {
  AnnCollectionView::{ coll, obs_idx, var_idx, }
}

///|
/// view 中的 obs 数量。
pub fn AnnCollectionView::n_obs(self : AnnCollectionView) -> Int32 {
  if self.obs_idx.length() > 0 {
    self.obs_idx.length()
  } else {
    self.coll.n_obs()
  }
}

///|
/// view 中的 var 数量。
pub fn AnnCollectionView::n_vars(self : AnnCollectionView) -> Int32 {
  if self.var_idx.length() > 0 {
    self.var_idx.length()
  } else {
    self.coll.n_vars()
  }
}

///|
/// 取 view 中第 `i` 个 obs 的完整表达向量。
pub fn AnnCollectionView::get_obs(
  self : AnnCollectionView,
  i : Int32,
) -> Array[Double]? {
  // 映射 i → global obs idx
  let global_i = if self.obs_idx.length() > 0 {
    if i >= 0 && i < self.obs_idx.length() {
      self.obs_idx[i]
    } else {
      return None
    }
  } else {
    i
  }
  let full_row = match self.coll.get_obs(global_i) {
    None => return None
    Some(r) => r
  }
  // 再按 var_idx 过滤
  if self.var_idx.length() == 0 {
    return Some(full_row)
  }
  let out : Array[Double] = Array(capacity=self.var_idx.length())
  for j = 0; j < self.var_idx.length(); j = j + 1 {
    let vj = self.var_idx[j]
    if vj >= 0 && vj < full_row.length() {
      out.push(full_row[vj])
    } else {
      out.push(0.0)
    }
  }
  Some(out)
}

///|
/// 取 view 中第 i 个 obs、第 j 个 var 的值。
pub fn AnnCollectionView::get_x(
  self : AnnCollectionView,
  i : Int32,
  j : Int32,
) -> Double? {
  // obs 映射 + 边界
  let global_i = if self.obs_idx.length() > 0 {
    if i >= 0 && i < self.obs_idx.length() {
      self.obs_idx[i]
    } else {
      return None
    }
  } else {
    if i < 0 || i >= self.coll.n_obs() {
      return None
    }
    i
  }
  // var 映射 + 边界
  let global_j = if self.var_idx.length() > 0 {
    if j >= 0 && j < self.var_idx.length() {
      self.var_idx[j]
    } else {
      return None
    }
  } else {
    if j < 0 || j >= self.coll.n_vars() {
      return None
    }
    j
  }
  self.coll.get_x(global_i, global_j)
}

///|
/// Materialize view → AnnData（带过滤后的 obs/var_names）。
pub fn AnnCollectionView::to_adata(self : AnnCollectionView) -> AnnData {
  let n_o = self.n_obs()
  let n_v = self.n_vars()
  let x : Array[Array[Double]] = Array(capacity=n_o)
  let obs_names : Array[String] = Array(capacity=n_o)
  // 取原 var_names（需要先拿到 full 源数据）
  // 先从 coll.to_adata 取一次完整的，方便拿 var_names + obs_names
  let full = self.coll.to_adata()
  let full_var = full.var_names
  let full_obs = full.obs_names
  // 构造过滤后的 var_names
  let var_names : Array[String] = Array(capacity=n_v)
  if self.var_idx.length() > 0 {
    for j in self.var_idx {
      if j >= 0 && j < full_var.length() {
        var_names.push(full_var[j])
      } else {
        var_names.push("var_" + j.to_string())
      }
    }
  } else {
    for v in full_var {
      var_names.push(v)
    }
  }
  // 构造过滤后的 x + obs_names
  for i = 0; i < n_o; i = i + 1 {
    let global_i = if self.obs_idx.length() > 0 { self.obs_idx[i] } else { i }
    // obs name
    if global_i >= 0 && global_i < full_obs.length() {
      obs_names.push(full_obs[global_i])
    } else {
      obs_names.push("obs_" + global_i.to_string())
    }
    // 表达向量
    match self.get_obs(i) {
      Some(row) => x.push(row)
      None => x.push(Array(capacity=n_v))
    }
  }
  AnnData::new(x, obs_names, var_names)
}
