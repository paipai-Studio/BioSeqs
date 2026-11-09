/// AnnData memory estimation.
///
/// Estimates memory usage of each AnnData component.

///|
/// Memory breakdown.
pub struct AnnDataMemory {
  x_bytes : Int64
  obs_bytes : Int64
  var_bytes : Int64
  layers_bytes : Int64
  obsm_bytes : Int64
  varm_bytes : Int64
  uns_bytes : Int64
  total_bytes : Int64
}

///|
pub fn AnnDataMemory::total_bytes(self : AnnDataMemory) -> Int64 {
  self.total_bytes
}

///|
pub fn AnnDataMemory::x_bytes(self : AnnDataMemory) -> Int64 {
  self.x_bytes
}

///|
/// Estimate memory usage of an AnnData object.
/// Each Double is 8 bytes, each String is ~50 bytes (overhead + chars).
pub fn anndata_estimate_memory(adata : AnnData) -> AnnDataMemory {
  let double_size = 8L
  let string_overhead = 50L
  // X: n_obs * n_vars * 8
  let x_bytes = adata.n_obs().to_int64() *
    adata.n_vars().to_int64() *
    double_size
  // obs_names: n_obs * string_overhead
  let obs_names_bytes = adata.n_obs().to_int64() * string_overhead
  // var_names: n_vars * string_overhead
  let var_names_bytes = adata.n_vars().to_int64() * string_overhead
  // obs annotation columns
  let mut obs_anno_bytes = 0L
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs_column(k) {
      Some(col) =>
        obs_anno_bytes = obs_anno_bytes +
          col.length().to_int64() * string_overhead
      None => ()
    }
  }
  let obs_bytes = obs_names_bytes + obs_anno_bytes
  // var annotation columns
  let mut var_anno_bytes = 0L
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.var_column(k) {
      Some(col) =>
        var_anno_bytes = var_anno_bytes +
          col.length().to_int64() * string_overhead
      None => ()
    }
  }
  let var_bytes = var_names_bytes + var_anno_bytes
  // layers
  let mut layers_bytes = 0L
  let layer_keys = adata.layers.keys().collect()
  for k in layer_keys {
    match adata.get_layer(k) {
      Some(mat) => {
        let mut n = 0L
        for row in mat {
          n = n + row.length().to_int64()
        }
        layers_bytes = layers_bytes + n * double_size
      }
      None => ()
    }
  }
  // obsm
  let mut obsm_bytes = 0L
  let obsm_keys = adata.obsm.keys().collect()
  for k in obsm_keys {
    match adata.get_obsm(k) {
      Some(mat) => {
        let mut n = 0L
        for row in mat {
          n = n + row.length().to_int64()
        }
        obsm_bytes = obsm_bytes + n * double_size
      }
      None => ()
    }
  }
  // varm
  let mut varm_bytes = 0L
  let varm_keys = adata.varm.keys().collect()
  for k in varm_keys {
    match adata.get_varm(k) {
      Some(mat) => {
        let mut n = 0L
        for row in mat {
          n = n + row.length().to_int64()
        }
        varm_bytes = varm_bytes + n * double_size
      }
      None => ()
    }
  }
  // uns (rough estimate: count keys * 100 bytes)
  let uns_keys = adata.uns.keys().collect()
  let uns_bytes = uns_keys.length().to_int64() * 100L
  let total = x_bytes +
    obs_bytes +
    var_bytes +
    layers_bytes +
    obsm_bytes +
    varm_bytes +
    uns_bytes
  {
    x_bytes,
    obs_bytes,
    var_bytes,
    layers_bytes,
    obsm_bytes,
    varm_bytes,
    uns_bytes,
    total_bytes: total,
  }
}

///|
/// Format bytes as human-readable string.
pub fn anndata_format_bytes(bytes : Int64) -> String {
  let b = bytes
  if b >= 1073741824L {
    (b.to_double() / 1073741824.0).to_string() + " GB"
  } else if b >= 1048576L {
    (b.to_double() / 1048576.0).to_string() + " MB"
  } else if b >= 1024L {
    (b.to_double() / 1024.0).to_string() + " KB"
  } else {
    b.to_string() + " B"
  }
}

///|
/// Print memory breakdown as a report string.
pub fn anndata_memory_report(adata : AnnData) -> String {
  let m = anndata_estimate_memory(adata)
  let mut out = "AnnData Memory Report\n"
  out = out + "  X:       " + anndata_format_bytes(m.x_bytes) + "\n"
  out = out + "  obs:     " + anndata_format_bytes(m.obs_bytes) + "\n"
  out = out + "  var:     " + anndata_format_bytes(m.var_bytes) + "\n"
  out = out + "  layers:  " + anndata_format_bytes(m.layers_bytes) + "\n"
  out = out + "  obsm:    " + anndata_format_bytes(m.obsm_bytes) + "\n"
  out = out + "  varm:    " + anndata_format_bytes(m.varm_bytes) + "\n"
  out = out + "  uns:     " + anndata_format_bytes(m.uns_bytes) + "\n"
  out = out + "  -------\n"
  out = out + "  total:   " + anndata_format_bytes(m.total_bytes) + "\n"
  out
}
