/// AnnData backed mode — memory-mapped file reading.
///
/// `anndata.read_h5ad(backed="r")` equivalent.
///
/// Keeps data on disk and loads chunks on demand.
/// This is a simplified simulation since MoonBit doesn't have
/// native mmap; we track file state and load on access.

///|
/// Backed AnnData: references a file and loads data lazily.
pub struct AnnDataBacked {
  filename : String
  n_obs : Int32
  n_vars : Int32
  obs_names : Array[String]
  var_names : Array[String]
  mut loaded : Bool // whether full data has been loaded
  cached : AnnData? // cached full load
  mode : String // "r" (read-only) or "r+" (writable)
  chunk_size : Int32 // chunk granularity for pad-aligned reads
}

///|
/// Open a backed AnnData from an h5ad file (read-only by default).
/// Only reads metadata (shape, names); X is loaded on demand.
/// `mode="r+"` enables chunk write-back (see `set_chunk` / `flush`).
pub fn anndata_backed_open(
  filename : String,
  content : String,
  mode~ : String = "r",
  chunk_size~ : Int32 = 1000,
) -> AnnDataBacked {
  let effective_chunk = if chunk_size <= 0 { 1000 } else { chunk_size }
  // Parse metadata from h5ad content
  match scanpy_read_h5ad(content) {
    Some(adata) =>
      {
        filename,
        n_obs: adata.n_obs(),
        n_vars: adata.n_vars(),
        obs_names: adata.get_obs_names(),
        var_names: adata.get_var_names(),
        loaded: false,
        cached: Some(adata),
        mode,
        chunk_size: effective_chunk,
      }
    None =>
      {
        filename,
        n_obs: 0,
        n_vars: 0,
        obs_names: [],
        var_names: [],
        loaded: false,
        cached: None,
        mode,
        chunk_size: effective_chunk,
      }
  }
}

///|
/// Number of observations (available without loading X).
pub fn AnnDataBacked::n_obs(self : AnnDataBacked) -> Int32 {
  self.n_obs
}

///|
/// Number of variables (available without loading X).
pub fn AnnDataBacked::n_vars(self : AnnDataBacked) -> Int32 {
  self.n_vars
}

///|
/// Shape (available without loading X).
pub fn AnnDataBacked::shape(self : AnnDataBacked) -> (Int32, Int32) {
  (self.n_obs, self.n_vars)
}

///|
/// Obs names (available without loading X).
pub fn AnnDataBacked::obs_names(self : AnnDataBacked) -> Array[String] {
  self.obs_names
}

///|
/// Var names (available without loading X).
pub fn AnnDataBacked::var_names(self : AnnDataBacked) -> Array[String] {
  self.var_names
}

///|
/// Load a single row from X (on demand).
pub fn AnnDataBacked::get_row(self : AnnDataBacked, i : Int32) -> Array[Double] {
  match self.cached {
    Some(adata) => {
      let row : Array[Double] = Array(capacity=self.n_vars)
      for j = 0; j < self.n_vars; j = j + 1 {
        row.push(adata.get_x(i, j))
      }
      row
    }
    None => Array::make(self.n_vars, 0.0)
  }
}

///|
/// Load a single column from X (on demand).
pub fn AnnDataBacked::get_col(self : AnnDataBacked, j : Int32) -> Array[Double] {
  match self.cached {
    Some(adata) => {
      let col : Array[Double] = Array(capacity=self.n_obs)
      for i = 0; i < self.n_obs; i = i + 1 {
        col.push(adata.get_x(i, j))
      }
      col
    }
    None => Array::make(self.n_obs, 0.0)
  }
}

///|
/// Load a single element from X.
pub fn AnnDataBacked::get_x(
  self : AnnDataBacked,
  i : Int32,
  j : Int32,
) -> Double {
  match self.cached {
    Some(adata) => adata.get_x(i, j)
    None => 0.0
  }
}

///|
/// Load a chunk of rows from X.
pub fn AnnDataBacked::get_chunk(
  self : AnnDataBacked,
  start : Int32,
  end_ : Int32,
) -> Array[Array[Double]] {
  match self.cached {
    Some(adata) => {
      let out : Array[Array[Double]] = Array(capacity=0)
      for i = start; i < end_ && i < self.n_obs; i = i + 1 {
        let row : Array[Double] = Array(capacity=self.n_vars)
        for j = 0; j < self.n_vars; j = j + 1 {
          row.push(adata.get_x(i, j))
        }
        out.push(row)
      }
      out
    }
    None => []
  }
}

///|
/// Load full data into memory (converts backed to in-memory).
pub fn AnnDataBacked::to_memory(self : AnnDataBacked) -> AnnData? {
  self.cached
}

///|
/// Check if data is fully loaded.
pub fn AnnDataBacked::is_loaded(self : AnnDataBacked) -> Bool {
  self.loaded
}

///|
/// Filename of the backing file.
pub fn AnnDataBacked::filename(self : AnnDataBacked) -> String {
  self.filename
}

///|
/// Compute row sums without loading full X (streaming).
pub fn AnnDataBacked::row_sums(self : AnnDataBacked) -> Array[Double] {
  match self.cached {
    Some(adata) => adata.row_sums()
    None => Array::make(self.n_obs, 0.0)
  }
}

///|
/// Compute col sums without loading full X (streaming).
pub fn AnnDataBacked::col_sums(self : AnnDataBacked) -> Array[Double] {
  match self.cached {
    Some(adata) => adata.col_sums()
    None => Array::make(self.n_vars, 0.0)
  }
}

// === r+ write-back (anndata backed="r+") ===

///|
/// True when the backing file was opened in `r+` mode (writable).
pub fn AnnDataBacked::is_writable(self : AnnDataBacked) -> Bool {
  self.mode == "r+"
}

///|
/// The open mode: "r" or "r+".
pub fn AnnDataBacked::mode(self : AnnDataBacked) -> String {
  self.mode
}

///|
/// Inplace write a chunk of rows back into the backing store (`r+` only).
/// `rows[i]` replaces source row `start + i`. Silently ignored in "r" mode
/// or when the data has not been loaded.
pub fn AnnDataBacked::set_chunk(
  self : AnnDataBacked,
  start : Int32,
  rows : Array[Array[Double]],
) -> Unit {
  if not(self.is_writable()) {
    return
  }
  match self.cached {
    Some(adata) => {
      if adata.is_sparse() {
        return
      }
      for i = 0; i < rows.length(); i = i + 1 {
        let r = start + i
        if r >= 0 && r < adata.x.length() {
          adata.x[r] = rows[i]
        }
      }
      self.loaded = true
    }
    None => ()
  }
}

///|
/// Flush in-place modifications back to the backing file (`r+` only).
/// Serializes the mutated in-memory object; the returned string is the new
/// file content to persist, or `None` when read-only / not loaded.
pub fn AnnDataBacked::flush(self : AnnDataBacked) -> String? {
  if not(self.is_writable()) {
    return None
  }
  match self.cached {
    Some(adata) => Some(scanpy_write_h5ad(adata))
    None => None
  }
}

// === pad-backed chunk-aligned reads ===

///|
/// Round `v` up to the next multiple of `chunk` (at least 1).
fn backed_pad_up(v : Int32, chunk : Int32) -> Int32 {
  if chunk <= 0 {
    return v
  }
  let rem = v % chunk
  if rem == 0 {
    v
  } else {
    v + (chunk - rem)
  }
}

///|
/// Chunk-aligned (padded) row read, mirroring `anndata._io.pad_backed`:
/// the range [start, end) is expanded outward to `chunk_size` boundaries
/// and the full padded chunk is returned (may extend past `n_obs` rows
/// are clamped to the data).
pub fn AnnDataBacked::pad_read(
  self : AnnDataBacked,
  start : Int32,
  end_ : Int32,
) -> Array[Array[Double]] {
  let chunk = if self.chunk_size <= 0 { 1000 } else { self.chunk_size }
  let s = if start < 0 { 0 } else { start }
  let e = if end_ > self.n_obs { self.n_obs } else { end_ }
  if s >= e {
    return []
  }
  let cs = s / chunk * chunk
  let ce = backed_pad_up(e, chunk)
  let out : Array[Array[Double]] = Array(capacity=0)
  match self.cached {
    Some(adata) =>
      for i = cs; i < ce && i < adata.n_obs(); i = i + 1 {
        let row : Array[Double] = Array(capacity=self.n_vars)
        for j = 0; j < self.n_vars; j = j + 1 {
          row.push(adata.get_x(i, j))
        }
        out.push(row)
      }
    None => ()
  }
  out
}

///|
/// Chunk size used for padded reads.
pub fn AnnDataBacked::chunk_size(self : AnnDataBacked) -> Int32 {
  self.chunk_size
}

// === backed subsetting ===

///|
/// Subset a backed AnnData along the obs axis. Mirrors anndata: subsetting
/// a backed object yields a new backed handle whose cached metadata (names,
/// shape) reflects the slice; X stays backed (rows resolved on demand).
pub fn AnnDataBacked::subset(
  self : AnnDataBacked,
  obs_idx : Array[Int32],
) -> AnnDataBacked {
  let new_obs_names : Array[String] = Array(capacity=obs_idx.length())
  for i in obs_idx {
    if i >= 0 && i < self.obs_names.length() {
      new_obs_names.push(self.obs_names[i])
    } else {
      new_obs_names.push("")
    }
  }
  let new_cached = match self.cached {
    Some(adata) => {
      let all_vars : Array[Int32] = Array(capacity=adata.n_vars())
      for j = 0; j < adata.n_vars(); j = j + 1 {
        all_vars.push(j)
      }
      Some(adata.subset(obs_idx, all_vars))
    }
    None => None
  }
  {
    filename: self.filename,
    n_obs: new_obs_names.length(),
    n_vars: self.n_vars,
    obs_names: new_obs_names,
    var_names: self.var_names.copy(),
    loaded: false,
    cached: new_cached,
    mode: self.mode,
    chunk_size: self.chunk_size,
  }
}

///|
/// True when the backed data has been materialized into memory.
/// `to_memory()` + any access marks the handle as loaded.
pub fn AnnDataBacked::mark_loaded(self : AnnDataBacked) -> Unit {
  self.loaded = true
}
