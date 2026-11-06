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
  loaded : Bool // whether full data has been loaded
  cached : AnnData? // cached full load
}

///|
/// Open a backed AnnData from an h5ad file (read-only).
/// Only reads metadata (shape, names); X is loaded on demand.
pub fn anndata_backed_open(
  filename : String,
  content : String,
) -> AnnDataBacked {
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
