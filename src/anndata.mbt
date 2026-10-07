///|
/// AnnData container (scanpy/anndata) — the single-cell data object that
/// scanpy's preprocessing, tools and plotting layers operate on. Mirrors
/// anndata's field set: `x` (n_obs × n_vars), per-observation annotation `obs`,
/// per-variable annotation `var`, named alternative matrices `layers`, cell
/// embeddings `obsm` / `varm` and pairwise matrices `obsp` / `varp`, plus
/// `obs_names` / `var_names` and unstructured annotation `uns` (string-valued
/// simplification of anndata's arbitrary nested dict).

///|
pub struct AnnData {
  x : Array[Array[Double]]
  x_sparse : SparseArray?            // sparse storage (when set, x is empty)
  obs : Map[String, Array[String]]
  variables : Map[String, Array[String]]
  layers : Map[String, Array[Array[Double]]]
  obsm : Map[String, Array[Array[Double]]]
  varm : Map[String, Array[Array[Double]]]
  obsp : Map[String, Array[Array[Double]]]
  varp : Map[String, Array[Array[Double]]]
  obs_names : Array[String]
  var_names : Array[String]
  uns : Map[String, String]
}

///|
/// Construct an AnnData from a dense `n_obs × n_vars` matrix `x` (rows are
/// observations/cells, columns are variables/genes) with the observation and
/// variable names. All metadata slots start empty.
pub fn AnnData::new(
  x : Array[Array[Double]],
  obs_names : Array[String],
  var_names : Array[String],
) -> AnnData {
  AnnData::{
    x,
    x_sparse: None,
    obs: Map([], capacity=4),
    variables: Map([], capacity=4),
    layers: Map([], capacity=4),
    obsm: Map([], capacity=4),
    varm: Map([], capacity=4),
    obsp: Map([], capacity=4),
    varp: Map([], capacity=4),
    obs_names,
    var_names,
    uns: Map([], capacity=4),
  }
}

///|
pub fn AnnData::n_obs(self : AnnData) -> Int {
  match self.x_sparse {
    Some(s) => {
      let dims = s.dim()
      if dims.length() >= 1 { dims[0] } else { 0 }
    }
    None => self.x.length()
  }
}

///|
pub fn AnnData::n_vars(self : AnnData) -> Int {
  match self.x_sparse {
    Some(s) => {
      let dims = s.dim()
      if dims.length() >= 2 { dims[1] } else { 0 }
    }
    None => {
      if self.x.length() == 0 {
        0
      } else {
        self.x[0].length()
      }
    }
  }
}

///|
/// (n_obs, n_vars) shape, matching `AnnData.shape`.
pub fn AnnData::shape(self : AnnData) -> (Int, Int) {
  (self.n_obs(), self.n_vars())
}

///|
/// Set an observation (cell) annotation column. Returns `self` for chaining.
pub fn AnnData::set_obs(self : AnnData, key : String, values : Array[String]) -> AnnData {
  self.obs[key] = values
  self
}

///|
/// Set a variable (gene) annotation column. Returns `self` for chaining.
pub fn AnnData::set_var(self : AnnData, key : String, values : Array[String]) -> AnnData {
  self.variables[key] = values
  self
}

///|
/// Set a named alternative matrix layer (`n_obs × n_vars`). Returns `self`.
pub fn AnnData::set_layer(
  self : AnnData,
  name : String,
  mat : Array[Array[Double]],
) -> AnnData {
  self.layers[name] = mat
  self
}

///|
/// Set a cell-embedding matrix under `name` (`n_obs × d`). Returns `self`.
pub fn AnnData::set_obsm(
  self : AnnData,
  name : String,
  mat : Array[Array[Double]],
) -> AnnData {
  self.obsm[name] = mat
  self
}

///|
/// Set a cell × cell pairwise matrix under `name`. Returns `self`.
pub fn AnnData::set_obsp(
  self : AnnData,
  name : String,
  mat : Array[Array[Double]],
) -> AnnData {
  self.obsp[name] = mat
  self
}

///|
/// Obs (cell) annotation column, or `None` if absent.
pub fn AnnData::obs_column(self : AnnData, key : String) -> Array[String]? {
  match self.obs.get(key) {
    Some(v) => Some(v)
    None => None
  }
}

///|
/// Var (gene) annotation column, or `None` if absent.
pub fn AnnData::var_column(self : AnnData, key : String) -> Array[String]? {
  match self.variables.get(key) {
    Some(v) => Some(v)
    None => None
  }
}

///|
fn pick_strs(a : Array[String], idx : Array[Int]) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for i in idx {
    if i < a.length() {
      out.push(a[i])
    } else {
      out.push("")
    }
  }
  out
}

///|
fn pick_rows(mat : Array[Array[Double]], idx : Array[Int]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for i in idx {
    if i < mat.length() {
      out.push(mat[i].copy())
    } else {
      out.push(Array(capacity=0))
    }
  }
  out
}

///|
fn pick_rc(
  mat : Array[Array[Double]],
  r_idx : Array[Int],
  c_idx : Array[Int],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for r in r_idx {
    let row : Array[Double] = Array(capacity=0)
    if r < mat.length() {
      for c in c_idx {
        if c < mat[r].length() {
          row.push(mat[r][c])
        } else {
          row.push(0.0)
        }
      }
    }
    out.push(row)
  }
  out
}

///|
fn pick_str_map_rows(
  m : Map[String, Array[String]],
  idx : Array[Int],
) -> Map[String, Array[String]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[String]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = pick_strs(m[k], idx)
  }
  out
}

///|
fn pick_rows_map(
  m : Map[String, Array[Array[Double]]],
  idx : Array[Int],
) -> Map[String, Array[Array[Double]]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = pick_rows(m[k], idx)
  }
  out
}

///|
fn pick_rc_map(
  m : Map[String, Array[Array[Double]]],
  r_idx : Array[Int],
  c_idx : Array[Int],
) -> Map[String, Array[Array[Double]]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = pick_rc(m[k], r_idx, c_idx)
  }
  out
}

///|
fn pick_square_map(
  m : Map[String, Array[Array[Double]]],
  idx : Array[Int],
) -> Map[String, Array[Array[Double]]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = pick_rc(m[k], idx, idx)
  }
  out
}

///|
/// Subset observations and variables by index (`AnnData[...]` slicing). Returns
/// a new AnnData whose `x`/`layers` are sliced on both axes, `obs`/`obsm`/`obsp`
/// along observations, `var`/`varm`/`varp` along variables, and `obs_names`/
/// `var_names` reindexed. `uns` is carried over unchanged.
pub fn AnnData::subset(
  self : AnnData,
  obs_idx : Array[Int],
  var_idx : Array[Int],
) -> AnnData {
  let x_sparse_new = match self.x_sparse {
    Some(s) => {
      Some(s.subset([obs_idx, var_idx]) catch { _ => s })
    }
    None => None
  }
  AnnData::{
    x: match x_sparse_new {
      Some(_) => []
      None => pick_rc(self.x, obs_idx, var_idx)
    },
    x_sparse: x_sparse_new,
    obs: pick_str_map_rows(self.obs, obs_idx),
    variables: pick_str_map_rows(self.variables, var_idx),
    layers: pick_rc_map(self.layers, obs_idx, var_idx),
    obsm: pick_rows_map(self.obsm, obs_idx),
    varm: pick_rows_map(self.varm, var_idx),
    obsp: pick_square_map(self.obsp, obs_idx),
    varp: pick_square_map(self.varp, var_idx),
    obs_names: pick_strs(self.obs_names, obs_idx),
    var_names: pick_strs(self.var_names, var_idx),
    uns: self.uns,
  }
}

///|
fn copy_str_map(m : Map[String, Array[String]]) -> Map[String, Array[String]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[String]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = m[k].copy()
  }
  out
}

///|
fn copy_mat_map(m : Map[String, Array[Array[Double]]]) -> Map[String, Array[Array[Double]]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = pick_rows(m[k], range(m[k].length()))
  }
  out
}

///|
/// Deep copy of the AnnData object.
pub fn AnnData::copy(self : AnnData) -> AnnData {
  AnnData::{
    x: match self.x_sparse {
      Some(_) => []
      None => pick_rows(self.x, range(self.x.length()))
    },
    x_sparse: self.x_sparse,
    obs: copy_str_map(self.obs),
    variables: copy_str_map(self.variables),
    layers: copy_mat_map(self.layers),
    obsm: copy_mat_map(self.obsm),
    varm: copy_mat_map(self.varm),
    obsp: copy_mat_map(self.obsp),
    varp: copy_mat_map(self.varp),
    obs_names: self.obs_names.copy(),
    var_names: self.var_names.copy(),
    uns: copy_string_map(self.uns),
  }
}
///|
/// Indices `0..n-1`.
fn range(n : Int) -> Array[Int] {
  let out : Array[Int] = Array(capacity=0)
  let mut i = 0
  while i < n {
    out.push(i)
    i = i + 1
  }
  out
}

///|
fn copy_string_map(m : Map[String, String]) -> Map[String, String] {
  let keys = m.keys().collect()
  let out : Map[String, String] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = m[k]
  }
  out
}

///|
/// Return the dense `n_obs × n_vars` expression matrix.
/// If sparse storage is used, materializes the dense matrix on demand.
pub fn AnnData::x_matrix(self : AnnData) -> Array[Array[Double]] {
  match self.x_sparse {
    Some(s) => {
      s.to_dense_matrix() catch { _ => [] }
    }
    None => self.x
  }
}

///|
/// Observation (cell) names.
pub fn AnnData::get_obs_names(self : AnnData) -> Array[String] {
  self.obs_names
}

///|
/// Variable (gene) names.
pub fn AnnData::get_var_names(self : AnnData) -> Array[String] {
  self.var_names
}


// ============================================================================
// Sparse storage support (P3 optimization)
//
// Single-cell expression matrices are ~90% zeros.  Storing them as a dense
// Array[Array[Double]] wastes memory: 10K cells × 20K genes ≈ 1.6 GB dense
// vs ~50 MB sparse.  When `x_sparse` is `Some(s)`, the dense `x` field is
// empty and all shape/matrix queries are served from the SparseArray.
// ============================================================================

///|
/// Construct an AnnData backed by a SparseArray (CSR-like COO internally).
/// The dense `x` field is left empty; `x_matrix()` materializes on demand.
pub fn AnnData::new_sparse(
  x_sparse : SparseArray,
  obs_names : Array[String],
  var_names : Array[String],
) -> AnnData {
  AnnData::{
    x: [],
    x_sparse: Some(x_sparse),
    obs: Map([], capacity=4),
    variables: Map([], capacity=4),
    layers: Map([], capacity=4),
    obsm: Map([], capacity=4),
    varm: Map([], capacity=4),
    obsp: Map([], capacity=4),
    varp: Map([], capacity=4),
    obs_names,
    var_names,
    uns: Map([], capacity=4),
  }
}

///|
/// True when the expression matrix is stored sparsely.
pub fn AnnData::is_sparse(self : AnnData) -> Bool {
  match self.x_sparse {
    Some(_) => true
    None => false
  }
}

///|
/// Underlying SparseArray, or `None` if dense storage is used.
pub fn AnnData::get_sparse(self : AnnData) -> SparseArray? {
  self.x_sparse
}

///|
/// Number of non-zero entries in the expression matrix.
/// For dense storage, counts non-zero elements in `x`.
pub fn AnnData::nnz(self : AnnData) -> Int {
  match self.x_sparse {
    Some(s) => s.nnzero()
    None => {
      let mut count = 0
      for row in self.x {
        for v in row {
          if v != 0.0 {
            count = count + 1
          }
        }
      }
      count
    }
  }
}

///|
/// Fill density (nnz / total).  Returns 0.0 for empty matrices.
pub fn AnnData::density(self : AnnData) -> Double {
  let n = self.n_obs()
  let m = self.n_vars()
  if n == 0 || m == 0 {
    return 0.0
  }
  let total = n * m
  self.nnz().to_double() / total.to_double()
}

///|
/// Row (per-cell) sums of the expression matrix.
/// Uses sparse-native operation when available.
pub fn AnnData::row_sums(self : AnnData) -> Array[Double] {
  match self.x_sparse {
    Some(s) => s.row_sums() catch { _ => [] }
    None => {
      let out : Array[Double] = Array(capacity=self.x.length())
      for row in self.x {
        let mut sum = 0.0
        for v in row {
          sum = sum + v
        }
        out.push(sum)
      }
      out
    }
  }
}

///|
/// Column (per-gene) sums of the expression matrix.
/// Uses sparse-native operation when available.
pub fn AnnData::col_sums(self : AnnData) -> Array[Double] {
  match self.x_sparse {
    Some(s) => s.column_sums() catch { _ => [] }
    None => {
      let n = self.n_vars()
      let out : Array[Double] = Array::make(n, 0.0)
      for row in self.x {
        for j in 0..<row.length() {
          out[j] = out[j] + row[j]
        }
      }
      out
    }
  }
}

///|
/// Row (per-cell) means of the expression matrix.
pub fn AnnData::row_means(self : AnnData) -> Array[Double] {
  let n_vars = self.n_vars().to_double()
  let sums = self.row_sums()
  if n_vars == 0.0 {
    return Array::make(sums.length(), 0.0)
  }
  let out : Array[Double] = Array(capacity=sums.length())
  for s in sums {
    out.push(s / n_vars)
  }
  out
}

///|
/// Column (per-gene) means of the expression matrix.
pub fn AnnData::col_means(self : AnnData) -> Array[Double] {
  let n_obs = self.n_obs().to_double()
  let sums = self.col_sums()
  if n_obs == 0.0 {
    return Array::make(sums.length(), 0.0)
  }
  let out : Array[Double] = Array(capacity=sums.length())
  for s in sums {
    out.push(s / n_obs)
  }
  out
}

///|
/// Convert dense storage to sparse.  Returns a new AnnData whose `x_sparse`
/// is populated and `x` is empty.  If already sparse, returns a copy.
pub fn AnnData::to_sparse(self : AnnData) -> AnnData {
  match self.x_sparse {
    Some(_) => self.copy()
    None => {
      let sparse = SparseArray::from_matrix(self.x) catch {
        _ => return self.copy()
      }
      AnnData::{
        x: [],
        x_sparse: Some(sparse),
        obs: copy_str_map(self.obs),
        variables: copy_str_map(self.variables),
        layers: copy_mat_map(self.layers),
        obsm: copy_mat_map(self.obsm),
        varm: copy_mat_map(self.varm),
        obsp: copy_mat_map(self.obsp),
        varp: copy_mat_map(self.varp),
        obs_names: self.obs_names.copy(),
        var_names: self.var_names.copy(),
        uns: copy_string_map(self.uns),
      }
    }
  }
}

///|
/// Convert sparse storage to dense.  Returns a new AnnData whose `x` is
/// populated and `x_sparse` is `None`.  If already dense, returns a copy.
pub fn AnnData::to_dense(self : AnnData) -> AnnData {
  match self.x_sparse {
    None => self.copy()
    Some(s) => {
      let dense = s.to_dense_matrix() catch { _ => [] }
      AnnData::{
        x: dense,
        x_sparse: None,
        obs: copy_str_map(self.obs),
        variables: copy_str_map(self.variables),
        layers: copy_mat_map(self.layers),
        obsm: copy_mat_map(self.obsm),
        varm: copy_mat_map(self.varm),
        obsp: copy_mat_map(self.obsp),
        varp: copy_mat_map(self.varp),
        obs_names: self.obs_names.copy(),
        var_names: self.var_names.copy(),
        uns: copy_string_map(self.uns),
      }
    }
  }
}

///|
/// Get a single element `x[i, j]` from the expression matrix.
/// Works for both sparse and dense storage.
pub fn AnnData::get_x(self : AnnData, i : Int, j : Int) -> Double {
  match self.x_sparse {
    Some(s) => s.get([i, j]) catch { _ => 0.0 }
    None => {
      if i < self.x.length() && j < self.x[i].length() {
        self.x[i][j]
      } else {
        0.0
      }
    }
  }
}

///|
/// Total sum of all elements in the expression matrix.
pub fn AnnData::x_sum(self : AnnData) -> Double {
  match self.x_sparse {
    Some(s) => s.sum()
    None => {
      let mut total = 0.0
      for row in self.x {
        for v in row {
          total = total + v
        }
      }
      total
    }
  }
}
