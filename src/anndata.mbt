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
  self.x.length()
}

///|
pub fn AnnData::n_vars(self : AnnData) -> Int {
  if self.x.length() == 0 {
    0
  } else {
    self.x[0].length()
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
  let out : Array[String] = Array::new()
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
  let out : Array[Array[Double]] = Array::new()
  for i in idx {
    if i < mat.length() {
      out.push(mat[i].copy())
    } else {
      out.push(Array::new())
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
  let out : Array[Array[Double]] = Array::new()
  for r in r_idx {
    let row : Array[Double] = Array::new()
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
  AnnData::{
    x: pick_rc(self.x, obs_idx, var_idx),
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
    x: pick_rows(self.x, range(self.x.length())),
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
  let out : Array[Int] = Array::new()
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
pub fn AnnData::x_matrix(self : AnnData) -> Array[Array[Double]] {
  self.x
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
