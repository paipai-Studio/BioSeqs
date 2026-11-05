/// AnnData raw attribute — `adata.raw` equivalent.
///
/// Frozen snapshot of pre-processing X and var annotations.
/// Allows downstream tools to access original counts after
/// normalization, filtering, or transformation.

///|
/// Raw AnnData snapshot
pub struct AnnDataRaw {
  x : Array[Array[Double]] // frozen expression matrix
  var_names : Array[String] // frozen variable names
  var_anno : Map[String, Array[String]] // frozen var annotations
}

///|
/// Capture a raw snapshot of an AnnData's current state.
/// Call this BEFORE preprocessing (normalize, filter, transform).
pub fn anndata_raw_capture(adata : AnnData) -> AnnDataRaw {
  let x_copy : Array[Array[Double]] = Array(capacity=adata.n_obs())
  for i = 0; i < adata.n_obs(); i = i + 1 {
    let row : Array[Double] = Array(capacity=adata.n_vars())
    for j = 0; j < adata.n_vars(); j = j + 1 {
      row.push(adata.get_x(i, j))
    }
    x_copy.push(row)
  }
  let var_keys = adata.variables.keys().collect()
  let var_anno : Map[String, Array[String]] = Map(
    [],
    capacity=var_keys.length(),
  )
  for k in var_keys {
    match adata.var_column(k) {
      Some(v) => var_anno[k] = v.copy()
      None => ()
    }
  }
  { x: x_copy, var_names: adata.get_var_names().copy(), var_anno, }
}

///|
/// Get raw expression value at (i, j).
pub fn AnnDataRaw::get_x(self : AnnDataRaw, i : Int32, j : Int32) -> Double {
  if i < self.x.length() && j < self.x[i].length() {
    self.x[i][j]
  } else {
    0.0
  }
}

///|
/// Number of observations in raw snapshot.
pub fn AnnDataRaw::n_obs(self : AnnDataRaw) -> Int32 {
  self.x.length()
}

///|
/// Number of variables in raw snapshot.
pub fn AnnDataRaw::n_vars(self : AnnDataRaw) -> Int32 {
  if self.x.length() > 0 {
    self.x[0].length()
  } else {
    0
  }
}

///|
/// Get raw var_names.
pub fn AnnDataRaw::var_names(self : AnnDataRaw) -> Array[String] {
  self.var_names
}

///|
/// Get a raw var annotation column.
pub fn AnnDataRaw::var_column(
  self : AnnDataRaw,
  key : String,
) -> Array[String]? {
  match self.var_anno.get(key) {
    Some(v) => Some(v)
    None => None
  }
}

///|
/// Get a full row from raw X (all genes for one cell).
pub fn AnnDataRaw::get_row(self : AnnDataRaw, i : Int32) -> Array[Double] {
  if i < self.x.length() {
    self.x[i].copy()
  } else {
    []
  }
}

///|
/// Get a full column from raw X (all cells for one gene).
pub fn AnnDataRaw::get_col(self : AnnDataRaw, j : Int32) -> Array[Double] {
  let out : Array[Double] = Array(capacity=self.x.length())
  for i = 0; i < self.x.length(); i = i + 1 {
    if j < self.x[i].length() {
      out.push(self.x[i][j])
    } else {
      out.push(0.0)
    }
  }
  out
}

///|
/// Subset raw along obs axis (e.g. after filtering cells).
pub fn AnnDataRaw::subset_obs(
  self : AnnDataRaw,
  idx : Array[Int32],
) -> AnnDataRaw {
  let x_new : Array[Array[Double]] = Array(capacity=idx.length())
  for i in idx {
    if i < self.x.length() {
      x_new.push(self.x[i].copy())
    } else {
      x_new.push([])
    }
  }
  { x: x_new, var_names: self.var_names.copy(), var_anno: self.var_anno, }
}

///|
/// Total sum of raw X.
pub fn AnnDataRaw::sum(self : AnnDataRaw) -> Double {
  let mut total = 0.0
  for row in self.x {
    for v in row {
      total = total + v
    }
  }
  total
}
