/// anndata_writeback2.mbt — Enhanced writeback mode for views
/// Item 14: writeback mode完善

///|
/// Writeback view: modifications to the view propagate to the source AnnData.
pub struct AnnDataWriteback {
  source : AnnData
  obs_indices : Array[Int32]
  var_indices : Array[Int32]
  is_obs_view : Bool
  is_var_view : Bool
}

///|
/// Create a writeback view.
pub fn AnnDataWriteback::create(
  source : AnnData,
  obs_indices~ : Array[Int32] = Array(capacity=0),
  var_indices~ : Array[Int32] = Array(capacity=0),
) -> AnnDataWriteback {
  let is_obs = obs_indices.length() > 0
  let is_var = var_indices.length() > 0
  AnnDataWriteback::{
    source,
    obs_indices,
    var_indices,
    is_obs_view: is_obs,
    is_var_view: is_var,
  }
}

///|
/// Read a value from the writeback view.
pub fn AnnDataWriteback::get_x(
  self : AnnDataWriteback,
  i : Int32,
  j : Int32,
) -> Double {
  let src_i = if self.is_obs_view && i < self.obs_indices.length() {
    self.obs_indices[i]
  } else {
    i
  }
  let src_j = if self.is_var_view && j < self.var_indices.length() {
    self.var_indices[j]
  } else {
    j
  }
  self.source.get_x(src_i, src_j)
}

///|
/// Write a value back to the source AnnData.
pub fn AnnDataWriteback::set_x(
  self : AnnDataWriteback,
  i : Int32,
  j : Int32,
  value : Double,
) -> Unit {
  let src_i = if self.is_obs_view && i < self.obs_indices.length() {
    self.obs_indices[i]
  } else {
    i
  }
  let src_j = if self.is_var_view && j < self.var_indices.length() {
    self.var_indices[j]
  } else {
    j
  }
  // Write to source X matrix
  let x = self.source.x_matrix()
  let si = src_i
  let sj = src_j
  if si < x.length() && sj < x[si].length() {
    x[si][sj] = value
  }
}

///|
/// Apply a function to all elements in the view, writing back.
pub fn AnnDataWriteback::apply(self : AnnDataWriteback, fn_id : Int32) -> Unit {
  let n_obs = if self.is_obs_view {
    self.obs_indices.length()
  } else {
    self.source.n_obs()
  }
  let n_vars = if self.is_var_view {
    self.var_indices.length()
  } else {
    self.source.n_vars()
  }
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let val = self.get_x(i, j)
      let new_val = ad_wb_apply_fn(fn_id, val)
      self.set_x(i, j, new_val)
    }
  }
}

///|
/// Commit all changes (no-op since we write directly, but provided for API parity).
pub fn AnnDataWriteback::commit(self : AnnDataWriteback) -> Unit {
  ignore(self)
}

///|
/// Get the underlying source AnnData.
pub fn AnnDataWriteback::get_source(self : AnnDataWriteback) -> AnnData {
  self.source
}

///|
fn ad_wb_apply_fn(fn_id : Int32, x : Double) -> Double {
  match fn_id {
    0 => @math.ln(1.0 + x) // log1p
    1 => if x > 0.0 { @math.ln(x) } else { 0.0 } // log
    2 => @math.log2(x) // log2
    3 => @math.pow(x, 0.5) // sqrt
    _ => x
  }
}

///|
pub fn anndata_writeback2_dummy() -> Unit {
  ignore(AnnDataWriteback::{
    source: AnnData::new(
      Array(capacity=0),
      Array(capacity=0),
      Array(capacity=0),
    ),
    obs_indices: Array(capacity=0),
    var_indices: Array(capacity=0),
    is_obs_view: false,
    is_var_view: false,
  })
}
