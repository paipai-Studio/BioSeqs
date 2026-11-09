/// AnnData views vs copies — inplace vs copy semantics.
///
/// In Python AnnData, `adata[mask]` returns a view (lazy).
/// `adata[mask].copy()` materializes it.
/// This module provides explicit view/copy control.

///|
/// View mode: which axis is subsetted.
pub enum AnnDataViewMode {
  ViewObs // subset along obs axis
  ViewVar // subset along var axis
  ViewBoth // subset along both axes
  ViewNone // no subsetting (full view)
}

///|
/// A lazy view descriptor. Does not copy data; stores indices.
pub struct AnnDataView {
  source : AnnData // reference to source
  obs_idx : Array[Int32] // obs indices (or all)
  var_idx : Array[Int32] // var indices (or all)
  mode : AnnDataViewMode
  mut modified : Bool // set once X is mutated through the view
}

///|
/// Create a view of an AnnData (no data copied).
pub fn anndata_view(adata : AnnData) -> AnnDataView {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let obs_idx : Array[Int32] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    obs_idx.push(i)
  }
  let var_idx : Array[Int32] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    var_idx.push(j)
  }
  { source: adata, obs_idx, var_idx, mode: ViewNone, modified: false, }
}

///|
/// Create a view with obs subsetting.
pub fn anndata_view_obs(adata : AnnData, idx : Array[Int32]) -> AnnDataView {
  let n_vars = adata.n_vars()
  let var_idx : Array[Int32] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    var_idx.push(j)
  }
  { source: adata, obs_idx: idx, var_idx, mode: ViewObs, modified: false, }
}

///|
/// Create a view with var subsetting.
pub fn anndata_view_var(adata : AnnData, idx : Array[Int32]) -> AnnDataView {
  let n_obs = adata.n_obs()
  let obs_idx : Array[Int32] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    obs_idx.push(i)
  }
  { source: adata, obs_idx, var_idx: idx, mode: ViewVar, modified: false, }
}

///|
/// Create a view with both axes subsetted.
pub fn anndata_view_both(
  adata : AnnData,
  obs_idx : Array[Int32],
  var_idx : Array[Int32],
) -> AnnDataView {
  { source: adata, obs_idx, var_idx, mode: ViewBoth, modified: false, }
}

///|
/// Number of obs in the view.
pub fn AnnDataView::n_obs(self : AnnDataView) -> Int32 {
  self.obs_idx.length()
}

///|
/// Number of var in the view.
pub fn AnnDataView::n_vars(self : AnnDataView) -> Int32 {
  self.var_idx.length()
}

///|
/// Get element (i, j) from the view.
pub fn AnnDataView::get_x(self : AnnDataView, i : Int32, j : Int32) -> Double {
  if i < self.obs_idx.length() && j < self.var_idx.length() {
    self.source.get_x(self.obs_idx[i], self.var_idx[j])
  } else {
    0.0
  }
}

///|
/// Materialize the view into a new AnnData (copy data).
pub fn AnnDataView::to_adata(self : AnnDataView) -> AnnData {
  self.source.subset(self.obs_idx, self.var_idx)
}

///|
/// Copy alias (same as to_adata).
pub fn AnnDataView::copy(self : AnnDataView) -> AnnData {
  self.to_adata()
}

///|
/// Get obs_names from the view (without materializing).
pub fn AnnDataView::obs_names(self : AnnDataView) -> Array[String] {
  let src_names = self.source.get_obs_names()
  let out : Array[String] = Array(capacity=self.obs_idx.length())
  for i in self.obs_idx {
    if i >= 0 && i < src_names.length() {
      out.push(src_names[i])
    } else {
      out.push("")
    }
  }
  out
}

///|
/// Get var_names from the view (without materializing).
pub fn AnnDataView::var_names(self : AnnDataView) -> Array[String] {
  let src_names = self.source.get_var_names()
  let out : Array[String] = Array(capacity=self.var_idx.length())
  for j in self.var_idx {
    if j >= 0 && j < src_names.length() {
      out.push(src_names[j])
    } else {
      out.push("")
    }
  }
  out
}

///|
/// Further subset an existing view.
pub fn AnnDataView::subset(
  self : AnnDataView,
  new_obs_idx : Array[Int32],
  new_var_idx : Array[Int32],
) -> AnnDataView {
  // Compose indices: map through the view's index arrays
  let composed_obs : Array[Int32] = Array(capacity=new_obs_idx.length())
  for i in new_obs_idx {
    if i >= 0 && i < self.obs_idx.length() {
      composed_obs.push(self.obs_idx[i])
    } else {
      composed_obs.push(-1)
    }
  }
  let composed_var : Array[Int32] = Array(capacity=new_var_idx.length())
  for j in new_var_idx {
    if j >= 0 && j < self.var_idx.length() {
      composed_var.push(self.var_idx[j])
    } else {
      composed_var.push(-1)
    }
  }
  // Filter out -1
  let valid_obs : Array[Int32] = Array(capacity=0)
  for i in composed_obs {
    if i >= 0 {
      valid_obs.push(i)
    }
  }
  let valid_var : Array[Int32] = Array(capacity=0)
  for j in composed_var {
    if j >= 0 {
      valid_var.push(j)
    }
  }
  {
    source: self.source,
    obs_idx: valid_obs,
    var_idx: valid_var,
    mode: ViewBoth,
    modified: false,
  }
}

///|
/// Check if the view is the full data (no subsetting).
pub fn AnnDataView::is_full(self : AnnDataView) -> Bool {
  match self.mode {
    ViewNone => true
    _ => false
  }
}

///|
/// Inplace modify X from a view's values.
/// Writes the view's data back to the source AnnData.
pub fn AnnDataView::writeback(
  self : AnnDataView,
  new_x : Array[Array[Double]],
) -> Unit {
  for i = 0; i < self.obs_idx.length() && i < new_x.length(); i = i + 1 {
    for j = 0; j < self.var_idx.length() && j < new_x[i].length(); j = j + 1 {
      let oi = self.obs_idx[i]
      let oj = self.var_idx[j]
      if oi < self.source.x.length() && oj < self.source.x[oi].length() {
        self.source.x[oi][oj] = new_x[i][j]
      }
    }
  }
  self.modified = true
}

// === SettingWithCopyWarning semantics ===
//
// In Python AnnData, mutating a view emits
// `ImplicitModificationWarning` (pandas' SettingWithCopyWarning analogue):
// the write goes through to the parent object, but the user is reminded
// that they are modifying a copy-like slice and that `.copy()` should be
// used instead. We mirror that: writes always reach the source, the view
// records `modified = true`, and `set_x`/`warning` surface the message.

///|
/// The implicit-modification warning message.
pub fn anndata_view_warning_text() -> String {
  "ImplicitModificationWarning: a value is trying to be set on a copy of " +
  "a slice from an AnnData. Try using .copy() instead."
}

///|
/// Set a single X element through the view (inplace write to the source).
/// Returns the implicit-modification warning text when the view is a real
/// subset (not full), or `None` when writing on a full view / out of range.
pub fn AnnDataView::set_x(
  self : AnnDataView,
  i : Int32,
  j : Int32,
  value : Double,
) -> String? {
  if i < 0 || i >= self.obs_idx.length() || j < 0 || j >= self.var_idx.length() {
    return None
  }
  let oi = self.obs_idx[i]
  let oj = self.var_idx[j]
  if oi < self.source.x.length() && oj < self.source.x[oi].length() {
    self.source.x[oi][oj] = value
  }
  self.modified = true
  if self.is_full() {
    None
  } else {
    Some(anndata_view_warning_text())
  }
}

///|
/// True once X has been mutated through this view.
pub fn AnnDataView::is_modified(self : AnnDataView) -> Bool {
  self.modified
}

///|
/// True when a subsequent inplace write should warn (real subset view).
pub fn AnnDataView::needs_copy(self : AnnDataView) -> Bool {
  self.modified && not(self.is_full())
}

///|
/// The pending implicit-modification warning for this view, or `None` if
/// the view has not been mutated through (or is a full view).
pub fn AnnDataView::warning(self : AnnDataView) -> String? {
  if self.needs_copy() {
    Some(anndata_view_warning_text())
  } else {
    None
  }
}
