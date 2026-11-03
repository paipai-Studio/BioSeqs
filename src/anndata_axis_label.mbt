/// AnnData axis label management.
///
/// `adata.obs_names.name = "cell_id"` and
/// `adata.var_names.name = "gene_id"` equivalents.
///
/// Manages axis labels and index names.

///|
/// Axis label info.
pub struct AnnDataAxisLabel {
  mut obs_axis_name : String // name of obs index (e.g. "cell_id")
  mut var_axis_name : String // name of var index (e.g. "gene_id")
  obs_col_names : Map[String, String] // obs column -> display name
  var_col_names : Map[String, String] // var column -> display name
}

///|
/// Create default axis labels.
pub fn anndata_axis_labels() -> AnnDataAxisLabel {
  {
    obs_axis_name: "obs",
    var_axis_name: "var",
    obs_col_names: Map([], capacity=4),
    var_col_names: Map([], capacity=4),
  }
}

///|
/// Set the obs axis name (index name).
pub fn AnnDataAxisLabel::set_obs_axis_name(
  self : AnnDataAxisLabel,
  name : String,
) -> Unit {
  self.obs_axis_name = name
}

///|
/// Set the var axis name (index name).
pub fn AnnDataAxisLabel::set_var_axis_name(
  self : AnnDataAxisLabel,
  name : String,
) -> Unit {
  self.var_axis_name = name
}

///|
/// Get the obs axis name.
pub fn AnnDataAxisLabel::obs_axis_name(self : AnnDataAxisLabel) -> String {
  self.obs_axis_name
}

///|
/// Get the var axis name.
pub fn AnnDataAxisLabel::var_axis_name(self : AnnDataAxisLabel) -> String {
  self.var_axis_name
}

///|
/// Set a display name for an obs column.
pub fn AnnDataAxisLabel::set_obs_col_display(
  self : AnnDataAxisLabel,
  col : String,
  display : String,
) -> Unit {
  self.obs_col_names[col] = display
}

///|
/// Set a display name for a var column.
pub fn AnnDataAxisLabel::set_var_col_display(
  self : AnnDataAxisLabel,
  col : String,
  display : String,
) -> Unit {
  self.var_col_names[col] = display
}

///|
/// Get display name for an obs column (falls back to col name).
pub fn AnnDataAxisLabel::obs_col_display(
  self : AnnDataAxisLabel,
  col : String,
) -> String {
  match self.obs_col_names.get(col) {
    Some(d) => d
    None => col
  }
}

///|
/// Get display name for a var column.
pub fn AnnDataAxisLabel::var_col_display(
  self : AnnDataAxisLabel,
  col : String,
) -> String {
  match self.var_col_names.get(col) {
    Some(d) => d
    None => col
  }
}

///|
/// Remove an obs column display name.
pub fn AnnDataAxisLabel::remove_obs_col_display(
  self : AnnDataAxisLabel,
  col : String,
) -> Unit {
  self.obs_col_names.remove(col)
}

///|
/// Remove a var column display name.
pub fn AnnDataAxisLabel::remove_var_col_display(
  self : AnnDataAxisLabel,
  col : String,
) -> Unit {
  self.var_col_names.remove(col)
}

///|
/// Generate a summary string of axis labels.
pub fn AnnDataAxisLabel::summary(self : AnnDataAxisLabel) -> String {
  let mut out = "Axis Labels:\n"
  out = out + "  obs axis: '" + self.obs_axis_name + "'\n"
  out = out + "  var axis: '" + self.var_axis_name + "'\n"
  let obs_keys = self.obs_col_names.keys().collect()
  if obs_keys.length() > 0 {
    out = out + "  obs column displays:\n"
    for k in obs_keys {
      out = out + "    " + k + " -> " + self.obs_col_display(k) + "\n"
    }
  }
  let var_keys = self.var_col_names.keys().collect()
  if var_keys.length() > 0 {
    out = out + "  var column displays:\n"
    for k in var_keys {
      out = out + "    " + k + " -> " + self.var_col_display(k) + "\n"
    }
  }
  out
}

///|
/// Copy axis labels.
pub fn AnnDataAxisLabel::copy(self : AnnDataAxisLabel) -> AnnDataAxisLabel {
  let result = anndata_axis_labels()
  result.set_obs_axis_name(self.obs_axis_name)
  result.set_var_axis_name(self.var_axis_name)
  let obs_keys = self.obs_col_names.keys().collect()
  for k in obs_keys {
    result.set_obs_col_display(k, self.obs_col_display(k))
  }
  let var_keys = self.var_col_names.keys().collect()
  for k in var_keys {
    result.set_var_col_display(k, self.var_col_display(k))
  }
  result
}

///|
/// Transfer axis labels from one AnnData context to another.
/// Useful after operations that preserve column structure.
pub fn anndata_transfer_axis_labels(
  src : AnnDataAxisLabel,
  dst : AnnDataAxisLabel,
) -> Unit {
  dst.set_obs_axis_name(src.obs_axis_name)
  dst.set_var_axis_name(src.var_axis_name)
  let obs_keys = src.obs_col_names.keys().collect()
  for k in obs_keys {
    dst.set_obs_col_display(k, src.obs_col_display(k))
  }
  let var_keys = src.var_col_names.keys().collect()
  for k in var_keys {
    dst.set_var_col_display(k, src.var_col_display(k))
  }
}
