/// anndata_to_df.mbt — Convert AnnData to DataFrame representation
/// Item 24: to_df / X_to_df

///|
/// DataFrame representation (simplified: column-major)
pub struct AnnDataFrame {
  columns : Map[String, Array[Double]]
  row_names : Array[String]
  col_names : Array[String]
}

///|
/// Convert AnnData X matrix to DataFrame.
/// Row names = obs_names, column names = var_names.
pub fn anndata_to_df(adata : AnnData, layer~ : String? = None) -> AnnDataFrame {
  let x = if layer is Some(key) && adata.layers.contains(key) {
    adata.layers[key]
  } else {
    adata.x_matrix()
  }
  let columns : Map[String, Array[Double]] = Map([], capacity=adata.n_vars())
  for j = 0; j < adata.n_vars(); j = j + 1 {
    let col : Array[Double] = Array(capacity=x.length())
    for i = 0; i < x.length(); i = i + 1 {
      if j < x[i].length() {
        col.push(x[i][j])
      } else {
        col.push(0.0)
      }
    }
    columns[adata.var_names[j]] = col
  }
  AnnDataFrame::{
    columns,
    row_names: adata.obs_names,
    col_names: adata.var_names,
  }
}

///|
/// Convert X to a simple string table (TSV).
pub fn anndata_x_to_df(adata : AnnData) -> String {
  let sb = StringBuilder()
  // Header
  let _ = sb.write("")
  for n in adata.var_names {
    let _ = sb.write("\t")
    let _ = sb.write(n)
  }
  let _ = sb.write("\n")
  // Rows
  let x = adata.x_matrix()
  for i = 0; i < x.length(); i = i + 1 {
    let _ = sb.write(adata.obs_names[i])
    for j = 0; j < x[i].length(); j = j + 1 {
      let _ = sb.write("\t")
      let _ = sb.write(x[i][j].to_string())
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
/// Get a column from the DataFrame.
pub fn AnnDataFrame::get_column(
  self : AnnDataFrame,
  name : String,
) -> Array[Double]? {
  if self.columns.contains(name) {
    Some(self.columns[name])
  } else {
    None
  }
}

///|
/// Get a row from the DataFrame by index.
pub fn AnnDataFrame::get_row(
  self : AnnDataFrame,
  idx : Int32,
) -> Array[Double]? {
  if idx < 0 || idx >= self.row_names.length() {
    return None
  }
  let row : Array[Double] = Array(capacity=self.col_names.length())
  for name in self.col_names {
    if self.columns.contains(name) {
      let col = self.columns[name]
      if idx < col.length() {
        row.push(col[idx])
      } else {
        row.push(0.0)
      }
    } else {
      row.push(0.0)
    }
  }
  Some(row)
}

///|
/// Get a row by name.
pub fn AnnDataFrame::get_row_by_name(
  self : AnnDataFrame,
  name : String,
) -> Array[Double]? {
  let mut idx = -1
  for i = 0; i < self.row_names.length(); i = i + 1 {
    if self.row_names[i] == name {
      idx = i
      break
    }
  }
  if idx >= 0 {
    self.get_row(idx)
  } else {
    None
  }
}

///|
/// Number of rows.
pub fn AnnDataFrame::n_rows(self : AnnDataFrame) -> Int32 {
  self.row_names.length()
}

///|
/// Number of columns.
pub fn AnnDataFrame::n_cols(self : AnnDataFrame) -> Int32 {
  self.col_names.length()
}

///|
/// Convert DataFrame to string representation.
pub fn AnnDataFrame::to_string(self : AnnDataFrame) -> String {
  let sb = StringBuilder()
  for name in self.col_names {
    let _ = sb.write(name)
    let _ = sb.write("\t")
  }
  let _ = sb.write("\n")
  for i = 0; i < self.row_names.length(); i = i + 1 {
    let _ = sb.write(self.row_names[i])
    let _ = sb.write("\t")
    for name in self.col_names {
      if self.columns.contains(name) {
        let col = self.columns[name]
        if i < col.length() {
          let _ = sb.write(col[i].to_string())
        }
      }
      let _ = sb.write("\t")
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
pub fn anndata_to_df_dummy() -> Unit {
  ignore(AnnDataFrame::{
    columns: Map([], capacity=0),
    row_names: Array(capacity=0),
    col_names: Array(capacity=0),
  })
}
