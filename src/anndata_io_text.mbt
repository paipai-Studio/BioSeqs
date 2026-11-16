/// anndata_io_text.mbt — Read AnnData from text/CSV files
/// Item 4: read_text / read_csv

///|
/// Read AnnData from a plain text matrix (whitespace or tab delimited).
/// First row may optionally be a header (gene names).
pub fn anndata_read_text(
  content : String,
  delimiter~ : String = "\t",
  has_header~ : Bool = true,
  first_col_names~ : Bool = true,
) -> AnnData? {
  let lines = content.trim().split("\n").to_array()
  if lines.length() == 0 {
    return None
  }
  let var_names : Array[String] = Array(capacity=0)
  let mut start_row = 0
  if has_header {
    let header_fields = lines[0].split(delimiter).to_array()
    let start = if first_col_names { 1 } else { 0 }
    for i = start; i < header_fields.length(); i = i + 1 {
      var_names.push(header_fields[i].to_owned())
    }
    start_row = 1
  }
  let x : Array[Array[Double]] = Array(capacity=lines.length() - start_row)
  let obs_names : Array[String] = Array(capacity=lines.length() - start_row)
  for r = start_row; r < lines.length(); r = r + 1 {
    let trimmed = lines[r].trim()
    if trimmed.length() == 0 {
      continue
    }
    let fields = trimmed.split(delimiter).to_array()
    let start = if first_col_names { 1 } else { 0 }
    if first_col_names && fields.length() > 0 {
      obs_names.push(fields[0].to_owned())
    }
    let row : Array[Double] = Array(capacity=fields.length() - start)
    for c = start; c < fields.length(); c = c + 1 {
      row.push(parse_double(fields[c].to_owned()).unwrap_or(0.0))
    }
    x.push(row)
  }
  // Fill in names if not provided
  if var_names.length() == 0 && x.length() > 0 {
    let n_vars = x[0].length()
    for j = 0; j < n_vars; j = j + 1 {
      var_names.push("gene_" + j.to_string())
    }
  }
  if obs_names.length() == 0 {
    for i = 0; i < x.length(); i = i + 1 {
      obs_names.push("cell_" + i.to_string())
    }
  }
  Some(AnnData::new(x, obs_names, var_names))
}

///|
/// Read AnnData from CSV (comma-delimited).
pub fn anndata_read_csv_text(
  content : String,
  has_header~ : Bool = true,
  first_col_names~ : Bool = true,
) -> AnnData? {
  anndata_read_text(content, delimiter=",", has_header~, first_col_names~)
}

///|
/// Write AnnData X matrix to plain text (tab-delimited).
pub fn anndata_write_text(adata : AnnData) -> String {
  let sb = StringBuilder()
  // Header: var_names
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

// ============================================================================
// Matrix Market (.mtx) read/write — AnnData 官方 API 风格入口
// ============================================================================

///|
/// Read a Matrix Market file into an `AnnData`.
///
/// 调用 scanpy_parse_mtx 解析 MTX 文本 → dense matrix → AnnData。
/// MTX 格式是 10x Genomics / Cell Ranger 官方输出格式。
///
/// 成功返回 `Some(adata)`，解析失败返回 `None`。
pub fn anndata_read_mtx(content : String) -> AnnData? {
  if content.trim().length() == 0 { return None }
  let (matrix, n_rows, n_cols) = scanpy_parse_mtx(content)
  if n_rows == 0 || n_cols == 0 { return None }
  // 构造 obs / var 名字（默认 "Cell_1", "Gene_1" ...）
  let obs_names : Array[String] = Array(capacity=n_rows)
  for i = 1; i <= n_rows; i = i + 1 {
    obs_names.push("Cell_" + i.to_string())
  }
  let var_names : Array[String] = Array(capacity=n_cols)
  for j = 1; j <= n_cols; j = j + 1 {
    var_names.push("Gene_" + j.to_string())
  }
  Some(AnnData::new(matrix, obs_names, var_names))
}

///|
/// Write an `AnnData` to Matrix Market coordinate format.
///
/// 输出标准 MTX 头 + 所有非零元素（1-based row/col index）。
/// dtype = real, storage = general。
pub fn anndata_write_mtx(adata : AnnData) -> String {
  let x = adata.x_matrix()
  let sb = StringBuilder()
  let _ = sb.write("%%MatrixMarket matrix coordinate real general\n")
  let _ = sb.write("%\n")
  let n_rows = x.length()
  let n_cols = if n_rows > 0 { x[0].length() } else { 0 }
  // 数非零元素
  let mut nnz : Int32 = 0
  for row in x {
    for v in row {
      if v != 0.0 { nnz = nnz + 1 }
    }
  }
  let _ = sb.write(n_rows.to_string())
  let _ = sb.write(" ")
  let _ = sb.write(n_cols.to_string())
  let _ = sb.write(" ")
  let _ = sb.write(nnz.to_string())
  let _ = sb.write("\n")
  // 非零元素（MTX 是 1-based）
  for i = 0; i < n_rows; i = i + 1 {
    let row = x[i]
    for j = 0; j < n_cols; j = j + 1 {
      let v = if j < row.length() { row[j] } else { 0.0 }
      if v != 0.0 {
        let _ = sb.write((i + 1).to_string())
        let _ = sb.write(" ")
        let _ = sb.write((j + 1).to_string())
        let _ = sb.write(" ")
        let _ = sb.write(v.to_string())
        let _ = sb.write("\n")
      }
    }
  }
  sb.to_string()
}
