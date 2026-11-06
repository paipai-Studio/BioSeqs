/// Scanpy-inspired: read_mtx — read Matrix Market format
/// Parses .mtx files (sparse matrix format).

///|
/// Read a Matrix Market file and return a dense matrix.
/// MTX format:
///   %%MatrixMarket matrix coordinate real general
///   n_rows n_cols n_nonzero
///   row col value
///   ...
///
/// Parameters:
/// - `content`: MTX file content as string
/// - `transpose` : if true, transpose the result (default false)
pub fn scanpy_read_mtx(
  content : String,
  transpose~ : Bool = false,
) -> Array[Array[Double]] {
  let lines = content.split("\n").to_array()
  if lines.length() == 0 {
    return []
  }
  let mut idx = 0
  let mut is_sparse = true
  let mut is_symmetric = false
  for line in lines {
    let trimmed = line.to_owned().trim().to_owned()
    if trimmed.length() == 0 {
      idx += 1
      continue
    }
    if trimmed[0:1].to_owned() == "%" {
      if trimmed.has_prefix("%%MatrixMarket") {
        if trimmed.contains("coordinate") {
          is_sparse = true
        } else if trimmed.contains("array") {
          is_sparse = false
        }
        if trimmed.contains("symmetric") {
          is_symmetric = true
        }
      }
      idx += 1
      continue
    }
    break
  }
  if idx >= lines.length() {
    return []
  }
  let dim_parts = lines[idx].to_owned().trim().to_owned().split(" ").to_array()
  if dim_parts.length() < 2 {
    return []
  }
  let n_rows = parse_int(dim_parts[0].to_owned())
  let n_cols = parse_int(dim_parts[1].to_owned())
  idx += 1
  let x : Array[Array[Double]] = []
  for _ in 0..<n_rows {
    let row : Array[Double] = []
    for _ in 0..<n_cols {
      row.push(0.0)
    }
    x.push(row)
  }
  if is_sparse {
    for line_idx in idx..<lines.length() {
      let trimmed = lines[line_idx].to_owned().trim().to_owned()
      if trimmed.length() == 0 {
        continue
      }
      let parts = trimmed.split(" ").to_array()
      if parts.length() >= 3 {
        let row = parse_int(parts[0].to_owned()) - 1
        let col = parse_int(parts[1].to_owned()) - 1
        let val = match parse_double(parts[2].to_owned()) {
          Some(v) => v
          None => 0.0
        }
        if row >= 0 && row < n_rows && col >= 0 && col < n_cols {
          x[row][col] = val
          if is_symmetric && row != col {
            x[col][row] = val
          }
        }
      }
    }
  } else {
    let mut col = 0
    let mut row = 0
    for line_idx in idx..<lines.length() {
      let trimmed = lines[line_idx].to_owned().trim().to_owned()
      if trimmed.length() == 0 {
        continue
      }
      let val = match parse_double(trimmed) {
        Some(v) => v
        None => 0.0
      }
      if row < n_rows && col < n_cols {
        x[row][col] = val
      }
      row += 1
      if row >= n_rows {
        row = 0
        col += 1
      }
    }
  }
  if transpose {
    let transposed : Array[Array[Double]] = []
    for j in 0..<n_cols {
      let row : Array[Double] = []
      for i in 0..<n_rows {
        row.push(x[i][j])
      }
      transposed.push(row)
    }
    return transposed
  }
  x
}

///|
/// Read MTX along with barcodes and features files (10x-style).
///
/// Parameters:
/// - `mtx_content`: MTX file content
/// - `barcodes`: cell barcodes (one per line)
/// - `features`: gene names (one per line)
pub fn scanpy_read_mtx_with_names(
  mtx_content : String,
  barcodes : Array[String],
  features : Array[String],
) -> AnnData {
  let x = scanpy_read_mtx(mtx_content, transpose=true)
  AnnData::new(x, barcodes, features)
}
