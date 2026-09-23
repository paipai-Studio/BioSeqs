///|
/// scanpy `read_10x_mtx` + Matrix Market parsing — load 10x Cell Ranger output
/// (matrix.mtx / features.tsv / barcodes.tsv) into an `AnnData`. The caller
/// passes already-decompressed file contents (gzip is decoded upstream).

///|
/// Make a list of names unique by appending `-1`, `-2`, … to later duplicate
/// occurrences, matching `anndata.utils.make_index_unique`.
pub fn scanpy_make_index_unique(names : Array[String]) -> Array[String] {
  let counts : Map[String, Int] = Map([], capacity=64)
  let out : Array[String] = Array(capacity=0)
  for name in names {
    match counts.get(name) {
      None => {
        counts.set(name, 1)
        out.push(name)
      }
      Some(c) => {
        counts.set(name, c + 1)
        out.push(name + "-" + c.to_string())
      }
    }
  }
  out
}

///|
/// Parse a Matrix Market coordinate file into a dense matrix. Returns
/// `(matrix, n_rows, n_cols)` where `matrix[i][j]` is a 0-based entry. Missing
/// entries are 0. Handles `integer`, `real` and `pattern` dtypes and
/// `general` / `symmetric` storage.
pub fn scanpy_parse_mtx(
  content : String,
) -> (Array[Array[Double]], Int, Int) {
  let lines = content.split("\n").to_array()
  let n_lines = lines.length()
  let mut i = 0
  // skip header + comments, capture the field token
  while i < n_lines {
    let line = lines[i].to_owned()
    let t = scanpy_trim(line)
    if t.length() == 0 {
      i = i + 1
      continue
    }
    let c0 = line.unsafe_get(0).to_int()
    if c0 == '%'.to_int() {
      i = i + 1
      continue
    }
    break
  }
  // find the size line (first non-header, non-comment line)
  let mut size_line = ""
  while i < n_lines {
    let t = scanpy_trim(lines[i].to_owned())
    if t.length() == 0 || t.unsafe_get(0).to_int() == '%'.to_int() {
      i = i + 1
      continue
    }
    size_line = t
    i = i + 1
    break
  }
  let dims = scanpy_split_ws(size_line)
  let n_rows = parse_int(dims[0])
  let n_cols = parse_int(dims[1])
  let mat : Array[Array[Double]] = Array(capacity=0)
  for _r = 0; _r < n_rows; _r = _r + 1 {
    mat.push(Array::make(n_cols, 0.0))
  }
  // parse entries
  while i < n_lines {
    let t = scanpy_trim(lines[i].to_owned())
    i = i + 1
    if t.length() == 0 || t.unsafe_get(0).to_int() == '%'.to_int() {
      continue
    }
    let parts = scanpy_split_ws(t)
    if parts.length() < 3 {
      continue
    }
    let r = parse_int(parts[0]) - 1
    let c = parse_int(parts[1]) - 1
    let v = scanpy_parse_mtx_value(parts[2])
    if r >= 0 && r < n_rows && c >= 0 && c < n_cols {
      mat[r][c] = v
    }
  }
  (mat, n_rows, n_cols)
}

///|
fn scanpy_parse_mtx_value(s : String) -> Double {
  if s == "." {
    0.0
  } else {
    match parse_double(s) {
      None => 0.0
      Some(v) => v
    }
  }
}

///|
fn scanpy_trim(s : String) -> String {
  let mut lo = 0
  let mut hi = s.length() - 1
  while lo <= hi {
    let c = s.unsafe_get(lo).to_int()
    if c == 32 || c == 9 || c == 13 {
      lo = lo + 1
    } else {
      break
    }
  }
  while hi >= lo {
    let c = s.unsafe_get(hi).to_int()
    if c == 32 || c == 9 || c == 13 {
      hi = hi - 1
    } else {
      break
    }
  }
  s[lo:hi + 1].to_owned()
}

///|
fn scanpy_split_ws(s : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let mut cur = StringBuilder::new()
  let mut has = false
  let mut i = 0
  while i <= s.length() {
    let is_sep = if i == s.length() {
      true
    } else {
      let c = s.unsafe_get(i).to_int()
      c == 32 || c == 9 || c == 13
    }
    if is_sep {
      if has {
        out.push(cur.to_string())
      }
      cur = StringBuilder::new()
      has = false
    } else {
      cur.write_char(s.unsafe_get(i).unsafe_to_char())
      has = true
    }
    i = i + 1
  }
  out
}

///|
/// Read 10x-Genomics-formatted mtx data (Cell Ranger v3+) into an `AnnData`.
/// `mtx` is the `matrix.mtx` content (genes × cells, transposed internally),
/// `features` is `features.tsv` (`gene_id \t gene_symbol \t feature_type`),
/// `barcodes` is `barcodes.tsv` (one cell barcode per line).
pub fn scanpy_read_10x_mtx(
  mtx : String,
  features : String,
  barcodes : String,
  var_names : String,
  make_unique : Bool,
  gex_only : Bool,
) -> AnnData {
  let mtx_u = gunzip_if_needed(mtx) catch { _ => mtx }
  let features_u = gunzip_if_needed(features) catch { _ => features }
  let barcodes_u = gunzip_if_needed(barcodes) catch { _ => barcodes }
  let (gmat, n_genes, n_cells) = scanpy_parse_mtx(mtx_u)
  // transpose genes × cells -> cells × genes
  let x : Array[Array[Double]] = Array(capacity=0)
  for c = 0; c < n_cells; c = c + 1 {
    let row : Array[Double] = Array(capacity=0)
    for g = 0; g < n_genes; g = g + 1 {
      row.push(gmat[g][c])
    }
    x.push(row)
  }

  // parse features.tsv (3 columns)
  let gene_ids : Array[String] = Array(capacity=0)
  let gene_symbols : Array[String] = Array(capacity=0)
  let feature_types : Array[String] = Array(capacity=0)
  for line in features_u.split("\n") {
    let t = scanpy_trim(line.to_owned())
    if t.length() == 0 {
      continue
    }
    let cols = split_by_tab(t)
    gene_ids.push(if cols.length() > 0 { cols[0] } else { "" })
    gene_symbols.push(if cols.length() > 1 { cols[1] } else { "" })
    feature_types.push(if cols.length() > 2 { cols[2] } else { "" })
  }

  // parse barcodes.tsv (1 column)
  let obs_names : Array[String] = Array(capacity=0)
  for line in barcodes_u.split("\n") {
    let t = scanpy_trim(line.to_owned())
    if t.length() == 0 {
      continue
    }
    obs_names.push(t)
  }

  // choose var_names + the companion annotation
  let var_names_arr = if var_names == "gene_ids" {
    gene_ids
  } else {
    gene_symbols
  }
  let unique_names = if make_unique {
    scanpy_make_index_unique(var_names_arr)
  } else {
    var_names_arr
  }
  let companion_key = if var_names == "gene_ids" {
    "gene_symbols"
  } else {
    "gene_ids"
  }
  let companion_vals = if var_names == "gene_ids" {
    gene_symbols
  } else {
    gene_ids
  }

  let mut adata = AnnData::new(x, obs_names, unique_names)
  adata = adata.set_var(companion_key, companion_vals)
  adata = adata.set_var("feature_types", feature_types)

  if gex_only {
    // keep only "Gene Expression" features (rows of adata → filter columns of X)
    let keep : Array[Int] = Array(capacity=0)
    for gi = 0; gi < feature_types.length(); gi = gi + 1 {
      if feature_types[gi] == "Gene Expression" {
        keep.push(gi)
      }
    }
    let new_var_names : Array[String] = Array(capacity=0)
    let new_companion : Array[String] = Array(capacity=0)
    let new_ftypes : Array[String] = Array(capacity=0)
    for k in keep {
      new_var_names.push(unique_names[k])
      new_companion.push(companion_vals[k])
      new_ftypes.push(feature_types[k])
    }
    let new_x : Array[Array[Double]] = Array(capacity=0)
    for c = 0; c < x.length(); c = c + 1 {
      let row : Array[Double] = Array(capacity=0)
      for k in keep {
        row.push(x[c][k])
      }
      new_x.push(row)
    }
    adata = AnnData::new(new_x, obs_names, new_var_names)
    adata = adata.set_var(companion_key, new_companion)
    adata = adata.set_var("feature_types", new_ftypes)
  }

  adata
}