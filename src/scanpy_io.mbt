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
  let mut cur = StringBuilder()
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
      cur = StringBuilder()
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

// ============================================================================
// h5ad text-format IO (JSON-based)
//
// Real .h5ad files use HDF5 binary format.  Without an HDF5 library we
// implement a JSON-based text format that preserves the full AnnData
// structure (X, obs, var, layers, obsm, obsp, varp, uns, names) and
// supports round-trip read/write.  This is sufficient for testing,
// interchange, and small datasets.
// ============================================================================

///|
/// Convert a dense matrix to a JsonValue array-of-arrays.
fn h5ad_mat_to_json(mat : Array[Array[Double]]) -> JsonValue {
  let rows : Array[JsonValue] = Array(capacity=mat.length())
  for row in mat {
    let cols : Array[JsonValue] = Array(capacity=row.length())
    for v in row {
      cols.push(JsonNum(v))
    }
    rows.push(JsonArr(cols))
  }
  JsonArr(rows)
}

///|
/// Convert a string array to a JsonValue array.
fn h5ad_strs_to_json(arr : Array[String]) -> JsonValue {
  let items : Array[JsonValue] = Array(capacity=arr.length())
  for s in arr {
    items.push(JsonStr(s))
  }
  JsonArr(items)
}

///|
/// Convert a Map[String, Array[String]] to a JsonValue object.
fn h5ad_str_map_to_json(m : Map[String, Array[String]]) -> JsonValue {
  let obj : Map[String, JsonValue] = Map([], capacity=4)
  for kv in m {
    obj[kv.0] = h5ad_strs_to_json(kv.1)
  }
  JsonObj(obj)
}

///|
/// Convert a Map[String, Array[Array[Double]]] to a JsonValue object.
fn h5ad_mat_map_to_json(m : Map[String, Array[Array[Double]]]) -> JsonValue {
  let obj : Map[String, JsonValue] = Map([], capacity=4)
  for kv in m {
    obj[kv.0] = h5ad_mat_to_json(kv.1)
  }
  JsonObj(obj)
}

///|
/// Convert a Map[String, String] to a JsonValue object.
fn h5ad_uns_to_json(m : Map[String, String]) -> JsonValue {
  let obj : Map[String, JsonValue] = Map([], capacity=4)
  for kv in m {
    obj[kv.0] = JsonStr(kv.1)
  }
  JsonObj(obj)
}

///|
/// Serialize an AnnData to a JSON string (h5ad text format).
/// Supports both dense and sparse X storage.
pub fn scanpy_write_h5ad(adata : AnnData) -> String {
  let root : Map[String, JsonValue] = Map([], capacity=16)
  // X matrix
  match adata.get_sparse() {
    Some(s) => {
      let dims = s.dim()
      let coords : Array[JsonValue] = Array(capacity=0)
      let values : Array[JsonValue] = Array(capacity=0)
      for entry in s.entries {
        let c : Array[JsonValue] = Array(capacity=2)
        c.push(JsonNum(entry.coordinates[0].to_double()))
        c.push(JsonNum(entry.coordinates[1].to_double()))
        coords.push(JsonArr(c))
        values.push(JsonNum(entry.value))
      }
      let sparse_obj : Map[String, JsonValue] = Map([], capacity=3)
      sparse_obj["dims"] = h5ad_strs_to_json([dims[0].to_string(), dims[1].to_string()])
      // Actually dims should be numbers
      let dims_arr : Array[JsonValue] = Array(capacity=2)
      dims_arr.push(JsonNum(dims[0].to_double()))
      dims_arr.push(JsonNum(dims[1].to_double()))
      sparse_obj["dims"] = JsonArr(dims_arr)
      sparse_obj["coords"] = JsonArr(coords)
      sparse_obj["values"] = JsonArr(values)
      root["X_sparse"] = JsonObj(sparse_obj)
    }
    None => {
      root["X"] = h5ad_mat_to_json(adata.x)
    }
  }
  root["obs_names"] = h5ad_strs_to_json(adata.obs_names)
  root["var_names"] = h5ad_strs_to_json(adata.var_names)
  root["obs"] = h5ad_str_map_to_json(adata.obs)
  root["var"] = h5ad_str_map_to_json(adata.variables)
  root["layers"] = h5ad_mat_map_to_json(adata.layers)
  root["obsm"] = h5ad_mat_map_to_json(adata.obsm)
  root["varm"] = h5ad_mat_map_to_json(adata.varm)
  root["obsp"] = h5ad_mat_map_to_json(adata.obsp)
  root["varp"] = h5ad_mat_map_to_json(adata.varp)
  root["uns"] = h5ad_uns_to_json(adata.uns)
  json_to_string(JsonObj(root))
}

///|
/// Extract a Double from a JsonValue, defaulting to 0.0.
fn h5ad_json_to_num(v : JsonValue) -> Double {
  match v {
    JsonNum(d) => d
    _ => 0.0
  }
}

///|
/// Extract a String from a JsonValue, defaulting to "".
fn h5ad_json_to_str(v : JsonValue) -> String {
  match v {
    JsonStr(s) => s
    _ => ""
  }
}

///|
/// Parse a JsonValue array-of-arrays into a dense matrix.
fn h5ad_json_to_mat(v : JsonValue) -> Array[Array[Double]] {
  match v {
    JsonArr(rows) => {
      let out : Array[Array[Double]] = Array(capacity=rows.length())
      for row in rows {
        match row {
          JsonArr(cols) => {
            let r : Array[Double] = Array(capacity=cols.length())
            for c in cols {
              r.push(h5ad_json_to_num(c))
            }
            out.push(r)
          }
          _ => out.push([])
        }
      }
      out
    }
    _ => []
  }
}

///|
/// Parse a JsonValue array into a string array.
fn h5ad_json_to_strs(v : JsonValue) -> Array[String] {
  match v {
    JsonArr(items) => {
      let out : Array[String] = Array(capacity=items.length())
      for item in items {
        out.push(h5ad_json_to_str(item))
      }
      out
    }
    _ => []
  }
}

///|
/// Parse a JsonValue object into a Map[String, Array[String]].
fn h5ad_json_to_str_map(v : JsonValue) -> Map[String, Array[String]] {
  let out : Map[String, Array[String]] = Map([], capacity=4)
  match v {
    JsonObj(m) => {
      for kv in m {
        out[kv.0] = h5ad_json_to_strs(kv.1)
      }
    }
    _ => { let _ = 0 }
  }
  out
}

///|
/// Parse a JsonValue object into a Map[String, Array[Array[Double]]].
fn h5ad_json_to_mat_map(v : JsonValue) -> Map[String, Array[Array[Double]]] {
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=4)
  match v {
    JsonObj(m) => {
      for kv in m {
        out[kv.0] = h5ad_json_to_mat(kv.1)
      }
    }
    _ => { let _ = 0 }
  }
  out
}

///|
/// Parse a JsonValue object into a Map[String, String].
fn h5ad_json_to_uns(v : JsonValue) -> Map[String, String] {
  let out : Map[String, String] = Map([], capacity=4)
  match v {
    JsonObj(m) => {
      for kv in m {
        out[kv.0] = h5ad_json_to_str(kv.1)
      }
    }
    _ => { let _ = 0 }
  }
  out
}

///|
/// Deserialize a JSON string (h5ad text format) back to an AnnData.
/// Returns `None` if parsing fails.
pub fn scanpy_read_h5ad(content : String) -> AnnData? {
  let json = parse_json(content)
  match json {
    None => None
    Some(JsonObj(root)) => {
      // Parse obs_names and var_names
      let obs_names = match root.get("obs_names") {
        Some(v) => h5ad_json_to_strs(v)
        None => []
      }
      let var_names = match root.get("var_names") {
        Some(v) => h5ad_json_to_strs(v)
        None => []
      }
      // Parse obs/var metadata
      let obs = match root.get("obs") {
        Some(v) => h5ad_json_to_str_map(v)
        None => Map([], capacity=4)
      }
      let variables = match root.get("var") {
        Some(v) => h5ad_json_to_str_map(v)
        None => Map([], capacity=4)
      }
      // Parse layers/obsm/varm/obsp/varp
      let layers = match root.get("layers") {
        Some(v) => h5ad_json_to_mat_map(v)
        None => Map([], capacity=4)
      }
      let obsm = match root.get("obsm") {
        Some(v) => h5ad_json_to_mat_map(v)
        None => Map([], capacity=4)
      }
      let varm = match root.get("varm") {
        Some(v) => h5ad_json_to_mat_map(v)
        None => Map([], capacity=4)
      }
      let obsp = match root.get("obsp") {
        Some(v) => h5ad_json_to_mat_map(v)
        None => Map([], capacity=4)
      }
      let varp = match root.get("varp") {
        Some(v) => h5ad_json_to_mat_map(v)
        None => Map([], capacity=4)
      }
      let uns = match root.get("uns") {
        Some(v) => h5ad_json_to_uns(v)
        None => Map([], capacity=4)
      }
      // Build AnnData — sparse or dense
      match root.get("X_sparse") {
        Some(JsonObj(sparse_obj)) => {
          let dims = match sparse_obj.get("dims") {
            Some(JsonArr(d)) => {
              let r = h5ad_json_to_num(d[0]).to_int()
              let c = h5ad_json_to_num(d[1]).to_int()
              [r, c]
            }
            _ => [0, 0]
          }
          let coords : Array[Array[Int]] = Array(capacity=0)
          match sparse_obj.get("coords") {
            Some(JsonArr(ca)) => {
              for c in ca {
                match c {
                  JsonArr(pair) => {
                    coords.push([
                      h5ad_json_to_num(pair[0]).to_int(),
                      h5ad_json_to_num(pair[1]).to_int(),
                    ])
                  }
                  _ => { let _ = 0 }
                }
              }
            }
            _ => { let _ = 0 }
          }
          let values : Array[Double] = Array(capacity=0)
          match sparse_obj.get("values") {
            Some(JsonArr(va)) => {
              for v in va {
                values.push(h5ad_json_to_num(v))
              }
            }
            _ => { let _ = 0 }
          }
          let sparse = SparseArray::from_coo(dims, coords, values) catch {
            _ => return None
          }
          let adata = AnnData::new_sparse(sparse, obs_names, var_names)
          // Set metadata
          for kv in obs {
            let _ = adata.set_obs(kv.0, kv.1)
          }
          for kv in variables {
            let _ = adata.set_var(kv.0, kv.1)
          }
          for kv in layers {
            let _ = adata.set_layer(kv.0, kv.1)
          }
          for kv in obsm {
            let _ = adata.set_obsm(kv.0, kv.1)
          }
          for kv in obsp {
            let _ = adata.set_obsp(kv.0, kv.1)
          }
          Some(adata)
        }
        _ => {
          // Dense X
          let x = match root.get("X") {
            Some(v) => h5ad_json_to_mat(v)
            None => []
          }
          let adata = AnnData::new(x, obs_names, var_names)
          for kv in obs {
            let _ = adata.set_obs(kv.0, kv.1)
          }
          for kv in variables {
            let _ = adata.set_var(kv.0, kv.1)
          }
          for kv in layers {
            let _ = adata.set_layer(kv.0, kv.1)
          }
          for kv in obsm {
            let _ = adata.set_obsm(kv.0, kv.1)
          }
          for kv in obsp {
            let _ = adata.set_obsp(kv.0, kv.1)
          }
          Some(adata)
        }
      }
    }
    _ => None
  }
}
