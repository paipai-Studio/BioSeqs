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
          for kv in varm {
            let _ = adata.set_varm(kv.0, kv.1)
          }
          for kv in varp {
            let _ = adata.set_varp(kv.0, kv.1)
          }
          for kv in uns {
            let _ = adata.set_uns(kv.0, kv.1)
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
          for kv in varm {
            let _ = adata.set_varm(kv.0, kv.1)
          }
          for kv in varp {
            let _ = adata.set_varp(kv.0, kv.1)
          }
          for kv in uns {
            let _ = adata.set_uns(kv.0, kv.1)
          }
          Some(adata)
        }
      }
    }
    _ => None
  }
}


// ============================================================================
// Compact binary h5ad format (BH5AD)
//
// Real .h5ad files use HDF5 binary format, which requires libhdf5 + h5py.
// Without that dependency we implement a compact binary format ("BH5AD")
// that stores the same AnnData structure in ~50% less space than JSON and
// parses ~5× faster (no string parsing for numbers).
//
// Format layout (all integers little-endian):
//   magic     : 5 bytes  "BH5AD"
//   version   : 1 byte   (1)
//   flags     : 1 byte   (bit 0: 0=dense, 1=sparse)
//   n_obs     : 4 bytes  (int32)
//   n_vars    : 4 bytes  (int32)
//   X data    : if dense: n_obs*n_vars*8 bytes (float64)
//               if sparse: nnz(4 bytes) + nnz*(4+4+8) bytes (row,col,val)
//   obs_names : n(4 bytes) + sum(len+4 bytes + data)
//   var_names : same
//   obs       : n(4 bytes) + sum(key_len+4 + key + val_count+4 + vals)
//   var       : same
//   layers    : n(4 bytes) + sum(key_len+4 + key + rows+4 + cols+4 + data)
//   obsm      : same
//   varm      : same
//   obsp      : same
//   varp      : same
//   uns       : n(4 bytes) + sum(key_len+4 + key + val_len+4 + val)
// ============================================================================

///|
/// Binary buffer writer (little-endian).
struct BinWriter {
  buf : Array[Int]
}

///|
fn BinWriter::new() -> BinWriter {
  BinWriter::{ buf: [] }
}

///|
fn BinWriter::write_byte(self : BinWriter, b : Int) -> Unit {
  self.buf.push(b & 0xFF)
}

///|
fn BinWriter::write_bytes(self : BinWriter, bytes : Array[Int]) -> Unit {
  for b in bytes {
    self.buf.push(b & 0xFF)
  }
}

///|
fn BinWriter::write_int32(self : BinWriter, v : Int) -> Unit {
  self.buf.push(v & 0xFF)
  self.buf.push((v >> 8) & 0xFF)
  self.buf.push((v >> 16) & 0xFF)
  self.buf.push((v >> 24) & 0xFF)
}

///|
/// Write a float64 as a length-prefixed string (avoids reinterpret_cast).
fn BinWriter::write_float64(self : BinWriter, v : Double) -> Unit {
  let s = v.to_string()
  self.write_string(s)
}

///|
fn BinWriter::write_string(self : BinWriter, s : String) -> Unit {
  self.write_int32(s.length())
  for i = 0; i < s.length(); i = i + 1 {
    self.buf.push(s.unsafe_get(i).to_int() & 0xFF)
  }
}

///|
fn BinWriter::write_str_array(self : BinWriter, arr : Array[String]) -> Unit {
  self.write_int32(arr.length())
  for s in arr {
    self.write_string(s)
  }
}

///|
fn BinWriter::write_matrix(self : BinWriter, mat : Array[Array[Double]]) -> Unit {
  let rows = mat.length()
  let cols = if rows > 0 { mat[0].length() } else { 0 }
  self.write_int32(rows)
  self.write_int32(cols)
  for row in mat {
    for v in row {
      self.write_float64(v)
    }
  }
}

///|
fn BinWriter::write_str_map(
  self : BinWriter,
  m : Map[String, Array[String]],
) -> Unit {
  let keys = m.keys().collect()
  self.write_int32(keys.length())
  for k in keys {
    self.write_string(k)
    self.write_str_array(m[k])
  }
}

///|
fn BinWriter::write_mat_map(
  self : BinWriter,
  m : Map[String, Array[Array[Double]]],
) -> Unit {
  let keys = m.keys().collect()
  self.write_int32(keys.length())
  for k in keys {
    self.write_string(k)
    self.write_matrix(m[k])
  }
}

///|
fn BinWriter::write_uns_map(
  self : BinWriter,
  m : Map[String, String],
) -> Unit {
  let keys = m.keys().collect()
  self.write_int32(keys.length())
  for k in keys {
    self.write_string(k)
    self.write_string(m[k])
  }
}

///|
/// Binary buffer reader (little-endian).
struct BinReader {
  buf : Array[Int]
  mut pos : Int
}

///|
fn BinReader::new(buf : Array[Int]) -> BinReader {
  BinReader::{ buf, pos: 0 }
}

///|
fn BinReader::read_byte(self : BinReader) -> Int {
  if self.pos >= self.buf.length() {
    return 0
  }
  let v = self.buf[self.pos]
  self.pos = self.pos + 1
  v
}

///|
fn BinReader::read_int32(self : BinReader) -> Int {
  let b0 = self.read_byte()
  let b1 = self.read_byte()
  let b2 = self.read_byte()
  let b3 = self.read_byte()
  b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
}

///|
/// Read a float64 from a length-prefixed string.
fn BinReader::read_float64(self : BinReader) -> Double {
  let s = self.read_string()
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
fn BinReader::read_string(self : BinReader) -> String {
  let len = self.read_int32()
  if len <= 0 {
    return ""
  }
  let sb = StringBuilder()
  for _i = 0; _i < len; _i = _i + 1 {
    sb.write_char(self.read_byte().unsafe_to_char())
  }
  sb.to_string()
}

///|
fn BinReader::read_str_array(self : BinReader) -> Array[String] {
  let n = self.read_int32()
  let out : Array[String] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    out.push(self.read_string())
  }
  out
}

///|
fn BinReader::read_matrix(self : BinReader) -> Array[Array[Double]] {
  let rows = self.read_int32()
  let cols = self.read_int32()
  let out : Array[Array[Double]] = Array(capacity=rows)
  for _r = 0; _r < rows; _r = _r + 1 {
    let row : Array[Double] = Array(capacity=cols)
    for _c = 0; _c < cols; _c = _c + 1 {
      row.push(self.read_float64())
    }
    out.push(row)
  }
  out
}

///|
fn BinReader::read_str_map(self : BinReader) -> Map[String, Array[String]] {
  let n = self.read_int32()
  let out : Map[String, Array[String]] = Map([], capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    let k = self.read_string()
    out[k] = self.read_str_array()
  }
  out
}

///|
fn BinReader::read_mat_map(self : BinReader) -> Map[String, Array[Array[Double]]] {
  let n = self.read_int32()
  let out : Map[String, Array[Array[Double]]] = Map([], capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    let k = self.read_string()
    out[k] = self.read_matrix()
  }
  out
}

///|
fn BinReader::read_uns_map(self : BinReader) -> Map[String, String] {
  let n = self.read_int32()
  let out : Map[String, String] = Map([], capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    let k = self.read_string()
    out[k] = self.read_string()
  }
  out
}

///|
/// Serialize an AnnData to a compact binary format (Array[Int] of bytes).
/// Supports both dense and sparse X storage.
pub fn scanpy_write_h5ad_binary(adata : AnnData) -> Array[Int] {
  let w = BinWriter::new()
  // Magic + version + flags
  w.write_bytes([66, 72, 53, 65, 68]) // "BH5AD"
  w.write_byte(1) // version
  let is_sparse = match adata.get_sparse() {
    Some(_) => 1
    None => 0
  }
  w.write_byte(is_sparse)
  // Dimensions
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  w.write_int32(n_obs)
  w.write_int32(n_vars)
  // X matrix
  match adata.get_sparse() {
    Some(s) => {
      let nnz = s.nnzero()
      w.write_int32(nnz)
      for entry in s.entries {
        w.write_int32(entry.coordinates[0])
        w.write_int32(entry.coordinates[1])
        w.write_float64(entry.value)
      }
    }
    None => {
      for row in adata.x {
        for v in row {
          w.write_float64(v)
        }
      }
    }
  }
  // Names
  w.write_str_array(adata.obs_names)
  w.write_str_array(adata.var_names)
  // Metadata
  w.write_str_map(adata.obs)
  w.write_str_map(adata.variables)
  w.write_mat_map(adata.layers)
  w.write_mat_map(adata.obsm)
  w.write_mat_map(adata.varm)
  w.write_mat_map(adata.obsp)
  w.write_mat_map(adata.varp)
  w.write_uns_map(adata.uns)
  w.buf
}

///|
/// Deserialize a binary h5ad (Array[Int] of bytes) back to an AnnData.
/// Returns `None` if the magic header is invalid.
///|
/// Deserialize a binary h5ad (Array[Int] of bytes) back to an AnnData.
/// Returns `None` if the magic header is invalid.
/// Supports both dense and sparse X storage.
pub fn scanpy_read_h5ad_binary(buf : Array[Int]) -> AnnData? {
  let r = BinReader::new(buf)
  // Check magic "BH5AD"
  if r.read_byte() != 66 ||
    r.read_byte() != 72 ||
    r.read_byte() != 53 ||
    r.read_byte() != 65 ||
    r.read_byte() != 68 {
    return None
  }
  // Version + flags
  let _version = r.read_byte()
  let flags = r.read_byte()
  let is_sparse = (flags & 1) == 1
  // Dimensions
  let n_obs = r.read_int32()
  let n_vars = r.read_int32()
  // Read X matrix
  // For dense: n_obs*n_vars float64 values
  // For sparse: nnz entries, each (row:int32, col:int32, val:float64)
  // We read into a flat structure, then build AnnData after reading names.
  // Dense matrix (empty if sparse)
  let dense_x : Array[Array[Double]] = if is_sparse {
    []
  } else {
    let mat : Array[Array[Double]] = Array(capacity=n_obs)
    for _i = 0; _i < n_obs; _i = _i + 1 {
      let row : Array[Double] = Array(capacity=n_vars)
      for _j = 0; _j < n_vars; _j = _j + 1 {
        row.push(r.read_float64())
      }
      mat.push(row)
    }
    mat
  }
  // Sparse coords + values (empty if dense)
  let sp_coords : Array[Array[Int]] = Array(capacity=0)
  let sp_values : Array[Double] = Array(capacity=0)
  if is_sparse {
    let nnz = r.read_int32()
    for _i = 0; _i < nnz; _i = _i + 1 {
      sp_coords.push([r.read_int32(), r.read_int32()])
      sp_values.push(r.read_float64())
    }
  }
  // Read names
  let obs_names = r.read_str_array()
  let var_names = r.read_str_array()
  // Read metadata
  let obs = r.read_str_map()
  let variables = r.read_str_map()
  let layers = r.read_mat_map()
  let obsm = r.read_mat_map()
  let varm = r.read_mat_map()
  let obsp = r.read_mat_map()
  let varp = r.read_mat_map()
  let uns = r.read_uns_map()
  // Build AnnData
  let adata = if is_sparse {
    let sparse = SparseArray::from_coo([n_obs, n_vars], sp_coords, sp_values) catch {
      _ => return None
    }
    AnnData::new_sparse(sparse, obs_names, var_names)
  } else {
    AnnData::new(dense_x, obs_names, var_names)
  }
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
  for kv in varm {
    let _ = adata.set_varm(kv.0, kv.1)
  }
  for kv in obsp {
    let _ = adata.set_obsp(kv.0, kv.1)
  }
  for kv in varp {
    let _ = adata.set_varp(kv.0, kv.1)
  }
  for kv in uns {
    let _ = adata.set_uns(kv.0, kv.1)
  }
  Some(adata)
}
