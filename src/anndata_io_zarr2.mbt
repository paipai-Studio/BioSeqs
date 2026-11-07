/// anndata_io_zarr2.mbt — chunked, compressed Zarr v2 store for AnnData.
///
/// Upgrades the text proxy in anndata_io_zarr to a real Zarr v2 group
/// hierarchy held in a virtual file system (`ZarrStore.files`):
/// `.zgroup` / `.zarray` metadata with real chunk grids (`X/<br>.<bc>`
/// chunk keys), a line-level RLE compressor, and obs/var/uns sidecars —
/// mirroring `adata.write_zarr()` chunking semantics.

///|
/// A Zarr v2 store: virtual directory of text files plus chunk parameters.
pub struct ZarrStore {
  dir : String
  files : Map[String, String] // virtual FS: path -> content
  chunk_rows : Int32
  chunk_cols : Int32
  compressor : String // "none" or "rle"
}

///|
/// Create an empty store rooted at `dir`.
pub fn ZarrStore::new(
  dir : String,
  chunk_rows~ : Int32 = 1000,
  chunk_cols~ : Int32 = 1000,
  compressor~ : String = "none",
) -> ZarrStore {
  ZarrStore::{
    dir,
    files: Map([], capacity=16),
    chunk_rows: if chunk_rows <= 0 {
      1000
    } else {
      chunk_rows
    },
    chunk_cols: if chunk_cols <= 0 {
      1000
    } else {
      chunk_cols
    },
    compressor: if compressor == "rle" {
      "rle"
    } else {
      "none"
    },
  }
}

///|
/// Write a virtual file.
pub fn ZarrStore::put(
  self : ZarrStore,
  path : String,
  content : String,
) -> Unit {
  self.files[path] = content
}

///|
/// Read a virtual file (or `None`).
pub fn ZarrStore::get(self : ZarrStore, path : String) -> String? {
  self.files.get(path)
}

///|
/// True when the path exists in the store.
pub fn ZarrStore::has(self : ZarrStore, path : String) -> Bool {
  self.files.contains(path)
}

///|
/// All paths in the store.
pub fn ZarrStore::paths(self : ZarrStore) -> Array[String] {
  self.files.keys().collect()
}

///|
/// Number of virtual files.
pub fn ZarrStore::n_files(self : ZarrStore) -> Int32 {
  self.files.keys().collect().length()
}

///|
/// Chunk grid dimension along an axis: ceil(n / c).
pub fn anndata_zarr_chunk_count(n : Int32, c : Int32) -> Int32 {
  if c <= 0 || n <= 0 {
    0
  } else {
    (n + c - 1) / c
  }
}

///|
/// Shape (rows, cols) of the chunk at block position (br, bc), accounting
/// for ragged edge chunks. `None` when the block is outside the grid.
pub fn anndata_zarr_chunk_shape(
  store : ZarrStore,
  n_rows : Int32,
  n_cols : Int32,
  br : Int32,
  bc : Int32,
) -> (Int32, Int32)? {
  if br < 0 ||
    bc < 0 ||
    br >= anndata_zarr_chunk_count(n_rows, store.chunk_rows) ||
    bc >= anndata_zarr_chunk_count(n_cols, store.chunk_cols) {
    return None
  }
  let r0 = br * store.chunk_rows
  let c0 = bc * store.chunk_cols
  let rows = if r0 + store.chunk_rows > n_rows {
    n_rows - r0
  } else {
    store.chunk_rows
  }
  let cols = if c0 + store.chunk_cols > n_cols {
    n_cols - c0
  } else {
    store.chunk_cols
  }
  Some((rows, cols))
}

///|
/// Line-level RLE: a line equal to its predecessor repeated k times is
/// replaced by a marker line "R*k". Data lines are numeric CSV, so the
/// "R*" prefix is unambiguous. Encoded output ends with a newline.
pub fn zarr_rle_encode_lines(lines : Array[String]) -> String {
  let sb = StringBuilder()
  let mut i = 0
  while i < lines.length() {
    let mut j = i + 1
    while j < lines.length() && lines[j] == lines[i] {
      j = j + 1
    }
    let k = j - i
    let _ = sb.write(lines[i])
    let _ = sb.write("\n")
    if k > 1 {
      let _ = sb.write("R*")
      let _ = sb.write(k.to_string())
      let _ = sb.write("\n")
    }
    i = j
  }
  sb.to_string()
}

///|
/// Decode the line-level RLE produced by `zarr_rle_encode_lines`.
pub fn zarr_rle_decode_lines(content : String) -> Array[String] {
  let raw = content.split("\n").to_array()
  let out : Array[String] = Array(capacity=0)
  for line in raw {
    let s = line.to_owned()
    if s.has_prefix("R*") {
      // "R*<k>" — the marker plus its predecessor represent k identical
      // lines, so emit k - 1 additional copies of the last line
      let parts = s.split("*").to_array()
      if parts.length() >= 2 {
        let k = parse_int(parts[1].to_owned())
        if out.length() > 0 {
          let last = out[out.length() - 1]
          for _r = 0; _r < k - 1; _r = _r + 1 {
            out.push(last)
          }
        }
      }
    } else if s.length() > 0 {
      out.push(s)
    }
  }
  out
}

///|
/// Serialize one chunk of X (rows r0..r0+rows, cols c0..c0+cols) into
/// CSV lines, applying the store's compressor.
fn zarr_encode_chunk(
  store : ZarrStore,
  x : Array[Array[Double]],
  r0 : Int32,
  c0 : Int32,
  rows : Int32,
  cols : Int32,
) -> String {
  let lines : Array[String] = Array(capacity=rows)
  for i = 0; i < rows; i = i + 1 {
    let sb = StringBuilder()
    for j = 0; j < cols; j = j + 1 {
      if j > 0 {
        let _ = sb.write(",")
      }
      let v = if r0 + i < x.length() && c0 + j < x[r0 + i].length() {
        x[r0 + i][c0 + j]
      } else {
        0.0
      }
      let _ = sb.write(v.to_string())
    }
    lines.push(sb.to_string())
  }
  if store.compressor == "rle" {
    zarr_rle_encode_lines(lines)
  } else {
    let sb = StringBuilder()
    for l in lines {
      let _ = sb.write(l)
      let _ = sb.write("\n")
    }
    sb.to_string()
  }
}

///|
/// Write AnnData into a real chunked Zarr v2 group:
/// `.zgroup`, `X/.zarray` (shape/chunks/dtype/compressor), one
/// `X/<br>.<bc>` file per chunk, obs/var column sidecars, `obs_names` /
/// `var_names` and root `.zattrs` holding `uns`.
pub fn anndata_zarr_write(
  adata : AnnData,
  dir : String,
  chunk_rows~ : Int32 = 1000,
  chunk_cols~ : Int32 = 1000,
  compressor~ : String = "none",
) -> ZarrStore {
  let store = ZarrStore::new(dir, chunk_rows~, chunk_cols~, compressor~)
  store.put(dir + "/.zgroup", "{\"zarr_format\":2}")
  // X/.zarray metadata
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let zarray : Map[String, JsonValue] = Map([], capacity=8)
  zarray["zarr_format"] = JsonNum(2.0)
  let shape : Array[JsonValue] = Array(capacity=2)
  shape.push(JsonNum(n_obs.to_double()))
  shape.push(JsonNum(n_vars.to_double()))
  zarray["shape"] = JsonArr(shape)
  let chunks : Array[JsonValue] = Array(capacity=2)
  chunks.push(JsonNum(store.chunk_rows.to_double()))
  chunks.push(JsonNum(store.chunk_cols.to_double()))
  zarray["chunks"] = JsonArr(chunks)
  zarray["dtype"] = JsonStr("<f8")
  match store.compressor {
    "rle" => {
      let comp : Map[String, JsonValue] = Map([], capacity=2)
      comp["id"] = JsonStr("zlib-rle")
      comp["level"] = JsonNum(1.0)
      zarray["compressor"] = JsonObj(comp)
    }
    _ => zarray["compressor"] = JsonNull
  }
  zarray["fill_value"] = JsonNum(0.0)
  zarray["order"] = JsonStr("C")
  zarray["filters"] = JsonNull
  store.put(dir + "/X/.zarray", json_to_string(JsonObj(zarray)))
  // X chunks
  let x = adata.x_matrix()
  let br_max = anndata_zarr_chunk_count(n_obs, store.chunk_rows)
  let bc_max = anndata_zarr_chunk_count(n_vars, store.chunk_cols)
  for br = 0; br < br_max; br = br + 1 {
    for bc = 0; bc < bc_max; bc = bc + 1 {
      let r0 = br * store.chunk_rows
      let c0 = bc * store.chunk_cols
      let rows = if r0 + store.chunk_rows > n_obs {
        n_obs - r0
      } else {
        store.chunk_rows
      }
      let cols = if c0 + store.chunk_cols > n_vars {
        n_vars - c0
      } else {
        store.chunk_cols
      }
      store.put(
        dir + "/X/" + br.to_string() + "." + bc.to_string(),
        zarr_encode_chunk(store, x, r0, c0, rows, cols),
      )
    }
  }
  // obs / var column sidecars
  for k in adata.obs.keys().collect() {
    let sb = StringBuilder()
    for v in adata.obs[k] {
      let _ = sb.write(v)
      let _ = sb.write("\n")
    }
    store.put(dir + "/obs/" + k, sb.to_string())
  }
  for k in adata.variables.keys().collect() {
    let sb = StringBuilder()
    for v in adata.variables[k] {
      let _ = sb.write(v)
      let _ = sb.write("\n")
    }
    store.put(dir + "/var/" + k, sb.to_string())
  }
  // names
  let sb_obs = StringBuilder()
  for n in adata.obs_names {
    let _ = sb_obs.write(n)
    let _ = sb_obs.write("\n")
  }
  store.put(dir + "/obs_names", sb_obs.to_string())
  let sb_var = StringBuilder()
  for n in adata.var_names {
    let _ = sb_var.write(n)
    let _ = sb_var.write("\n")
  }
  store.put(dir + "/var_names", sb_var.to_string())
  // root .zattrs = uns
  let uns_obj : Map[String, JsonValue] = Map(
    [],
    capacity=adata.uns.keys().collect().length(),
  )
  for k in adata.uns.keys().collect() {
    uns_obj[k] = adata.uns[k]
  }
  store.put(dir + "/.zattrs", json_to_string(JsonObj(uns_obj)))
  store
}

///|
/// Split a virtual text file into non-empty lines.
fn zarr_read_lines(store : ZarrStore, path : String) -> Array[String] {
  match store.get(path) {
    Some(c) => {
      let out : Array[String] = Array(capacity=0)
      for line in c.split("\n").to_array() {
        let s = line.to_owned()
        if s.length() > 0 {
          out.push(s)
        }
      }
      out
    }
    None => []
  }
}

///|
/// Last "/"-separated segment of a path (column key extraction).
fn zarr_path_leaf(path : String) -> String {
  let parts = path.split("/").to_array()
  if parts.length() > 0 {
    parts[parts.length() - 1].to_owned()
  } else {
    ""
  }
}

///|
/// Read a chunked Zarr store back into an AnnData (inverse of
/// `anndata_zarr_write`). Returns `None` when metadata is missing/malformed.
pub fn anndata_zarr_read(store : ZarrStore) -> AnnData? {
  let dir = store.dir
  // .zarray -> shape
  let zarray_content = match store.get(dir + "/X/.zarray") {
    Some(c) => c
    None => return None
  }
  let zarray = match parse_json(zarray_content) {
    Some(JsonObj(m)) => m
    _ => return None
  }
  let dims : Array[Int32] = Array(capacity=2)
  match zarray.get("shape") {
    Some(JsonArr(arr)) => {
      if arr.length() < 2 {
        return None
      }
      for v in arr {
        match v {
          JsonNum(d) => dims.push(d.to_int())
          _ => dims.push(0)
        }
      }
    }
    _ => return None
  }
  let n_obs = dims[0]
  let n_vars = dims[1]
  // reassemble X from chunks
  let x : Array[Array[Double]] = Array(capacity=n_obs)
  for _i = 0; _i < n_obs; _i = _i + 1 {
    let row : Array[Double] = Array(capacity=n_vars)
    for _j = 0; _j < n_vars; _j = _j + 1 {
      row.push(0.0)
    }
    x.push(row)
  }
  let br_max = anndata_zarr_chunk_count(n_obs, store.chunk_rows)
  let bc_max = anndata_zarr_chunk_count(n_vars, store.chunk_cols)
  for br = 0; br < br_max; br = br + 1 {
    for bc = 0; bc < bc_max; bc = bc + 1 {
      let path = dir + "/X/" + br.to_string() + "." + bc.to_string()
      match store.get(path) {
        Some(content) => {
          let lines = match store.compressor {
            "rle" => zarr_rle_decode_lines(content)
            _ => zarr_read_lines(store, path)
          }
          let r0 = br * store.chunk_rows
          let c0 = bc * store.chunk_cols
          for i = 0; i < lines.length() && r0 + i < n_obs; i = i + 1 {
            let fields = lines[i].split(",").to_array()
            for j = 0; j < fields.length() && c0 + j < n_vars; j = j + 1 {
              x[r0 + i][c0 + j] = parse_double(fields[j].to_owned()).unwrap_or(
                0.0,
              )
            }
          }
        }
        None => ()
      }
    }
  }
  // names
  let obs_names = zarr_read_lines(store, dir + "/obs_names")
  let var_names = zarr_read_lines(store, dir + "/var_names")
  let adata = AnnData::new(x, obs_names, var_names)
  // obs / var columns
  for path in store.paths() {
    if path.has_prefix(dir + "/obs/") && not(path.has_suffix(".zarray")) {
      adata.obs[zarr_path_leaf(path)] = zarr_read_lines(store, path)
    } else if path.has_prefix(dir + "/var/") && not(path.has_suffix(".zarray")) {
      adata.variables[zarr_path_leaf(path)] = zarr_read_lines(store, path)
    }
  }
  // uns
  match store.get(dir + "/.zattrs") {
    Some(c) =>
      match parse_json(c) {
        Some(JsonObj(m)) =>
          for k in m.keys().collect() {
            adata.uns[k] = m[k]
          }
        _ => ()
      }
    None => ()
  }
  Some(adata)
}

///|
/// Sorted virtual-path listing (zarr hierarchy dump).
pub fn anndata_zarr_hierarchy(store : ZarrStore) -> Array[String] {
  let out = store.paths()
  out.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  out
}
