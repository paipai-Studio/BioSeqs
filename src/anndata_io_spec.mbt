/// anndata_io_spec.mbt — full h5ad encoding-type / encoding-version spec.
///
/// Mirrors `anndata.io.specs`: every element written into an h5ad store is
/// tagged with an `encoding-type` and `encoding-version`. This module
/// implements the complete registry (array, dataframe, dict, categorical,
/// categorical-masked, string-array, nullable-*, csr/csc, masked-array) and
/// uses it to serialize a whole AnnData to / from a spec-conformant JSON
/// representation of an h5ad file (anndata_h5ad_encode / decode).

///|
/// Full h5ad encoding registry (anndata.io.specs). Names are disjoint from
/// the legacy AnnDataEncoding enum in anndata_elem.
pub(all) enum AnnDataEnc {
  SpecArray // "array" (dense n-d array)
  SpecDataFrame // "dataframe"
  SpecDict // "dict"
  SpecCategorical // "categorical"
  SpecCategoricalMasked // "categorical-masked"
  SpecStringArray // "string-array"
  SpecNullableInteger // "nullable-integer"
  SpecNullableBoolean // "nullable-boolean"
  SpecNullableStringArray // "nullable-string-array"
  SpecCSR // "csr_matrix"
  SpecCSC // "csc_matrix"
  SpecScalar // "scalar"
  SpecMaskedArray // "masked-array" (deprecated in anndata, kept for reads)
}

///|
/// Canonical encoding-type name.
pub fn anndata_enc_name(enc : AnnDataEnc) -> String {
  match enc {
    SpecArray => "array"
    SpecDataFrame => "dataframe"
    SpecDict => "dict"
    SpecCategorical => "categorical"
    SpecCategoricalMasked => "categorical-masked"
    SpecStringArray => "string-array"
    SpecNullableInteger => "nullable-integer"
    SpecNullableBoolean => "nullable-boolean"
    SpecNullableStringArray => "nullable-string-array"
    SpecCSR => "csr_matrix"
    SpecCSC => "csc_matrix"
    SpecScalar => "scalar"
    SpecMaskedArray => "masked-array"
  }
}

///|
/// Canonical encoding-version for each encoding-type (matches anndata
/// 0.11's io specs registry).
pub fn anndata_enc_version(enc : AnnDataEnc) -> String {
  match enc {
    SpecArray => "0.2.0"
    SpecDataFrame => "0.2.0"
    SpecDict => "0.1.0"
    SpecCategorical => "0.2.0"
    SpecCategoricalMasked => "0.2.0"
    SpecStringArray => "0.2.0"
    SpecNullableInteger => "0.1.0"
    SpecNullableBoolean => "0.1.0"
    SpecNullableStringArray => "0.1.0"
    SpecCSR => "0.1.0"
    SpecCSC => "0.1.0"
    SpecScalar => "0.1.0"
    SpecMaskedArray => "0.1.0"
  }
}

///|
/// Look up an encoding by its canonical name.
pub fn anndata_enc_from_name(name : String) -> AnnDataEnc? {
  if name == "array" {
    Some(SpecArray)
  } else if name == "dataframe" {
    Some(SpecDataFrame)
  } else if name == "dict" {
    Some(SpecDict)
  } else if name == "categorical" {
    Some(SpecCategorical)
  } else if name == "categorical-masked" {
    Some(SpecCategoricalMasked)
  } else if name == "string-array" {
    Some(SpecStringArray)
  } else if name == "nullable-integer" {
    Some(SpecNullableInteger)
  } else if name == "nullable-boolean" {
    Some(SpecNullableBoolean)
  } else if name == "nullable-string-array" {
    Some(SpecNullableStringArray)
  } else if name == "csr_matrix" {
    Some(SpecCSR)
  } else if name == "csc_matrix" {
    Some(SpecCSC)
  } else if name == "scalar" {
    Some(SpecScalar)
  } else if name == "masked-array" {
    Some(SpecMaskedArray)
  } else {
    None
  }
}

///|
/// The {"encoding-type", "encoding-version"} attribute pair for an encoding.
pub fn anndata_enc_attrs(enc : AnnDataEnc) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=2)
  m["encoding-type"] = JsonStr(anndata_enc_name(enc))
  m["encoding-version"] = JsonStr(anndata_enc_version(enc))
  JsonObj(m)
}

///|
/// Read the encoding declared by an attrs object (or `None`).
pub fn anndata_enc_from_attrs(attrs : JsonValue) -> AnnDataEnc? {
  match attrs {
    JsonObj(m) =>
      match m.get("encoding-type") {
        Some(JsonStr(s)) => anndata_enc_from_name(s)
        _ => None
      }
    _ => None
  }
}

///|
/// Validate a declared version against the registry. Returns an error
/// message on mismatch, or `None` when the pair is spec-conformant.
pub fn anndata_enc_check(name : String, version : String) -> String? {
  match anndata_enc_from_name(name) {
    None => Some("unknown encoding-type: " + name)
    Some(enc) =>
      if anndata_enc_version(enc) != version {
        Some(
          "encoding-version mismatch for " +
          name +
          ": expected " +
          anndata_enc_version(enc) +
          ", got " +
          version,
        )
      } else {
        None
      }
  }
}

///|
/// Structural validation of a CSR/CSC element (data/indices/indptr/shape).
/// Returns an error message or `None` when the element is well-formed.
pub fn anndata_sparse_spec_check(
  data_len : Int32,
  indices_len : Int32,
  indptr_len : Int32,
  n_rows : Int32,
) -> String? {
  if data_len != indices_len {
    Some("csr/csc: data and indices must have the same length")
  } else if indptr_len != n_rows + 1 {
    Some("csr/csc: indptr must have length n_rows + 1")
  } else {
    None
  }
}

// === spec-conformant JSON representation of a full h5ad file ===

///|
/// Wrap a payload with its encoding attributes.
fn io_spec_wrap(enc : AnnDataEnc, value : JsonValue) -> JsonValue {
  match anndata_enc_attrs(enc) {
    JsonObj(attrs) => {
      attrs["value"] = value
      JsonObj(attrs)
    }
    _ => JsonObj(Map([], capacity=0))
  }
}

///|
/// Array[String] -> JSON array of strings.
fn io_spec_json_strs(arr : Array[String]) -> JsonValue {
  let out : Array[JsonValue] = Array(capacity=arr.length())
  for s in arr {
    out.push(JsonStr(s))
  }
  JsonArr(out)
}

///|
/// Dense matrix -> nested JSON array.
fn io_spec_json_mat(m : Array[Array[Double]]) -> JsonValue {
  let rows : Array[JsonValue] = Array(capacity=m.length())
  for row in m {
    let cols : Array[JsonValue] = Array(capacity=row.length())
    for v in row {
      cols.push(JsonNum(v))
    }
    rows.push(JsonArr(cols))
  }
  JsonArr(rows)
}

///|
/// JSON array of strings -> Array[String].
fn io_spec_strs(v : JsonValue) -> Array[String] {
  match v {
    JsonArr(arr) => {
      let out : Array[String] = Array(capacity=arr.length())
      for item in arr {
        match item {
          JsonStr(s) => out.push(s)
          _ => out.push("")
        }
      }
      out
    }
    _ => []
  }
}

///|
/// Nested JSON array -> dense matrix.
fn io_spec_mat(v : JsonValue) -> Array[Array[Double]] {
  match v {
    JsonArr(rows) => {
      let out : Array[Array[Double]] = Array(capacity=rows.length())
      for row in rows {
        match row {
          JsonArr(cols) => {
            let r : Array[Double] = Array(capacity=cols.length())
            for c in cols {
              match c {
                JsonNum(d) => r.push(d)
                _ => r.push(0.0)
              }
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
/// JSON array of numbers -> Array[Double].
fn io_spec_nums(v : JsonValue) -> Array[Double] {
  match v {
    JsonArr(arr) => {
      let out : Array[Double] = Array(capacity=arr.length())
      for item in arr {
        match item {
          JsonNum(d) => out.push(d)
          _ => out.push(0.0)
        }
      }
      out
    }
    _ => []
  }
}

///|
/// Does the uns metadata mark this obs/var column as categorical?
fn io_spec_is_cat(
  uns : Map[String, JsonValue],
  axis : String,
  col : String,
) -> Bool {
  match uns.get("__typed__/" + axis + "/" + col) {
    Some(JsonObj(m)) =>
      match m.get("kind") {
        Some(JsonStr(s)) => s == "categorical"
        _ => false
      }
    _ => false
  }
}

///|
/// Encode one obs/var column value with its spec-conformant encoding.
fn io_spec_col_value(
  uns : Map[String, JsonValue],
  axis : String,
  col : String,
  values : Array[String],
) -> JsonValue {
  if io_spec_is_cat(uns, axis, col) {
    // categorical: categories + codes
    let cats = anndata_obs_var_unique_values(values)
    let cat_arr : Array[JsonValue] = Array(capacity=cats.length())
    for c in cats {
      cat_arr.push(JsonStr(c))
    }
    let codes : Array[JsonValue] = Array(capacity=values.length())
    for v in values {
      let mut code = -1
      for k = 0; k < cats.length(); k = k + 1 {
        if cats[k] == v {
          code = k
          break
        }
      }
      codes.push(JsonNum(code.to_double()))
    }
    let inner : Map[String, JsonValue] = Map([], capacity=3)
    inner["categories"] = JsonArr(cat_arr)
    inner["codes"] = JsonArr(codes)
    inner["ordered"] = JsonBool(false)
    io_spec_wrap(SpecCategorical, JsonObj(inner))
  } else {
    io_spec_wrap(SpecStringArray, io_spec_json_strs(values))
  }
}

///|
/// Unique values of a string column (order of first appearance).
fn anndata_obs_var_unique_values(values : Array[String]) -> Array[String] {
  let seen : Set[String] = Set([])
  let out : Array[String] = Array(capacity=0)
  for v in values {
    if not(seen.contains(v)) {
      seen.add(v)
      out.push(v)
    }
  }
  out
}

///|
/// Encode an obs/var dataframe with the dataframe encoding.
fn io_spec_df_encode(
  names : Array[String],
  cols : Map[String, Array[String]],
  uns : Map[String, JsonValue],
  axis : String,
) -> JsonValue {
  let df : Map[String, JsonValue] = Map([], capacity=4)
  df["_index"] = io_spec_wrap(SpecStringArray, io_spec_json_strs(names))
  let keys = cols.keys().collect()
  let order : Array[JsonValue] = Array(capacity=keys.length())
  for k in keys {
    order.push(JsonStr(k))
    df[k] = io_spec_col_value(uns, axis, k, cols[k])
  }
  df["column-order"] = JsonArr(order)
  io_spec_wrap(SpecDataFrame, JsonObj(df))
}

///|
/// Decode a dataframe element: returns (index, columns).
fn io_spec_df_decode(
  v : JsonValue,
) -> (Array[String], Map[String, Array[String]])? {
  // the dataframe payload sits under the "value" key of the wrapping attrs
  let inner = io_spec_unwrap(v)
  match inner {
    JsonObj(m) => {
      let names = match m.get("_index") {
        Some(elem) => io_spec_strs(io_spec_unwrap(elem))
        None => []
      }
      let cols : Map[String, Array[String]] = Map([], capacity=4)
      match m.get("column-order") {
        Some(JsonArr(order)) =>
          for k in order {
            match k {
              JsonStr(key) =>
                match m.get(key) {
                  Some(col_elem) =>
                    match io_spec_col_decode(col_elem) {
                      Some((_kind, vals, _meta)) => cols[key] = vals
                      None => cols[key] = []
                    }
                  _ => cols[key] = []
                }
              _ => ()
            }
          }
        _ => ()
      }
      Some((names, cols))
    }
    _ => None
  }
}

///|
/// Unwrap a {"encoding-type", "encoding-version", "value": v} element.
fn io_spec_unwrap(v : JsonValue) -> JsonValue {
  match v {
    JsonObj(m) =>
      match m.get("value") {
        Some(inner) => inner
        None => v
      }
    _ => v
  }
}

///|
/// Decode a single column element. Returns
/// (encoding-name, values, raw-inner-map-for-metadata).
fn io_spec_col_decode(
  v : JsonValue,
) -> (String, Array[String], Map[String, JsonValue])? {
  match v {
    JsonObj(m) => {
      let et = match m.get("encoding-type") {
        Some(JsonStr(s)) => s
        _ => "string-array"
      }
      let payload = match m.get("value") {
        Some(inner) => inner
        None => JsonObj(m)
      }
      if et == "categorical" {
        match payload {
          JsonObj(inner) => {
            let cats = match inner.get("categories") {
              Some(c) => io_spec_strs(c)
              None => []
            }
            let codes = match inner.get("codes") {
              Some(c) => io_spec_nums(c)
              None => []
            }
            let out : Array[String] = Array(capacity=codes.length())
            for code in codes {
              let k = code.to_int()
              if k >= 0 && k < cats.length() {
                out.push(cats[k])
              } else {
                out.push("")
              }
            }
            Some((et, out, inner))
          }
          _ => Some((et, [], Map([], capacity=0)))
        }
      } else {
        Some((et, io_spec_strs(payload), Map([], capacity=0)))
      }
    }
    _ => None
  }
}

///|
/// Encode X (dense or sparse) with its spec element.
fn io_spec_x_encode(adata : AnnData) -> JsonValue {
  match adata.get_sparse() {
    Some(s) => {
      let dims = s.dim()
      let coords = s.nzcoordinates()
      let values = s.nzvalues()
      let data : Array[JsonValue] = Array(capacity=values.length())
      let indices : Array[JsonValue] = Array(capacity=values.length())
      for v in values {
        data.push(JsonNum(v))
      }
      for c in coords {
        indices.push(JsonNum(c[1].to_double()))
      }
      let indptr : Array[JsonValue] = Array(capacity=dims[0] + 1)
      let counts : Array[Int32] = Array::make(dims[0] + 1, 0)
      for c in coords {
        if c[0] >= 0 && c[0] < dims[0] {
          counts[c[0] + 1] = counts[c[0] + 1] + 1
        }
      }
      let mut cum : Int32 = 0
      indptr.push(JsonNum(0.0))
      for r = 0; r < dims[0]; r = r + 1 {
        cum = cum + counts[r + 1]
        indptr.push(JsonNum(cum.to_double()))
      }
      indptr.push(JsonNum(cum.to_double()))
      let inner : Map[String, JsonValue] = Map([], capacity=4)
      let shape : Array[JsonValue] = Array(capacity=2)
      shape.push(JsonNum(dims[0].to_double()))
      shape.push(JsonNum(dims[1].to_double()))
      inner["shape"] = JsonArr(shape)
      inner["data"] = JsonArr(data)
      inner["indices"] = JsonArr(indices)
      inner["indptr"] = JsonArr(indptr)
      io_spec_wrap(SpecCSR, JsonObj(inner))
    }
    None => io_spec_wrap(SpecArray, io_spec_json_mat(adata.x))
  }
}

///|
/// Full spec-conformant encoding of an AnnData as a JSON value.
pub fn anndata_h5ad_encode(adata : AnnData) -> JsonValue {
  let root : Map[String, JsonValue] = Map([], capacity=16)
  root["encoding-type"] = JsonStr("anndata")
  root["encoding-version"] = JsonStr("0.1.0")
  root["X"] = io_spec_x_encode(adata)
  root["obs"] = io_spec_df_encode(adata.obs_names, adata.obs, adata.uns, "obs")
  root["var"] = io_spec_df_encode(
    adata.var_names,
    adata.variables,
    adata.uns,
    "var",
  )
  root["obs_names"] = io_spec_wrap(
    SpecStringArray,
    io_spec_json_strs(adata.obs_names),
  )
  root["var_names"] = io_spec_wrap(
    SpecStringArray,
    io_spec_json_strs(adata.var_names),
  )
  // dict-typed slots
  let dict_slot = fn(m : Map[String, Array[Array[Double]]]) -> JsonValue {
    let inner : Map[String, JsonValue] = Map([], capacity=4)
    for k in m.keys().collect() {
      inner[k] = io_spec_wrap(SpecArray, io_spec_json_mat(m[k]))
    }
    io_spec_wrap(SpecDict, JsonObj(inner))
  }
  root["layers"] = dict_slot(adata.layers)
  root["obsm"] = dict_slot(adata.obsm)
  root["varm"] = dict_slot(adata.varm)
  root["obsp"] = dict_slot(adata.obsp)
  root["varp"] = dict_slot(adata.varp)
  root["uns"] = io_spec_wrap(SpecDict, JsonObj(adata.uns))
  // raw
  match adata.raw {
    Some(r) => {
      let raw : Map[String, JsonValue] = Map([], capacity=4)
      raw["X"] = io_spec_wrap(SpecArray, io_spec_json_mat(r.x))
      raw["var_names"] = io_spec_wrap(
        SpecStringArray,
        io_spec_json_strs(r.var_names),
      )
      raw["var"] = io_spec_df_encode(r.var_names, r.var_anno, adata.uns, "var")
      root["raw"] = JsonObj(raw)
    }
    None => root["raw"] = JsonNull
  }
  JsonObj(root)
}

///|
/// Serialize the spec-conformant h5ad JSON to a string.
pub fn anndata_h5ad_encode_string(adata : AnnData) -> String {
  json_to_string(anndata_h5ad_encode(adata))
}

///|
/// Decode a spec-conformant h5ad JSON value back into an AnnData.
/// Elements missing required encodings degrade to empty slots.
pub fn anndata_h5ad_decode(value : JsonValue) -> AnnData? {
  match value {
    JsonObj(root) => {
      // X
      let x_elem = match root.get("X") {
        Some(x) => x
        None => return None
      }
      let mut x_mat : Array[Array[Double]] = []
      let mut x_sparse : SparseArray? = None
      let et = match x_elem {
        JsonObj(m) =>
          match m.get("encoding-type") {
            Some(JsonStr(s)) => s
            _ => "array"
          }
        _ => "array"
      }
      if et == "array" {
        x_mat = io_spec_mat(io_spec_unwrap(x_elem))
      } else {
        // csr/csc: rebuild sparse from data/indices/indptr/shape
        let payload = io_spec_unwrap(x_elem)
        match payload {
          JsonObj(m) => {
            let shape = match m.get("shape") {
              Some(s) => io_spec_nums(s)
              None => []
            }
            let data = io_spec_nums(
              match m.get("data") {
                Some(d) => d
                None => JsonArr([])
              },
            )
            let indices = io_spec_nums(
              match m.get("indices") {
                Some(d) => d
                None => JsonArr([])
              },
            )
            let indptr = io_spec_nums(
              match m.get("indptr") {
                Some(d) => d
                None => JsonArr([])
              },
            )
            if shape.length() == 2 {
              let n_rows = shape[0].to_int()
              let n_cols = shape[1].to_int()
              let coords : Array[Array[Int32]] = Array(capacity=data.length())
              let vals : Array[Double] = Array(capacity=data.length())
              for r = 0; r < n_rows; r = r + 1 {
                if r + 1 < indptr.length() {
                  let start = indptr[r].to_int()
                  let end_ = if r + 1 < indptr.length() {
                    indptr[r + 1].to_int()
                  } else {
                    start
                  }
                  for k = start; k < end_ && k < indices.length(); k = k + 1 {
                    let c : Array[Int32] = Array(capacity=2)
                    c.push(r)
                    c.push(indices[k].to_int())
                    coords.push(c)
                    vals.push(if k < data.length() { data[k] } else { 0.0 })
                  }
                }
              }
              let arr = SparseArray::from_coo([n_rows, n_cols], coords, vals) catch {
                _ => return None
              }
              x_sparse = Some(arr)
            }
          }
          _ => ()
        }
      }
      // obs / var dataframes (decoded once; names reused when the
      // dedicated obs_names/var_names elements are absent)
      let obs_names : Array[String] = []
      let obs_cols : Map[String, Array[String]] = Map([], capacity=4)
      match root.get("obs") {
        Some(o) =>
          match io_spec_df_decode(o) {
            Some((names, cols)) => {
              for k in cols.keys().collect() {
                obs_cols[k] = cols[k]
              }
              if names.length() > 0 {
                for s in names {
                  obs_names.push(s)
                }
              }
            }
            None => ()
          }
        None => ()
      }
      let obs_names_final = match root.get("obs_names") {
        Some(v) => io_spec_strs(io_spec_unwrap(v))
        None => obs_names
      }
      let var_names : Array[String] = []
      let var_cols : Map[String, Array[String]] = Map([], capacity=4)
      match root.get("var") {
        Some(v) =>
          match io_spec_df_decode(v) {
            Some((names, cols)) => {
              for k in cols.keys().collect() {
                var_cols[k] = cols[k]
              }
              if names.length() > 0 {
                for s in names {
                  var_names.push(s)
                }
              }
            }
            None => ()
          }
        None => ()
      }
      let var_names_final = match root.get("var_names") {
        Some(v) => io_spec_strs(io_spec_unwrap(v))
        None => var_names
      }
      let adata = match x_sparse {
        Some(s) => AnnData::new_sparse(s, obs_names_final, var_names_final)
        None => AnnData::new(x_mat, obs_names_final, var_names_final)
      }
      for k in obs_cols.keys().collect() {
        adata.obs[k] = obs_cols[k]
      }
      for k in var_cols.keys().collect() {
        adata.variables[k] = var_cols[k]
      }
      // dict slots
      let dict_decode = fn(
        v : JsonValue,
        target : Map[String, Array[Array[Double]]],
      ) -> Unit {
        match v {
          JsonObj(m) =>
            match m.get("value") {
              Some(JsonObj(inner)) =>
                for k in inner.keys().collect() {
                  target[k] = io_spec_mat(io_spec_unwrap(inner[k]))
                }
              _ => ()
            }
          _ => ()
        }
      }
      match root.get("layers") {
        Some(v) => dict_decode(v, adata.layers)
        None => ()
      }
      match root.get("obsm") {
        Some(v) => dict_decode(v, adata.obsm)
        None => ()
      }
      match root.get("varm") {
        Some(v) => dict_decode(v, adata.varm)
        None => ()
      }
      match root.get("obsp") {
        Some(v) => dict_decode(v, adata.obsp)
        None => ()
      }
      match root.get("varp") {
        Some(v) => dict_decode(v, adata.varp)
        None => ()
      }
      // uns
      match root.get("uns") {
        Some(JsonObj(m)) =>
          match m.get("value") {
            Some(JsonObj(inner)) =>
              for k in inner.keys().collect() {
                adata.uns[k] = inner[k]
              }
            _ => ()
          }
        _ => ()
      }
      // raw
      match root.get("raw") {
        Some(JsonObj(raw)) => {
          let raw_x = match raw.get("X") {
            Some(x) => io_spec_mat(io_spec_unwrap(x))
            None => []
          }
          let raw_var_names = match raw.get("var_names") {
            Some(v) => io_spec_strs(io_spec_unwrap(v))
            None => []
          }
          let raw_var : Map[String, Array[String]] = Map([], capacity=4)
          match raw.get("var") {
            Some(v) =>
              match io_spec_df_decode(v) {
                Some((_names, cols)) =>
                  for k in cols.keys().collect() {
                    raw_var[k] = cols[k]
                  }
                None => ()
              }
            None => ()
          }
          adata.raw = Some({
            x: raw_x,
            var_names: raw_var_names,
            var_anno: raw_var,
          })
        }
        _ => ()
      }
      Some(adata)
    }
    _ => None
  }
}

///|
/// Parse + decode a spec-conformant h5ad JSON string.
pub fn anndata_h5ad_decode_string(content : String) -> AnnData? {
  match parse_json(content) {
    Some(v) => anndata_h5ad_decode(v)
    None => None
  }
}
