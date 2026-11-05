/// AnnData uns dict operations.
///
/// `adata.uns["key"] = value` and nested dict access.
///
/// The AnnData `uns` field stores arbitrary unstructured annotations.
/// This module provides typed accessors for common patterns.

///|
/// Get a string value from uns.
pub fn anndata_uns_get_str(adata : AnnData, key : String) -> String? {
  match adata.get_uns(key) {
    Some(JsonStr(s)) => Some(s)
    _ => None
  }
}

///|
/// Get a numeric value from uns.
pub fn anndata_uns_get_num(adata : AnnData, key : String) -> Double? {
  match adata.get_uns(key) {
    Some(JsonNum(n)) => Some(n)
    _ => None
  }
}

///|
/// Get a string array from uns.
pub fn anndata_uns_get_str_array(
  adata : AnnData,
  key : String,
) -> Array[String]? {
  match adata.get_uns(key) {
    Some(JsonArr(arr)) => {
      let out : Array[String] = Array(capacity=arr.length())
      for v in arr {
        match v {
          JsonStr(s) => out.push(s)
          _ => out.push("")
        }
      }
      Some(out)
    }
    _ => None
  }
}

///|
/// Get a numeric array from uns.
pub fn anndata_uns_get_num_array(
  adata : AnnData,
  key : String,
) -> Array[Double]? {
  match adata.get_uns(key) {
    Some(JsonArr(arr)) => {
      let out : Array[Double] = Array(capacity=arr.length())
      for v in arr {
        match v {
          JsonNum(n) => out.push(n)
          _ => out.push(0.0)
        }
      }
      Some(out)
    }
    _ => None
  }
}

///|
/// Get a 2D numeric array from uns (e.g. PCA components).
pub fn anndata_uns_get_matrix(
  adata : AnnData,
  key : String,
) -> Array[Array[Double]]? {
  match adata.get_uns(key) {
    Some(JsonArr(rows)) => {
      let out : Array[Array[Double]] = Array(capacity=0)
      for row in rows {
        match row {
          JsonArr(cols) => {
            let r : Array[Double] = Array(capacity=cols.length())
            for v in cols {
              match v {
                JsonNum(n) => r.push(n)
                _ => r.push(0.0)
              }
            }
            out.push(r)
          }
          _ => out.push([])
        }
      }
      Some(out)
    }
    _ => None
  }
}

///|
/// Set a string in uns.
pub fn anndata_uns_set_str(
  adata : AnnData,
  key : String,
  value : String,
) -> Unit {
  ignore(adata.set_uns(key, JsonStr(value)))
}

///|
/// Set a number in uns.
pub fn anndata_uns_set_num(
  adata : AnnData,
  key : String,
  value : Double,
) -> Unit {
  ignore(adata.set_uns(key, JsonNum(value)))
}

///|
/// Set a string array in uns.
pub fn anndata_uns_set_str_array(
  adata : AnnData,
  key : String,
  values : Array[String],
) -> Unit {
  let arr : Array[JsonValue] = Array(capacity=values.length())
  for s in values {
    arr.push(JsonStr(s))
  }
  ignore(adata.set_uns(key, JsonArr(arr)))
}

///|
/// Set a numeric array in uns.
pub fn anndata_uns_set_num_array(
  adata : AnnData,
  key : String,
  values : Array[Double],
) -> Unit {
  let arr : Array[JsonValue] = Array(capacity=values.length())
  for v in values {
    arr.push(JsonNum(v))
  }
  ignore(adata.set_uns(key, JsonArr(arr)))
}

///|
/// Set a 2D matrix in uns.
pub fn anndata_uns_set_matrix(
  adata : AnnData,
  key : String,
  values : Array[Array[Double]],
) -> Unit {
  let rows : Array[JsonValue] = Array(capacity=values.length())
  for row in values {
    let r : Array[JsonValue] = Array(capacity=row.length())
    for v in row {
      r.push(JsonNum(v))
    }
    rows.push(JsonArr(r))
  }
  ignore(adata.set_uns(key, JsonArr(rows)))
}

///|
/// Delete a key from uns.
pub fn anndata_uns_del(adata : AnnData, key : String) -> Unit {
  adata.uns.remove(key)
}

///|
/// List all uns keys.
pub fn anndata_uns_keys(adata : AnnData) -> Array[String] {
  adata.uns.keys().collect()
}

///|
/// Check if uns has a key.
pub fn anndata_uns_has(adata : AnnData, key : String) -> Bool {
  match adata.get_uns(key) {
    Some(_) => true
    None => false
  }
}

///|
/// Copy uns from one AnnData to another.
pub fn anndata_uns_copy(src : AnnData, dst : AnnData) -> Unit {
  let keys = src.uns.keys().collect()
  for k in keys {
    match src.get_uns(k) {
      Some(v) => ignore(dst.set_uns(k, v))
      None => ()
    }
  }
}

///|
/// Pop a key from uns (get and remove).
pub fn anndata_uns_pop(adata : AnnData, key : String) -> JsonValue? {
  match adata.get_uns(key) {
    Some(v) => {
      adata.uns.remove(key)
      Some(v)
    }
    None => None
  }
}
