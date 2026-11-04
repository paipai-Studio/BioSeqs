/// Scanpy-inspired: read_zarr / write_zarr — Zarr format I/O
/// Zarr is a chunked, compressed array storage format.
/// Since MoonBit doesn't have native Zarr, we use JSON as intermediate.

///|
/// Read AnnData from Zarr-like JSON representation.
/// Zarr AnnData structure: "X", "obs", "var", "obsm", "varm", "layers", "uns".
///
/// Parameters:
/// - `json_str`: JSON string with Zarr-AnnData data
pub fn scanpy_read_zarr(json_str : String) -> AnnData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_zarr_json(j)
    None => AnnData::new([], [], [])
  }
}

///|
fn scanpy_parse_zarr_json(j : JsonValue) -> AnnData {
  let x : Array[Array[Double]] = []
  let obs_names : Array[String] = []
  let var_names : Array[String] = []
  let obs : Map[String, Array[String]] = Map([], capacity=0)
  let var_ : Map[String, Array[String]] = Map([], capacity=0)
  let obsm : Map[String, Array[Array[Double]]] = Map([], capacity=0)
  let layers : Map[String, Array[Array[Double]]] = Map([], capacity=0)
  match j {
    JsonObj(m) => {
      match m.get("X") {
        Some(JsonArr(rows)) =>
          for row in rows {
            match row {
              JsonArr(cols) => {
                let r : Array[Double] = []
                for v in cols {
                  match v {
                    JsonNum(n) => r.push(n)
                    _ => r.push(0.0)
                  }
                }
                x.push(r)
              }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("obs") {
        Some(JsonObj(om)) =>
          for obs_key in om.keys() {
            match om.get(obs_key) {
              Some(JsonObj(col_map)) =>
                match col_map.get("values") {
                  Some(JsonArr(arr)) => {
                    let vals : Array[String] = []
                    for item in arr {
                      match item {
                        JsonStr(s) => vals.push(s)
                        JsonNum(n) => vals.push(n.to_string())
                        _ => vals.push("")
                      }
                    }
                    obs[obs_key] = vals
                  }
                  _ => ()
                }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("var") {
        Some(JsonObj(vm)) =>
          for var_key in vm.keys() {
            match vm.get(var_key) {
              Some(JsonObj(col_map)) =>
                match col_map.get("values") {
                  Some(JsonArr(arr)) => {
                    let vals : Array[String] = []
                    for item in arr {
                      match item {
                        JsonStr(s) => vals.push(s)
                        JsonNum(n) => vals.push(n.to_string())
                        _ => vals.push("")
                      }
                    }
                    var_[var_key] = vals
                  }
                  _ => ()
                }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("obsm") {
        Some(JsonObj(om)) =>
          for obsm_key in om.keys() {
            match om.get(obsm_key) {
              Some(JsonArr(arr)) => {
                let mat : Array[Array[Double]] = []
                for row in arr {
                  match row {
                    JsonArr(cols) => {
                      let r : Array[Double] = []
                      for v in cols {
                        match v {
                          JsonNum(n) => r.push(n)
                          _ => r.push(0.0)
                        }
                      }
                      mat.push(r)
                    }
                    _ => ()
                  }
                }
                obsm[obsm_key] = mat
              }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("layers") {
        Some(JsonObj(lm)) =>
          for layer_key in lm.keys() {
            match lm.get(layer_key) {
              Some(JsonArr(arr)) => {
                let mat : Array[Array[Double]] = []
                for row in arr {
                  match row {
                    JsonArr(cols) => {
                      let r : Array[Double] = []
                      for v in cols {
                        match v {
                          JsonNum(n) => r.push(n)
                          _ => r.push(0.0)
                        }
                      }
                      mat.push(r)
                    }
                    _ => ()
                  }
                }
                layers[layer_key] = mat
              }
              _ => ()
            }
          }
        _ => ()
      }
    }
    _ => ()
  }
  match obs.get("_index") {
    Some(v) =>
      for s in v {
        obs_names.push(s)
      }
    None => ()
  }
  match var_.get("_index") {
    Some(v) =>
      for s in v {
        var_names.push(s)
      }
    None => ()
  }
  if x.length() > 0 {
    if obs_names.length() == 0 {
      for i in 0..<x.length() {
        obs_names.push("cell_" + i.to_string())
      }
    }
    if var_names.length() == 0 {
      for j in 0..<x[0].length() {
        var_names.push("gene_" + j.to_string())
      }
    }
  }
  let adata = AnnData::new(x, obs_names, var_names)
  for obs_key in obs.keys() {
    if obs_key != "_index" {
      match obs.get(obs_key) {
        Some(v) => ignore(adata.set_obs(obs_key, v))
        None => ()
      }
    }
  }
  for var_key in var_.keys() {
    if var_key != "_index" {
      match var_.get(var_key) {
        Some(v) => ignore(adata.set_var(var_key, v))
        None => ()
      }
    }
  }
  for obsm_key in obsm.keys() {
    match obsm.get(obsm_key) {
      Some(v) => ignore(adata.set_obsm(obsm_key, v))
      None => ()
    }
  }
  for layer_key in layers.keys() {
    match layers.get(layer_key) {
      Some(v) => ignore(adata.set_layer(layer_key, v))
      None => ()
    }
  }
  adata
}

///|
/// Write AnnData to Zarr-like JSON format.
pub fn scanpy_write_zarr(adata : AnnData) -> String {
  let root : Map[String, JsonValue] = Map([], capacity=0)
  let x_rows : Array[JsonValue] = []
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  for i in 0..<n_obs {
    let row : Array[JsonValue] = []
    for j in 0..<n_vars {
      row.push(JsonNum(adata.get_x(i, j)))
    }
    x_rows.push(JsonArr(row))
  }
  root["X"] = JsonArr(x_rows)
  let obs_json : Map[String, JsonValue] = Map([], capacity=0)
  let obs_idx_map : Map[String, JsonValue] = Map([], capacity=0)
  obs_idx_map["values"] = JsonArr(
    Array::from_iter(adata.get_obs_names().iter().map(fn(s) { JsonStr(s) })),
  )
  obs_json["_index"] = JsonObj(obs_idx_map)
  root["obs"] = JsonObj(obs_json)
  let var_json : Map[String, JsonValue] = Map([], capacity=0)
  let var_idx_map : Map[String, JsonValue] = Map([], capacity=0)
  var_idx_map["values"] = JsonArr(
    Array::from_iter(adata.get_var_names().iter().map(fn(s) { JsonStr(s) })),
  )
  var_json["_index"] = JsonObj(var_idx_map)
  root["var"] = JsonObj(var_json)
  root["shape"] = JsonArr([
    JsonNum(adata.n_obs().to_double()),
    JsonNum(adata.n_vars().to_double()),
  ])
  json_to_string(JsonObj(root))
}
