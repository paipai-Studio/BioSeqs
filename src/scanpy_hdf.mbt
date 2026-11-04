/// Scanpy-inspired: read_hdf — generic HDF5 file reading
/// Since MoonBit doesn't have native HDF5, this reads a JSON representation
/// of HDF5 structure with groups and datasets.

///|
/// Read a generic HDF5-like JSON structure as AnnData.
/// Expected structure: "/X" (matrix), "/obs" (cell metadata),
/// "/var" (gene metadata), "/obs_names", "/var_names".
///
/// Parameters:
/// - `json_str`: JSON string with HDF5-like data
/// - `key`: key within HDF5 for the AnnData group (default "")
pub fn scanpy_read_hdf(json_str : String, key~ : String = "") -> AnnData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_hdf_json(j, key)
    None => AnnData::new([], [], [])
  }
}

///|
fn scanpy_parse_hdf_json(j : JsonValue, key : String) -> AnnData {
  let root = if key.length() > 0 {
    match j {
      JsonObj(m) =>
        match m.get(key) {
          Some(v) => v
          None => j
        }
      _ => j
    }
  } else {
    j
  }
  let x : Array[Array[Double]] = []
  let obs_names : Array[String] = []
  let var_names : Array[String] = []
  let obs : Map[String, Array[String]] = Map([], capacity=0)
  let var_ : Map[String, Array[String]] = Map([], capacity=0)
  let uns : Map[String, JsonValue] = Map([], capacity=0)
  match root {
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
      match m.get("obs_names") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => obs_names.push(s)
              _ => obs_names.push("")
            }
          }
        _ => ()
      }
      match m.get("var_names") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => var_names.push(s)
              _ => var_names.push("")
            }
          }
        _ => ()
      }
      match m.get("obs") {
        Some(JsonObj(om)) =>
          for obs_key in om.keys() {
            match om.get(obs_key) {
              Some(JsonArr(arr)) => {
                let vals : Array[String] = []
                for item in arr {
                  match item {
                    JsonStr(s) => vals.push(s)
                    JsonNum(n) => vals.push(n.to_string())
                    JsonBool(b) => vals.push(if b { "true" } else { "false" })
                    _ => vals.push("")
                  }
                }
                if vals.length() > 0 {
                  obs[obs_key] = vals
                }
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
              Some(JsonArr(arr)) => {
                let vals : Array[String] = []
                for item in arr {
                  match item {
                    JsonStr(s) => vals.push(s)
                    JsonNum(n) => vals.push(n.to_string())
                    JsonBool(b) => vals.push(if b { "true" } else { "false" })
                    _ => vals.push("")
                  }
                }
                if vals.length() > 0 {
                  var_[var_key] = vals
                }
              }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("uns") {
        Some(JsonObj(um)) =>
          for uns_key in um.keys() {
            match um.get(uns_key) {
              Some(v) => uns[uns_key] = v
              None => ()
            }
          }
        _ => ()
      }
    }
    _ => ()
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
    match obs.get(obs_key) {
      Some(v) => ignore(adata.set_obs(obs_key, v))
      None => ()
    }
  }
  for var_key in var_.keys() {
    match var_.get(var_key) {
      Some(v) => ignore(adata.set_var(var_key, v))
      None => ()
    }
  }
  for uns_key in uns.keys() {
    match uns.get(uns_key) {
      Some(v) => ignore(adata.set_uns(uns_key, v))
      None => ()
    }
  }
  adata
}

///|
/// Write AnnData to HDF5-like JSON format.
pub fn scanpy_write_hdf(adata : AnnData) -> String {
  let root : Map[String, JsonValue] = Map([], capacity=0)
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let x_rows : Array[JsonValue] = []
  for i in 0..<n_obs {
    let row : Array[JsonValue] = []
    for j in 0..<n_vars {
      row.push(JsonNum(adata.get_x(i, j)))
    }
    x_rows.push(JsonArr(row))
  }
  root["X"] = JsonArr(x_rows)
  root["obs_names"] = JsonArr(
    Array::from_iter(adata.get_obs_names().iter().map(fn(s) { JsonStr(s) })),
  )
  root["var_names"] = JsonArr(
    Array::from_iter(adata.get_var_names().iter().map(fn(s) { JsonStr(s) })),
  )
  root["shape"] = JsonArr([
    JsonNum(adata.n_obs().to_double()),
    JsonNum(adata.n_vars().to_double()),
  ])
  json_to_string(JsonObj(root))
}
