/// Scanpy-inspired: read_10x_h5 — read 10x HDF5 format
/// Parses the standard 10x Genomics HDF5 output structure.
/// Since MoonBit doesn't have native HDF5 support, this reads a simplified
/// JSON-based representation of the 10x HDF5 structure.

///|
/// Read 10x data from a JSON representation of HDF5 structure.
/// The JSON should contain: "matrix/data", "matrix/indices", "matrix/indptr",
/// "matrix/shape", "matrix/barcodes", "matrix/features".
///
/// Parameters:
/// - `json_str`: JSON string with 10x data
/// - `genome`: genome key (default "matrix")
pub fn scanpy_read_10x_h5(
  json_str : String,
  _genome~ : String = "matrix",
) -> AnnData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_10x_json(j)
    None => AnnData::new([], [], [])
  }
}

///|
fn scanpy_parse_10x_json(j : JsonValue) -> AnnData {
  // Extract shape
  let (n_obs, n_vars) = match j {
    JsonObj(m) => {
      let shape = match m.get("shape") {
        Some(JsonArr(s)) => {
          let obs = match s[0] {
            JsonNum(n) => n.floor().to_int()
            _ => 0
          }
          let vars = match s[1] {
            JsonNum(n) => n.floor().to_int()
            _ => 0
          }
          (obs, vars)
        }
        _ => (0, 0)
      }
      shape
    }
    _ => (0, 0)
  }
  // Extract barcodes (obs_names)
  let obs_names = Array(capacity=n_obs)
  match j {
    JsonObj(m) =>
      match m.get("barcodes") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => obs_names.push(s)
              _ => obs_names.push("")
            }
          }
        _ => ()
      }
    _ => ()
  }
  // Fill if empty
  if obs_names.length() == 0 {
    for i in 0..<n_obs {
      obs_names.push("cell_" + i.to_string())
    }
  }
  // Extract features (var_names)
  let var_names = Array(capacity=n_vars)
  match j {
    JsonObj(m) =>
      match m.get("features") {
        Some(JsonArr(arr)) =>
          for item in arr {
            match item {
              JsonStr(s) => var_names.push(s)
              _ => var_names.push("")
            }
          }
        _ => ()
      }
    _ => ()
  }
  if var_names.length() == 0 {
    for j in 0..<n_vars {
      var_names.push("gene_" + j.to_string())
    }
  }
  // Extract data (sparse or dense)
  let x = Array(capacity=n_obs)
  for _ in 0..<n_obs {
    let row = Array(capacity=n_vars)
    for _ in 0..<n_vars {
      row.push(0.0)
    }
    x.push(row)
  }
  match j {
    JsonObj(m) =>
      // Try dense "data" first
      match m.get("data") {
        Some(JsonArr(arr)) =>
          // Dense: data is row-major
          for i in 0..<n_obs {
            for j in 0..<n_vars {
              let idx = i * n_vars + j
              if idx < arr.length() {
                match arr[idx] {
                  JsonNum(v) => x[i][j] = v
                  _ => ()
                }
              }
            }
          }
        _ => {
          // Try sparse: data, indices, indptr
          let data = Array(capacity=0)
          let indices = Array(capacity=0)
          let indptr = Array(capacity=0)
          match m.get("sparse_data") {
            Some(JsonArr(arr)) =>
              for item in arr {
                match item {
                  JsonNum(v) => data.push(v)
                  _ => data.push(0.0)
                }
              }
            _ => ()
          }
          match m.get("sparse_indices") {
            Some(JsonArr(arr)) =>
              for item in arr {
                match item {
                  JsonNum(v) => indices.push(v.floor().to_int())
                  _ => indices.push(0)
                }
              }
            _ => ()
          }
          match m.get("sparse_indptr") {
            Some(JsonArr(arr)) =>
              for item in arr {
                match item {
                  JsonNum(v) => indptr.push(v.floor().to_int())
                  _ => indptr.push(0)
                }
              }
            _ => ()
          }
          // Fill sparse data (CSR format)
          if indptr.length() > 0 {
            for i in 0..<n_obs {
              if i + 1 < indptr.length() {
                for k in indptr[i]..<indptr[i + 1] {
                  if k < indices.length() && k < data.length() {
                    let j = indices[k]
                    if j < n_vars {
                      x[i][j] = data[k]
                    }
                  }
                }
              }
            }
          }
        }
      }
    _ => ()
  }
  AnnData::new(x, obs_names, var_names)
}
