/// Scanpy-inspired: read_loom / write_loom — Loom file format I/O
/// Loom is a HDF5-based format for genomic data.
/// Since MoonBit doesn't have native HDF5, we use JSON as intermediate format.

///|
/// Read a Loom file from JSON representation.
/// Loom structure: "matrix" (n_obs × n_vars), "row_attrs" (gene metadata),
/// "col_attrs" (cell metadata), "row_graphs", "col_graphs".
///
/// Parameters:
/// - `json_str`: JSON string with Loom data
/// - `transpose` : if true, swap rows/cols (default false; loom is gene × cell)
pub fn scanpy_read_loom(
  json_str : String,
  transpose~ : Bool = false,
) -> AnnData {
  let json = parse_json(json_str)
  match json {
    Some(j) => scanpy_parse_loom_json(j, transpose)
    None => AnnData::new([], [], [])
  }
}

///|
fn scanpy_parse_loom_json(j : JsonValue, transpose : Bool) -> AnnData {
  let x : Array[Array[Double]] = []
  let obs_names : Array[String] = []
  let var_names : Array[String] = []
  let obs : Map[String, Array[String]] = Map([], capacity=0)
  let var_ : Map[String, Array[String]] = Map([], capacity=0)
  match j {
    JsonObj(m) => {
      match m.get("matrix") {
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
      match m.get("row_attrs") {
        Some(JsonObj(ra)) =>
          for ra_key in ra.keys() {
            match ra.get(ra_key) {
              Some(JsonArr(arr)) => {
                let vals : Array[String] = []
                for item in arr {
                  match item {
                    JsonStr(s) => vals.push(s)
                    JsonNum(n) => vals.push(n.to_string())
                    _ => vals.push("")
                  }
                }
                var_[ra_key] = vals
              }
              _ => ()
            }
          }
        _ => ()
      }
      match m.get("col_attrs") {
        Some(JsonObj(ca)) =>
          for ca_key in ca.keys() {
            match ca.get(ca_key) {
              Some(JsonArr(arr)) => {
                let vals : Array[String] = []
                for item in arr {
                  match item {
                    JsonStr(s) => vals.push(s)
                    JsonNum(n) => vals.push(n.to_string())
                    _ => vals.push("")
                  }
                }
                obs[ca_key] = vals
              }
              _ => ()
            }
          }
        _ => ()
      }
    }
    _ => ()
  }
  match var_.get("Accession") {
    Some(v) =>
      for s in v {
        var_names.push(s)
      }
    None =>
      match var_.get("Gene") {
        Some(v) =>
          for s in v {
            var_names.push(s)
          }
        None => ()
      }
  }
  match obs.get("CellID") {
    Some(v) =>
      for s in v {
        obs_names.push(s)
      }
    None => ()
  }
  if x.length() > 0 {
    if var_names.length() == 0 {
      for j in 0..<x[0].length() {
        var_names.push("gene_" + j.to_string())
      }
    }
    if obs_names.length() == 0 {
      for i in 0..<x.length() {
        obs_names.push("cell_" + i.to_string())
      }
    }
  }
  if !transpose && x.length() > 0 && x[0].length() > 0 {
    let n_genes = x.length()
    let n_cells = x[0].length()
    let transposed : Array[Array[Double]] = []
    for i in 0..<n_cells {
      let row : Array[Double] = []
      for j in 0..<n_genes {
        row.push(x[j][i])
      }
      transposed.push(row)
    }
    let adata = AnnData::new(transposed, obs_names, var_names)
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
    return adata
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
  adata
}

///|
/// Write AnnData to Loom JSON format.
/// Loom stores as gene × cell (transposed from AnnData's cell × gene).
pub fn scanpy_write_loom(adata : AnnData) -> String {
  let n_genes = adata.n_vars()
  let n_cells = adata.n_obs()
  let rows : Array[JsonValue] = []
  for j in 0..<n_genes {
    let row : Array[JsonValue] = []
    for i in 0..<n_cells {
      row.push(JsonNum(adata.get_x(i, j)))
    }
    rows.push(JsonArr(row))
  }
  let row_attrs : Map[String, JsonValue] = Map([], capacity=0)
  row_attrs["Accession"] = JsonArr(
    Array::from_iter(adata.get_var_names().iter().map(fn(s) { JsonStr(s) })),
  )
  let col_attrs : Map[String, JsonValue] = Map([], capacity=0)
  col_attrs["CellID"] = JsonArr(
    Array::from_iter(adata.get_obs_names().iter().map(fn(s) { JsonStr(s) })),
  )
  let root : Map[String, JsonValue] = Map([], capacity=0)
  root["matrix"] = JsonArr(rows)
  root["row_attrs"] = JsonObj(row_attrs)
  root["col_attrs"] = JsonObj(col_attrs)
  root["shape"] = JsonArr([
    JsonNum(n_genes.to_double()),
    JsonNum(n_cells.to_double()),
  ])
  json_to_string(JsonObj(root))
}
