///|
/// BIOM 1.0 (Biological Observation Matrix, JSON) table — parse and write.
/// Mirrors scikit-bio's `skbio.io.biom` `Table`. BIOM 2.1 (HDF5) is out of
/// scope (no HDF5 dependency).

///|
pub struct BiomTable {
  table_id : String
  matrix_type : String
  matrix_element_type : String
  observation_ids : Array[String]
  sample_ids : Array[String]
  observation_metadata : Array[Map[String, JsonValue]]
  sample_metadata : Array[Map[String, JsonValue]]
  data : Array[Array[Double]]
}

///|
fn json_as_str(v : JsonValue) -> String? {
  match v {
    JsonStr(s) => Some(s)
    _ => None
  }
}

///|
fn json_as_num(v : JsonValue) -> Double? {
  match v {
    JsonNum(d) => Some(d)
    _ => None
  }
}

///|
fn json_as_obj(v : JsonValue) -> Map[String, JsonValue]? {
  match v {
    JsonObj(m) => Some(m)
    _ => None
  }
}

///|
fn json_as_arr(v : JsonValue) -> Array[JsonValue]? {
  match v {
    JsonArr(a) => Some(a)
    _ => None
  }
}

///|
/// Extract (id, metadata) from a BIOM row/column entry.
fn biom_id_metadata(
  entry : JsonValue,
) -> (String, Map[String, JsonValue])? {
  let m = match json_as_obj(entry) {
    Some(x) => x
    None => return None
  }
  let id = match m.get("id") {
    Some(x) => {
      match json_as_str(x) {
        Some(s) => s
        None => return None
      }
    }
    None => return None
  }
  let meta = match m.get("metadata") {
    Some(JsonObj(mm)) => mm
    _ => Map([], capacity=0)
  }
  Some((id, meta))
}

///|
/// Parse a BIOM 1.0 JSON document into a `BiomTable`.
pub fn parse_biom(text : String) -> BiomTable? {
  let root = match parse_json(text) {
    Some(x) => x
    None => return None
  }
  let obj = match json_as_obj(root) {
    Some(x) => x
    None => return None
  }

  let table_id = match obj.get("id") {
    Some(x) => json_as_str(x).unwrap_or("")
    None => ""
  }
  let matrix_type = match obj.get("matrix_type") {
    Some(x) => json_as_str(x).unwrap_or("dense")
    None => "dense"
  }
  let matrix_element_type = match obj.get("matrix_element_type") {
    Some(x) => json_as_str(x).unwrap_or("float")
    None => "float"
  }

  let shape = match json_as_arr(obj.get("shape").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }
  let n_rows = json_as_num(shape[0]).unwrap_or(0.0).to_int()
  let n_cols = json_as_num(shape[1]).unwrap_or(0.0).to_int()

  let rows_arr = match json_as_arr(obj.get("rows").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }
  let cols_arr = match json_as_arr(obj.get("columns").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }

  let observation_ids : Array[String] = []
  let observation_metadata : Array[Map[String, JsonValue]] = []
  for r in rows_arr {
    match biom_id_metadata(r) {
      Some((id, meta)) => {
        observation_ids.push(id)
        observation_metadata.push(meta)
      }
      None => return None
    }
  }
  let sample_ids : Array[String] = []
  let sample_metadata : Array[Map[String, JsonValue]] = []
  for c in cols_arr {
    match biom_id_metadata(c) {
      Some((id, meta)) => {
        sample_ids.push(id)
        sample_metadata.push(meta)
      }
      None => return None
    }
  }

  // data: 稠密或稀疏
  let dens : Array[Array[Double]] = Array(capacity=n_rows)
  for i = 0; i < n_rows; i = i + 1 {
    dens.push(Array::make(n_cols, 0.0))
  }
  let data_node = obj.get("data").unwrap_or(JsonNull)
  if matrix_type == "sparse" {
    let entries = match json_as_arr(data_node) {
      Some(x) => x
      None => return None
    }
    for e in entries {
      let t = match json_as_arr(e) {
        Some(x) => x
        None => return None
      }
      let r = json_as_num(t[0]).unwrap_or(0.0).to_int()
      let c = json_as_num(t[1]).unwrap_or(0.0).to_int()
      let v = json_as_num(t[2]).unwrap_or(0.0)
      if r >= 0 && r < n_rows && c >= 0 && c < n_cols {
        dens[r][c] = v
      }
    }
  } else {
    let rows_data = match json_as_arr(data_node) {
      Some(x) => x
      None => return None
    }
    for i = 0; i < n_rows; i = i + 1 {
      let vals = match json_as_arr(rows_data[i]) {
        Some(x) => x
        None => return None
      }
      for j = 0; j < n_cols; j = j + 1 {
        dens[i][j] = json_as_num(vals[j]).unwrap_or(0.0)
      }
    }
  }

  Some(BiomTable::{
    table_id,
    matrix_type,
    matrix_element_type,
    observation_ids,
    sample_ids,
    observation_metadata,
    sample_metadata,
    data: dens,
  })
}

///|
/// Serialize a `BiomTable` back to BIOM 1.0 JSON (compact).
pub fn BiomTable::to_json(self : BiomTable) -> String {
  let n_rows = self.observation_ids.length()
  let n_cols = self.sample_ids.length()

  let rows : Array[JsonValue] = []
  for i = 0; i < n_rows; i = i + 1 {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m.set("id", JsonStr(self.observation_ids[i]))
    let meta = self.observation_metadata[i]
    if meta.length() > 0 {
      m.set("metadata", JsonObj(meta))
    }
    rows.push(JsonObj(m))
  }
  let cols : Array[JsonValue] = []
  for j = 0; j < n_cols; j = j + 1 {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m.set("id", JsonStr(self.sample_ids[j]))
    let meta = self.sample_metadata[j]
    if meta.length() > 0 {
      m.set("metadata", JsonObj(meta))
    }
    cols.push(JsonObj(m))
  }

  let data_node : JsonValue =
    if self.matrix_type == "sparse" {
      let entries : Array[JsonValue] = []
      for i = 0; i < n_rows; i = i + 1 {
        for j = 0; j < n_cols; j = j + 1 {
          if self.data[i][j] != 0.0 {
            entries.push(
              JsonArr([
                JsonNum(i.to_double()),
                JsonNum(j.to_double()),
                JsonNum(self.data[i][j]),
              ]),
            )
          }
        }
      }
      JsonArr(entries)
    } else {
      let rows_d : Array[JsonValue] = []
      for i = 0; i < n_rows; i = i + 1 {
        let vals : Array[JsonValue] = []
        for j = 0; j < n_cols; j = j + 1 {
          vals.push(JsonNum(self.data[i][j]))
        }
        rows_d.push(JsonArr(vals))
      }
      JsonArr(rows_d)
    }

  let root : Map[String, JsonValue] = Map([], capacity=8)
  root.set("id", JsonStr(self.table_id))
  root.set("format", JsonStr("Biological Observation Matrix 1.0.0"))
  root.set("matrix_type", JsonStr(self.matrix_type))
  root.set("matrix_element_type", JsonStr(self.matrix_element_type))
  root.set(
    "shape",
    JsonArr([JsonNum(n_rows.to_double()), JsonNum(n_cols.to_double())]),
  )
  root.set("rows", JsonArr(rows))
  root.set("columns", JsonArr(cols))
  root.set("data", data_node)

  json_to_string(JsonObj(root))
}