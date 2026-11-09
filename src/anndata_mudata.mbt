/// anndata_mudata.mbt — multimodal container (MuData / mudata analogue).
///
/// A MuData bundles several AnnData modalities (e.g. RNA + ADT + ATAC)
/// that share (a subset of) observations, plus global obsm / varm / uns
/// slots — mirroring the mudata ecosystem's `MuData` object and its
/// `.h5mu` storage.

///|
/// Multimodal container over several AnnData modalities.
pub struct MuData {
  mods : Map[String, AnnData]
  obs_names : Array[String] // global (joined) observation names
  obsm : Map[String, Array[Array[Double]]] // global embeddings
  varm : Map[String, Array[Array[Double]]] // global feature embeddings
  uns : Map[String, JsonValue] // global unstructured annotation
}

///|
/// Build a MuData from modality AnnDatas. `join` selects the global
/// observation set: "inner" (intersection, order follows the first mod)
/// or "outer" (union, first-seen order).
pub fn mudata_new(
  mods : Map[String, AnnData],
  join~ : String = "inner",
) -> MuData {
  let keys = mods.keys().collect()
  let outer : Array[String] = Array(capacity=0)
  let seen : Set[String] = Set([])
  for k in keys {
    let m = mods[k]
    for n in m.obs_names {
      if not(seen.contains(n)) {
        seen.add(n)
        outer.push(n)
      }
    }
  }
  let global : Array[String] = if join == "outer" || keys.length() == 0 {
    outer
  } else {
    // inner: keep names present in every mod, in first-mod order
    let sets : Array[Set[String]] = Array(capacity=keys.length())
    for k in keys {
      let s : Set[String] = Set([])
      for n in mods[k].obs_names {
        s.add(n)
      }
      sets.push(s)
    }
    let out : Array[String] = Array(capacity=0)
    for n in outer {
      let mut in_all = true
      for s in sets {
        if not(s.contains(n)) {
          in_all = false
          break
        }
      }
      if in_all {
        out.push(n)
      }
    }
    out
  }
  MuData::{
    mods,
    obs_names: global,
    obsm: Map([], capacity=4),
    varm: Map([], capacity=4),
    uns: Map([], capacity=4),
  }
}

///|
/// Number of global observations.
pub fn MuData::n_obs(self : MuData) -> Int32 {
  self.obs_names.length()
}

///|
/// Number of modalities.
pub fn MuData::n_mods(self : MuData) -> Int32 {
  self.mods.keys().collect().length()
}

///|
/// Modality names (unsorted map order).
pub fn MuData::mod_names(self : MuData) -> Array[String] {
  self.mods.keys().collect()
}

///|
/// A modality by name, or `None`.
pub fn MuData::mod(self : MuData, name : String) -> AnnData? {
  self.mods.get(name)
}

///|
/// Variable names of a modality (empty when the mod is absent).
pub fn MuData::mod_var_names(self : MuData, name : String) -> Array[String] {
  match self.mods.get(name) {
    Some(m) => m.var_names
    None => []
  }
}

///|
/// Variable count of a modality.
pub fn MuData::mod_n_vars(self : MuData, name : String) -> Int32 {
  match self.mods.get(name) {
    Some(m) => m.n_vars()
    None => 0
  }
}

///|
/// Set a global embedding matrix. Returns `self` for chaining.
pub fn MuData::set_global_obsm(
  self : MuData,
  name : String,
  mat : Array[Array[Double]],
) -> MuData {
  self.obsm[name] = mat
  self
}

///|
/// Global embedding matrix, or `None`.
pub fn MuData::get_global_obsm(
  self : MuData,
  name : String,
) -> Array[Array[Double]]? {
  self.obsm.get(name)
}

///|
/// Global obsm keys.
pub fn MuData::global_obsm_keys(self : MuData) -> Array[String] {
  self.obsm.keys().collect()
}

///|
/// Set a global var embedding matrix. Returns `self` for chaining.
pub fn MuData::set_global_varm(
  self : MuData,
  name : String,
  mat : Array[Array[Double]],
) -> MuData {
  self.varm[name] = mat
  self
}

///|
/// Global var embedding matrix, or `None`.
pub fn MuData::get_global_varm(
  self : MuData,
  name : String,
) -> Array[Array[Double]]? {
  self.varm.get(name)
}

///|
/// Set a global unstructured annotation. Returns `self` for chaining.
pub fn MuData::set_global_uns(
  self : MuData,
  key : String,
  value : JsonValue,
) -> MuData {
  self.uns[key] = value
  self
}

///|
/// Global unstructured annotation value, or `None`.
pub fn MuData::get_global_uns(self : MuData, key : String) -> JsonValue? {
  self.uns.get(key)
}

///|
/// Global uns keys.
pub fn MuData::uns_keys(self : MuData) -> Array[String] {
  self.uns.keys().collect()
}

///|
/// Materialize the MuData into a single AnnData by column-block
/// concatenation of the modality X matrices (rows aligned by the global
/// obs names; modalities missing an observation contribute zeros, and
/// `mod:<name>` is prefixed to each var name like mudata does).
pub fn MuData::to_adata(self : MuData) -> AnnData {
  let keys = self.mods.keys().collect()
  // build per-mod row lookup
  let lookups : Array[Map[String, Int32]] = Array(capacity=keys.length())
  for k in keys {
    let m : Map[String, Int32] = Map([], capacity=self.mods[k].n_obs())
    for i = 0; i < self.mods[k].n_obs(); i = i + 1 {
      m[self.mods[k].obs_names[i]] = i
    }
    lookups.push(m)
  }
  // global var list: mod:name blocks in map order
  let var_names : Array[String] = Array(capacity=0)
  for k in keys {
    for v in self.mods[k].var_names {
      var_names.push("mod:" + k + ":" + v)
    }
  }
  let x : Array[Array[Double]] = Array(capacity=self.obs_names.length())
  for i = 0; i < self.obs_names.length(); i = i + 1 {
    let row : Array[Double] = Array(capacity=var_names.length())
    for mk = 0; mk < keys.length(); mk = mk + 1 {
      let mod_data = self.mods[keys[mk]]
      let ax = mod_data.x_matrix()
      let src_row = match lookups[mk].get(self.obs_names[i]) {
        Some(r) => r
        None => -1
      }
      for j = 0; j < mod_data.n_vars(); j = j + 1 {
        if src_row >= 0 && src_row < ax.length() && j < ax[src_row].length() {
          row.push(ax[src_row][j])
        } else {
          row.push(0.0)
        }
      }
    }
    x.push(row)
  }
  let out = AnnData::new(x, self.obs_names.copy(), var_names)
  // global slots carried over
  for k in self.obsm.keys().collect() {
    out.obsm[k] = self.obsm[k]
  }
  for k in self.varm.keys().collect() {
    out.varm[k] = self.varm[k]
  }
  for k in self.uns.keys().collect() {
    out.uns[k] = self.uns[k]
  }
  out
}

///|
/// Serialize a MuData to a text (.h5mu analogue): header line, global
/// uns/obsm/varm as JSON lines, mod count, then per mod a name line
/// followed by its single-line h5ad JSON.
pub fn mudata_write_json(mdata : MuData) -> String {
  let keys = mdata.mods.keys().collect()
  let sb = StringBuilder()
  let _ = sb.write("mudata 0.1.0\n")
  // global slots as JSON objects
  let uns_obj : Map[String, JsonValue] = Map(
    [],
    capacity=mdata.uns.keys().collect().length(),
  )
  for k in mdata.uns.keys().collect() {
    uns_obj[k] = mdata.uns[k]
  }
  let _ = sb.write("uns ")
  let _ = sb.write(json_to_string(JsonObj(uns_obj)))
  let _ = sb.write("\n")
  let dict_to_json = fn(m : Map[String, Array[Array[Double]]]) -> JsonValue {
    let obj : Map[String, JsonValue] = Map([], capacity=4)
    for k in m.keys().collect() {
      let rows : Array[JsonValue] = Array(capacity=m[k].length())
      for row in m[k] {
        let cols : Array[JsonValue] = Array(capacity=row.length())
        for v in row {
          cols.push(JsonNum(v))
        }
        rows.push(JsonArr(cols))
      }
      obj[k] = JsonArr(rows)
    }
    JsonObj(obj)
  }
  let _ = sb.write("obsm ")
  let _ = sb.write(json_to_string(dict_to_json(mdata.obsm)))
  let _ = sb.write("\n")
  let _ = sb.write("varm ")
  let _ = sb.write(json_to_string(dict_to_json(mdata.varm)))
  let _ = sb.write("\n")
  let _ = sb.write("mods ")
  let _ = sb.write(keys.length().to_string())
  let _ = sb.write("\n")
  for k in keys {
    let _ = sb.write("mod ")
    let _ = sb.write(k)
    let _ = sb.write("\n")
    let _ = sb.write(scanpy_write_h5ad(mdata.mods[k]))
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
/// Parse a MuData written by `mudata_write_json`.
pub fn mudata_read_json(content : String) -> MuData? {
  let lines_raw = content.split("\n").to_array()
  if lines_raw.length() < 5 || lines_raw[0].to_owned() != "mudata 0.1.0" {
    return None
  }
  // global uns
  let uns_line = lines_raw[1].to_owned().split(" ").to_array()
  if uns_line.length() < 2 || uns_line[0].to_owned() != "uns" {
    return None
  }
  let uns_tmp : Map[String, JsonValue] = Map([], capacity=4)
  match parse_json(uns_line[1].to_owned()) {
    Some(JsonObj(m)) =>
      for k in m.keys().collect() {
        uns_tmp[k] = m[k]
      }
    _ => ()
  }
  // global obsm / varm
  let read_dict_line = fn(
    line : String,
    target : Map[String, Array[Array[Double]]],
  ) -> Bool {
    let parts = line.split(" ").to_array()
    if parts.length() < 2 ||
      (parts[0].to_owned() != "obsm" && parts[0].to_owned() != "varm") {
      return false
    }
    match parse_json(parts[1].to_owned()) {
      Some(JsonObj(m)) => {
        for k in m.keys().collect() {
          match m[k] {
            JsonArr(rows) => {
              let mat : Array[Array[Double]] = Array(capacity=rows.length())
              for row in rows {
                match row {
                  JsonArr(cols) => {
                    let r : Array[Double] = Array(capacity=cols.length())
                    for c in cols {
                      match c {
                        JsonNum(v) => r.push(v)
                        _ => r.push(0.0)
                      }
                    }
                    mat.push(r)
                  }
                  _ => mat.push([])
                }
              }
              target[k] = mat
            }
            _ => ()
          }
        }
        true
      }
      _ => false
    }
  }
  let obsm_tmp : Map[String, Array[Array[Double]]] = Map([], capacity=4)
  if not(read_dict_line(lines_raw[2].to_owned(), obsm_tmp)) {
    return None
  }
  let varm_tmp : Map[String, Array[Array[Double]]] = Map([], capacity=4)
  if not(read_dict_line(lines_raw[3].to_owned(), varm_tmp)) {
    return None
  }
  // mods
  let header = lines_raw[4].to_owned().split(" ").to_array()
  if header.length() < 2 || header[0].to_owned() != "mods" {
    return None
  }
  let mods : Map[String, AnnData] = Map([], capacity=4)
  let mut i = 5
  while i + 1 < lines_raw.length() {
    let tag = lines_raw[i].to_owned().split(" ").to_array()
    if tag.length() < 2 || tag[0].to_owned() != "mod" {
      return None
    }
    let name = tag[1].to_owned()
    match scanpy_read_h5ad(lines_raw[i + 1].to_owned()) {
      Some(a) => mods[name] = a
      None => return None
    }
    i = i + 2
  }
  let mdata = mudata_new(mods)
  for k in uns_tmp.keys().collect() {
    mdata.uns[k] = uns_tmp[k]
  }
  for k in obsm_tmp.keys().collect() {
    mdata.obsm[k] = obsm_tmp[k]
  }
  for k in varm_tmp.keys().collect() {
    mdata.varm[k] = varm_tmp[k]
  }
  Some(mdata)
}
