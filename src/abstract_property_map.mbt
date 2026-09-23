///|
/// abstract_property_map - Unified residue property access interface.
/// Inspired by Biopython Bio.PDB.AbstractPropertyMap.
/// Provides a generic container for per-residue properties (DSSP, HSExposure,
/// ResidueDepth, NACCESS, SASA) indexed by (chain_id, resseq, icode).

///|
/// Key: (chain_id, resseq, icode) uniquely identifies a residue.
pub struct ResidueKey {
  chain_id : String
  resseq : Int
  icode : String
}

///|
/// A single residue property record.
pub struct PropertyRecord {
  key : ResidueKey
  props : Map[String, String]
}

///|
/// Abstract property map: dictionary of residue properties.
pub struct AbstractPropertyMap {
  records : Map[String, PropertyRecord]
  property_names : Array[String]
  mut n_residues : Int
}

///|
/// Construct a residue key.
pub fn ResidueKey::new(chain_id : String, resseq : Int, icode~ : String = "") -> ResidueKey {
  ResidueKey::{ chain_id, resseq, icode }
}

///|
/// String key for hashing: "chain:resseq:icode".
pub fn residue_key_to_string(key : ResidueKey) -> String {
  key.chain_id + ":" + key.resseq.to_string() + ":" + key.icode
}

///|
/// Construct an empty property map with known property names.
pub fn AbstractPropertyMap::new(property_names : Array[String]) -> AbstractPropertyMap {
  AbstractPropertyMap::{ records: Map([], capacity=64), property_names, n_residues: 0 }
}

///|
/// Add a property record for a residue.
pub fn AbstractPropertyMap::add(
  self : AbstractPropertyMap,
  chain_id : String,
  resseq : Int,
  props : Map[String, String],
  icode~ : String = ""
) -> Unit {
  let key = ResidueKey::new(chain_id, resseq, icode~)
  let skey = residue_key_to_string(key)
  let record = PropertyRecord::{ key, props }
  if !(self.records.contains(skey)) {
    self.n_residues = self.n_residues + 1
  }
  self.records[skey] = record
}

///|
/// Look up a single property value for a residue.
pub fn AbstractPropertyMap::get_property(
  self : AbstractPropertyMap,
  chain_id : String,
  resseq : Int,
  prop_name : String,
  icode~ : String = ""
) -> String? {
  let key = ResidueKey::new(chain_id, resseq, icode~)
  let skey = residue_key_to_string(key)
  match self.records.get(skey) {
    Some(record) => record.props.get(prop_name)
    None => None
  }
}

///|
/// Look up all properties for a residue.
pub fn AbstractPropertyMap::get_residue(
  self : AbstractPropertyMap,
  chain_id : String,
  resseq : Int,
  icode~ : String = ""
) -> Map[String, String]? {
  let key = ResidueKey::new(chain_id, resseq, icode~)
  let skey = residue_key_to_string(key)
  match self.records.get(skey) {
    Some(record) => Some(record.props)
    None => None
  }
}

///|
/// Check if a residue exists in the map.
pub fn AbstractPropertyMap::has_residue(
  self : AbstractPropertyMap,
  chain_id : String,
  resseq : Int,
  icode~ : String = ""
) -> Bool {
  let key = ResidueKey::new(chain_id, resseq, icode~)
  self.records.contains(residue_key_to_string(key))
}

///|
/// Iterate all residue keys in the map.
pub fn AbstractPropertyMap::keys(self : AbstractPropertyMap) -> Array[ResidueKey] {
  let result : Array[ResidueKey] = []
  for m in self.records.values() {
    result.push(m.key)
  }
  result
}

///|
/// Get all records as an array.
pub fn AbstractPropertyMap::all_records(self : AbstractPropertyMap) -> Array[PropertyRecord] {
  let result : Array[PropertyRecord] = []
  for r in self.records.values() {
    result.push(r)
  }
  result
}

///|
/// Filter residues by a predicate on a specific property.
pub fn AbstractPropertyMap::filter_by_property(
  self : AbstractPropertyMap,
  prop_name : String,
  predicate : (String) -> Bool
) -> AbstractPropertyMap {
  let result = AbstractPropertyMap::new(self.property_names)
  for record in self.records.values() {
    match record.props.get(prop_name) {
      Some(val) => {
        if predicate(val) {
          let props_copy : Map[String, String] = Map([], capacity=record.props.length())
          for k in record.props.keys() {
            match record.props.get(k) {
              Some(v) => props_copy[k] = v
              None => continue
            }
          }
          result.add(record.key.chain_id, record.key.resseq, props_copy, icode=record.key.icode)
        }
      }
      None => continue
    }
  }
  result
}

///|
/// Compute summary statistics for a numeric property.
///   Returns (count, min, max, mean, sum).
pub fn AbstractPropertyMap::stats_property(
  self : AbstractPropertyMap,
  prop_name : String
) -> (Int, Double, Double, Double, Double) {
  let mut count = 0
  let mut sum = 0.0
  let mut min_v = 1.0e30
  let mut max_v = -1.0e30
  for record in self.records.values() {
    match record.props.get(prop_name) {
      Some(val) => {
        let v = abstract_property_atof(val)
        count = count + 1
        sum = sum + v
        if v < min_v {
          min_v = v
        }
        if v > max_v {
          max_v = v
        }
      }
      None => continue
    }
  }
  let mean = if count > 0 { sum / count.to_double() } else { 0.0 }
  (count, min_v, max_v, mean, sum)
}

///|
/// Export to TSV (tab-separated values) format.
pub fn AbstractPropertyMap::to_tsv(self : AbstractPropertyMap) -> String {
  let sb = StringBuilder()
  sb.write_string("chain_id\tresseq\ticode")
  let mut i = 0
  while i < self.property_names.length() {
    sb.write_string("\t")
    sb.write_string(self.property_names[i])
    i = i + 1
  }
  sb.write_string("\n")
  for record in self.records.values() {
    sb.write_string(record.key.chain_id)
    sb.write_string("\t")
    sb.write_string(record.key.resseq.to_string())
    sb.write_string("\t")
    sb.write_string(record.key.icode)
    let mut j = 0
    while j < self.property_names.length() {
      sb.write_string("\t")
      match record.props.get(self.property_names[j]) {
        Some(v) => sb.write_string(v)
        None => sb.write_string("NA")
      }
      j = j + 1
    }
    sb.write_string("\n")
  }
  sb.to_string()
}

///|
/// Merge two property maps (second overrides first on conflict).
pub fn AbstractPropertyMap::merge(
  a : AbstractPropertyMap,
  b : AbstractPropertyMap
) -> AbstractPropertyMap {
  let merged_names : Array[String] = []
  let seen : Map[String, Bool] = Map([], capacity=32)
  let mut i = 0
  while i < a.property_names.length() {
    let n = a.property_names[i]
    if !(seen.contains(n)) {
      merged_names.push(n)
      seen[n] = true
    }
    i = i + 1
  }
  let mut j = 0
  while j < b.property_names.length() {
    let n = b.property_names[j]
    if !(seen.contains(n)) {
      merged_names.push(n)
      seen[n] = true
    }
    j = j + 1
  }
  let result = AbstractPropertyMap::new(merged_names)
  for record in a.records.values() {
    let props_copy : Map[String, String] = Map([], capacity=record.props.length())
    for k in record.props.keys() {
      match record.props.get(k) {
        Some(v) => props_copy[k] = v
        None => continue
      }
    }
    result.add(record.key.chain_id, record.key.resseq, props_copy, icode=record.key.icode)
  }
  for record in b.records.values() {
    let skey = residue_key_to_string(record.key)
    match result.records.get(skey) {
      Some(existing) => {
        for k in record.props.keys() {
          match record.props.get(k) {
            Some(v) => existing.props[k] = v
            None => continue
          }
        }
      }
      None => {
        let props_copy : Map[String, String] = Map([], capacity=record.props.length())
        for k in record.props.keys() {
          match record.props.get(k) {
            Some(v) => props_copy[k] = v
            None => continue
          }
        }
        result.add(record.key.chain_id, record.key.resseq, props_copy, icode=record.key.icode)
      }
    }
  }
  result
}

///|
/// Simple string to double.
fn abstract_property_atof(s : String) -> Double {
  let s_trim = s.trim().to_string()
  if s_trim.length() == 0 {
    return 0.0
  }
  let neg = if s_trim[0] == 45 { true } else { false }
  let start = if neg { 1 } else { 0 }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut in_frac = false
  let mut frac_div = 10.0
  let mut k = start
  while k < s_trim.length() {
    let c = s_trim[k]
    if c >= 48 && c <= 57 {
      let digit = (c - 48).to_int().to_double()
      if in_frac {
        frac_part = frac_part + digit / frac_div
        frac_div = frac_div * 10.0
      } else {
        int_part = int_part * 10.0 + digit
      }
    } else if c == 46 {
      in_frac = true
    }
    k = k + 1
  }
  if neg {
    -(int_part + frac_part)
  } else {
    int_part + frac_part
  }
}

///|
/// Example: build a small DSSP-like property map.
pub fn abstract_property_map_example() -> AbstractPropertyMap {
  let pm = AbstractPropertyMap::new(["SS", "ASA", "NHCO"])
  let r1 : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=3)
    m["SS"] = "H"
    m["ASA"] = "75.5"
    m["NHCO"] = "1"
    m
  }
  let r2 : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=3)
    m["SS"] = "E"
    m["ASA"] = "120.0"
    m["NHCO"] = "0"
    m
  }
  let r3 : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=3)
    m["SS"] = "C"
    m["ASA"] = "45.3"
    m["NHCO"] = "1"
    m
  }
  pm.add("A", 1, r1)
  pm.add("A", 2, r2)
  pm.add("A", 3, r3)
  pm
}