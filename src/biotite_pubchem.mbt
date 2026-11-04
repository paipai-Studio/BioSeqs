/// biotite_pubchem.mbt — PubChem compound query URL builder
/// Inspired by biotite.database.pubchem

///|
/// PubChem search type.
pub enum BiotitePubChemSearchType {
  BiotitePubChemByName(String)
  BiotitePubChemBySmiles(String)
  BiotitePubChemByInchiKey(String)
  BiotitePubChemByInchi
  BiotitePubChemByFormula(String)
  BiotitePubChemByCid(Int32)
}

///|
/// Build a PubChem REST API URL for compound search.
pub fn biotite_pubchem_query_url(
  search : BiotitePubChemSearchType,
  properties : Array[String],
) -> String {
  let base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
  let search_part = match search {
    BiotitePubChemByName(name) => "name/" + name + "/JSON"
    BiotitePubChemBySmiles(smiles) => "smiles/" + smiles + "/JSON"
    BiotitePubChemByInchiKey(key) => "inchikey/" + key + "/JSON"
    BiotitePubChemByInchi => "inchi/JSON"
    BiotitePubChemByFormula(formula) => "formula/" + formula + "/JSON"
    BiotitePubChemByCid(cid) => "cid/" + cid.to_string() + "/JSON"
  }
  let url = base + "/" + search_part
  if properties.length() > 0 {
    let props = biotite_str_join(properties, ",")
    url + "?properties=" + props
  } else {
    url
  }
}

///|
/// Build a PubChem URL for downloading compound properties by CID.
pub fn biotite_pubchem_properties_url(
  cid : Int32,
  properties : Array[String],
) -> String {
  let base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  let props = biotite_str_join(properties, ",")
  base + cid.to_string() + "/property/" + props + "/JSON"
}

///|
/// Build a PubChem URL for downloading compound SDF by CID.
pub fn biotite_pubchem_sdf_url(cid : Int32) -> String {
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" +
  cid.to_string() +
  "/SDF?record_type=3d"
}

///|
/// Parse a simple JSON response to extract a string property value.
/// This is a minimal parser for flat JSON objects.
pub fn biotite_pubchem_parse_json_string(
  json_str : String,
  key : String,
) -> String? {
  // Search for "key":"value" pattern
  let search = "\"" + key + "\":\""
  let pos = json_str.find(search.to_string())
  if pos is Some(start) {
    let val_start = start + search.length()
    // Find closing quote
    let rest = json_str[val_start:]
    let end = rest.find("\"".to_string())
    if end is Some(e) {
      return Some(rest[0:e].to_owned())
    }
  }
  None
}

///|
/// Parse a simple JSON response to extract a numeric property value.
pub fn biotite_pubchem_parse_json_number(
  json_str : String,
  key : String,
) -> Double? {
  // Search for "key":number pattern
  let search = "\"" + key + "\":"
  let pos = json_str.find(search.to_string())
  if pos is Some(start) {
    let val_start = start + search.length()
    let rest = json_str[val_start:]
    // Find end of number (comma, }, or whitespace)
    let mut end = 0
    for i = 0; i < rest.length(); i = i + 1 {
      let ch = rest[i].to_int()
      if ch == 44 || ch == 125 || ch == 32 || ch == 10 { // comma, }, space, newline
        end = i
        break
      }
      end = i + 1
    }
    if end > 0 {
      let num_str = rest[0:end].to_owned()
      let parsed = biotite_parse_dbl(num_str)
      if parsed is Some(v) {
        return Some(v)
      }
    }
  }
  None
}

///|
fn biotite_str_join(arr : Array[String], sep : String) -> String {
  let sb = StringBuilder()
  for i = 0; i < arr.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write(sep)
    }
    let _ = sb.write(arr[i])
  }
  sb.to_string()
}

///|
fn biotite_parse_dbl(s : String) -> Double? {
  let mut result = 0.0
  let mut decimal = false
  let mut divisor = 1.0
  let mut negative = false
  let mut started = false
  for i = 0; i < s.length(); i = i + 1 {
    let ch = s[i].to_int()
    if ch == 45 { // '-'
      negative = true
      continue
    }
    if ch == 46 { // '.'
      decimal = true
      continue
    }
    if ch >= 48 && ch <= 57 {
      let digit = (ch - 48).to_double()
      if decimal {
        divisor = divisor * 10.0
        result = result + digit / divisor
      } else {
        result = result * 10.0 + digit
      }
      started = true
    }
  }
  if !started {
    None
  } else {
    Some(if negative { -result } else { result })
  }
}

///|
/// Use unused constructors.
fn _biotite_pubchem_use_constructors() -> Unit {
  ignore(BiotitePubChemByName(""))
  ignore(BiotitePubChemBySmiles(""))
  ignore(BiotitePubChemByInchiKey(""))
  ignore(BiotitePubChemByInchi)
  ignore(BiotitePubChemByFormula(""))
  ignore(BiotitePubChemByCid(0))
}
