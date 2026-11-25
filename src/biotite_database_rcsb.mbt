///|
/// biotite_database_rcsb.mbt — RCSB download/API URL builders and JSON
/// response extraction (no network access).

///|
pub fn biotite_rcsb_fetch_pdb_url(pdb_id : String) -> String {
  "https://files.rcsb.org/download/" + pdb_id + ".pdb"
}

///|
pub fn biotite_rcsb_fetch_cif_url(pdb_id : String) -> String {
  "https://files.rcsb.org/download/" + pdb_id + ".cif"
}

///|
pub fn biotite_rcsb_entry_api_url(pdb_id : String) -> String {
  "https://data.rcsb.org/rest/v1/core/entry/" + pdb_id
}

///|
pub fn biotite_rcsb_assembly_api_url(pdb_id : String) -> String {
  "https://data.rcsb.org/rest/v1/core/assembly/" + pdb_id
}

///|
/// Extract every `"assembly_id":"..."` value from an RCSB assembly JSON.
pub fn biotite_rcsb_parse_assembly_ids(json : String) -> Array[String] {
  let needle = "\"assembly_id\":\""
  let result : Array[String] = []
  let mut start = 0
  while start < json.length() {
    let rest = json[start:]
    match rest.find(needle.to_string()) {
      Some(pos) => {
        let value_start = start + pos + needle.length()
        let tail = json[value_start:]
        match tail.find("\"".to_string()) {
          Some(end) => {
            result.push(tail[0:end].to_owned())
            start = value_start + end + 1
          }
          None => break
        }
      }
      None => break
    }
  }
  result
}

///|
/// Extract the first `"title":"..."` value from an RCSB entry JSON.
pub fn biotite_rcsb_parse_title(json : String) -> String? {
  let needle = "\"title\":\""
  match json.find(needle.to_string()) {
    Some(pos) => {
      let tail = json[pos + needle.length():]
      match tail.find("\"".to_string()) {
        Some(end) => Some(tail[0:end].to_owned())
        None => None
      }
    }
    None => None
  }
}
