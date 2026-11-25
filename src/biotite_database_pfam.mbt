///|
/// biotite_database_pfam.mbt — Pfam / InterPro URL builders and parsing.

///|
pub fn biotite_pfam_entry_url(accession : String) -> String {
  "https://www.ebi.ac.uk/interpro/api/entry/pfam/" + accession + "/?format=json"
}

///|
pub fn biotite_pfam_alignment_url(accession : String) -> String {
  "https://pfam.xfam.org/family/" + accession + "/alignment/stockholm"
}

///|
/// Extract the first `"description":"..."` value from a Pfam JSON response.
pub fn biotite_pfam_parse_description(json : String) -> String? {
  let needle = "\"description\":\""
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

///|
/// Extract the first `"accession":"..."` value.
pub fn biotite_pfam_parse_accession(json : String) -> String? {
  let needle = "\"accession\":\""
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
