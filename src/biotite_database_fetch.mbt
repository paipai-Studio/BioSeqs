///|
/// biotite_database_fetch.mbt — unified URL dispatch and cache paths.

///|
/// Resolve a database identifier to a download URL.
/// Databases without a public stable URL template return "".
pub fn biotite_database_fetch_url(
  database : String,
  identifier : String,
) -> String {
  match database {
    "rcsb" | "pdb" => biotite_rcsb_fetch_pdb_url(identifier)
    "pfam" => biotite_pfam_entry_url(identifier)
    _ => ""
  }
}

///|
/// Resolve a database identifier using an explicit base URL (for DOB and
/// other configurable services).
pub fn biotite_database_fetch_url_with_base(
  database : String,
  identifier : String,
  base_url : String,
) -> String {
  if database == "dob" {
    biotite_dob_entry_url(biotite_dob_client(base_url), identifier)
  } else {
    biotite_database_fetch_url(database, identifier)
  }
}

///|
/// Deterministic cache path for a URL under `root`.
pub fn biotite_database_cache_path(root : String, url : String) -> String {
  let sb = StringBuilder()
  for ch in url {
    let code = ch.to_int()
    let safe = (code >= 65 && code <= 90) ||
      (code >= 97 && code <= 122) ||
      (code >= 48 && code <= 57) ||
      ch == '-' ||
      ch == '_' ||
      ch == '.'
    if safe {
      let _ = sb.write_char(ch)
    } else {
      let _ = sb.write_char('_')
    }
  }
  root + "/" + sb.to_string() + ".dat"
}
