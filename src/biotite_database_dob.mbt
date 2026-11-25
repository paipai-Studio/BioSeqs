///|
/// biotite_database_dob.mbt — configurable Database-of-Biological-Data client.
///
/// The base URL is supplied by the caller so no external host is hard-coded.

///|
pub(all) struct BiotiteDobClient {
  base_url : String
}

///|
pub fn biotite_dob_client(base_url : String) -> BiotiteDobClient {
  { base_url, }
}

///|
pub fn biotite_dob_entry_url(
  client : BiotiteDobClient,
  entry_id : String,
) -> String {
  client.base_url + "/entry/" + entry_id
}

///|
pub fn biotite_dob_search_url(
  client : BiotiteDobClient,
  query : String,
) -> String {
  client.base_url + "/search?q=" + biotite_dob_encode(query)
}

///|
fn biotite_dob_encode(text : String) -> String {
  let sb = StringBuilder()
  let hex = [
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F",
  ]
  for ch in text {
    let code = ch.to_int()
    let keep = (code >= 65 && code <= 90) ||
      (code >= 97 && code <= 122) ||
      (code >= 48 && code <= 57) ||
      ch == '-' ||
      ch == '_' ||
      ch == '.' ||
      ch == '~'
    if keep {
      let _ = sb.write_char(ch)
    } else {
      let hi = code / 16 % 16
      let lo = code % 16
      let _ = sb.write("%" + hex[hi] + hex[lo])
    }
  }
  sb.to_string()
}
