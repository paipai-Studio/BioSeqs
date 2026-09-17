///|
/// AlphaFold Protein Structure Database client.
///
/// Port of Biopython's Bio.PDB.alphafold_db. Queries
/// `https://alphafold.com/api/prediction/{accession}` for the list of
/// predicted-structure entries and exposes their `entryId` / `cifUrl` /
/// `pdbUrl` / `pdbId` metadata. The HTTP transfer is injected via `fetch`, so
/// the module is fully offline-testable; file I/O is deliberately left to the
/// host (this library only returns URLs and file names).

///|
pub suberror AlphafoldError {
  AlphafoldError(String)
}

///|
fn[T] af_fail(msg : String) -> T raise AlphafoldError {
  raise AlphafoldError(msg)
}

///|
/// A single AlphaFold prediction record. Known fields are decoded into typed
/// accessors; every other field is preserved in `fields`.
pub struct AlphafoldPrediction {
  entry_id : String
  cif_url : String
  pdb_url : String
  pdb_id : String
  uniprot_accession : String
  uniprot_name : String
  gene : String
  uniprot_start : Int
  uniprot_end : Int
  latest_version : Int
  fields : Map[String, String]
}

///|
/// The predictions URL for a UniProt accession.
pub fn alphafold_predictions_url(qualifier : String) -> String {
  "https://alphafold.com/api/prediction/" + qualifier
}

///|
/// The file name (last path segment) of an AlphaFold `cifUrl`.
pub fn alphafold_cif_file_name(cif_url : String) -> String {
  let parts = split_by_char(cif_url, '/'.to_int())
  let mut i = parts.length() - 1
  while i >= 0 {
    if parts[i].length() > 0 {
      return parts[i]
    }
    i = i - 1
  }
  ""
}

// ---------------------------------------------------------------------------
// Minimal JSON array-of-objects parser
// ---------------------------------------------------------------------------

///|
fn af_is_json_ws(cp : Int) -> Bool {
  cp == 32 || cp == 9 || cp == 10 || cp == 13
}

///|
fn af_skip_ws(s : String, i : Int) -> Int {
  let mut j = i
  while j < s.length() && af_is_json_ws(s.unsafe_get(j).to_int()) {
    j = j + 1
  }
  j
}

///|
fn af_hex4(s : String, i : Int) -> Int raise AlphafoldError {
  let mut v = 0
  for k = 0; k < 4; k = k + 1 {
    let c = s.unsafe_get(i + k).to_int()
    let d = if c >= 48 && c <= 57 {
      c - 48
    } else if c >= 97 && c <= 102 {
      c - 97 + 10
    } else if c >= 65 && c <= 70 {
      c - 65 + 10
    } else {
      af_fail("Invalid \\u escape in JSON string")
    }
    v = v * 16 + d
  }
  v
}

///|
/// Parse a JSON string starting at the opening quote (handles common escapes
/// and \uXXXX). Returns (string, next index).
fn af_parse_string(s : String, i : Int) -> (String, Int) raise AlphafoldError {
  let buf = StringBuilder::new()
  let mut j = i + 1
  while j < s.length() {
    let c = s.unsafe_get(j).to_int()
    if c == 34 {
      return (buf.to_string(), j + 1)
    }
    if c == 92 {
      let e = s.unsafe_get(j + 1).to_int()
      if e == 117 {
        let code = af_hex4(s, j + 2)
        buf.write_char(Char::from_int(code))
        j = j + 6
        continue
      }
      if e == 34 {
        buf.write_char('"')
      } else if e == 92 {
        buf.write_char('\\')
      } else if e == 47 {
        buf.write_char('/')
      } else if e == 110 {
        buf.write_char('\n')
      } else if e == 116 {
        buf.write_char('\t')
      } else if e == 114 {
        buf.write_char('\r')
      } else if e == 98 {
        buf.write_char(Char::from_int(8))
      } else if e == 102 {
        buf.write_char(Char::from_int(12))
      } else {
        af_fail("Invalid JSON escape character")
      }
      j = j + 2
      continue
    }
    buf.write_char(s.unsafe_get(j).unsafe_to_char())
    j = j + 1
  }
  af_fail("Unterminated JSON string")
}

///|
/// Parse a JSON value (string, number, true/false, null) into its text form.
/// Numbers are kept verbatim; `null` becomes the empty string.
fn af_parse_value(s : String, i : Int) -> (String, Int) raise AlphafoldError {
  if s.unsafe_get(i).to_int() == 34 {
    return af_parse_string(s, i)
  }
  let mut j = i
  while j < s.length() {
    let d = s.unsafe_get(j).to_int()
    if d == 44 || d == 125 || d == 93 {
      break
    }
    j = j + 1
  }
  let raw = s[i:j].to_string()
  let raw = trim(raw)
  if raw == "null" {
    ("", j)
  } else {
    (raw, j)
  }
}

///|
/// Parse a JSON object into a string map. Assumes `s[i] == '{'`.
fn af_parse_object(
  s : String,
  i : Int,
) -> (Map[String, String], Int) raise AlphafoldError {
  let m : Map[String, String] = Map([], capacity=16)
  let mut j = af_skip_ws(s, i + 1)
  if s.unsafe_get(j).to_int() == 125 {
    return (m, j + 1)
  }
  while j < s.length() {
    let (key, k1) = af_parse_string(s, j)
    let k2 = af_skip_ws(s, k1)
    if s.unsafe_get(k2).to_int() != 58 {
      af_fail("Malformed JSON object: expected ':'")
    }
    let k3 = af_skip_ws(s, k2 + 1)
    let (value, k4) = af_parse_value(s, k3)
    m.set(key, value)
    let k5 = af_skip_ws(s, k4)
    let d = s.unsafe_get(k5).to_int()
    if d == 44 {
      j = af_skip_ws(s, k5 + 1)
    } else if d == 125 {
      return (m, k5 + 1)
    } else {
      af_fail("Malformed JSON object: expected ',' or '}'")
    }
  }
  af_fail("Unterminated JSON object")
}

///|
/// Parse a JSON array of objects into a list of string maps.
fn af_parse_array(
  s : String,
  i : Int,
) -> (Array[Map[String, String]], Int) raise AlphafoldError {
  let arr : Array[Map[String, String]] = []
  let mut j = af_skip_ws(s, i + 1)
  if s.unsafe_get(j).to_int() == 93 {
    return (arr, j + 1)
  }
  while j < s.length() {
    if s.unsafe_get(j).to_int() != 123 {
      af_fail("Expected a JSON object in the array")
    }
    let (obj, k1) = af_parse_object(s, j)
    arr.push(obj)
    let k2 = af_skip_ws(s, k1)
    let d = s.unsafe_get(k2).to_int()
    if d == 44 {
      j = af_skip_ws(s, k2 + 1)
    } else if d == 93 {
      return (arr, k2 + 1)
    } else {
      af_fail("Malformed JSON array: expected ',' or ']'")
    }
  }
  af_fail("Unterminated JSON array")
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

///|
/// Parse an AlphaFold predictions response (a JSON array of objects).
pub fn alphafold_parse_predictions(
  json : String,
) -> Array[AlphafoldPrediction] raise AlphafoldError {
  let i = af_skip_ws(json, 0)
  if i >= json.length() || json.unsafe_get(i).to_int() != 91 {
    af_fail("Expected a JSON array")
  }
  let (arr, _) = af_parse_array(json, i)
  let out : Array[AlphafoldPrediction] = []
  for m in arr {
    out.push(AlphafoldPrediction::{
      entry_id: m.get("entryId").unwrap_or(""),
      cif_url: m.get("cifUrl").unwrap_or(""),
      pdb_url: m.get("pdbUrl").unwrap_or(""),
      pdb_id: m.get("pdbId").unwrap_or(""),
      uniprot_accession: m.get("uniprotAccession").unwrap_or(""),
      uniprot_name: m.get("uniprotName").unwrap_or(""),
      gene: m.get("gene").unwrap_or(""),
      uniprot_start: parse_int(m.get("uniprotStart").unwrap_or("0")),
      uniprot_end: parse_int(m.get("uniprotEnd").unwrap_or("0")),
      latest_version: parse_int(m.get("latestVersion").unwrap_or("0")),
      fields: m,
    })
  }
  out
}

///|
/// Fetch and parse the AlphaFold predictions for a UniProt accession.
pub fn alphafold_get_predictions(
  fetch~ : (String, String?) -> String raise AlphafoldError,
  qualifier~ : String,
) -> Array[AlphafoldPrediction] raise AlphafoldError {
  let url = alphafold_predictions_url(qualifier)
  let json = (fetch)(url, None)
  alphafold_parse_predictions(json)
}

///|
/// A sample predictions response for tests and demos.
pub fn alphafold_sample_json() -> String {
  "[{\"entryId\":\"P00520\",\"cifUrl\":\"https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.cif\",\"pdbUrl\":\"https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.pdb\",\"pdbId\":\"AF-P00520-F1\",\"uniprotAccession\":\"P00520\",\"uniprotName\":\"Tyrosine-protein kinase ABL1\",\"uniprotStart\":1,\"uniprotEnd\":284,\"latestVersion\":4,\"gene\":\"ABL1\"},{\"entryId\":\"P69905\",\"cifUrl\":\"https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.cif\",\"pdbUrl\":\"https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb\",\"pdbId\":\"AF-P69905-F1\",\"uniprotAccession\":\"P69905\",\"uniprotName\":\"Hemoglobin subunit alpha\",\"uniprotStart\":1,\"uniprotEnd\":142,\"latestVersion\":4,\"gene\":\"HBA1\"}]"
}