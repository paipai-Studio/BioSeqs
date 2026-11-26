/// Biotite `structure.repair` — functionalities for repairing
/// malformed structures.
///
/// Faithful port of the official implementation
/// (`src/biotite/structure/repair.py` in biotite-dev/biotite):
/// - `create_continuous_res_ids` rebuilds residue IDs that increment
///   by 1 per residue (official `cumsum` over residue starts),
///   optionally restarting at 1 for each chain,
/// - `infer_elements` guesses the element of each atom from its atom
///   name: digits are removed (`1H` -> `H`), names starting with one
///   of the biomolecule elements C/N/O/S/H contribute their first
///   character, otherwise the name is matched exactly against the
///   element table (two characters first, then one),
/// - `create_atom_names` enumerates atoms per element as
///   `{element}{n}`; for a single atom the element itself is used as
///   name,
/// - unresolvable element names yield an empty string (official
///   implementation emits a `UserWarning` first; this port has no
///   warnings).

///|
/// The official element table of `_guess_element` (all known
/// elements, uppercase abbreviations).
let biotite_repair_elements : Array[String] = [
  "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P",
  "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU",
  "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC",
  "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA",
  "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
  "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT",
  "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES",
  "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH",
  "FL", "MC", "LV", "TS", "OG",
]

///|
/// Guess the element of an atom from its atom name
/// (official `_guess_element`): digits are stripped, biomolecule
/// elements (C/N/O/S/H) win by their first character, otherwise the
/// remaining name is matched exactly against the element table (two
/// characters first, then one). An empty string is returned when no
/// element can be inferred.
fn biotite_repair_guess_element(atom_name : String) -> String {
  let buf = StringBuilder::new()
  for c in atom_name {
    if !(c >= '0' && c <= '9') {
      buf.write_char(c)
    }
  }
  let elem = buf.to_string().to_upper()
  if elem == "" {
    return ""
  }
  let echars : Array[Char] = []
  for c in elem {
    echars.push(c)
  }
  let first = echars[0]
  if first == 'C' ||
    first == 'N' ||
    first == 'O' ||
    first == 'S' ||
    first == 'H' {
    return first.to_string()
  }
  // Exact match of the (up to) two leading characters, then of the
  // first character alone.
  if echars.length() >= 2 {
    let two = echars[0].to_string() + echars[1].to_string()
    for e in biotite_repair_elements {
      if e == two {
        return e
      }
    }
  }
  let one = first.to_string()
  for e in biotite_repair_elements {
    if e == one {
      return e
    }
  }
  ""
}

///|
/// Infer the elements of atoms based on their atom names
/// (official `infer_elements`; pass `arr.names` when working on an
/// AtomArray).
pub fn biotite_infer_elements(atom_names : Array[String]) -> Array[String] {
  let out : Array[String] = Array(capacity=atom_names.length())
  for name in atom_names {
    out.push(biotite_repair_guess_element(name))
  }
  out
}

///|
/// Create atom names for a single residue based on elements: each
/// element is enumerated separately as `{element}{n}`; for a single
/// atom the element itself is the name (official `create_atom_names`).
/// Only use on single residues — the enumeration would continue into
/// the next residue otherwise.
pub fn biotite_create_atom_names(elements : Array[String]) -> Array[String] {
  if elements.length() == 1 {
    return [elements[0]]
  }
  let counter : Map[String, Int32] = Map([], capacity=16)
  let out : Array[String] = Array(capacity=elements.length())
  for e in elements {
    let cur = match counter.get(e) {
      Some(v) => v
      None => 0
    }
    counter[e] = cur + 1
    out.push(e + (cur + 1).to_string())
  }
  out
}

///|
/// Create an array of continuous residue IDs: incremented by 1 for
/// each residue (official `create_continuous_res_ids`). With
/// `restart_each_chain` the IDs are reset to 1 for each chain.
pub fn biotite_create_continuous_res_ids(
  arr : AtomArray,
  restart_each_chain~ : Bool = true,
) -> Array[Int32] {
  let n = arr.len()
  // Official construction: a 1 at each residue start, cumsum over all
  // atoms yields per-atom continuous residue IDs.
  let diff = Array::make(n, 0)
  for s in arr.get_residue_starts() {
    diff[s] = 1
  }
  let res_ids : Array[Int32] = Array(capacity=n)
  let mut run = 0
  for i = 0; i < n; i = i + 1 {
    run = run + diff[i]
    res_ids.push(run)
  }
  if restart_each_chain {
    for s in arr.get_chain_starts() {
      let shift = res_ids[s] - 1
      for i in s..<n {
        res_ids[i] = res_ids[i] - shift
      }
    }
  }
  res_ids
}

///|
pub fn _biotite_structure_repair_dummy() -> Unit {
  let arr = AtomArray::new()
  let _ = biotite_create_continuous_res_ids(arr)
  let _ = biotite_infer_elements([])
  let _ = biotite_create_atom_names([])
  let _ = biotite_repair_guess_element("CA")
  ignore(biotite_repair_elements)
}
