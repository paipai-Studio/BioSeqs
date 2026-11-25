///|
/// biotite_structure_deduplicate.mbt — duplicate atom/structure removal.

///|
fn biotite_full_atom_signature(arr : AtomArray, i : Int32) -> String {
  let c = arr.get_coord(i)
  arr.get_name(i) +
  "|" +
  arr.get_res_name(i) +
  "|" +
  arr.get_chain_id(i).to_string() +
  "|" +
  arr.res_ids[i].to_string() +
  "|" +
  arr.elements[i] +
  "|" +
  c.x.to_string() +
  "," +
  c.y.to_string() +
  "," +
  c.z.to_string()
}

///|
/// Remove atoms whose full signature (identity + coordinates) was already seen.
pub fn biotite_deduplicate_structure(arr : AtomArray) -> AtomArray {
  let seen : Map[String, Bool] = Map([])
  let keep : Array[Bool] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    let key = biotite_full_atom_signature(arr, i)
    if seen.contains(key) {
      keep.push(false)
    } else {
      seen[key] = true
      keep.push(true)
    }
  }
  arr.filter(keep)
}

///|
/// Mask of atoms in `a` that also occur in `b` (identity + coordinates).
pub fn biotite_filter_identical(a : AtomArray, b : AtomArray) -> Array[Bool] {
  let present : Map[String, Bool] = Map([])
  for j in 0..<b.len() {
    present[biotite_full_atom_signature(b, j)] = true
  }
  let mask : Array[Bool] = Array(capacity=a.len())
  for i in 0..<a.len() {
    mask.push(present.contains(biotite_full_atom_signature(a, i)))
  }
  mask
}

///|
/// Remove duplicate structures, preserving first occurrence order.
pub fn biotite_deduplicate_structures(
  structures : Array[AtomArray],
) -> Array[AtomArray] {
  let seen : Map[String, Bool] = Map([])
  let result : Array[AtomArray] = []
  for arr in structures {
    let key = arr.to_pdb_string()
    if !seen.contains(key) {
      seen[key] = true
      result.push(arr)
    }
  }
  result
}

///|
/// Indices where a boolean mask is true.
pub fn biotite_filter_indices(mask : Array[Bool]) -> Array[Int32] {
  let result : Array[Int32] = []
  for i in 0..<mask.length() {
    if mask[i] {
      result.push(i)
    }
  }
  result
}
