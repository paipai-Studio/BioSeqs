///|
/// biotite_structure_compare.mbt — atom-by-atom structure comparison.

///|
pub(all) struct BiotiteCompareResult {
  identical : Array[Bool]
  n_identical : Int32
  n_different : Int32
  rmsd : Double
}

///|
fn biotite_atom_identity_equal(
  a : AtomArray,
  i : Int32,
  b : AtomArray,
  j : Int32,
) -> Bool {
  a.get_name(i) == b.get_name(j) &&
  a.get_res_name(i) == b.get_res_name(j) &&
  a.get_chain_id(i) == b.get_chain_id(j) &&
  a.res_ids[i] == b.res_ids[j] &&
  a.elements[i] == b.elements[j]
}

///|
/// Compare two structures position by position.
///
/// An atom pair is considered identical when annotation and coordinates
/// agree within `tolerance`. Length differences count as different.
pub fn biotite_compare(
  a : AtomArray,
  b : AtomArray,
  tolerance~ : Double = 1.0e-6,
) -> BiotiteCompareResult {
  let n = if a.len() > b.len() { a.len() } else { b.len() }
  let identical : Array[Bool] = Array(capacity=n)
  let mut n_identical = 0
  let mut sum_sq = 0.0
  for i in 0..<n {
    if i < a.len() &&
      i < b.len() &&
      biotite_atom_identity_equal(a, i, b, i) &&
      a.get_coord(i).distance(b.get_coord(i)) <= tolerance {
      identical.push(true)
      n_identical = n_identical + 1
    } else {
      identical.push(false)
    }
    if i < a.len() && i < b.len() {
      let d = a.get_coord(i).distance(b.get_coord(i))
      sum_sq = sum_sq + d * d
    }
  }
  let compared = if a.len() < b.len() { a.len() } else { b.len() }
  let rmsd = if compared > 0 {
    @math.pow(sum_sq / compared.to_double(), 0.5)
  } else {
    0.0
  }
  { identical, n_identical, n_different: n - n_identical, rmsd, }
}

///|
/// All (index in a, index in b) pairs that are identical within tolerance.
pub fn biotite_get_identical_atoms(
  a : AtomArray,
  b : AtomArray,
  tolerance~ : Double = 1.0e-6,
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = []
  for i in 0..<a.len() {
    for j in 0..<b.len() {
      if biotite_atom_identity_equal(a, i, b, j) &&
        a.get_coord(i).distance(b.get_coord(j)) <= tolerance {
        result.push((i, j))
        break
      }
    }
  }
  result
}
