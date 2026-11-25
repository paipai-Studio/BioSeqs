///|
/// biotite_structure_residue.mbt — residue iteration helpers.

///|
/// Residue descriptors: (chain ID, residue ID, residue name, start atom index).
pub fn biotite_get_residues(
  arr : AtomArray,
) -> Array[(String, Int32, String, Int32)] {
  let starts = arr.get_residue_starts()
  let result : Array[(String, Int32, String, Int32)] = Array(
    capacity=starts.length(),
  )
  for start in starts {
    result.push(
      (
        arr.get_chain_id(start).to_string(),
        arr.res_ids[start],
        arr.get_res_name(start),
        start,
      ),
    )
  }
  result
}

///|
/// Start atom index of every residue.
pub fn biotite_residue_iter(arr : AtomArray) -> Array[Int32] {
  arr.get_residue_starts()
}

///|
/// Average residue mass for every residue in the structure.
pub fn biotite_residue_masses(arr : AtomArray) -> Array[Double] {
  let residues = biotite_get_residues(arr)
  let result : Array[Double] = Array(capacity=residues.length())
  for (_chain, _res_id, res_name, _start) in residues {
    result.push(biotite_residue_mass(res_name))
  }
  result
}

///|
/// Number of residues in the structure.
pub fn biotite_get_residue_count(arr : AtomArray) -> Int32 {
  arr.get_residue_count()
}
