///|
/// biotite_structure_integrity.mbt — residue, bond and charge checks.

///|
/// Report residues missing any required atom name.
///
/// Each entry is "CHAIN:RESID:RESNAME missing ATOM".
pub fn biotite_check_residues(
  arr : AtomArray,
  required_atoms : Array[String],
) -> Array[String] {
  let starts = arr.get_residue_starts()
  let result : Array[String] = []
  for ri in 0..<starts.length() {
    let start = starts[ri]
    let end = if ri + 1 < starts.length() { starts[ri + 1] } else { arr.len() }
    for atom_name in required_atoms {
      let mut found = false
      for ai in start..<end {
        if arr.get_name(ai) == atom_name {
          found = true
          break
        }
      }
      if !found {
        result.push(
          arr.get_chain_id(start).to_string() +
          ":" +
          arr.res_ids[start].to_string() +
          ":" +
          arr.get_res_name(start) +
          " missing " +
          atom_name,
        )
      }
    }
  }
  result
}

///|
/// Report double bonds whose length exceeds `max_length`.
pub fn biotite_check_double_bonds(
  arr : AtomArray,
  max_length~ : Double = 1.45,
) -> Array[(Int32, Int32)] {
  let bonds = arr.compute_bonds()
  let result : Array[(Int32, Int32)] = []
  for bond in bonds.bonds {
    if bond.bond_type is BondType::Double {
      let d = arr.get_coord(bond.atom1).distance(arr.get_coord(bond.atom2))
      if d > max_length {
        result.push((bond.atom1, bond.atom2))
      }
    }
  }
  result
}

///|
/// Verify that formal charges sum to the expected total.
pub fn biotite_check_charges(
  atoms : Array[BiotiteFormalChargeAtom],
  expected : Int32,
) -> Bool {
  biotite_verify_total_charge(atoms, expected)
}
