/// biotite_formal_charges.mbt — Formal charge inference from element and bonds
/// Inspired by biotite.structure.formal_charges

///|
/// Typical valence electrons for common elements.
fn biotite_valence_electrons(element : String) -> Int32 {
  match element {
    "H" => 1
    "C" => 4
    "N" => 5
    "O" => 6
    "F" => 7
    "P" => 5
    "S" => 6
    "Cl" => 7
    "Br" => 7
    "I" => 7
    "B" => 3
    "Li" => 1
    "Na" => 1
    "K" => 1
    "Mg" => 2
    "Ca" => 2
    "Zn" => 2
    "Fe" => 2
    "Cu" => 1
    "Mn" => 2
    "Al" => 3
    "Si" => 4
    "Se" => 6
    "As" => 5
    _ => 0
  }
}

///|
/// Typical number of lone pairs for common elements in neutral state.
fn biotite_typical_lone_pairs(
  element : String,
  bond_order_sum : Double,
) -> Int32 {
  // For neutral atoms, lone pairs = (valence - bonds) / 2
  let ve = biotite_valence_electrons(element)
  if ve == 0 {
    return 0
  }
  let lp = (ve.to_double() - bond_order_sum) / 2.0
  lp.floor().to_int()
}

///|
/// Calculate formal charge for an atom.
///
/// Formal charge = valence_electrons - (non_bonding_electrons + bonding_electrons/2)
///              = valence_electrons - (2*lone_pairs + bond_order_sum)
///
/// Parameters:
/// - element: element symbol (e.g. "C", "N", "O")
/// - bond_order_sum: sum of bond orders to all neighbors
/// - n_lone_pairs: number of lone pairs (-1 = auto-compute)
pub fn biotite_formal_charge(
  element : String,
  bond_order_sum : Double,
  n_lone_pairs~ : Int32 = -1,
) -> Int32 {
  let ve = biotite_valence_electrons(element)
  if ve == 0 {
    return 0
  }
  let lp = if n_lone_pairs >= 0 {
    n_lone_pairs
  } else {
    biotite_typical_lone_pairs(element, bond_order_sum)
  }
  // Formal charge = VE - (2*LP + BOS)
  ve - (2 * lp + bond_order_sum.floor().to_int())
}

///|
/// Atom info for batch formal charge calculation.
pub struct BiotiteFormalChargeAtom {
  element : String
  bond_order_sum : Double
  n_lone_pairs : Int32 // -1 = auto
}

///|
/// Calculate formal charges for all atoms in a molecule.
///
/// atoms: array of (element, bond_order_sum, n_lone_pairs)
/// Returns: array of formal charges
pub fn biotite_formal_charges(
  atoms : Array[BiotiteFormalChargeAtom],
) -> Array[Int32] {
  let results : Array[Int32] = Array::make(atoms.length(), 0)
  for i = 0; i < atoms.length(); i = i + 1 {
    let atom = atoms[i]
    results[i] = biotite_formal_charge(
      atom.element,
      atom.bond_order_sum,
      n_lone_pairs=atom.n_lone_pairs,
    )
  }
  results
}

///|
/// Verify that formal charges sum to the expected total charge.
pub fn biotite_verify_total_charge(
  atoms : Array[BiotiteFormalChargeAtom],
  expected_charge : Int32,
) -> Bool {
  let charges = biotite_formal_charges(atoms)
  let mut total = 0
  for i = 0; i < charges.length(); i = i + 1 {
    total = total + charges[i]
  }
  total == expected_charge
}

///|
/// Estimate partial charges using electronegativity-based method.
/// This is a simplified model: assigns partial charges based on
/// electronegativity differences between bonded atoms.
pub fn biotite_partial_charges(
  elements : Array[String],
  bond_orders : Array[(Int32, Int32, Double)], // (atom_i, atom_j, order)
) -> Array[Double] {
  let n = elements.length()
  let charges : Array[Double] = Array::make(n, 0.0)
  // Electronegativity (Pauling scale)
  let en = fn(elem : String) -> Double {
    match elem {
      "H" => 2.20
      "C" => 2.55
      "N" => 3.04
      "O" => 3.44
      "F" => 3.98
      "P" => 2.19
      "S" => 2.58
      "Cl" => 3.16
      "Br" => 2.96
      "I" => 2.66
      _ => 2.50
    }
  }
  for bond_idx = 0; bond_idx < bond_orders.length(); bond_idx = bond_idx + 1 {
    let (i, j, order) = bond_orders[bond_idx]
    if i < n && j < n {
      let en_i = en(elements[i])
      let en_j = en(elements[j])
      let diff = en_j - en_i
      // Partial charge proportional to EN difference * bond order
      let delta = diff * order * 0.25
      charges[i] = charges[i] + delta
      charges[j] = charges[j] - delta
    }
  }
  charges
}

///|
/// Default constructor for BiotiteFormalChargeAtom.
pub fn BiotiteFormalChargeAtom::default() -> BiotiteFormalChargeAtom {
  { element: "C", bond_order_sum: 4.0, n_lone_pairs: 0, }
}
