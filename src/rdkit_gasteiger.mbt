/// Gasteiger partial atomic charges.
///
/// RDKit equivalent: `rdkit.Chem.PartialCharges.ComputeGasteigerCharges`
///
/// The Gasteiger-Marsili method iteratively computes partial charges
/// based on electronegativity equalization. Each atom's charge is
/// refined through a damping parameter.

///|
/// Element electronegativity values (Pauling scale, simplified)
fn gasteiger_en(symbol : String) -> Double {
  match symbol {
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
    "B" => 2.04
    "Si" => 1.90
    "Se" => 2.55
    "As" => 2.18
    _ => 2.5 // default
  }
}

///|
/// Element-specific charge damping parameter (higher = faster convergence)
fn gasteiger_damping(symbol : String) -> Double {
  match symbol {
    "H" => 0.5
    "C" => 1.0
    "N" => 1.5
    "O" => 2.0
    "F" => 2.5
    "P" => 1.0
    "S" => 1.2
    "Cl" => 1.8
    "Br" => 1.5
    "I" => 1.2
    _ => 1.0
  }
}

///|
/// Compute Gasteiger partial charges for all atoms.
/// Returns array of partial charges (one per atom).
/// Iterative method: q_i = q_i + alpha * (EN_avg - EN_i) / damping_i
pub fn rdkit_gasteiger_charges(
  mol : Molecule,
  max_iter~ : Int32 = 20,
  alpha~ : Double = 0.5,
) -> Array[Double] {
  let n = mol.num_atoms()
  if n == 0 {
    return []
  }
  // Initialize charges to formal charge / n
  let charges : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    charges[i] = mol.atoms[i].charge.to_double()
  }
  // Build adjacency list
  let adj : Array[Array[Int32]] = Array::make(n, [])
  for bond in mol.bonds {
    adj[bond.from].push(bond.to)
    adj[bond.to].push(bond.from)
  }
  // Iterative refinement
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let new_charges = charges.copy()
    for i = 0; i < n; i = i + 1 {
      let en_i = gasteiger_en(mol.atoms[i].symbol)
      let damping_i = gasteiger_damping(mol.atoms[i].symbol)
      let neighbors = adj[i]
      let mut delta = 0.0
      for j in neighbors {
        let en_j = gasteiger_en(mol.atoms[j].symbol)
        // Charge flows from lower EN to higher EN
        let diff = en_j - en_i
        let q_diff = charges[i] - charges[j]
        delta = delta + alpha * diff * (1.0 - q_diff * 0.5) / damping_i
      }
      new_charges[i] = charges[i] + delta / neighbors.length().to_double()
    }
    // Copy back
    for i = 0; i < n; i = i + 1 {
      charges[i] = new_charges[i]
    }
  }
  // Normalize: total charge should equal formal charge
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    total = total + charges[i]
  }
  let formal = mol_formal_charge(mol).to_double()
  if n > 0 {
    let correction = (formal - total) / n.to_double()
    for i = 0; i < n; i = i + 1 {
      charges[i] = charges[i] + correction
    }
  }
  charges
}

///|
/// Get the charge of a specific atom (by index).
pub fn rdkit_atom_charge(mol : Molecule, atom_idx : Int32) -> Double {
  let charges = rdkit_gasteiger_charges(mol)
  if atom_idx >= 0 && atom_idx < charges.length() {
    charges[atom_idx]
  } else {
    0.0
  }
}

///|
/// Compute the most positive and most negative charges.
pub fn rdkit_charge_extrema(mol : Molecule) -> (Double, Double) {
  let charges = rdkit_gasteiger_charges(mol)
  if charges.length() == 0 {
    return (0.0, 0.0)
  }
  let mut max_q = charges[0]
  let mut min_q = charges[0]
  for q in charges {
    if q > max_q {
      max_q = q
    }
    if q < min_q {
      min_q = q
    }
  }
  (max_q, min_q)
}
