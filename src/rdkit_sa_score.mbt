/// Synthetic Accessibility (SA) score.
///
/// RDKit equivalent: `rdkit.Chem.SA_Score.sascorer`
///
/// Estimates how difficult a molecule is to synthesize.
/// Based on fragment complexity and stereochemistry penalty.
/// Score range: 1 (easy) to 10 (very difficult).
///
/// Reference: Ertl & Schuffenhauer, J Cheminf 2009.

///|
/// Compute synthetic accessibility score.
/// Lower = easier to synthesize.
pub fn rdkit_sa_score(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  if n == 0 {
    return 0.0
  }
  // Component 1: fragment complexity (based on size and rings)
  let fragment_score = rdkit_sa_fragment_score(mol)
  // Component 2: stereochemistry penalty
  let stereo_score = rdkit_sa_stereo_penalty(mol)
  // Component 3: macrocycle penalty
  let macro_score = rdkit_sa_macrocycle_penalty(mol)
  // Component 4: complexity penalty (large molecules harder)
  let size_score = rdkit_sa_size_penalty(mol)
  // Total SA score (clamped to [1, 10])
  let total = fragment_score + stereo_score + macro_score + size_score
  if total < 1.0 {
    1.0
  } else if total > 10.0 {
    10.0
  } else {
    total
  }
}

///|
/// Fragment complexity: based on ring count, sp3 fraction, and atom count.
fn rdkit_sa_fragment_score(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let rings = mol_ring_count(mol)
  let aromatic_rings = mol_aromatic_ring_count(mol)
  let fsp3 = mol_fraction_sp3(mol)
  // Base score: 1.0 for simple molecules
  let mut score = 1.0
  // Ring complexity: each ring adds ~0.3
  score = score + rings.to_double() * 0.3
  // Aromatic rings slightly easier than aliphatic
  score = score - aromatic_rings.to_double() * 0.1
  // High sp3 content = more complex stereochemistry
  score = score + fsp3 * 0.5
  // Very small or very large molecules are harder
  if n < 5 {
    score = score + 0.5
  } else if n > 30 {
    score = score + (n - 30).to_double() * 0.1
  }
  score
}

///|
/// Stereochemistry penalty: each stereocenter adds difficulty.
fn rdkit_sa_stereo_penalty(mol : Molecule) -> Double {
  let centers = find_stereocenters(mol)
  centers.length().to_double() * 0.5
}

///|
/// Macrocycle penalty: rings with > 8 atoms are hard to synthesize.
fn rdkit_sa_macrocycle_penalty(mol : Molecule) -> Double {
  let rings = mol_find_rings(mol)
  let mut penalty = 0.0
  for ring in rings {
    if ring.length() > 8 {
      penalty = penalty + (ring.length() - 8).to_double() * 0.2
    }
  }
  penalty
}

///|
/// Size penalty: very large molecules are harder.
fn rdkit_sa_size_penalty(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let heavy = mol_heavy_atom_count(mol)
  // Logarithmic size penalty
  if heavy > 0 {
    @math.ln(heavy.to_double()) * 0.3
  } else {
    0.0
  }
}

///|
/// Categorize SA score into difficulty levels.
pub fn rdkit_sa_category(mol : Molecule) -> String {
  let score = rdkit_sa_score(mol)
  if score < 3.0 {
    "Easy"
  } else if score < 5.0 {
    "Moderate"
  } else if score < 7.0 {
    "Difficult"
  } else {
    "Very_Difficult"
  }
}

///|
/// Check if molecule is easily synthesizable (SA < 4).
pub fn rdkit_sa_is_easy(mol : Molecule) -> Bool {
  rdkit_sa_score(mol) < 4.0
}
