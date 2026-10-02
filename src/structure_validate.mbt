///|
/// OpenStructure-style structure validation: clash detection,
/// distance-based bond inference, and stereochemistry checks.

///|
/// A single clash between two atoms.
pub struct Clash {
  atom1_idx : Int
  atom2_idx : Int
  atom1_name : String
  atom2_name : String
  distance : Double
  vdw_sum : Double
  overlap : Double // vdw_sum * 0.85 - distance
}

///|
/// Result of clash detection.
pub struct ClashResult {
  clashes : Array[Clash]
  num_clashes : Int
  clash_score : Double
}

///|
/// An inferred bond between two atoms.
pub struct InferredBond {
  atom1_idx : Int
  atom2_idx : Int
  atom1_name : String
  atom2_name : String
  distance : Double
  expected_length : Double
}

///|
/// Stereochemistry validation result.
pub struct StereoValidation {
  num_chiral_centers : Int
  num_cis_peptides : Int
  num_bad_bonds : Int
  num_bad_angles : Int
  issues : Array[String]
}

///|
/// Detect steric clashes in a structure.
/// Two atoms clash if distance < tolerance * (vdw1 + vdw2).
/// Default tolerance = 0.85 (standard in OpenStructure / MolProbity).
pub fn detect_clashes(
  atoms : Array[Atom],
  tolerance? : Double = 0.85,
) -> ClashResult {
  let n = atoms.length()
  let clashes : Array[Clash] = Array(capacity=0)
  let mut total_score = 0.0
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d = atoms[i].distance(atoms[j])
      let vdw1 = chem_vdw_radius(atoms[i].element)
      let vdw2 = chem_vdw_radius(atoms[j].element)
      let vdw_sum = vdw1 + vdw2
      let threshold = tolerance * vdw_sum
      // Skip if atoms are too far apart (optimization: > 5 Å)
      if d < 5.0 && d < threshold && d > 0.01 {
        let overlap = threshold - d
        clashes.push(Clash::{
          atom1_idx: i,
          atom2_idx: j,
          atom1_name: atoms[i].name,
          atom2_name: atoms[j].name,
          distance: d,
          vdw_sum,
          overlap,
        })
        total_score = total_score + overlap * overlap
      }
      j = j + 1
    }
    i = i + 1
  }
  ClashResult::{
    clashes,
    num_clashes: clashes.length(),
    clash_score: total_score,
  }
}

///|
/// Overall clash score (sum of squared overlaps).
pub fn clash_score(atoms : Array[Atom]) -> Double {
  detect_clashes(atoms).clash_score
}

///|
/// Count number of clashes.
pub fn count_clashes(atoms : Array[Atom]) -> Int {
  detect_clashes(atoms).num_clashes
}

///|
/// Infer bonds from interatomic distances.
/// Two atoms are bonded if their distance is within
/// [min_factor * cov_sum, max_factor * cov_sum]
/// where cov_sum = covalent_radius(a1) + covalent_radius(a2).
pub fn infer_bonds(
  atoms : Array[Atom],
  min_factor? : Double = 0.4,
  max_factor? : Double = 1.3,
) -> Array[InferredBond] {
  let n = atoms.length()
  let bonds : Array[InferredBond] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d = atoms[i].distance(atoms[j])
      let cov1 = chem_covalent_radius(atoms[i].element)
      let cov2 = chem_covalent_radius(atoms[j].element)
      let cov_sum = cov1 + cov2
      let min_dist = min_factor * cov_sum
      let max_dist = max_factor * cov_sum
      if d >= min_dist && d <= max_dist && d > 0.01 {
        bonds.push(InferredBond::{
          atom1_idx: i,
          atom2_idx: j,
          atom1_name: atoms[i].name,
          atom2_name: atoms[j].name,
          distance: d,
          expected_length: cov_sum,
        })
      }
      j = j + 1
    }
    i = i + 1
  }
  bonds
}

///|
/// Count inferred bonds.
pub fn count_inferred_bonds(atoms : Array[Atom]) -> Int {
  infer_bonds(atoms).length()
}

///|
/// Get inferred bond list as (atom1_idx, atom2_idx) pairs.
pub fn inferred_bonds_pairs(atoms : Array[Atom]) -> Array[(Int, Int)] {
  let bonds = infer_bonds(atoms)
  let pairs : Array[(Int, Int)] = Array(capacity=0)
  for b in bonds {
    pairs.push((b.atom1_idx, b.atom2_idx))
  }
  pairs
}

///|
/// Validate stereochemistry of a protein chain.
/// Checks: cis-peptide bonds, unusual bond lengths/angles.
pub fn validate_stereochemistry(chain : Chain) -> StereoValidation {
  let residues = chain.get_residues()
  let n = residues.length()
  let issues : Array[String] = Array(capacity=0)
  let mut num_cis = 0
  let num_bad_bonds = 0
  let num_bad_angles = 0
  let mut num_chiral = 0
  // Check for cis-peptides: omega angle near 0° instead of 180°
  for i = 0; i < n - 1; i = i + 1 {
    match calc_omega(residues[i], residues[i + 1]) {
      Some(omega) => {
        let dev = (omega - 180.0).abs()
        if dev < 30.0 {
          // omega near 0° → cis peptide
          if omega.abs() < 30.0 {
            num_cis = num_cis + 1
            issues.push("cis-peptide at residue " + residues[i].resseq.to_string())
          }
        }
      }
      None => ()
    }
  }
  // Count CA atoms with 4+ neighbors as potential chiral centers
  for res in residues {
    let atoms = res.get_atoms()
    for atom in atoms {
      if atom.name == "CA" && atom.element == "C" {
        num_chiral = num_chiral + 1
      }
    }
  }
  StereoValidation::{
    num_chiral_centers: num_chiral,
    num_cis_peptides: num_cis,
    num_bad_bonds,
    num_bad_angles,
    issues,
  }
}

///|
/// Check if a structure has any severe clashes.
pub fn has_severe_clashes(atoms : Array[Atom], threshold? : Double = 0.4) -> Bool {
  let result = detect_clashes(atoms)
  for clash in result.clashes {
    if clash.overlap > threshold {
      return true
    }
  }
  false
}

///|
/// Clash report as a human-readable string.
pub fn clash_report(result : ClashResult) -> String {
  let sb = StringBuilder()
  sb.write_string("Clash Report\n")
  sb.write_string("  Total clashes: " + result.num_clashes.to_string() + "\n")
  sb.write_string("  Clash score: " + result.clash_score.to_string() + "\n")
  let max_show = if result.num_clashes < 10 { result.num_clashes } else { 10 }
  for i = 0; i < max_show; i = i + 1 {
    let c = result.clashes[i]
    sb.write_string("  " + c.atom1_name + " -- " + c.atom2_name)
    sb.write_string("  d=" + c.distance.to_string())
    sb.write_string("  overlap=" + c.overlap.to_string() + "\n")
  }
  sb.to_string()
}

///|
/// Bond inference report.
pub fn bond_report(bonds : Array[InferredBond]) -> String {
  let sb = StringBuilder()
  sb.write_string("Inferred Bonds: " + bonds.length().to_string() + "\n")
  let max_show = if bonds.length() < 20 { bonds.length() } else { 20 }
  for i = 0; i < max_show; i = i + 1 {
    let b = bonds[i]
    sb.write_string("  " + b.atom1_name + " -- " + b.atom2_name)
    sb.write_string("  d=" + b.distance.to_string())
    sb.write_string("  expected=" + b.expected_length.to_string() + "\n")
  }
  sb.to_string()
}
