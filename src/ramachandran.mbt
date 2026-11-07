///|
/// Ramachandran plot analysis — backbone dihedral angle regions.
///
/// Biotite provides Ramachandran plot analysis. This module classifies
/// residues into favored/allowed/disallowed regions.

///|
/// Ramachandran region classification.
pub enum RamaRegion {
  Favored // most populated region
  Allowed // less populated but sterically allowed
  Outlier // disallowed region
}

///|
/// Get region as string.
pub fn RamaRegion::to_string(self : RamaRegion) -> String {
  match self {
    Favored => "favored"
    Allowed => "allowed"
    Outlier => "outlier"
  }
}

///|
/// Classify a residue's (phi, psi) pair into Ramachandran regions.
///
/// Uses a simplified boundary-based classification.
/// General case (non-Gly, non-Pro, non-prePro):
/// - Favored: |phi| in [60, 120] and |psi| in [30, 120] (beta)
///            or phi in [-120, -60] and psi in [-60, 0] (alpha-R)
/// - Allowed: broader boundaries
/// - Outlier: everything else
pub fn classify_rama(
  phi : Double,
  psi : Double,
  residue_type~ : String = "general",
) -> RamaRegion {
  match residue_type {
    "GLY" => classify_rama_gly(phi, psi)
    "PRO" => classify_rama_pro(phi, psi)
    _ => classify_rama_general(phi, psi)
  }
}

///|
/// General case Ramachandran classification.
fn classify_rama_general(phi : Double, psi : Double) -> RamaRegion {
  // Alpha helix region (favored)
  if phi >= -120.0 && phi <= -30.0 && psi >= -90.0 && psi <= -10.0 {
    return Favored
  }
  // Beta sheet region (favored)
  if phi >= -150.0 && phi <= -60.0 && (psi >= 90.0 || psi <= -150.0) {
    let p = if psi < 0.0 { psi + 360.0 } else { psi }
    if p >= 90.0 && p <= 180.0 {
      return Favored
    }
    if p >= 180.0 && p <= 240.0 {
      return Favored
    }
  }
  // Left-handed helix (allowed but rare)
  if phi >= 30.0 && phi <= 120.0 && psi >= 0.0 && psi <= 60.0 {
    return Allowed
  }
  // Extended beta region (allowed)
  if phi >= -180.0 && phi <= -60.0 && psi >= 100.0 && psi <= 180.0 {
    return Favored
  }
  // Check broader allowed region
  if phi >= -180.0 && phi <= 0.0 {
    if psi >= -180.0 && psi <= 180.0 {
      return Allowed
    }
  }
  Outlier
}

///|
/// Glycine Ramachandran classification (broader allowed regions).
fn classify_rama_gly(phi : Double, psi : Double) -> RamaRegion {
  // Glycine can access more conformational space
  if phi >= -180.0 && phi <= 180.0 && psi >= -180.0 && psi <= 180.0 {
    // Check if in highly populated regions
    if (phi >= -120.0 && phi <= -30.0 && psi >= -90.0 && psi <= -10.0) ||
      (phi >= -150.0 && phi <= -60.0 && psi >= 100.0 && psi <= 180.0) ||
      (phi >= 30.0 && phi <= 120.0 && psi >= 0.0 && psi <= 60.0) {
      return Favored
    }
    return Allowed
  }
  Outlier
}

///|
/// Proline Ramachandran classification (restricted phi).
fn classify_rama_pro(phi : Double, psi : Double) -> RamaRegion {
  // Proline has restricted phi around -60°
  if phi >= -90.0 && phi <= -30.0 {
    if psi >= -90.0 && psi <= -10.0 {
      return Favored
    }
    if psi >= 100.0 && psi <= 180.0 {
      return Allowed
    }
  }
  if phi >= -100.0 && phi <= -20.0 {
    return Allowed
  }
  Outlier
}

///|
/// Classify all residues in a structure.
///
/// Returns an array of (residue_index, region) pairs.
pub fn classify_all_ramachandran(
  phi_psi_pairs : Array[(Double, Double)],
  residue_types : Array[String],
) -> Array[(Int32, RamaRegion)] {
  let result : Array[(Int32, RamaRegion)] = []
  for i in 0..<phi_psi_pairs.length() {
    let phi = phi_psi_pairs[i].0
    let psi = phi_psi_pairs[i].1
    let rtype = if i < residue_types.length() {
      residue_types[i]
    } else {
      "general"
    }
    let region = classify_rama(phi, psi, residue_type=rtype)
    result.push((i, region))
  }
  result
}

///|
/// Count residues in each Ramachandran region.
pub fn count_rama_regions(
  classifications : Array[(Int32, RamaRegion)],
) -> (Int32, Int32, Int32) {
  let mut favored = 0
  let mut allowed = 0
  let mut outlier = 0
  for pair in classifications {
    match pair.1 {
      Favored => favored = favored + 1
      Allowed => allowed = allowed + 1
      Outlier => outlier = outlier + 1
    }
  }
  (favored, allowed, outlier)
}

///|
/// Calculate Ramachandran outlier percentage.
pub fn rama_outlier_percentage(
  classifications : Array[(Int32, RamaRegion)],
) -> Double {
  let total = classifications.length()
  if total == 0 {
    return 0.0
  }
  let (_, _, outlier) = count_rama_regions(classifications)
  outlier.to_double() / total.to_double() * 100.0
}
