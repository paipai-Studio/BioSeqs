/// biopython_torsion_angles.mbt — Protein backbone torsion angles
/// Inspired by Biopython Bio.PDB.Polypeptide: phi, psi, omega calculation

///|
/// Calculate phi (φ) angle for residue i.
/// phi = dihedral(C_{i-1}, N_i, CA_i, C_i)
/// Requires 4 atoms: C of previous residue, N, CA, C of current residue.
/// Each atom is a 3D coordinate [x, y, z].
pub fn biopython_phi_angle(
  prev_c : Array[Double],
  n : Array[Double],
  ca : Array[Double],
  c : Array[Double],
) -> Double {
  biopython_dihedral_angle(prev_c, n, ca, c)
}

///|
/// Calculate psi (ψ) angle for residue i.
/// psi = dihedral(N_i, CA_i, C_i, N_{i+1})
pub fn biopython_psi_angle(
  n : Array[Double],
  ca : Array[Double],
  c : Array[Double],
  next_n : Array[Double],
) -> Double {
  biopython_dihedral_angle(n, ca, c, next_n)
}

///|
/// Calculate omega (ω) angle for residue i.
/// omega = dihedral(CA_i, C_i, N_{i+1}, CA_{i+1})
pub fn biopython_omega_angle(
  ca : Array[Double],
  c : Array[Double],
  next_n : Array[Double],
  next_ca : Array[Double],
) -> Double {
  biopython_dihedral_angle(ca, c, next_n, next_ca)
}

///|
/// Calculate dihedral angle from 4 points.
/// Returns angle in degrees.
fn biopython_dihedral_angle(
  p1 : Array[Double],
  p2 : Array[Double],
  p3 : Array[Double],
  p4 : Array[Double],
) -> Double {
  // Vectors
  let b1 = biopython_vec_sub(p2, p1)
  let b2 = biopython_vec_sub(p3, p2)
  let b3 = biopython_vec_sub(p4, p3)
  // Normal vectors to planes
  let n1 = biopython_vec_cross(b1, b2)
  let n2 = biopython_vec_cross(b2, b3)
  // Normalize
  let n1_norm = biopython_vec_norm(n1)
  let n2_norm = biopython_vec_norm(n2)
  if n1_norm < 0.0000000001 || n2_norm < 0.0000000001 {
    return 0.0
  }
  let n1_unit = [n1[0] / n1_norm, n1[1] / n1_norm, n1[2] / n1_norm]
  let n2_unit = [n2[0] / n2_norm, n2[1] / n2_norm, n2[2] / n2_norm]
  // cos(angle) = n1 · n2
  let cos_angle = biopython_vec_dot(n1_unit, n2_unit)
  // Clamp to [-1, 1]
  let clamped = if cos_angle > 1.0 {
    1.0
  } else if cos_angle < -1.0 {
    -1.0
  } else {
    cos_angle
  }
  let angle = @math.acos(clamped)
  // Sign: (n1 × n2) · b2
  let cross = biopython_vec_cross(n1_unit, n2_unit)
  let sign = biopython_vec_dot(cross, b2)
  let pi = 3.14159265358979323846
  let result = if sign < 0.0 { -angle } else { angle }
  // Convert to degrees and shift to [-180, 180)
  result * 180.0 / pi
}

///|
/// Calculate all backbone torsion angles for a protein.
/// Input: array of (N, CA, C) coordinates for each residue.
/// Returns: array of (phi, psi, omega) for each residue.
pub fn biopython_backbone_torsion_angles(
  residues : Array[(Array[Double], Array[Double], Array[Double])],
) -> Array[(Double, Double, Double)] {
  let n = residues.length()
  if n == 0 {
    return []
  }
  let result : Array[(Double, Double, Double)] = Array::make(n, (0.0, 0.0, 0.0))
  for i = 0; i < n; i = i + 1 {
    let (n_i, ca_i, c_i) = residues[i]
    // Phi: needs C of previous residue
    let phi = if i > 0 {
      let (_, _, c_prev) = residues[i - 1]
      biopython_phi_angle(c_prev, n_i, ca_i, c_i)
    } else {
      0.0
    }
    // Psi: needs N of next residue
    let psi = if i < n - 1 {
      let (n_next, _, _) = residues[i + 1]
      biopython_psi_angle(n_i, ca_i, c_i, n_next)
    } else {
      0.0
    }
    // Omega: needs N and CA of next residue
    let omega = if i < n - 1 {
      let (n_next, ca_next, _) = residues[i + 1]
      biopython_omega_angle(ca_i, c_i, n_next, ca_next)
    } else {
      0.0
    }
    result[i] = (phi, psi, omega)
  }
  result
}

///|
/// Classify secondary structure from phi/psi angles.
/// Returns: "H" (alpha helix), "E" (beta sheet), "C" (coil).
pub fn biopython_ss_from_torsion(phi : Double, psi : Double) -> String {
  // Alpha helix: phi ≈ -60, psi ≈ -45
  if phi > -120.0 && phi < -30.0 && psi > -90.0 && psi < 0.0 {
    "H"
    // Beta sheet: phi ≈ -120, psi ≈ 120
  } else if phi > -180.0 && phi < -60.0 && psi > 60.0 && psi < 180.0 {
    "E"
    // Left-handed helix: phi ≈ 60, psi ≈ 45
  } else if phi > 30.0 && phi < 120.0 && psi > 0.0 && psi < 90.0 {
    "H"
  } else {
    "C"
  }
}

// Vector operations

///|
fn biopython_vec_sub(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

///|
fn biopython_vec_cross(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ]
}

///|
fn biopython_vec_dot(a : Array[Double], b : Array[Double]) -> Double {
  a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

///|
fn biopython_vec_norm(a : Array[Double]) -> Double {
  biopython_torsion_sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
}

///|
fn biopython_torsion_sqrt(x : Double) -> Double {
  @math.pow(x, 0.5)
}
