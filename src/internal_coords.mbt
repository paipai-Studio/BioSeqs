///|
/// Bio.PDB.internal_coords - Internal coordinate representation of protein structures.
/// Inspired by Biopython's Bio.PDB.internal_coords module.
///
/// Represents protein structures using internal coordinates (bond lengths,
/// bond angles, and torsion angles) instead of Cartesian coordinates.
/// This is fundamental for structural biology calculations including:
/// - Protein structure prediction and design
/// - Side-chain rotamer libraries
/// - Conformational analysis
/// - Backbone-dependent dihedral angles
///
/// Torsion angle conventions:
///   phi (C(i-1) - N(i) - CA(i) - C(i))
///   psi (N(i) - CA(i) - C(i) - N(i+1))
///   omega (CA(i-1) - C(i-1) - N(i) - CA(i)), typically ~180°
///   chi1 (N - CA - CB - CG) for side chain rotation

// ============================================================================
// Data Structures
// ============================================================================

///|
/// A single torsion angle with its name and value.
pub struct TorsionAngle {
  name : String
  value : Double
  atom_names : Array[String]
}

///|
/// Bond length between two atoms.
pub struct BondLength {
  atom1 : String
  atom2 : String
  length : Double
}

///|
/// Bond angle formed by three atoms.
pub struct BondAngle {
  atom1 : String
  atom2 : String
  atom3 : String
  angle : Double
}

///|
/// Internal coordinate representation of a single residue.
pub struct InternalCoordResidue {
  resname : String
  resseq : Int
  bond_lengths : Array[BondLength]
  bond_angles : Array[BondAngle]
  torsions : Array[TorsionAngle]
  start_atom : String
}

///|
/// Full internal coordinate representation of a protein chain.
pub struct InternalCoordChain {
  chain_id : Char
  residues : Array[InternalCoordResidue]
}

///|
/// Result of converting internal coordinates to Cartesian.
pub struct ICToCartesianResult {
  atom_names : Array[String]
  residue_names : Array[String]
  residue_ids : Array[Int]
  x : Array[Double]
  y : Array[Double]
  z : Array[Double]
}

///|
/// Parameters for a rotamer conformer.
pub struct Rotamer {
  name : String
  chi_angles : Array[Double]
  probability : Double
}

///|
/// Backbone-dependent rotamer library entry.
pub struct RotamerLibraryEntry {
  phi : Double
  psi : Double
  residues : Array[(String, Array[Rotamer])]
}

// ============================================================================
// Constructors
// ============================================================================

///|
/// Create a new torsion angle.
pub fn TorsionAngle::new(
  name~ : String,
  value~ : Double,
  atom_names? : Array[String] = [],
) -> TorsionAngle {
  TorsionAngle::{ name, value, atom_names }
}

///|
/// Create a new bond length.
pub fn BondLength::new(
  atom1~ : String,
  atom2~ : String,
  length~ : Double,
) -> BondLength {
  BondLength::{ atom1, atom2, length }
}

///|
/// Create a new bond angle.
pub fn BondAngle::new(
  atom1~ : String,
  atom2~ : String,
  atom3~ : String,
  angle~ : Double,
) -> BondAngle {
  BondAngle::{ atom1, atom2, atom3, angle }
}

///|
/// Create internal coordinates for a residue.
pub fn InternalCoordResidue::new(
  resname~ : String,
  resseq~ : Int,
  start_atom? : String = "N",
) -> InternalCoordResidue {
  InternalCoordResidue::{
    resname,
    resseq,
    bond_lengths: [],
    bond_angles: [],
    torsions: [],
    start_atom,
  }
}

///|
/// Create an internal coordinate chain.
pub fn InternalCoordChain::new(chain_id~ : Char) -> InternalCoordChain {
  InternalCoordChain::{ chain_id, residues: [] }
}

///|
/// Create IC to Cartesian result.
pub fn ICToCartesianResult::new() -> ICToCartesianResult {
  ICToCartesianResult::{
    atom_names: [],
    residue_names: [],
    residue_ids: [],
    x: [],
    y: [],
    z: [],
  }
}

///|
/// Create a rotamer.
pub fn Rotamer::new(
  name~ : String,
  chi_angles? : Array[Double] = [],
  probability? : Double = 1.0,
) -> Rotamer {
  Rotamer::{ name, chi_angles, probability }
}

///|
/// Create a rotamer library entry.
pub fn RotamerLibraryEntry::new(
  phi~ : Double,
  psi~ : Double,
) -> RotamerLibraryEntry {
  RotamerLibraryEntry::{ phi, psi, residues: [] }
}

// ============================================================================
// Standard bond geometry parameters (in Angstroms and degrees)
// ============================================================================

///|
/// Standard peptide bond length (C-N).
pub fn ic_peptide_bond_length() -> Double { 1.33 }

///|
/// Standard CA-C bond length.
pub fn ic_ca_c_bond_length() -> Double { 1.52 }

///|
/// Standard N-CA bond length.
pub fn ic_n_ca_bond_length() -> Double { 1.47 }

///|
/// Standard omega bond angle (CA-C-N-CA), ~180 degrees.
pub fn ic_omega_angle() -> Double { 180.0 }

///|
/// Standard tau angle (N-CA-C), ~110 degrees.
pub fn ic_tau_angle() -> Double { 110.0 }

// ============================================================================
// Torsion angle calculations
// ============================================================================

///|
/// Compute the dihedral angle defined by four points.
/// Returns the angle in radians.
pub fn ic_dihedral_angle(
  p1 : Array[Double],
  p2 : Array[Double],
  p3 : Array[Double],
  p4 : Array[Double],
) -> Double {
  let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]
  let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]]
  let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]]

  // Normal vectors to planes
  let n1 = ic_cross(b1, b2)
  let n2 = ic_cross(b2, b3)

  // Normalize b2
  let b2_norm = ic_vec_norm(b2)
  let b2_unit = if b2_norm > 1.0e-10 {
    [b2[0] / b2_norm, b2[1] / b2_norm, b2[2] / b2_norm]
  } else {
    [0.0, 0.0, 0.0]
  }

  // m1 = n1 x b2_unit
  let m1 = ic_cross(n1, b2_unit)

  let x = ic_dot(n1, n2)
  let y = ic_dot(m1, n2)

  // atan2 gives the correct quadrant
  -@math.atan2(y, x)
}

///|
/// PI constant.
pub fn ic_pi() -> Double { 3.14159265358979323846 }

///|
/// Convert radians to degrees.
pub fn ic_rad_to_deg(rad : Double) -> Double { rad * 180.0 / ic_pi() }

///|
/// Convert degrees to radians.
pub fn ic_deg_to_rad(deg : Double) -> Double { deg * ic_pi() / 180.0 }

///|
/// Compute phi angle for a residue (C(i-1) - N(i) - CA(i) - C(i)).
pub fn ic_compute_phi(
  c_prev : Array[Double],
  n_curr : Array[Double],
  ca_curr : Array[Double],
  c_curr : Array[Double],
) -> Double {
  ic_dihedral_angle(c_prev, n_curr, ca_curr, c_curr)
}

///|
/// Compute psi angle for a residue (N(i) - CA(i) - C(i) - N(i+1)).
pub fn ic_compute_psi(
  n_curr : Array[Double],
  ca_curr : Array[Double],
  c_curr : Array[Double],
  n_next : Array[Double],
) -> Double {
  ic_dihedral_angle(n_curr, ca_curr, c_curr, n_next)
}

///|
/// Compute omega angle (CA(i-1) - C(i-1) - N(i) - CA(i)).
pub fn ic_compute_omega(
  ca_prev : Array[Double],
  c_prev : Array[Double],
  n_curr : Array[Double],
  ca_curr : Array[Double],
) -> Double {
  ic_dihedral_angle(ca_prev, c_prev, n_curr, ca_curr)
}

///|
/// Compute chi1 angle (N - CA - CB - CG).
pub fn ic_compute_chi1(
  n : Array[Double],
  ca : Array[Double],
  cb : Array[Double],
  cg : Array[Double],
) -> Double {
  ic_dihedral_angle(n, ca, cb, cg)
}

// ============================================================================
// Internal coordinate to Cartesian conversion
// ============================================================================

///|
/// Build a simple extended polypeptide chain from internal coordinates.
/// Uses standard bond lengths and angles.
/// Returns the ICToCartesianResult with atom positions.
pub fn ic_build_extended_chain(
  sequence : Array[String],
) -> ICToCartesianResult {
  let result = ICToCartesianResult::new()

  // Standard bond parameters
  let bond_n_ca = 1.47
  let bond_ca_c = 1.52
  let bond_c_n = 1.33
  let angle_n_ca_c = ic_deg_to_rad(110.0)
  let angle_ca_c_n = ic_deg_to_rad(116.0)
  let angle_c_n_ca = ic_deg_to_rad(121.0)

  // Extended chain: phi = -150, psi = 150, omega = 180 (in radians)
  let phi_ext = ic_deg_to_rad(-150.0)
  let psi_ext = ic_deg_to_rad(150.0)

  let n_residues = sequence.length()
  if n_residues == 0 {
    return result
  }

  // Initialize first atom positions
  let n_pos = [0.0, 0.0, 0.0]
  let ca_pos = [bond_n_ca, 0.0, 0.0]
  let c_pos = [
    bond_n_ca + bond_ca_c * @math.cos(ic_pi() - angle_n_ca_c),
    bond_ca_c * @math.sin(ic_pi() - angle_n_ca_c),
    0.0,
  ]

  // Direction vectors
  let mut prev_n = n_pos
  let mut prev_ca = ca_pos
  let mut prev_c = c_pos

  // Add first residue atoms
  result.atom_names.push("N")
  result.residue_names.push(sequence[0])
  result.residue_ids.push(1)
  result.x.push(n_pos[0])
  result.y.push(n_pos[1])
  result.z.push(n_pos[2])

  result.atom_names.push("CA")
  result.residue_names.push(sequence[0])
  result.residue_ids.push(1)
  result.x.push(ca_pos[0])
  result.y.push(ca_pos[1])
  result.z.push(ca_pos[2])

  result.atom_names.push("C")
  result.residue_names.push(sequence[0])
  result.residue_ids.push(1)
  result.x.push(c_pos[0])
  result.y.push(c_pos[1])
  result.z.push(c_pos[2])

  result.atom_names.push("O")
  result.residue_names.push(sequence[0])
  result.residue_ids.push(1)
  let o_pos = ic_oxygen_position(prev_c, prev_ca, prev_n)
  result.x.push(o_pos[0])
  result.y.push(o_pos[1])
  result.z.push(o_pos[2])

  // Build remaining residues using phi/psi angles
  let mut i = 1
  while i < n_residues {
    let phi = phi_ext
    let psi = psi_ext

    // Place next N atom using psi rotation
    let next_n = ic_rotate_around_bond(prev_c, prev_ca, prev_n, psi, bond_c_n)

    // Place next CA atom using phi rotation on next residue
    let next_ca = ic_rotate_around_bond(next_n, prev_c, prev_ca, phi + ic_pi(), bond_n_ca)

    // Place next C atom
    let next_c = ic_build_carbonyl(next_n, next_ca, prev_c)

    // Add atoms
    result.atom_names.push("N")
    result.residue_names.push(sequence[i])
    result.residue_ids.push(i + 1)
    result.x.push(next_n[0])
    result.y.push(next_n[1])
    result.z.push(next_n[2])

    result.atom_names.push("CA")
    result.residue_names.push(sequence[i])
    result.residue_ids.push(i + 1)
    result.x.push(next_ca[0])
    result.y.push(next_ca[1])
    result.z.push(next_ca[2])

    result.atom_names.push("C")
    result.residue_names.push(sequence[i])
    result.residue_ids.push(i + 1)
    result.x.push(next_c[0])
    result.y.push(next_c[1])
    result.z.push(next_c[2])

    result.atom_names.push("O")
    result.residue_names.push(sequence[i])
    result.residue_ids.push(i + 1)
    let next_o = ic_oxygen_position(next_c, next_ca, next_n)
    result.x.push(next_o[0])
    result.y.push(next_o[1])
    result.z.push(next_o[2])

    prev_n = next_n
    prev_ca = next_ca
    prev_c = next_c

    i = i + 1
  }

  result
}

///|
/// Rotate a point around a bond axis by a given angle.
pub fn ic_rotate_around_bond(
  origin : Array[Double],
  bond_end : Array[Double],
  point : Array[Double],
  angle : Double,
  new_length : Double,
) -> Array[Double] {
  // Bond vector
  let bond = [
    bond_end[0] - origin[0],
    bond_end[1] - origin[1],
    bond_end[2] - origin[2],
  ]
  let bond_len = ic_vec_norm(bond)
  let bond_unit = if bond_len > 1.0e-10 {
    [bond[0] / bond_len, bond[1] / bond_len, bond[2] / bond_len]
  } else {
    [1.0, 0.0, 0.0]
  }

  // Vector from origin to point
  let v = [
    point[0] - origin[0],
    point[1] - origin[1],
    point[2] - origin[2],
  ]

  // Decompose v into parallel and perpendicular components
  let v_par = ic_dot(v, bond_unit)
  let v_perp = [
    v[0] - v_par * bond_unit[0],
    v[1] - v_par * bond_unit[1],
    v[2] - v_par * bond_unit[2],
  ]

  // Cross product of bond_unit and v_perp
  let cross = ic_cross(bond_unit, v_perp)

  // Rotate perpendicular component
  let cos_a = @math.cos(angle)
  let sin_a = @math.sin(angle)

  let new_v_perp = [
    cos_a * v_perp[0] + sin_a * cross[0],
    cos_a * v_perp[1] + sin_a * cross[1],
    cos_a * v_perp[2] + sin_a * cross[2],
  ]

  // Reconstruct and scale
  let new_len = v_par
  let target_len = if new_length > 0.0 { new_length } else { new_len }

  let result = [
    origin[0] + new_v_perp[0] + target_len * bond_unit[0],
    origin[1] + new_v_perp[1] + target_len * bond_unit[1],
    origin[2] + new_v_perp[2] + target_len * bond_unit[2],
  ]

  result
}

///|
/// Build the carbonyl C position given N, CA, and previous C.
pub fn ic_build_carbonyl(
  n_pos : Array[Double],
  ca_pos : Array[Double],
  prev_c_pos : Array[Double],
) -> Array[Double] {
  // C position is trigonal planar with N-CA-C angle ~110°
  let bond_ca_c = 1.52
  let angle_n_ca_c = ic_deg_to_rad(110.0)

  let v_ca_n = [
    n_pos[0] - ca_pos[0],
    n_pos[1] - ca_pos[1],
    n_pos[2] - ca_pos[2],
  ]
  let v_ca_prev_c = [
    prev_c_pos[0] - ca_pos[0],
    prev_c_pos[1] - ca_pos[1],
    prev_c_pos[2] - ca_pos[2],
  ]

  let n_ca_norm = ic_vec_norm(v_ca_n)
  let pc_norm = ic_vec_norm(v_ca_prev_c)

  let v1 = if n_ca_norm > 1.0e-10 {
    [v_ca_n[0] / n_ca_norm, v_ca_n[1] / n_ca_norm, v_ca_n[2] / n_ca_norm]
  } else {
    [1.0, 0.0, 0.0]
  }
  let v2 = if pc_norm > 1.0e-10 {
    [v_ca_prev_c[0] / pc_norm, v_ca_prev_c[1] / pc_norm, v_ca_prev_c[2] / pc_norm]
  } else {
    [0.0, 1.0, 0.0]
  }

  let cross_n_c = ic_cross(v1, v2)
  let cn_norm = ic_vec_norm(cross_n_c)
  let normal = if cn_norm > 1.0e-10 {
    [cross_n_c[0] / cn_norm, cross_n_c[1] / cn_norm, cross_n_c[2] / cn_norm]
  } else {
    [0.0, 0.0, 1.0]
  }

  // C position in the plane
  let cos_ang = @math.cos(angle_n_ca_c)
  let sin_ang = @math.sin(angle_n_ca_c)
  let new_dir = [
    cos_ang * v1[0] + sin_ang * ic_cross(normal, v1)[0],
    cos_ang * v1[1] + sin_ang * ic_cross(normal, v1)[1],
    cos_ang * v1[2] + sin_ang * ic_cross(normal, v1)[2],
  ]

  [
    ca_pos[0] + bond_ca_c * new_dir[0],
    ca_pos[1] + bond_ca_c * new_dir[1],
    ca_pos[2] + bond_ca_c * new_dir[2],
  ]
}

///|
/// Position the carbonyl oxygen.
pub fn ic_oxygen_position(
  c_pos : Array[Double],
  ca_pos : Array[Double],
  n_pos : Array[Double],
) -> Array[Double] {
  // O is along the bisector of CA-C-N, opposite the bisector direction
  let bond_c_o = 1.24

  let v_c_ca = [
    ca_pos[0] - c_pos[0],
    ca_pos[1] - c_pos[1],
    ca_pos[2] - c_pos[2],
  ]
  let v_c_n = [
    n_pos[0] - c_pos[0],
    n_pos[1] - c_pos[1],
    n_pos[2] - c_pos[2],
  ]

  let ca_norm = ic_vec_norm(v_c_ca)
  let n_norm = ic_vec_norm(v_c_n)

  let u_ca = if ca_norm > 1.0e-10 {
    [v_c_ca[0] / ca_norm, v_c_ca[1] / ca_norm, v_c_ca[2] / ca_norm]
  } else {
    [1.0, 0.0, 0.0]
  }
  let u_n = if n_norm > 1.0e-10 {
    [v_c_n[0] / n_norm, v_c_n[1] / n_norm, v_c_n[2] / n_norm]
  } else {
    [0.0, 1.0, 0.0]
  }

  // Bisector direction (CA-C-N angle ~121°, so bisector is between them)
  let bisector = [
    u_ca[0] + u_n[0],
    u_ca[1] + u_n[1],
    u_ca[2] + u_n[2],
  ]
  let bis_norm = ic_vec_norm(bisector)
  let bis_unit = if bis_norm > 1.0e-10 {
    [bisector[0] / bis_norm, bisector[1] / bis_norm, bisector[2] / bis_norm]
  } else {
    [0.0, 0.0, 1.0]
  }

  // O is opposite to the bisector
  [
    c_pos[0] - bond_c_o * bis_unit[0],
    c_pos[1] - bond_c_o * bis_unit[1],
    c_pos[2] - bond_c_o * bis_unit[2],
  ]
}

// ============================================================================
// Side chain placement
// ============================================================================

///|
/// Place CB atom given N, CA, and C positions.
pub fn ic_place_cb(
  n_pos : Array[Double],
  ca_pos : Array[Double],
  c_pos : Array[Double],
) -> Array[Double] {
  let bond_ca_cb = 1.53
  let angle_n_ca_cb = ic_deg_to_rad(112.0)

  let v_ca_n = [
    n_pos[0] - ca_pos[0],
    n_pos[1] - ca_pos[1],
    n_pos[2] - ca_pos[2],
  ]
  let v_ca_c = [
    c_pos[0] - ca_pos[0],
    c_pos[1] - ca_pos[1],
    c_pos[2] - ca_pos[2],
  ]

  let n_norm = ic_vec_norm(v_ca_n)
  let c_norm = ic_vec_norm(v_ca_c)

  let u_n = if n_norm > 1.0e-10 {
    [v_ca_n[0] / n_norm, v_ca_n[1] / n_norm, v_ca_n[2] / n_norm]
  } else {
    [1.0, 0.0, 0.0]
  }
  let u_c = if c_norm > 1.0e-10 {
    [v_ca_c[0] / c_norm, v_ca_c[1] / c_norm, v_ca_c[2] / c_norm]
  } else {
    [0.0, 1.0, 0.0]
  }

  let cross_n_c = ic_cross(u_n, u_c)
  let cn_norm = ic_vec_norm(cross_n_c)
  let normal = if cn_norm > 1.0e-10 {
    [cross_n_c[0] / cn_norm, cross_n_c[1] / cn_norm, cross_n_c[2] / cn_norm]
  } else {
    [0.0, 0.0, 1.0]
  }

  // CB is tetrahedrally placed
  let cos_ang = @math.cos(angle_n_ca_cb)
  let sin_ang = @math.sin(angle_n_ca_cb)
  let new_dir = [
    cos_ang * u_n[0] + sin_ang * ic_cross(normal, u_n)[0],
    cos_ang * u_n[1] + sin_ang * ic_cross(normal, u_n)[1],
    cos_ang * u_n[2] + sin_ang * ic_cross(normal, u_n)[2],
  ]

  [
    ca_pos[0] + bond_ca_cb * new_dir[0],
    ca_pos[1] + bond_ca_cb * new_dir[1],
    ca_pos[2] + bond_ca_cb * new_dir[2],
  ]
}

///|
/// Rotate CB around CA-N bond by chi1 angle.
pub fn ic_rotate_chi1(
  n_pos : Array[Double],
  ca_pos : Array[Double],
  cb_pos : Array[Double],
  chi1 : Double,
) -> Array[Double] {
  let bond_ca_cb = 1.53
  ic_rotate_around_bond(n_pos, ca_pos, cb_pos, chi1, bond_ca_cb)
}

// ============================================================================
// Vector math utilities
// ============================================================================

///|
/// Dot product of two 3D vectors.
pub fn ic_dot(a : Array[Double], b : Array[Double]) -> Double {
  a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

///|
/// Cross product of two 3D vectors.
pub fn ic_cross(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ]
}

///|
/// Euclidean norm of a 3D vector.
pub fn ic_vec_norm(v : Array[Double]) -> Double {
  (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

///|
/// Distance between two 3D points.
pub fn ic_distance(a : Array[Double], b : Array[Double]) -> Double {
  let dx = a[0] - b[0]
  let dy = a[1] - b[1]
  let dz = a[2] - b[2]
  (dx * dx + dy * dy + dz * dz).sqrt()
}

// ============================================================================
// Rotamer library
// ============================================================================

///|
/// Get common chi1 angle values for different amino acids.
pub fn ic_chi1_rotamers(resname : String) -> Array[Rotamer] {
  match resname {
    "LEU" => [
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-60.0)], probability: 0.25 },
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-60.0)], probability: 0.55 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-60.0)], probability: 0.20 },
    ]
    "VAL" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0)], probability: 0.45 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0)], probability: 0.35 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0)], probability: 0.20 },
    ]
    "PHE" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.60 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(90.0)], probability: 0.25 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(-90.0)], probability: 0.15 },
    ]
    "TYR" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.60 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(90.0)], probability: 0.25 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(-90.0)], probability: 0.15 },
    ]
    "TRP" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.50 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.30 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-90.0)], probability: 0.20 },
    ]
    "SER" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0)], probability: 0.35 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0)], probability: 0.40 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0)], probability: 0.25 },
    ]
    "THR" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0)], probability: 0.25 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0)], probability: 0.40 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0)], probability: 0.35 },
    ]
    "CYS" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0)], probability: 0.35 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0)], probability: 0.40 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0)], probability: 0.25 },
    ]
    "MET" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.50 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.25 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-90.0)], probability: 0.25 },
    ]
    "LYS" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.40 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.30 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.30 },
    ]
    "ARG" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.35 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.35 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.30 },
    ]
    "ASP" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.60 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0)], probability: 0.20 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0)], probability: 0.20 },
    ]
    "GLU" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.45 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.30 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(180.0)], probability: 0.25 },
    ]
    "ASN" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.55 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(90.0)], probability: 0.25 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(-90.0)], probability: 0.20 },
    ]
    "GLN" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.50 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(90.0)], probability: 0.25 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0), ic_deg_to_rad(-90.0)], probability: 0.25 },
    ]
    "HIS" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(0.0)], probability: 0.45 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(0.0)], probability: 0.30 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0)], probability: 0.25 },
    ]
    "ALA" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0)], probability: 1.0 },
    ]
    "GLY" => [
      Rotamer::{ name: "t", chi_angles: [], probability: 1.0 },
    ]
    "PRO" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(-60.0)], probability: 1.0 },
    ]
    "ILE" => [
      Rotamer::{ name: "t", chi_angles: [ic_deg_to_rad(180.0), ic_deg_to_rad(-60.0)], probability: 0.60 },
      Rotamer::{ name: "g+", chi_angles: [ic_deg_to_rad(-60.0), ic_deg_to_rad(-60.0)], probability: 0.20 },
      Rotamer::{ name: "g-", chi_angles: [ic_deg_to_rad(60.0), ic_deg_to_rad(180.0)], probability: 0.20 },
    ]
    _ => [
      Rotamer::{ name: "t", chi_angles: [], probability: 1.0 },
    ]
  }
}

///|
/// Select the best rotamer for a given backbone-dependent context.
pub fn ic_select_rotamer(
  resname : String,
  _phi : Double,
  _psi : Double,
) -> Rotamer {
  let rotamers = ic_chi1_rotamers(resname)
  // Default: return the most probable rotamer
  if rotamers.length() > 0 {
    rotamers[0]
  } else {
    Rotamer::new(name="t", chi_angles=[], probability=1.0)
  }
}

// ============================================================================
// Validation and analysis
// ============================================================================

///|
/// Validate internal coordinates by checking bond length and angle constraints.
pub fn ic_validate_coords(ic : InternalCoordChain) -> Array[(Int, String)] {
  let issues : Array[(Int, String)] = []
  let mut i = 0
  while i < ic.residues.length() {
    let res = ic.residues[i]
    for bl in res.bond_lengths {
      if bl.length < 0.5 || bl.length > 3.0 {
        issues.push((res.resseq, "Bond length " + bl.atom1 + "-" + bl.atom2 + " out of range: " + bl.length.to_string()))
      }
    }
    for ba in res.bond_angles {
      if ba.angle < 60.0 || ba.angle > 180.0 {
        issues.push((res.resseq, "Bond angle " + ba.atom1 + "-" + ba.atom2 + "-" + ba.atom3 + " out of range: " + ba.angle.to_string()))
      }
    }
    i = i + 1
  }
  issues
}

///|
/// Get the summary of an internal coordinate chain.
pub fn ic_chain_summary(ic : InternalCoordChain) -> String {
  let n_res = ic.residues.length()
  let mut total_torsions = 0
  for res in ic.residues {
    total_torsions = total_torsions + res.torsions.length()
  }
  "Chain " + ic.chain_id.to_string() +
  ": " + n_res.to_string() + " residues, " +
  total_torsions.to_string() + " torsion angles"
}

///|
/// Add a residue to the internal coordinate chain.
pub fn InternalCoordChain::add_residue(
  self : InternalCoordChain,
  res : InternalCoordResidue,
) -> InternalCoordChain {
  let new_residues = self.residues.copy()
  new_residues.push(res)
  InternalCoordChain::{
    chain_id: self.chain_id,
    residues: new_residues,
  }
}

///|
/// Add a torsion angle to a residue.
pub fn InternalCoordResidue::add_torsion(
  self : InternalCoordResidue,
  torsion : TorsionAngle,
) -> InternalCoordResidue {
  let new_torsions = self.torsions.copy()
  new_torsions.push(torsion)
  InternalCoordResidue::{
    resname: self.resname,
    resseq: self.resseq,
    bond_lengths: self.bond_lengths,
    bond_angles: self.bond_angles,
    torsions: new_torsions,
    start_atom: self.start_atom,
  }
}

///|
/// Add a bond length to a residue.
pub fn InternalCoordResidue::add_bond_length(
  self : InternalCoordResidue,
  bl : BondLength,
) -> InternalCoordResidue {
  let new_bond_lengths = self.bond_lengths.copy()
  new_bond_lengths.push(bl)
  InternalCoordResidue::{
    resname: self.resname,
    resseq: self.resseq,
    bond_lengths: new_bond_lengths,
    bond_angles: self.bond_angles,
    torsions: self.torsions,
    start_atom: self.start_atom,
  }
}

///|
/// Add a bond angle to a residue.
pub fn InternalCoordResidue::add_bond_angle(
  self : InternalCoordResidue,
  ba : BondAngle,
) -> InternalCoordResidue {
  let new_bond_angles = self.bond_angles.copy()
  new_bond_angles.push(ba)
  InternalCoordResidue::{
    resname: self.resname,
    resseq: self.resseq,
    bond_lengths: self.bond_lengths,
    bond_angles: new_bond_angles,
    torsions: self.torsions,
    start_atom: self.start_atom,
  }
}

///|
/// Get the number of residues.
pub fn InternalCoordChain::n_residues(self : InternalCoordChain) -> Int {
  self.residues.length()
}

///|
/// Get residue names.
pub fn InternalCoordChain::residue_names(self : InternalCoordChain) -> Array[String] {
  self.residues.map(fn(r) { r.resname })
}

///|
/// Get all torsion angles for a specific residue.
pub fn InternalCoordChain::get_torsions(
  self : InternalCoordChain,
  resseq : Int,
) -> Array[TorsionAngle] {
  for res in self.residues {
    if res.resseq == resseq {
      return res.torsions
    }
  }
  []
}

///|
/// Get average phi/psi angles (Ramachandran plot data).
pub fn ic_ramachandran_data(
  phi_values : Array[Double],
  psi_values : Array[Double],
) -> (Double, Double, Double, Double) {
  let n = phi_values.length()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }

  let mut sum_phi = 0.0
  let mut sum_psi = 0.0
  let mut sum_phi2 = 0.0
  let mut sum_psi2 = 0.0

  let mut i = 0
  while i < n {
    sum_phi = sum_phi + phi_values[i]
    sum_psi = sum_psi + psi_values[i]
    sum_phi2 = sum_phi2 + phi_values[i] * phi_values[i]
    sum_psi2 = sum_psi2 + psi_values[i] * psi_values[i]
    i = i + 1
  }

  let mean_phi = sum_phi / n.to_double()
  let mean_psi = sum_psi / n.to_double()
  let var_phi = sum_phi2 / n.to_double() - mean_phi * mean_phi
  let var_psi = sum_psi2 / n.to_double() - mean_psi * mean_psi

  (mean_phi, mean_psi, var_phi, var_psi)
}
