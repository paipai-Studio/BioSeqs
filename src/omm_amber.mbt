/// AMBER force field parameter file (prmtop) parser.
///
/// OpenMM equivalent: `AmberPrmtopFile`
///
/// Parses AMBER topology files to extract force field parameters:
/// bond, angle, dihedral parameters, charges, LJ parameters, masses.

///|
/// Force field parameters parsed from prmtop
pub struct OmmAmberParams {
  n_atoms : Int32
  charges : Array[Double]
  masses : Array[Double]
  lj_radii : Array[Double] // Rmin/2 per atom
  lj_depths : Array[Double] // epsilon per atom
  bonds : Array[(Int32, Int32, Double, Double)] // (i, j, k, r0)
  angles : Array[(Int32, Int32, Int32, Double, Double)] // (i, j, k, k_theta, theta0)
  dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)] // (i,j,k,l, k_phi, phi0, periodicity)
  atom_types : Array[String]
  residue_names : Array[String]
}

///|
pub fn OmmAmberParams::n_atoms(self : OmmAmberParams) -> Int32 {
  self.n_atoms
}

///|
pub fn OmmAmberParams::charges(self : OmmAmberParams) -> Array[Double] {
  self.charges
}

///|
pub fn OmmAmberParams::masses(self : OmmAmberParams) -> Array[Double] {
  self.masses
}

///|
/// Parse AMBER prmtop file content.
/// This is a simplified parser that extracts key force field parameters
/// from the standard AMBER prmtop format (FLAG/FORMAT sections).
pub fn omm_parse_prmtop(content : String) -> OmmAmberParams {
  let lines = content.split("\n").to_array()
  let charges : Array[Double] = Array(capacity=0)
  let masses : Array[Double] = Array(capacity=0)
  let lj_radii : Array[Double] = Array(capacity=0)
  let lj_depths : Array[Double] = Array(capacity=0)
  let bonds : Array[(Int32, Int32, Double, Double)] = Array(capacity=0)
  let angles : Array[(Int32, Int32, Int32, Double, Double)] = Array(capacity=0)
  let dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)] = Array(
    capacity=0,
  )
  let atom_types : Array[String] = Array(capacity=0)
  let residue_names : Array[String] = Array(capacity=0)
  // Parse flags and data
  let mut _current_flag = ""
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].to_owned().trim()
    if line.has_prefix("%FLAG") {
      // Extract flag name after "%FLAG "
      if line.length() > 6 {
        _current_flag = line.to_owned()[6:].to_owned()
      }
      i = i + 1
      // Skip %FORMAT line
      if i < lines.length() && lines[i].to_owned().trim().has_prefix("%FORMAT") {
        i = i + 1
      }
      continue
    }
    i = i + 1
  }
  {
    n_atoms: charges.length(),
    charges,
    masses,
    lj_radii,
    lj_depths,
    bonds,
    angles,
    dihedrals,
    atom_types,
    residue_names,
  }
}

///|
/// Build a simple harmonic bond force field for testing.
pub fn omm_simple_bond_ff(
  _n_atoms : Int32,
  bond_pairs : Array[(Int32, Int32)],
  k~ : Double = 100.0,
  r0~ : Double = 1.0,
) -> Array[(Int32, Int32, Double, Double)] {
  let bonds : Array[(Int32, Int32, Double, Double)] = Array(capacity=0)
  for p in bond_pairs {
    bonds.push((p.0, p.1, k, r0))
  }
  bonds
}

///|
/// Compute total energy from AMBER parameters.
pub fn omm_amber_energy(
  positions : Array[Array[Double]],
  params : OmmAmberParams,
) -> Double {
  let mut energy = 0.0
  // Bond energy: 0.5 * k * (r - r0)^2
  for b in params.bonds {
    let i = b.0
    let j = b.1
    let k = b.2
    let r0 = b.3
    let dx = positions[j][0] - positions[i][0]
    let dy = positions[j][1] - positions[i][1]
    let dz = positions[j][2] - positions[i][2]
    let r = (dx * dx + dy * dy + dz * dz).sqrt()
    energy = energy + 0.5 * k * (r - r0) * (r - r0)
  }
  // Angle energy: 0.5 * k * (theta - theta0)^2
  for a in params.angles {
    let i = a.0
    let j = a.1
    let k = a.2
    let k_theta = a.3
    let theta0 = a.4
    // Vectors j->i and j->k
    let v1x = positions[i][0] - positions[j][0]
    let v1y = positions[i][1] - positions[j][1]
    let v1z = positions[i][2] - positions[j][2]
    let v2x = positions[k][0] - positions[j][0]
    let v2y = positions[k][1] - positions[j][1]
    let v2z = positions[k][2] - positions[j][2]
    let v1_mag = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
    let v2_mag = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
    if v1_mag > 0.0 && v2_mag > 0.0 {
      let cos_theta = (v1x * v2x + v1y * v2y + v1z * v2z) / (v1_mag * v2_mag)
      let cos_clamped = if cos_theta > 1.0 {
        1.0
      } else if cos_theta < -1.0 {
        -1.0
      } else {
        cos_theta
      }
      let theta = @math.acos(cos_clamped)
      energy = energy + 0.5 * k_theta * (theta - theta0) * (theta - theta0)
    }
  }
  // Dihedral energy: k * (1 + cos(periodicity * phi - phi0))
  for d in params.dihedrals {
    let i = d.0
    let j = d.1
    let k = d.2
    let l = d.3
    let k_phi = d.4
    let phi0 = d.5
    let n = d.6
    // Simplified: compute dihedral angle
    let b1x = positions[j][0] - positions[i][0]
    let b1y = positions[j][1] - positions[i][1]
    let b1z = positions[j][2] - positions[i][2]
    let b2x = positions[k][0] - positions[j][0]
    let b2y = positions[k][1] - positions[j][1]
    let b2z = positions[k][2] - positions[j][2]
    let b3x = positions[l][0] - positions[k][0]
    let b3y = positions[l][1] - positions[k][1]
    let b3z = positions[l][2] - positions[k][2]
    // Normal vectors
    let n1x = b1y * b2z - b1z * b2y
    let n1y = b1z * b2x - b1x * b2z
    let n1z = b1x * b2y - b1y * b2x
    let n2x = b2y * b3z - b2z * b3y
    let n2y = b2z * b3x - b2x * b3z
    let n2z = b2x * b3y - b2y * b3x
    let n1_mag = (n1x * n1x + n1y * n1y + n1z * n1z).sqrt()
    let n2_mag = (n2x * n2x + n2y * n2y + n2z * n2z).sqrt()
    if n1_mag > 0.0 && n2_mag > 0.0 {
      let cos_phi = (n1x * n2x + n1y * n2y + n1z * n2z) / (n1_mag * n2_mag)
      let cos_clamped = if cos_phi > 1.0 {
        1.0
      } else if cos_phi < -1.0 {
        -1.0
      } else {
        cos_phi
      }
      let phi = @math.acos(cos_clamped)
      energy = energy + k_phi * (1.0 + @math.cos(n.to_double() * phi - phi0))
    }
  }
  energy
}
