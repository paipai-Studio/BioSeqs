/// omm_charmm.mbt — CHARMM PSF topology + parameter file parser and energy
/// Item 1: CHARMM force field support

///|
pub struct OmmCharmmmParams {
  n_atoms : Int32
  atom_types : Array[String]
  charges : Array[Double]
  masses : Array[Double]
  // Bond parameters: (atom_i, atom_j, k, r0)
  bonds : Array[(Int32, Int32, Double, Double)]
  // Angle parameters: (i, j, k, theta0, k_theta)
  angles : Array[(Int32, Int32, Int32, Double, Double)]
  // Dihedral parameters: (i, j, k, l, phi0, k_dih, n)
  dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)]
}

///|
/// Parse a simplified CHARMM PSF file content.
pub fn omm_parse_psf(content : String) -> OmmCharmmmParams {
  let lines = content.split("\n").to_array()
  let n_atoms = 0
  let atom_types : Array[String] = Array(capacity=0)
  let charges : Array[Double] = Array(capacity=0)
  let masses : Array[Double] = Array(capacity=0)
  let bonds : Array[(Int32, Int32, Double, Double)] = Array(capacity=0)
  let angles : Array[(Int32, Int32, Int32, Double, Double)] = Array(capacity=0)
  let dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)] = Array(
    capacity=0,
  )
  let mut i = 0
  let mut na = n_atoms
  while i < lines.length() {
    let line = lines[i].to_owned().trim().to_owned()
    if line.contains("NATOM") {
      let parts = line.split(" ").to_array()
      // Find the number after NATOM
      for p in parts {
        let pp = p.to_owned().trim().to_owned()
        let v = parse_int(pp)
        if v > 0 {
          na = v
          break
        }
      }
      i = i + 1
      // Read atom lines
      for a = 0; a < na; a = a + 1 {
        if i + a < lines.length() {
          let al = lines[i + a].to_owned().trim().to_owned()
          let ap = al.split(" ").to_array()
          // PSF: index segid resid resname atomname atom_type charge mass
          let mut atype = "C"
          let mut chg = 0.0
          let mut mass = 12.0
          if ap.length() >= 6 {
            atype = ap[5].to_owned().trim().to_owned()
          }
          if ap.length() >= 7 {
            chg = parse_double(ap[6].to_owned()).unwrap_or(0.0)
          }
          if ap.length() >= 8 {
            mass = parse_double(ap[7].to_owned()).unwrap_or(12.0)
          }
          atom_types.push(atype)
          charges.push(chg)
          masses.push(mass)
        }
      }
      i = i + na
      continue
    }
    if line.contains("NBOND") {
      i = i + 1
      // Bond lines: pairs of atom indices
      while i < lines.length() {
        let bl = lines[i].to_owned().trim().to_owned()
        if bl.length() == 0 || bl.contains("NTHETA") || bl.contains("NPHI") {
          break
        }
        let bp = bl.split(" ").to_array()
        let mut j = 0
        while j + 1 < bp.length() {
          let a1 = parse_int(bp[j].to_owned())
          let a2 = parse_int(bp[j + 1].to_owned())
          if a1 > 0 && a2 > 0 {
            bonds.push((a1 - 1, a2 - 1, 200.0, 1.5))
          }
          j = j + 2
        }
        i = i + 1
      }
      continue
    }
    i = i + 1
  }
  OmmCharmmmParams::{
    n_atoms: na,
    atom_types,
    charges,
    masses,
    bonds,
    angles,
    dihedrals,
  }
}

///|
/// Compute CHARMM bond energy from parameters and coordinates.
pub fn omm_charmm_bond_energy(
  params : OmmCharmmmParams,
  coords : Array[(Double, Double, Double)],
) -> Double {
  let mut energy = 0.0
  for b in params.bonds {
    let (i, j, k, r0) = b
    if i < coords.length() && j < coords.length() {
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = omm_ch_abs(dx * dx + dy * dy + dz * dz)
      let dr = r - r0
      energy = energy + k * dr * dr
    }
  }
  energy
}

///|
/// Compute CHARMM electrostatic energy.
pub fn omm_charmm_elec_energy(
  params : OmmCharmmmParams,
  coords : Array[(Double, Double, Double)],
  cutoff~ : Double = 12.0,
) -> Double {
  let mut energy = 0.0
  let n = params.n_atoms
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = omm_ch_abs(dx * dx + dy * dy + dz * dz)
      if r > cutoff || r < 0.0000000001 {
        continue
      }
      energy = energy + 138.935485 * params.charges[i] * params.charges[j] / r
    }
  }
  energy
}

///|
/// Parse CHARMM parameter file (simplified — bond/angle/dihedral params).
pub fn omm_parse_charmm_params(
  content : String,
) -> (
  Array[(String, String, Double, Double)],
  Array[(String, String, String, Double, Double)],
) {
  let lines = content.split("\n").to_array()
  let bond_params : Array[(String, String, Double, Double)] = Array(capacity=0)
  let angle_params : Array[(String, String, String, Double, Double)] = Array(
    capacity=0,
  )
  let mut section = "none"
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned().trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if line.contains("BONDS") {
      section = "bonds"
      continue
    }
    if line.contains("ANGLES") {
      section = "angles"
      continue
    }
    if line.contains("DIHEDRALS") {
      section = "dihedrals"
      continue
    }
    if section == "bonds" {
      let parts = line.split(" ").to_array()
      if parts.length() >= 4 {
        let t1 = parts[0].to_owned()
        let t2 = parts[1].to_owned()
        let k = parse_double(parts[2].to_owned()).unwrap_or(0.0)
        let r0 = parse_double(parts[3].to_owned()).unwrap_or(0.0)
        bond_params.push((t1, t2, k, r0))
      }
    }
    if section == "angles" {
      let parts = line.split(" ").to_array()
      if parts.length() >= 5 {
        let t1 = parts[0].to_owned()
        let t2 = parts[1].to_owned()
        let t3 = parts[2].to_owned()
        let k = parse_double(parts[3].to_owned()).unwrap_or(0.0)
        let th0 = parse_double(parts[4].to_owned()).unwrap_or(0.0)
        angle_params.push((t1, t2, t3, k, th0))
      }
    }
  }
  (bond_params, angle_params)
}

///|
fn omm_ch_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
pub fn omm_charmm_dummy() -> Unit {
  ignore(OmmCharmmmParams::{
    n_atoms: 0,
    atom_types: Array(capacity=0),
    charges: Array(capacity=0),
    masses: Array(capacity=0),
    bonds: Array(capacity=0),
    angles: Array(capacity=0),
    dihedrals: Array(capacity=0),
  })
}
