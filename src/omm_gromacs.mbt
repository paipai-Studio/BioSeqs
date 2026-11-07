/// omm_gromacs.mbt — GROMACS TOP/GRO parser and energy
/// Item 2: GROMACS force field support

///|
pub struct OmmGromacsTop {
  n_atoms : Int32
  atom_names : Array[String]
  atom_types : Array[String]
  charges : Array[Double]
  masses : Array[Double]
  residues : Array[String]
  // Bonded interactions: (i, j, func, k, r0)
  bonds : Array[(Int32, Int32, Int32, Double, Double)]
  // Angles: (i, j, k, func, theta0, k)
  angles : Array[(Int32, Int32, Int32, Int32, Double, Double)]
  // Dihedrals: (i, j, k, l, func, phi0, k, n)
  dihedrals : Array[(Int32, Int32, Int32, Int32, Int32, Double, Double, Int32)]
}

///|
pub struct OmmGromacsGro {
  n_atoms : Int32
  coords : Array[(Double, Double, Double)]
  velocities : Array[(Double, Double, Double)]
  box : (Double, Double, Double)
}

///|
/// Parse a simplified GROMACS .top file.
pub fn omm_parse_gromacs_top(content : String) -> OmmGromacsTop {
  let lines = content.split("\n").to_array()
  let atom_names : Array[String] = Array(capacity=0)
  let atom_types : Array[String] = Array(capacity=0)
  let charges : Array[Double] = Array(capacity=0)
  let masses : Array[Double] = Array(capacity=0)
  let residues : Array[String] = Array(capacity=0)
  let bonds : Array[(Int32, Int32, Int32, Double, Double)] = Array(capacity=0)
  let angles : Array[(Int32, Int32, Int32, Int32, Double, Double)] = Array(
    capacity=0,
  )
  let dihedrals : Array[
    (Int32, Int32, Int32, Int32, Int32, Double, Double, Int32),
  ] = Array(capacity=0)
  let mut section = "none"
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned().trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if line.has_prefix("[ atoms ]") {
      section = "atoms"
      continue
    }
    if line.has_prefix("[ bonds ]") {
      section = "bonds"
      continue
    }
    if line.has_prefix("[ angles ]") {
      section = "angles"
      continue
    }
    if line.has_prefix("[ dihedrals ]") {
      section = "dihedrals"
      continue
    }
    if line.has_prefix("[") {
      section = "other"
      continue
    }
    if line.has_prefix(";") {
      continue
    }
    if section == "atoms" {
      let parts = line.split(" ").to_array()
      // nr type resnr residue atom cgnr charge mass
      if parts.length() >= 7 {
        atom_types.push(parts[1].to_owned())
        residues.push(parts[3].to_owned())
        atom_names.push(parts[4].to_owned())
        charges.push(parse_double(parts[6].to_owned()).unwrap_or(0.0))
        if parts.length() >= 8 {
          masses.push(parse_double(parts[7].to_owned()).unwrap_or(12.0))
        } else {
          masses.push(12.0)
        }
      }
    }
    if section == "bonds" {
      let parts = line.split(" ").to_array()
      if parts.length() >= 5 {
        let a1 = parse_int(parts[0].to_owned()) - 1
        let a2 = parse_int(parts[1].to_owned()) - 1
        let func = parse_int(parts[2].to_owned())
        let k = parse_double(parts[3].to_owned()).unwrap_or(0.0)
        let r0 = parse_double(parts[4].to_owned()).unwrap_or(0.0)
        bonds.push((a1, a2, func, k, r0))
      }
    }
    if section == "angles" {
      let parts = line.split(" ").to_array()
      if parts.length() >= 6 {
        let a1 = parse_int(parts[0].to_owned()) - 1
        let a2 = parse_int(parts[1].to_owned()) - 1
        let a3 = parse_int(parts[2].to_owned()) - 1
        let func = parse_int(parts[3].to_owned())
        let th0 = parse_double(parts[4].to_owned()).unwrap_or(0.0)
        let k = parse_double(parts[5].to_owned()).unwrap_or(0.0)
        angles.push((a1, a2, a3, func, th0, k))
      }
    }
    if section == "dihedrals" {
      let parts = line.split(" ").to_array()
      if parts.length() >= 8 {
        let a1 = parse_int(parts[0].to_owned()) - 1
        let a2 = parse_int(parts[1].to_owned()) - 1
        let a3 = parse_int(parts[2].to_owned()) - 1
        let a4 = parse_int(parts[3].to_owned()) - 1
        let func = parse_int(parts[4].to_owned())
        let phi0 = parse_double(parts[5].to_owned()).unwrap_or(0.0)
        let k = parse_double(parts[6].to_owned()).unwrap_or(0.0)
        let n = parse_int(parts[7].to_owned())
        dihedrals.push((a1, a2, a3, a4, func, phi0, k, n))
      }
    }
  }
  OmmGromacsTop::{
    n_atoms: atom_names.length(),
    atom_names,
    atom_types,
    charges,
    masses,
    residues,
    bonds,
    angles,
    dihedrals,
  }
}

///|
/// Parse a GROMACS .gro coordinate file.
pub fn omm_parse_gro(content : String) -> OmmGromacsGro {
  let lines = content.split("\n").to_array()
  let coords : Array[(Double, Double, Double)] = Array(capacity=0)
  let velocities : Array[(Double, Double, Double)] = Array(capacity=0)
  let mut n_atoms = 0
  let mut box = (0.0, 0.0, 0.0)
  if lines.length() >= 2 {
    n_atoms = parse_int(lines[1].to_owned().trim().to_owned())
  }
  // Atom lines start at index 2
  for a = 0; a < n_atoms; a = a + 1 {
    let idx = 2 + a
    if idx < lines.length() {
      let line = lines[idx].to_owned()
      // GRO format: fixed columns, but we use split for simplicity
      // resnr(5) resname(5) atomname(5) atomnr(5) x(8) y(8) z(8) [vx vy vz]
      if line.length() >= 20 {
        // Try fixed-column parse: x at col 20-28, y at 28-36, z at 36-44
        let xs = if line.length() >= 28 {
          line[20:28].to_owned().trim().to_owned()
        } else {
          "0"
        }
        let ys = if line.length() >= 36 {
          line[28:36].to_owned().trim().to_owned()
        } else {
          "0"
        }
        let zs = if line.length() >= 44 {
          line[36:44].to_owned().trim().to_owned()
        } else {
          "0"
        }
        let x = parse_double(xs).unwrap_or(0.0)
        let y = parse_double(ys).unwrap_or(0.0)
        let z = parse_double(zs).unwrap_or(0.0)
        coords.push((x, y, z))
        // Velocities if present
        if line.length() >= 68 {
          let vx = parse_double(line[44:52].to_owned().trim().to_owned()).unwrap_or(
            0.0,
          )
          let vy = parse_double(line[52:60].to_owned().trim().to_owned()).unwrap_or(
            0.0,
          )
          let vz = parse_double(line[60:68].to_owned().trim().to_owned()).unwrap_or(
            0.0,
          )
          velocities.push((vx, vy, vz))
        } else {
          velocities.push((0.0, 0.0, 0.0))
        }
      }
    }
  }
  // Box line
  let box_idx = 2 + n_atoms
  if box_idx < lines.length() {
    let parts = lines[box_idx]
      .to_owned()
      .trim()
      .to_owned()
      .split(" ")
      .to_array()
    if parts.length() >= 3 {
      box = (
        parse_double(parts[0].to_owned()).unwrap_or(0.0),
        parse_double(parts[1].to_owned()).unwrap_or(0.0),
        parse_double(parts[2].to_owned()).unwrap_or(0.0),
      )
    }
  }
  OmmGromacsGro::{ n_atoms, coords, velocities, box, }
}

///|
/// Compute GROMACS bond energy (harmonic: 0.5 * k * (r-r0)^2).
pub fn omm_gromacs_bond_energy(
  top : OmmGromacsTop,
  coords : Array[(Double, Double, Double)],
) -> Double {
  let mut energy = 0.0
  for b in top.bonds {
    let (i, j, _func, k, r0) = b
    if i >= 0 && j >= 0 && i < coords.length() && j < coords.length() {
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = omm_gr_abs(dx * dx + dy * dy + dz * dz)
      let dr = r - r0
      energy = energy + 0.5 * k * dr * dr
    }
  }
  energy
}

///|
fn omm_gr_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
pub fn omm_gromacs_dummy() -> Unit {
  ignore(OmmGromacsTop::{
    n_atoms: 0,
    atom_names: Array(capacity=0),
    atom_types: Array(capacity=0),
    charges: Array(capacity=0),
    masses: Array(capacity=0),
    residues: Array(capacity=0),
    bonds: Array(capacity=0),
    angles: Array(capacity=0),
    dihedrals: Array(capacity=0),
  })
  ignore(OmmGromacsGro::{
    n_atoms: 0,
    coords: Array(capacity=0),
    velocities: Array(capacity=0),
    box: (0.0, 0.0, 0.0),
  })
}
