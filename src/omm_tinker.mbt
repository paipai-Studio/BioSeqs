/// omm_tinker.mbt — Tinker parameter/XYZ parser
/// Item 4: Tinker force field support

///|
pub struct OmmTinkerParams {
  // Atom type parameters: (type_id, element, atomic_mass, charge)
  atom_types : Array[(Int32, String, Double, Double)]
  // Bond parameters: (type_i, type_j, k, r0)
  bond_params : Array[(Int32, Int32, Double, Double)]
  // Angle parameters: (ti, tj, tk, k, theta0)
  angle_params : Array[(Int32, Int32, Int32, Double, Double)]
  // Torsion parameters: (ti, tj, tk, tl, k, phi0, n)
  torsion_params : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)]
  // van der Waals: (type_id, radius, epsilon)
  vdw_params : Array[(Int32, Double, Double)]
}

///|
pub struct OmmTinkerXYZ {
  n_atoms : Int32
  title : String
  atoms : Array[(String, Double, Double, Double, Int32)] // (element, x, y, z, type)
  bonds : Array[Array[Int32]] // connectivity list per atom
}

///|
/// Parse a Tinker .params file (simplified).
pub fn omm_parse_tinker_params(content : String) -> OmmTinkerParams {
  let lines = content.split("\n").to_array()
  let atom_types : Array[(Int32, String, Double, Double)] = Array(capacity=0)
  let bond_params : Array[(Int32, Int32, Double, Double)] = Array(capacity=0)
  let angle_params : Array[(Int32, Int32, Int32, Double, Double)] = Array(
    capacity=0,
  )
  let torsion_params : Array[
    (Int32, Int32, Int32, Int32, Double, Double, Int32),
  ] = Array(capacity=0)
  let vdw_params : Array[(Int32, Double, Double)] = Array(capacity=0)
  let mut section = "none"
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned().trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if line.has_prefix("atom") {
      section = "atom"
      continue
    }
    if line.has_prefix("bond") {
      section = "bond"
      continue
    }
    if line.has_prefix("angle") {
      section = "angle"
      continue
    }
    if line.has_prefix("torsion") {
      section = "torsion"
      continue
    }
    if line.has_prefix("vdw") {
      section = "vdw"
      continue
    }
    let parts = line.split(" ").to_array()
    if section == "atom" && parts.length() >= 4 {
      let tid = parse_int(parts[0].to_owned())
      let el = parts[1].to_owned()
      let mass = parse_double(parts[2].to_owned()).unwrap_or(12.0)
      let chg = if parts.length() >= 5 {
        parse_double(parts[3].to_owned()).unwrap_or(0.0)
      } else {
        0.0
      }
      atom_types.push((tid, el, mass, chg))
    }
    if section == "bond" && parts.length() >= 4 {
      let ti = parse_int(parts[0].to_owned())
      let tj = parse_int(parts[1].to_owned())
      let k = parse_double(parts[2].to_owned()).unwrap_or(0.0)
      let r0 = parse_double(parts[3].to_owned()).unwrap_or(0.0)
      bond_params.push((ti, tj, k, r0))
    }
    if section == "angle" && parts.length() >= 5 {
      let ti = parse_int(parts[0].to_owned())
      let tj = parse_int(parts[1].to_owned())
      let tk = parse_int(parts[2].to_owned())
      let k = parse_double(parts[3].to_owned()).unwrap_or(0.0)
      let th0 = parse_double(parts[4].to_owned()).unwrap_or(0.0)
      angle_params.push((ti, tj, tk, k, th0))
    }
    if section == "torsion" && parts.length() >= 7 {
      let ti = parse_int(parts[0].to_owned())
      let tj = parse_int(parts[1].to_owned())
      let tk = parse_int(parts[2].to_owned())
      let tl = parse_int(parts[3].to_owned())
      let k = parse_double(parts[4].to_owned()).unwrap_or(0.0)
      let phi0 = parse_double(parts[5].to_owned()).unwrap_or(0.0)
      let n = parse_int(parts[6].to_owned())
      torsion_params.push((ti, tj, tk, tl, k, phi0, n))
    }
    if section == "vdw" && parts.length() >= 3 {
      let tid = parse_int(parts[0].to_owned())
      let r = parse_double(parts[1].to_owned()).unwrap_or(0.0)
      let e = parse_double(parts[2].to_owned()).unwrap_or(0.0)
      vdw_params.push((tid, r, e))
    }
  }
  OmmTinkerParams::{
    atom_types,
    bond_params,
    angle_params,
    torsion_params,
    vdw_params,
  }
}

///|
/// Parse a Tinker .xyz coordinate file.
pub fn omm_parse_tinker_xyz(content : String) -> OmmTinkerXYZ {
  let lines = content.split("\n").to_array()
  let atoms : Array[(String, Double, Double, Double, Int32)] = Array(capacity=0)
  let bonds : Array[Array[Int32]] = Array(capacity=0)
  let mut title = ""
  let mut n_atoms = 0
  if lines.length() >= 1 {
    n_atoms = parse_int(lines[0].to_owned().trim().to_owned())
  }
  if lines.length() >= 2 {
    title = lines[1].to_owned().trim().to_owned()
  }
  for a = 0; a < n_atoms; a = a + 1 {
    let idx = 2 + a
    if idx < lines.length() {
      let parts = lines[idx].to_owned().trim().to_owned().split(" ").to_array()
      if parts.length() >= 6 {
        let el = parts[1].to_owned()
        let x = parse_double(parts[2].to_owned()).unwrap_or(0.0)
        let y = parse_double(parts[3].to_owned()).unwrap_or(0.0)
        let z = parse_double(parts[4].to_owned()).unwrap_or(0.0)
        let tid = parse_int(parts[5].to_owned())
        atoms.push((el, x, y, z, tid))
        // Bonded atoms: remaining fields
        let bond_list : Array[Int32] = Array(capacity=0)
        for b = 6; b < parts.length(); b = b + 1 {
          let bv = parse_int(parts[b].to_owned())
          if bv > 0 {
            bond_list.push(bv - 1)
          }
        }
        bonds.push(bond_list)
      }
    }
  }
  OmmTinkerXYZ::{ n_atoms, title, atoms, bonds, }
}

///|
/// Compute Tinker bond energy.
pub fn omm_tinker_bond_energy(
  params : OmmTinkerParams,
  xyz : OmmTinkerXYZ,
) -> Double {
  let mut energy = 0.0
  for a = 0; a < xyz.n_atoms; a = a + 1 {
    if a >= xyz.atoms.length() || a >= xyz.bonds.length() {
      continue
    }
    let (_, ax, ay, az, atype) = xyz.atoms[a]
    for b in xyz.bonds[a] {
      if b <= a || b >= xyz.atoms.length() {
        continue
      }
      let (_, bx, by, bz, btype) = xyz.atoms[b]
      let dx = ax - bx
      let dy = ay - by
      let dz = az - bz
      let r = omm_tk_abs(dx * dx + dy * dy + dz * dz)
      // Look up bond params
      for (ti, tj, k, r0) in params.bond_params {
        if (ti == atype && tj == btype) || (ti == btype && tj == atype) {
          let dr = r - r0
          energy = energy + k * dr * dr
          break
        }
      }
    }
  }
  energy
}

///|
fn omm_tk_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
pub fn omm_tinker_dummy() -> Unit {
  ignore(OmmTinkerParams::{
    atom_types: Array(capacity=0),
    bond_params: Array(capacity=0),
    angle_params: Array(capacity=0),
    torsion_params: Array(capacity=0),
    vdw_params: Array(capacity=0),
  })
  ignore(OmmTinkerXYZ::{
    n_atoms: 0,
    title: "",
    atoms: Array(capacity=0),
    bonds: Array(capacity=0),
  })
}
