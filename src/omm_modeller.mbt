/// omm_modeller.mbt — OpenMM `Modeller` equivalent: system construction and
/// editing tools.
///
/// Provides `add_solvent` (periodic water box filling + counter-ion
/// neutralization), `add_hydrogens` (valence-rule protonation),
/// `center_in_box` and `delete_water`, mirroring
/// `openmm.app.Modeller.addSolvent / addHydrogens / centerInBox`.
///
/// Water model: TIP3P (r_OH = 0.9572 A, angle = 104.52 deg, q_O = -0.834,
/// q_H = +0.417, sigma_O = 3.1507 A, epsilon_O = 0.1521 kcal/mol).
/// Units: A, kcal/mol, amu.

///|
const OMM_TIP3P_R_OH : Double = 0.9572

///|
const OMM_TIP3P_HALF_ANGLE : Double = 52.26 // degrees, half of 104.52

///|
const OMM_TIP3P_Q_O : Double = -0.834

///|
const OMM_TIP3P_Q_H : Double = 0.417

///|
const OMM_TIP3P_SIGMA_O : Double = 3.1507

///|
const OMM_TIP3P_EPS_O : Double = 0.1521

///|
const OMM_MASS_O : Double = 15.9994

///|
const OMM_MASS_H : Double = 1.008

///|
/// One atom in the model being built (per-atom metadata + LJ parameters).
pub(all) struct OmmModelAtom {
  name : String
  resname : String
  resid : Int32
  element : String
  mass : Double
  charge : Double
  sigma : Double
  epsilon : Double
}

///|
pub fn OmmModelAtom::new(
  name~ : String = "X",
  resname~ : String = "UNK",
  resid~ : Int32 = 1,
  element~ : String = "X",
  mass~ : Double = 0.0,
  charge~ : Double = 0.0,
  sigma~ : Double = 0.0,
  epsilon~ : Double = 0.0,
) -> OmmModelAtom {
  { name, resname, resid, element, mass, charge, sigma, epsilon, }
}

///|
/// Build a TIP3P water atom ("O", "H1" or "H2").
fn omm_water_atom(role : String, resid : Int32) -> OmmModelAtom {
  if role == "O" {
    OmmModelAtom::{
      name: "O",
      resname: "HOH",
      resid,
      element: "O",
      mass: OMM_MASS_O,
      charge: OMM_TIP3P_Q_O,
      sigma: OMM_TIP3P_SIGMA_O,
      epsilon: OMM_TIP3P_EPS_O,
    }
  } else {
    OmmModelAtom::{
      name: role,
      resname: "HOH",
      resid,
      element: "H",
      mass: OMM_MASS_H,
      charge: OMM_TIP3P_Q_H,
      sigma: 0.0,
      epsilon: 0.0,
    }
  }
}

///|
/// Modeller: atom metadata + positions + periodic box.
pub struct OmmModeller {
  mut atoms : Array[OmmModelAtom]
  mut positions : Array[Array[Double]]
  mut box : Array[Double]
  mut next_resid : Int32
}

///|
pub fn OmmModeller::new() -> OmmModeller {
  { atoms: [], positions: [], box: [0.0, 0.0, 0.0], next_resid: 1, }
}

///|
/// Number of atoms currently in the model.
pub fn OmmModeller::num_atoms(self : OmmModeller) -> Int32 {
  self.atoms.length()
}

///|
/// Atom metadata (read-only view).
pub fn OmmModeller::atom_list(self : OmmModeller) -> Array[OmmModelAtom] {
  self.atoms
}

///|
/// Current positions.
pub fn OmmModeller::positions(self : OmmModeller) -> Array[Array[Double]] {
  self.positions
}

///|
/// Periodic box lengths [lx, ly, lz].
pub fn OmmModeller::box_lengths(self : OmmModeller) -> Array[Double] {
  self.box
}

///|
/// Append an atom; returns its index.
pub fn OmmModeller::add_atom(
  self : OmmModeller,
  atom~ : OmmModelAtom,
  position~ : Array[Double],
) -> Int32 {
  let idx = self.atoms.length()
  self.atoms.push(atom)
  self.positions.push([position[0], position[1], position[2]])
  if atom.resid >= self.next_resid {
    self.next_resid = atom.resid + 1
  }
  idx
}

///|
/// Set the periodic box lengths.
pub fn OmmModeller::set_box(
  self : OmmModeller,
  lx : Double,
  ly : Double,
  lz : Double,
) -> Unit {
  self.box = [lx, ly, lz]
}

///|
/// Compute the net charge of the model.
pub fn OmmModeller::net_charge(self : OmmModeller) -> Double {
  let mut q = 0.0
  for a in self.atoms {
    q = q + a.charge
  }
  q
}

///|
/// True when the atom metadata marks a solvent water (HOH/WAT).
fn omm_is_water_res(resname : String) -> Bool {
  resname == "HOH" || resname == "WAT"
}

///|
/// OpenMM Modeller::centerInBox — translate all atoms so that the geometric
/// center of the model sits at the center of the periodic box.
/// No-op when the box is not set.
pub fn OmmModeller::center_in_box(self : OmmModeller) -> Unit {
  if self.box[0] <= 0.0 || self.box[1] <= 0.0 || self.box[2] <= 0.0 {
    return
  }
  let n = self.positions.length()
  if n == 0 {
    return
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for p in self.positions {
    cx = cx + p[0]
    cy = cy + p[1]
    cz = cz + p[2]
  }
  let dx = 0.5 * self.box[0] - cx / n.to_double()
  let dy = 0.5 * self.box[1] - cy / n.to_double()
  let dz = 0.5 * self.box[2] - cz / n.to_double()
  for p in self.positions {
    p[0] = p[0] + dx
    p[1] = p[1] + dy
    p[2] = p[2] + dz
  }
}

///|
/// OpenMM Modeller::deleteWater — remove all water atoms.
/// Returns the number of removed atoms.
pub fn OmmModeller::delete_water(self : OmmModeller) -> Int32 {
  let kept_atoms : Array[OmmModelAtom] = Array(capacity=0)
  let kept_pos : Array[Array[Double]] = Array(capacity=0)
  let mut removed = 0
  for i = 0; i < self.atoms.length(); i = i + 1 {
    if omm_is_water_res(self.atoms[i].resname) {
      removed = removed + 1
    } else {
      kept_atoms.push(self.atoms[i])
      kept_pos.push(self.positions[i])
    }
  }
  if removed > 0 {
    self.atoms = kept_atoms
    self.positions = kept_pos
  }
  removed
}

///|
/// OpenMM Modeller::addSolvent — fill the periodic box with TIP3P water.
///
/// If the box is not set, it is derived from the solute bounding box plus
/// `padding` on each side and the solute is centered. Water oxygens are
/// placed on a cubic lattice with spacing `spacing`, skipping lattice sites
/// closer than `min_dist` (minimum-image) to any existing atom.
///
/// Returns the number of water molecules added.
pub fn OmmModeller::add_solvent(
  self : OmmModeller,
  spacing~ : Double = 3.1,
  min_dist~ : Double = 2.4,
  padding~ : Double = 1.0,
) -> Int32 {
  let n_solute = self.positions.length()
  if n_solute == 0 {
    return 0
  }
  // Derive box from solute bounding box when not set.
  if self.box[0] <= 0.0 || self.box[1] <= 0.0 || self.box[2] <= 0.0 {
    let lo = [0.0, 0.0, 0.0]
    let hi = [0.0, 0.0, 0.0]
    for i = 0; i < n_solute; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        if i == 0 || self.positions[i][d] < lo[d] {
          lo[d] = self.positions[i][d]
        }
        if i == 0 || self.positions[i][d] > hi[d] {
          hi[d] = self.positions[i][d]
        }
      }
    }
    let lx = hi[0] - lo[0] + 2.0 * padding
    let ly = hi[1] - lo[1] + 2.0 * padding
    let lz = hi[2] - lo[2] + 2.0 * padding
    self.set_box(lx, ly, lz)
    self.center_in_box()
  }
  let lx = self.box[0]
  let ly = self.box[1]
  let lz = self.box[2]
  // TIP3P hydrogen offsets (H in the y-z plane).
  let hy = OMM_TIP3P_R_OH *
    @math.sin(OMM_TIP3P_HALF_ANGLE * 3.141592653589793 / 180.0)
  let hz = OMM_TIP3P_R_OH *
    @math.cos(OMM_TIP3P_HALF_ANGLE * 3.141592653589793 / 180.0)
  let mut n_water = 0
  let mut x = 0.5 * spacing
  while x < lx - 0.25 * spacing {
    let mut y = 0.5 * spacing
    while y < ly - 0.25 * spacing {
      let mut z = 0.5 * spacing
      while z < lz - 0.25 * spacing {
        let po = [x, y, z]
        // Reject lattice sites too close to any existing atom (min image).
        let mut ok = true
        for i = 0; i < self.positions.length(); i = i + 1 {
          let d = md_pbc_distance(po, self.positions[i], self.box)
          if d < min_dist {
            ok = false
            break
          }
        }
        if ok {
          let resid = self.next_resid
          self.next_resid = self.next_resid + 1
          let io = self.add_atom(atom=omm_water_atom("O", resid), position=po)
          ignore(io)
          ignore(
            self.add_atom(atom=omm_water_atom("H1", resid), position=[
              x,
              y + hy,
              z + hz,
            ]),
          )
          ignore(
            self.add_atom(atom=omm_water_atom("H2", resid), position=[
              x,
              y - hy,
              z + hz,
            ]),
          )
          n_water = n_water + 1
        }
        z = z + spacing
      }
      y = y + spacing
    }
    x = x + spacing
  }
  n_water
}

///|
/// OpenMM Modeller::addSolvent neutralization — replace water molecules with
/// monoatomic counter-ions (NA+ / CL-) until the net charge is neutralized.
/// Returns the number of ions added.
pub fn OmmModeller::neutralize(
  self : OmmModeller,
  positive~ : String = "NA",
  negative~ : String = "CL",
) -> Int32 {
  let net = self.net_charge()
  if net.abs() < 1.0e-9 {
    return 0
  }
  let want_positive = net < 0.0
  let n_ions = net.abs().ceil().to_int()
  // Collect water oxygen indices (replacement candidates).
  let oxygens : Array[Int32] = Array(capacity=0)
  for i = 0; i < self.atoms.length(); i = i + 1 {
    if omm_is_water_res(self.atoms[i].resname) && self.atoms[i].name == "O" {
      oxygens.push(i)
    }
  }
  let n_water = oxygens.length()
  if n_water == 0 {
    return 0
  }
  let n_add = if n_ions > n_water { n_water } else { n_ions }
  for k = 0; k < n_add; k = k + 1 {
    // Spread ion choices evenly over the water list.
    let oi = oxygens[k * n_water / n_add]
    let o = self.atoms[oi]
    let ion_name = if want_positive { positive } else { negative }
    let ion_mass = if want_positive { 22.9898 } else { 35.453 }
    let ion_sigma = if want_positive { 3.3284 } else { 4.40104 }
    let ion_eps = if want_positive { 0.0031 } else { 0.1 }
    let ion_charge = if want_positive { 1.0 } else { -1.0 }
    self.atoms[oi] = {
      name: ion_name,
      resname: ion_name,
      resid: o.resid,
      element: if want_positive {
        "Na"
      } else {
        "Cl"
      },
      mass: ion_mass,
      charge: ion_charge,
      sigma: ion_sigma,
      epsilon: ion_eps,
    }
    // Remove the two hydrogens of this water (search by resid, descending).
    let mut h_removed = 0
    let mut i = self.atoms.length() - 1
    while i >= 0 {
      if h_removed < 2 &&
        self.atoms[i].resname == "HOH" &&
        self.atoms[i].resid == o.resid {
        ignore(self.atoms.remove(i))
        ignore(self.positions.remove(i))
        h_removed = h_removed + 1
      }
      i = i - 1
    }
  }
  n_add
}

///|
/// Standard valence (number of bonds) for common organic elements.
fn omm_valence(element : String) -> Int32 {
  match element {
    "C" => 4
    "N" => 3
    "O" => 2
    "S" => 2
    "P" => 3
    _ => 0
  }
}

///|
/// Typical X-H bond length (A) for hydrogen placement.
fn omm_bond_length_h(element : String) -> Double {
  match element {
    "C" => 1.09
    "N" => 1.01
    "O" => 0.96
    "S" => 1.34
    "P" => 1.42
    _ => 1.1
  }
}

///|
/// Normalize a 3-vector; returns (1, 0, 0) for the null vector.
fn omm_unit(v : Array[Double]) -> Array[Double] {
  let m = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
  if m < 1.0e-12 {
    return [1.0, 0.0, 0.0]
  }
  [v[0] / m, v[1] / m, v[2] / m]
}

///|
/// Cross product of two 3-vectors.
fn omm_cross(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ]
}

///|
/// OpenMM Modeller::addHydrogens — add hydrogens to heavy atoms according to
/// standard valence rules (C:4, N:3, O:2, S:2, P:3). Bonds are detected
/// geometrically (heavy-heavy < 2.0 A, X-H < 1.6 A). pH-dependent
/// protonation states are not modeled (ph is accepted for API parity).
/// Returns the number of hydrogens added.
pub fn OmmModeller::add_hydrogens(
  self : OmmModeller,
  ph~ : Double = 7.0,
) -> Int32 {
  ignore(ph)
  let n = self.atoms.length()
  if n == 0 {
    return 0
  }
  let is_h : Array[Bool] = Array::make(n, false)
  for i = 0; i < n; i = i + 1 {
    is_h[i] = self.atoms[i].element == "H" || self.atoms[i].element == "D"
  }
  let bonded : Array[Array[Array[Double]]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    bonded[i] = []
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let r = md_dist(self.positions[i], self.positions[j])
      if r < 1.0e-9 {
        continue
      }
      let pair_heavy = !is_h[i] && !is_h[j] && r < 2.0
      let pair_h = is_h[i] != is_h[j] && r < 1.6
      if pair_heavy || pair_h {
        let u = omm_unit([
          self.positions[j][0] - self.positions[i][0],
          self.positions[j][1] - self.positions[i][1],
          self.positions[j][2] - self.positions[i][2],
        ])
        if !is_h[i] {
          bonded[i].push([u[0], u[1], u[2]])
        }
        if !is_h[j] {
          bonded[j].push([-u[0], -u[1], -u[2]])
        }
      }
    }
  }
  let mut added = 0
  for i = 0; i < n; i = i + 1 {
    if is_h[i] {
      continue
    }
    let val = omm_valence(self.atoms[i].element)
    if val == 0 {
      continue
    }
    let missing = val - bonded[i].length()
    if missing <= 0 {
      continue
    }
    let parent = self.atoms[i]
    let p = self.positions[i]
    // Reference direction: opposite the sum of existing bond vectors.
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    for u in bonded[i] {
      sx = sx + u[0]
      sy = sy + u[1]
      sz = sz + u[2]
    }
    let axis = omm_unit([-sx, -sy, -sz])
    // Orthogonal basis around the axis.
    let ref_up = if axis[0].abs() < 0.9 {
      [1.0, 0.0, 0.0]
    } else {
      [0.0, 1.0, 0.0]
    }
    let u_vec = omm_unit(omm_cross(axis, ref_up))
    let w_vec = omm_cross(axis, u_vec)
    let theta = 60.0 * 3.141592653589793 / 180.0
    let r_h = omm_bond_length_h(parent.element)
    for k = 0; k < missing; k = k + 1 {
      let phi = 2.0 * 3.141592653589793 * k.to_double() / missing.to_double()
      let ct = @math.cos(theta)
      let st = @math.sin(theta)
      let cp = @math.cos(phi)
      let sp = @math.sin(phi)
      let dir = [
        ct * axis[0] + st * (cp * u_vec[0] + sp * w_vec[0]),
        ct * axis[1] + st * (cp * u_vec[1] + sp * w_vec[1]),
        ct * axis[2] + st * (cp * u_vec[2] + sp * w_vec[2]),
      ]
      let hname = if missing == 1 { "H" } else { "H" + (k + 1).to_string() }
      ignore(
        self.add_atom(
          atom=OmmModelAtom::{
            name: hname,
            resname: parent.resname,
            resid: parent.resid,
            element: "H",
            mass: OMM_MASS_H,
            charge: 0.0,
            sigma: 0.0,
            epsilon: 0.0,
          },
          position=[
            p[0] + r_h * dir[0],
            p[1] + r_h * dir[1],
            p[2] + r_h * dir[2],
          ],
        ),
      )
      added = added + 1
    }
  }
  added
}
