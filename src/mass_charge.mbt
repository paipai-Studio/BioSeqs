///|
/// Mass properties and charge analysis for molecular structures.
///
/// Biotite provides mass center, moments of inertia, and charge calculations.

///|
/// 10^x helper.
fn mc_pow10(x : Double) -> Double {
  let ln10 = 2.302585092994046
  @math.exp(x * ln10)
}

///|
/// Amino acid monoisotopic masses (Da).
fn aa_mass(resname : String) -> Double {
  match resname {
    "ALA" => 71.03711
    "ARG" => 156.10111
    "ASN" => 114.04293
    "ASP" => 115.02694
    "CYS" => 103.00919
    "GLN" => 128.05858
    "GLU" => 129.04259
    "GLY" => 57.02146
    "HIS" => 137.05891
    "ILE" => 113.08406
    "LEU" => 113.08406
    "LYS" => 128.09496
    "MET" => 131.04049
    "PHE" => 147.06841
    "PRO" => 97.05276
    "SER" => 87.03203
    "THR" => 101.04768
    "TRP" => 186.07931
    "TYR" => 163.06333
    "VAL" => 99.06841
    _ => 110.0
  }
}

///|
/// Atom element masses (Da).
fn element_mass(element : String) -> Double {
  match element {
    "H" => 1.00794
    "C" => 12.0107
    "N" => 14.0067
    "O" => 15.9994
    "S" => 32.065
    "P" => 30.9738
    "FE" => 55.845
    "ZN" => 65.38
    "CU" => 63.546
    "MG" => 24.305
    "CA" => 40.078
    "NA" => 22.9898
    "K" => 39.0983
    "CL" => 35.453
    "F" => 18.9984
    "BR" => 79.904
    "I" => 126.904
    "SE" => 78.96
    _ => 12.0
  }
}

///|
/// Mass properties of a structure.
pub struct MassProperties {
  total_mass : Double
  center_of_mass : Vector3
  inertia : (Double, Double, Double) // principal moments
  radius_of_gyration : Double
}

///|
/// Get total mass.
pub fn MassProperties::get_total_mass(self : MassProperties) -> Double {
  self.total_mass
}

///|
/// Get center of mass.
pub fn MassProperties::get_center_of_mass(self : MassProperties) -> Vector3 {
  self.center_of_mass
}

///|
/// Get radius of gyration.
pub fn MassProperties::get_radius_of_gyration(self : MassProperties) -> Double {
  self.radius_of_gyration
}

///|
/// Calculate mass properties from atomic coordinates and elements.
pub fn calculate_mass_properties(
  coords : Array[Vector3],
  elements : Array[String],
) -> MassProperties {
  let n = coords.length()
  if n == 0 {
    return MassProperties::{
      total_mass: 0.0,
      center_of_mass: Vector3::new(0.0, 0.0, 0.0),
      inertia: (0.0, 0.0, 0.0),
      radius_of_gyration: 0.0,
    }
  }
  let mut total_mass = 0.0
  let mut wx = 0.0
  let mut wy = 0.0
  let mut wz = 0.0
  for i in 0..<n {
    let m = element_mass(elements[i])
    total_mass = total_mass + m
    wx = wx + m * coords[i].x
    wy = wy + m * coords[i].y
    wz = wz + m * coords[i].z
  }
  let com = if total_mass > 0.0 {
    Vector3::new(wx / total_mass, wy / total_mass, wz / total_mass)
  } else {
    Vector3::new(0.0, 0.0, 0.0)
  }
  let mut ixx = 0.0
  let mut iyy = 0.0
  let mut izz = 0.0
  let mut rg_sum = 0.0
  for i in 0..<n {
    let m = element_mass(elements[i])
    let dx = coords[i].x - com.x
    let dy = coords[i].y - com.y
    let dz = coords[i].z - com.z
    ixx = ixx + m * (dy * dy + dz * dz)
    iyy = iyy + m * (dx * dx + dz * dz)
    izz = izz + m * (dx * dx + dy * dy)
    rg_sum = rg_sum + m * (dx * dx + dy * dy + dz * dz)
  }
  let rg = if total_mass > 0.0 { (rg_sum / total_mass).sqrt() } else { 0.0 }
  MassProperties::{
    total_mass,
    center_of_mass: com,
    inertia: (ixx, iyy, izz),
    radius_of_gyration: rg,
  }
}

///|
/// Calculate protein molecular weight from residue sequence.
pub fn calculate_molecular_weight(residue_names : Array[String]) -> Double {
  let water_mass = 18.01056
  let mut total = water_mass
  for r in residue_names {
    total = total + aa_mass(r)
  }
  total
}

///|
/// Residue side chain charge at given pH.
fn residue_sidechain_charge(resname : String, ph : Double) -> Double {
  let pka_arg = 12.48
  let pka_lys = 10.53
  let pka_his = 6.00
  let pka_asp = 3.65
  let pka_glu = 4.25
  let pka_cys = 8.18
  let pka_tyr = 10.07
  match resname {
    "ARG" => 1.0 / (1.0 + mc_pow10(ph - pka_arg))
    "LYS" => 1.0 / (1.0 + mc_pow10(ph - pka_lys))
    "HIS" => 1.0 / (1.0 + mc_pow10(ph - pka_his))
    "ASP" => -1.0 / (1.0 + mc_pow10(pka_asp - ph))
    "GLU" => -1.0 / (1.0 + mc_pow10(pka_glu - ph))
    "CYS" => -1.0 / (1.0 + mc_pow10(pka_cys - ph))
    "TYR" => -1.0 / (1.0 + mc_pow10(pka_tyr - ph))
    _ => 0.0
  }
}

///|
/// Calculate net charge at a given pH.
pub fn calculate_net_charge(
  residue_names : Array[String],
  ph : Double,
) -> Double {
  let pka_n_term = 9.0
  let pka_c_term = 2.0
  let mut charge = 0.0
  charge = charge + 1.0 / (1.0 + mc_pow10(ph - pka_n_term))
  charge = charge + -1.0 / (1.0 + mc_pow10(pka_c_term - ph))
  for r in residue_names {
    charge = charge + residue_sidechain_charge(r, ph)
  }
  charge
}

///|
/// Calculate isoelectric point (pI) using binary search.
pub fn calculate_pI(residue_names : Array[String]) -> Double {
  let mut low = 0.0
  let mut high = 14.0
  for _ in 0..<100 {
    let mid = (low + high) / 2.0
    let charge = calculate_net_charge(residue_names, mid)
    if charge > 0.0 {
      low = mid
    } else {
      high = mid
    }
  }
  (low + high) / 2.0
}
