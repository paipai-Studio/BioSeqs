/// omm_water_models.mbt — SPC/SPC-E/TIP3P-FB/TIP4P-Ew/TIP5P parameter library
/// Item 20: Water model parameters and assignment

///|
/// Water model parameters
pub struct OmmWaterModel {
  name : String
  // Geometry
  dOH : Double // O-H bond length (Å)
  dHH : Double // H-H distance (Å)
  angle_HOH : Double // H-O-H angle (degrees)
  // Charges
  qO : Double
  qH : Double
  // Lennard-Jones (on oxygen)
  sigma_O : Double // Å
  epsilon_O : Double // kcal/mol
  // Lennard-Jones (on hydrogen, usually 0)
  sigma_H : Double
  epsilon_H : Double
  // Masses
  mass_O : Double
  mass_H : Double
  // Virtual site (for TIP4P/TIP5P)
  has_virtual_site : Bool
  qV : Double // virtual site charge
  dOV : Double // O-V distance
  // Model type
  n_sites : Int32 // 3, 4, or 5
}

///|
/// SPC water model (Simple Point Charge)
pub fn omm_water_spc() -> OmmWaterModel {
  OmmWaterModel::{
    name: "SPC",
    dOH: 1.0,
    dHH: 1.63299,
    angle_HOH: 109.47,
    qO: -0.82,
    qH: 0.41,
    sigma_O: 3.166,
    epsilon_O: 0.1553,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 15.9994,
    mass_H: 1.008,
    has_virtual_site: false,
    qV: 0.0,
    dOV: 0.0,
    n_sites: 3,
  }
}

///|
/// SPC/E water model (Simple Point Charge Extended)
pub fn omm_water_spce() -> OmmWaterModel {
  OmmWaterModel::{
    name: "SPC/E",
    dOH: 1.0,
    dHH: 1.63299,
    angle_HOH: 109.47,
    qO: -0.8476,
    qH: 0.4238,
    sigma_O: 3.166,
    epsilon_O: 0.1553,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 15.9994,
    mass_H: 1.008,
    has_virtual_site: false,
    qV: 0.0,
    dOV: 0.0,
    n_sites: 3,
  }
}

///|
/// TIP3P-FB water model (Force Balance optimized)
pub fn omm_water_tip3p_fb() -> OmmWaterModel {
  OmmWaterModel::{
    name: "TIP3P-FB",
    dOH: 0.9572,
    dHH: 1.5139,
    angle_HOH: 104.52,
    qO: -0.410,
    qH: 0.205,
    sigma_O: 3.18195,
    epsilon_O: 0.0901,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 15.9994,
    mass_H: 1.008,
    has_virtual_site: false,
    qV: 0.0,
    dOV: 0.0,
    n_sites: 3,
  }
}

///|
/// TIP4P-Ew water model (Ewald-optimized 4-site)
pub fn omm_water_tip4p_ew() -> OmmWaterModel {
  OmmWaterModel::{
    name: "TIP4P-Ew",
    dOH: 0.9572,
    dHH: 1.5139,
    angle_HOH: 104.52,
    qO: 0.0, // charge moved to virtual site
    qH: 0.52422,
    sigma_O: 3.16435,
    epsilon_O: 0.16275,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 15.9994,
    mass_H: 1.008,
    has_virtual_site: true,
    qV: -1.04844,
    dOV: 0.125, // M-site distance from O
    n_sites: 4,
  }
}

///|
/// TIP5P water model (5-site with lone pairs)
pub fn omm_water_tip5p() -> OmmWaterModel {
  OmmWaterModel::{
    name: "TIP5P",
    dOH: 0.9572,
    dHH: 1.5139,
    angle_HOH: 104.52,
    qO: 0.0,
    qH: 0.241,
    sigma_O: 3.097,
    epsilon_O: 0.0989,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 15.9994,
    mass_H: 1.008,
    has_virtual_site: true,
    qV: -0.241, // each lone pair
    dOV: 0.70, // lone pair distance from O
    n_sites: 5,
  }
}

///|
/// Get water model by name.
pub fn omm_water_model_by_name(name : String) -> OmmWaterModel? {
  let lower = name.to_lower()
  if lower == "spc" {
    Some(omm_water_spc())
  } else if lower == "spc/e" || lower == "spce" {
    Some(omm_water_spce())
  } else if lower == "tip3p" || lower == "tip3p-fb" {
    Some(omm_water_tip3p_fb())
  } else if lower == "tip4p" || lower == "tip4p-ew" || lower == "tip4pew" {
    Some(omm_water_tip4p_ew())
  } else if lower == "tip5p" {
    Some(omm_water_tip5p())
  } else {
    None
  }
}

///|
/// Generate initial water molecule coordinates.
/// O at origin, H1 and H2 in xz-plane.
pub fn omm_water_coords(
  model : OmmWaterModel,
) -> Array[(Double, Double, Double)] {
  let dOH = model.dOH
  let half_angle = model.angle_HOH * 3.141592653589793 / 360.0 // half of H-O-H angle in radians
  let h1x = dOH * @math.sin(half_angle)
  let h1z = dOH * @math.cos(half_angle)
  let h2x = -h1x
  let h2z = h1z
  let coords : Array[(Double, Double, Double)] = Array(capacity=model.n_sites)
  coords.push((0.0, 0.0, 0.0)) // O at origin
  coords.push((h1x, 0.0, h1z)) // H1
  coords.push((h2x, 0.0, h2z)) // H2
  if model.has_virtual_site && model.n_sites >= 4 {
    // Virtual site (M or lone pair) along negative z from O
    coords.push((0.0, 0.0, -model.dOV))
  }
  if model.n_sites >= 5 {
    // Second lone pair
    coords.push((0.0, 0.0, -model.dOV))
  }
  coords
}

///|
/// Get charges for water model.
pub fn omm_water_charges(model : OmmWaterModel) -> Array[Double] {
  let charges : Array[Double] = Array(capacity=model.n_sites)
  charges.push(model.qO)
  charges.push(model.qH)
  charges.push(model.qH)
  if model.has_virtual_site && model.n_sites >= 4 {
    charges.push(model.qV)
  }
  if model.n_sites >= 5 {
    charges.push(model.qV)
  }
  charges
}

///|
/// Get masses for water model.
pub fn omm_water_masses(model : OmmWaterModel) -> Array[Double] {
  let masses : Array[Double] = Array(capacity=model.n_sites)
  masses.push(model.mass_O)
  masses.push(model.mass_H)
  masses.push(model.mass_H)
  if model.has_virtual_site && model.n_sites >= 4 {
    masses.push(0.0) // virtual site has no mass
  }
  if model.n_sites >= 5 {
    masses.push(0.0)
  }
  masses
}

///|
/// List all available water model names.
pub fn omm_water_model_names() -> Array[String] {
  ["SPC", "SPC/E", "TIP3P-FB", "TIP4P-Ew", "TIP5P"]
}

///|
pub fn omm_water_models_dummy() -> Unit {
  ignore(OmmWaterModel::{
    name: "",
    dOH: 0.0,
    dHH: 0.0,
    angle_HOH: 0.0,
    qO: 0.0,
    qH: 0.0,
    sigma_O: 0.0,
    epsilon_O: 0.0,
    sigma_H: 0.0,
    epsilon_H: 0.0,
    mass_O: 0.0,
    mass_H: 0.0,
    has_virtual_site: false,
    qV: 0.0,
    dOV: 0.0,
    n_sites: 0,
  })
}
