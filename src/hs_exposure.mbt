///|
/// Bio.PDB.HSExposure - Hemisphere exposure calculation.

///|
fn hse_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = if x > 1.0 { x } else { 1.0 }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = if next > guess { next - guess } else { guess - next }
    guess = next
    if diff < 0.0000000000001 {
      return guess
    }
    i = i + 1
  }
  guess
}

///|
pub struct PDBAtom {
  name : String
  x : Double
  y : Double
  z : Double
}

///|
pub fn PDBAtom::new(name : String, x : Double, y : Double, z : Double) -> PDBAtom {
  PDBAtom::{ name, x, y, z }
}

///|
pub struct HSExposureResult {
  residue_name : String
  hse_up : Double
  hse_down : Double
  hse_total : Double
}

///|
pub fn HSExposureResult::new(residue_name : String, hse_up : Double, hse_down : Double, hse_total : Double) -> HSExposureResult {
  HSExposureResult::{ residue_name, hse_up, hse_down, hse_total }
}

///|
fn atom_distance(a1 : PDBAtom, a2 : PDBAtom) -> Double {
  let dx = a1.x - a2.x
  let dy = a1.y - a2.y
  let dz = a1.z - a2.z
  
  hse_sqrt(dx * dx + dy * dy + dz * dz)
}

///|
pub fn calculate_hse(ca_atom : PDBAtom, cb_atom : PDBAtom, n_atom : PDBAtom, all_ca_atoms : Array[PDBAtom]) -> HSExposureResult {
  let mut hse_up = 0
  let mut hse_down = 0
  
  let mut i = 0
  while i < all_ca_atoms.length() {
    let other_ca = all_ca_atoms[i]
    
    if other_ca.name == ca_atom.name {
      i = i + 1
      continue
    }
    
    let dist = atom_distance(ca_atom, other_ca)
    
    if dist <= 10.0 && dist > 0.0 {
      let dx = other_ca.x - ca_atom.x
      let dy = other_ca.y - ca_atom.y
      let dz = other_ca.z - ca_atom.z
      
      let cb_dx = cb_atom.x - ca_atom.x
      let cb_dy = cb_atom.y - ca_atom.y
      let cb_dz = cb_atom.z - ca_atom.z
      
      let n_dx = n_atom.x - ca_atom.x
      let n_dy = n_atom.y - ca_atom.y
      let n_dz = n_atom.z - ca_atom.z
      
      let nx = cb_dy * n_dz - cb_dz * n_dy
      let ny = cb_dz * n_dx - cb_dx * n_dz
      let nz = cb_dx * n_dy - cb_dy * n_dx
      
      let dot = dx * nx + dy * ny + dz * nz
      
      if dot > 0.0 {
        hse_up = hse_up + 1
      } else if dot < 0.0 {
        hse_down = hse_down + 1
      }
    }
    
    i = i + 1
  }
  
  HSExposureResult::new(ca_atom.name, hse_up.to_double(), hse_down.to_double(), (hse_up + hse_down).to_double())
}

///|
pub fn create_example_atoms() -> Array[PDBAtom] {
  let atoms : Array[PDBAtom] = Array::new()
  
  atoms.push(PDBAtom::new("CA1", 0.0, 0.0, 0.0))
  atoms.push(PDBAtom::new("CB1", 1.5, 0.0, 0.0))
  atoms.push(PDBAtom::new("N1", 0.0, 1.5, 0.0))
  
  atoms.push(PDBAtom::new("CA2", 3.0, 0.0, 0.0))
  atoms.push(PDBAtom::new("CB2", 4.5, 0.0, 0.0))
  atoms.push(PDBAtom::new("N2", 3.0, 1.5, 0.0))
  
  atoms.push(PDBAtom::new("CA3", 0.0, 3.0, 0.0))
  atoms.push(PDBAtom::new("CB3", 1.5, 3.0, 0.0))
  atoms.push(PDBAtom::new("N3", 0.0, 4.5, 0.0))
  
  atoms
}

///|
pub fn classify_exposure(hse_total : Double) -> String {
  if hse_total < 5.0 {
    "buried"
  } else if hse_total < 15.0 {
    "partially_exposed"
  } else {
    "exposed"
  }
}