///|
/// ResidueDepth module - Calculate residue depth in protein structures.
/// Inspired by Biopython Bio.PDB.ResidueDepth.
///
/// Residue depth is a measure of how deeply buried a residue is within
/// a protein structure. It is calculated as the distance from the residue's
/// center of mass to the nearest solvent-accessible surface point.

///|
/// Represents 3D coordinates for ResidueDepth.
pub struct RDPoint3D {
  x : Double
  y : Double
  z : Double
}

///|
/// Create a new RDPoint3D.
pub fn RDPoint3D::new(x : Double, y : Double, z : Double) -> RDPoint3D {
  RDPoint3D::{ x, y, z }
}

///|
/// Calculate distance between two points.
pub fn RDPoint3D::distance(self : RDPoint3D, other : RDPoint3D) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// Vector subtraction.
pub fn RDPoint3D::subtract(self : RDPoint3D, other : RDPoint3D) -> RDPoint3D {
  RDPoint3D::new(self.x - other.x, self.y - other.y, self.z - other.z)
}

///|
/// Vector addition.
pub fn RDPoint3D::add(self : RDPoint3D, other : RDPoint3D) -> RDPoint3D {
  RDPoint3D::new(self.x + other.x, self.y + other.y, self.z + other.z)
}

///|
/// Scalar multiplication.
pub fn RDPoint3D::multiply(self : RDPoint3D, scalar : Double) -> RDPoint3D {
  RDPoint3D::new(self.x * scalar, self.y * scalar, self.z * scalar)
}

///|
/// Calculate vector length.
pub fn RDPoint3D::length(self : RDPoint3D) -> Double {
  (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
}

///|
/// Normalize vector.
pub fn RDPoint3D::normalize(self : RDPoint3D) -> RDPoint3D {
  let len = self.length()
  if len < 0.000001 {
    RDPoint3D::new(0.0, 0.0, 0.0)
  } else {
    RDPoint3D::new(self.x / len, self.y / len, self.z / len)
  }
}

///|
/// An atom for ResidueDepth.
pub struct RDAtom {
  name : String
  coord : RDPoint3D
  element : String
}

///|
/// Create a new RDAtom.
pub fn RDAtom::new(name : String, coord : RDPoint3D, element : String) -> RDAtom {
  RDAtom::{ name, coord, element }
}

///|
/// A residue for ResidueDepth.
pub struct RDResidue {
  resname : String
  resseq : Int
  atoms : Array[RDAtom]
}

///|
/// Create a new RDResidue.
pub fn RDResidue::new(resname : String, resseq : Int, atoms : Array[RDAtom]) -> RDResidue {
  RDResidue::{ resname, resseq, atoms }
}

///|
/// Get the center of mass of a residue.
pub fn RDResidue::center_of_mass(self : RDResidue) -> RDPoint3D {
  if self.atoms.length() == 0 {
    return RDPoint3D::new(0.0, 0.0, 0.0)
  }
  
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_z = 0.0
  
  let mut i = 0
  while i < self.atoms.length() {
    let atom = self.atoms[i]
    sum_x = sum_x + atom.coord.x
    sum_y = sum_y + atom.coord.y
    sum_z = sum_z + atom.coord.z
    i = i + 1
  }
  
  let n = self.atoms.length().to_double()
  RDPoint3D::new(sum_x / n, sum_y / n, sum_z / n)
}

///|
/// Get the CA atom if it exists.
pub fn RDResidue::get_ca(self : RDResidue) -> RDAtom? {
  let mut i = 0
  while i < self.atoms.length() {
    let atom = self.atoms[i]
    if atom.name == "CA" {
      return Some(atom)
    }
    i = i + 1
  }
  None
}

///|
/// ResidueDepth result struct.
pub struct RDResult {
  resseq : Int
  resname : String
  ca_depth : Double
  com_depth : Double
  ca_sasa : Double
  com_sasa : Double
}

///|
/// Create a new RDResult.
pub fn RDResult::new(
  resseq : Int,
  resname : String,
  ca_depth : Double,
  com_depth : Double,
  ca_sasa : Double,
  com_sasa : Double
) -> RDResult {
  RDResult::{ resseq, resname, ca_depth, com_depth, ca_sasa, com_sasa }
}

///|
/// Calculate the solvent-accessible surface area (SASA) using the rolling ball algorithm.
pub fn calculate_sasa(atom : RDAtom, all_atoms : Array[RDAtom], probe_radius : Double) -> Double {
  let atom_radius = get_atom_radius(atom.element)
  let total_radius = atom_radius + probe_radius
  
  let mut sasa = 0.0
  let n_points = 92
  
  let mut i = 0
  while i < n_points {
    let theta = (i.to_double() / n_points.to_double()) * 2.0 * 3.141592653589793
    let mut j = 0
    while j < n_points {
      let phi = (j.to_double() / n_points.to_double()) * 3.141592653589793
      
      let x = total_radius * @math.cos(theta) * @math.sin(phi)
      let y = total_radius * @math.sin(theta) * @math.sin(phi)
      let z = total_radius * @math.cos(phi)
      
      let test_point = RDPoint3D::new(atom.coord.x + x, atom.coord.y + y, atom.coord.z + z)
      
      let mut is_accessible = true
      let mut k = 0
      while k < all_atoms.length() {
        let other_atom = all_atoms[k]
        if other_atom.name == atom.name && other_atom.coord.x == atom.coord.x && other_atom.coord.y == atom.coord.y && other_atom.coord.z == atom.coord.z {
          k = k + 1
          continue
        }
        
        let other_radius = get_atom_radius(other_atom.element) + probe_radius
        let dist = atom.coord.distance(other_atom.coord)
        
        if dist < other_radius + atom_radius {
          is_accessible = false
          break
        }
        k = k + 1
      }
      
      if is_accessible {
        let area_element = 4.0 * 3.141592653589793 * total_radius * total_radius / (n_points.to_double() * n_points.to_double())
        sasa = sasa + area_element
      }
      j = j + 1
    }
    i = i + 1
  }
  
  sasa
}

///|
/// Get the van der Waals radius for an atom element.
pub fn get_atom_radius(element : String) -> Double {
  if element == "H" || element == "h" { 1.20 }
  else if element == "C" || element == "c" { 1.70 }
  else if element == "N" || element == "n" { 1.55 }
  else if element == "O" || element == "o" { 1.52 }
  else if element == "S" || element == "s" { 1.80 }
  else if element == "P" || element == "p" { 1.80 }
  else { 1.50 }
}

///|
/// Calculate depth for a point using ray casting.
pub fn calculate_point_depth(point : RDPoint3D, all_atoms : Array[RDAtom], probe_radius : Double) -> Double {
  let mut min_distance = 10000.0
  
  let mut i = 0
  while i < all_atoms.length() {
    let atom = all_atoms[i]
    let atom_radius = get_atom_radius(atom.element) + probe_radius
    let dist = point.distance(atom.coord)
    
    let depth = dist - atom_radius
    if depth < min_distance {
      min_distance = depth
    }
    i = i + 1
  }
  
  if min_distance < 0.0 {
    0.0
  } else {
    min_distance
  }
}

///|
/// Calculate depth using grid-based approach.
pub fn calculate_residue_depth(
  residue : RDResidue,
  all_residues : Array[RDResidue],
  probe_radius : Double,
  grid_spacing : Double
) -> RDResult {
  let all_atoms = collect_all_atoms(all_residues)
  
  let com = residue.center_of_mass()
  
  let ca_coord = match residue.get_ca() {
    Some(a) => a.coord
    None => com
  }
  
  let com_depth = calculate_point_depth(com, all_atoms, probe_radius)
  let ca_depth = calculate_point_depth(ca_coord, all_atoms, probe_radius)
  
  let ca_atom = match residue.get_ca() {
    Some(a) => a
    None => RDAtom::new("CA", ca_coord, "C")
  }
  
  let com_atom = RDAtom::new("COM", com, "C")
  
  let ca_sasa = calculate_sasa(ca_atom, all_atoms, probe_radius)
  let com_sasa = calculate_sasa(com_atom, all_atoms, probe_radius)
  
  RDResult::new(residue.resseq, residue.resname, ca_depth, com_depth, ca_sasa, com_sasa)
}

///|
/// Calculate depth using grid-based approach (with defaults).
pub fn calculate_residue_depth_default(
  residue : RDResidue,
  all_residues : Array[RDResidue]
) -> RDResult {
  calculate_residue_depth(residue, all_residues, 1.4, 0.5)
}

///|
/// Collect all atoms from residues.
pub fn collect_all_atoms(residues : Array[RDResidue]) -> Array[RDAtom] {
  let atoms : Array[RDAtom] = Array::new()
  
  let mut i = 0
  while i < residues.length() {
    let residue = residues[i]
    let mut j = 0
    while j < residue.atoms.length() {
      atoms.push(residue.atoms[j])
      j = j + 1
    }
    i = i + 1
  }
  
  atoms
}

///|
/// Analyze all residues in a structure.
pub fn analyze_residue_depth(
  residues : Array[RDResidue],
  probe_radius : Double,
  grid_spacing : Double
) -> Array[RDResult] {
  let results : Array[RDResult] = Array::new()
  
  let mut i = 0
  while i < residues.length() {
    let result = calculate_residue_depth(residues[i], residues, probe_radius, grid_spacing)
    results.push(result)
    i = i + 1
  }
  
  results
}

///|
/// Analyze all residues in a structure (with defaults).
pub fn analyze_residue_depth_default(residues : Array[RDResidue]) -> Array[RDResult] {
  analyze_residue_depth(residues, 1.4, 0.5)
}

///|
/// Calculate average depth for a structure.
pub fn calculate_average_depth(results : Array[RDResult]) -> (Double, Double) {
  if results.length() == 0 {
    return (0.0, 0.0)
  }
  
  let mut sum_ca = 0.0
  let mut sum_com = 0.0
  
  let mut i = 0
  while i < results.length() {
    let result = results[i]
    sum_ca = sum_ca + result.ca_depth
    sum_com = sum_com + result.com_depth
    i = i + 1
  }
  
  let n = results.length().to_double()
  (sum_ca / n, sum_com / n)
}

///|
/// Find surface residues (depth below threshold).
pub fn find_surface_residues(results : Array[RDResult], depth_threshold : Double) -> Array[RDResult] {
  let surface : Array[RDResult] = Array::new()
  
  let mut i = 0
  while i < results.length() {
    let result = results[i]
    if result.ca_depth < depth_threshold {
      surface.push(result)
    }
    i = i + 1
  }
  
  surface
}

///|
/// Find surface residues (with default threshold).
pub fn find_surface_residues_default(results : Array[RDResult]) -> Array[RDResult] {
  find_surface_residues(results, 3.0)
}

///|
/// Find core residues (depth above threshold).
pub fn find_core_residues(results : Array[RDResult], depth_threshold : Double) -> Array[RDResult] {
  let core : Array[RDResult] = Array::new()
  
  let mut i = 0
  while i < results.length() {
    let result = results[i]
    if result.ca_depth > depth_threshold {
      core.push(result)
    }
    i = i + 1
  }
  
  core
}

///|
/// Find core residues (with default threshold).
pub fn find_core_residues_default(results : Array[RDResult]) -> Array[RDResult] {
  find_core_residues(results, 8.0)
}

///|
/// Bio API functions.

pub fn bio_residue_depth_analyze(residues : Array[RDResidue]) -> Array[RDResult] {
  analyze_residue_depth_default(residues)
}

pub fn bio_residue_depth_calculate(residue : RDResidue, all_residues : Array[RDResidue]) -> RDResult {
  calculate_residue_depth_default(residue, all_residues)
}

pub fn bio_residue_depth_find_surface(residues : Array[RDResidue], threshold : Double) -> Array[RDResult] {
  let results = analyze_residue_depth_default(residues)
  find_surface_residues(results, threshold)
}

pub fn bio_residue_depth_find_core(residues : Array[RDResidue], threshold : Double) -> Array[RDResult] {
  let results = analyze_residue_depth_default(residues)
  find_core_residues(results, threshold)
}

pub fn bio_residue_depth_average(results : Array[RDResult]) -> (Double, Double) {
  calculate_average_depth(results)
}