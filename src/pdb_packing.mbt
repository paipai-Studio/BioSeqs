///|
/// Bio.PDB.Packing - Protein packing density analysis.
/// Inspired by Biopython's Bio.PDB.Packing module.
///
/// Provides tools for calculating local packing density around protein atoms,
/// per-residue packing analysis, and solvent accessible surface area (SASA)
/// using the Lee & Richards algorithm.

///|
/// An atom for packing density analysis, with coordinates and residue info.
pub struct PackingAtom {
  x : Double
  y : Double
  z : Double
  residue_name : String
  residue_seq : Int
  chain_id : Char
  atom_name : String
  element : String
  vdw_radius : Double
}

///|
/// Construct a PackingAtom.
pub fn PackingAtom::new(
  x~ : Double,
  y~ : Double,
  z~ : Double,
  residue_name? : String = "UNK",
  residue_seq? : Int = 1,
  chain_id? : Char = ' ',
  atom_name? : String = "",
  element? : String = "",
  vdw_radius? : Double = 1.7,
) -> PackingAtom {
  PackingAtom::{
    x,
    y,
    z,
    residue_name,
    residue_seq,
    chain_id,
    atom_name,
    element,
    vdw_radius,
  }
}

///|
/// Get coordinates as an array.
pub fn PackingAtom::get_coords(self : PackingAtom) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Calculate distance to another PackingAtom.
pub fn PackingAtom::distance(self : PackingAtom, other : PackingAtom) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// A point on a unit sphere, used for shell-based calculations.
pub struct SpherePoint {
  x : Double
  y : Double
  z : Double
}

///|
/// Construct a SpherePoint using positional arguments.
pub fn SpherePoint::new(x : Double, y : Double, z : Double) -> SpherePoint {
  SpherePoint::{ x, y, z }
}

///|
/// Get coordinates as array.
pub fn SpherePoint::get_coords(self : SpherePoint) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Result of packing density calculation for a single residue.
pub struct PackingResult {
  residue_name : String
  residue_seq : Int
  chain_id : Char
  density : Double
  n_contacting_atoms : Int
  n_shell_points : Int
}

///|
/// Construct a PackingResult.
pub fn PackingResult::new(
  residue_name~ : String,
  residue_seq~ : Int,
  chain_id? : Char = ' ',
  density? : Double = 0.0,
  n_contacting_atoms? : Int = 0,
  n_shell_points? : Int = 0,
) -> PackingResult {
  PackingResult::{
    residue_name,
    residue_seq,
    chain_id,
    density,
    n_contacting_atoms,
    n_shell_points,
  }
}

///|
/// Get the density value.
pub fn PackingResult::get_density(self : PackingResult) -> Double {
  self.density
}

///|
/// Get the residue identifier string.
pub fn PackingResult::get_residue_id(self : PackingResult) -> String {
  self.chain_id.to_string() + self.residue_seq.to_string() + self.residue_name
}

///|
/// Result of comprehensive packing analysis.
pub struct PackingAnalysisResult {
  results : Array[PackingResult]
  mean_density : Double
  min_density : Double
  max_density : Double
  median_density : Double
  n_residues : Int
  n_buried : Int
  n_exposed : Int
}

///|
/// Construct a PackingAnalysisResult.
pub fn PackingAnalysisResult::new(
  results~ : Array[PackingResult],
  mean_density? : Double = 0.0,
  min_density? : Double = 0.0,
  max_density? : Double = 0.0,
  median_density? : Double = 0.0,
  n_residues? : Int = 0,
  n_buried? : Int = 0,
  n_exposed? : Int = 0,
) -> PackingAnalysisResult {
  PackingAnalysisResult::{
    results,
    mean_density,
    min_density,
    max_density,
    median_density,
    n_residues,
    n_buried,
    n_exposed,
  }
}

///|
/// Get results array.
pub fn PackingAnalysisResult::get_results(self : PackingAnalysisResult) -> Array[PackingResult] {
  self.results
}

///|
/// Get summary statistics as a formatted string.
pub fn PackingAnalysisResult::summary(self : PackingAnalysisResult) -> String {
  let mut s = "Packing Analysis Summary:\n"
  s = s + "  Residues analyzed: " + self.n_residues.to_string() + "\n"
  s = s + "  Mean density: " + self.mean_density.to_string() + "\n"
  s = s + "  Min density: " + self.min_density.to_string() + "\n"
  s = s + "  Max density: " + self.max_density.to_string() + "\n"
  s = s + "  Median density: " + self.median_density.to_string() + "\n"
  s = s + "  Buried residues: " + self.n_buried.to_string() + "\n"
  s = s + "  Exposed residues: " + self.n_exposed.to_string()
  s
}

///|
/// Generate approximately uniform points on a unit sphere
/// using the golden spiral (Fibonacci sphere) algorithm.
/// n_points: number of points to generate (default 100).
pub fn generate_packing_sphere_points(n_points? : Int = 100) -> Array[SpherePoint] {
  let points : Array[SpherePoint] = Array::new()
  let n = if n_points < 1 { 1 } else { n_points }
  let golden_angle = 3.14159265358979323846 * (3.0 - 5.0.sqrt())

  let mut i = 0
  while i < n {
    let y = 1.0 - (2.0 * i.to_double()) / (n.to_double() - 1.0)
    let radius_at_y = (1.0 - y * y).sqrt()
    let theta = golden_angle * i.to_double()

    let x = @math.cos(theta) * radius_at_y
    let z = @math.sin(theta) * radius_at_y

    points.push(SpherePoint::new(x, y, z))
    i = i + 1
  }
  points
}

///|
/// Check if a point is inside a sphere defined by center and radius.
pub fn atom_in_sphere(
  center : Array[Double],
  point : Array[Double],
  radius : Double,
) -> Bool {
  let dx = point[0] - center[0]
  let dy = point[1] - center[1]
  let dz = point[2] - center[2]
  (dx * dx + dy * dy + dz * dz) <= radius * radius
}

///|
/// Count how many atoms from an array are within a sphere
/// centered at a given point with a given radius.
/// Excludes self (atom at given index).
pub fn count_atoms_in_sphere(
  center : Array[Double],
  atoms : Array[PackingAtom],
  radius : Double,
  exclude_index? : Int = -1,
) -> Int {
  let mut count = 0
  let mut i = 0
  while i < atoms.length() {
    if i != exclude_index {
      let atom = atoms[i]
      if atom_in_sphere(center, [atom.x, atom.y, atom.z], radius) {
        count = count + 1
      }
    }
    i = i + 1
  }
  count
}

///|
/// Calculate local packing density for a set of atoms.
/// Uses Lee & Richards algorithm: places sphere points on
/// the atom's van der Waals+probe sphere and counts how many
/// are buried by other atoms.
///
/// radius: probe radius in Angstroms (default 1.4 Å for water).
/// n_points: number of sphere points for shell (default 100).
pub fn calculate_packing_density(
  atoms : Array[PackingAtom],
  radius? : Double = 1.4,
  n_points? : Int = 100,
) -> PackingResult {
  if atoms.length() == 0 {
    return PackingResult::new(residue_name="", residue_seq=0)
  }

  let sphere_points = generate_packing_sphere_points(n_points=n_points)
  let atom = atoms[0]
  let atom_radius = atom.vdw_radius + radius

  let mut n_buried = 0
  let n_total = sphere_points.length()

  for sp in sphere_points {
    let px = atom.x + sp.x * atom_radius
    let py = atom.y + sp.y * atom_radius
    let pz = atom.z + sp.z * atom_radius

    let mut buried = false
    let mut j = 1
    while j < atoms.length() {
      let other = atoms[j]
      let dx = px - other.x
      let dy = py - other.y
      let dz = pz - other.z
      let dist_sq = dx * dx + dy * dy + dz * dz
      let other_r = other.vdw_radius + radius
      if dist_sq < other_r * other_r {
        buried = true
        break
      }
      j = j + 1
    }

    if buried {
      n_buried = n_buried + 1
    }
  }

  let density = if n_total > 0 {
    n_buried.to_double() / n_total.to_double()
  } else {
    0.0
  }

  PackingResult::new(
    residue_name=atom.residue_name,
    residue_seq=atom.residue_seq,
    chain_id=atom.chain_id,
    density~,
    n_contacting_atoms=n_buried,
    n_shell_points=n_total,
  )
}

///|
/// Calculate packing density for each residue in the atom set.
/// Groups atoms by residue and calculates per-residue density.
///
/// radius: shell radius in Angstroms (default 7.0 Å).
/// n_points: number of sphere points (default 100).
pub fn packing_density_per_residue(
  atoms : Array[PackingAtom],
  radius? : Double = 1.4,
  n_points? : Int = 100,
) -> PackingAnalysisResult {
  if atoms.length() == 0 {
    return PackingAnalysisResult::new(results=[])
  }

  let residue_groups : Map[String, Array[PackingAtom]] = Map([], capacity=16)

  for atom in atoms {
    let key = atom.chain_id.to_string() + ":" + atom.residue_seq.to_string() + ":" + atom.residue_name
    match residue_groups.get(key) {
      Some(arr) => {
        arr.push(atom)
      }
      None => {
        let arr : Array[PackingAtom] = Array::new()
        arr.push(atom)
        residue_groups[key] = arr
      }
    }
  }

  let results : Array[PackingResult] = Array::new()
  let all_keys = residue_groups.keys().collect()
  let mut ki = 0
  while ki < all_keys.length() {
    let key = all_keys[ki]
    let res_atoms = residue_groups[key]
    let density_result = calculate_packing_density(res_atoms, radius~, n_points~)
    results.push(density_result)
    ki = ki + 1
  }

  let n = results.length()
  if n == 0 {
    return PackingAnalysisResult::new(results=[])
  }

  let mut sum = 0.0
  let mut min_d = results[0].density
  let mut max_d = results[0].density
  let mut i = 0
  while i < n {
    let d = results[i].density
    sum = sum + d
    if d < min_d {
      min_d = d
    }
    if d > max_d {
      max_d = d
    }
    i = i + 1
  }

  let mean_d = sum / n.to_double()

  let sorted_densities : Array[Double] = Array::new()
  for r in results {
    sorted_densities.push(r.density)
  }
  sorted_densities.sort()

  let median = if n % 2 == 0 {
    (sorted_densities[n / 2 - 1] + sorted_densities[n / 2]) / 2.0
  } else {
    sorted_densities[n / 2]
  }

  let mut n_buried = 0
  for r in results {
    if r.density > mean_d {
      n_buried = n_buried + 1
    }
  }

  let n_exposed = n - n_buried

  PackingAnalysisResult::new(
    results~,
    mean_density=mean_d,
    min_density=min_d,
    max_density=max_d,
    median_density=median,
    n_residues=n,
    n_buried~,
    n_exposed~,
  )
}

///|
/// Calculate approximate Solvent Accessible Surface Area (SASA)
/// using the Lee & Richards (Shrake-Rupley) algorithm.
///
/// atoms: array of PackingAtom with coordinates and vdw radii.
/// probe_radius: water probe radius in Angstroms (default 1.4 Å).
/// n_points: number of sphere points per atom (default 100).
pub fn calculate_packing_sasa(
  atoms : Array[PackingAtom],
  probe_radius? : Double = 1.4,
  n_points? : Int = 100,
) -> Double {
  let n = atoms.length()
  if n == 0 {
    return 0.0
  }

  let sphere_points = generate_packing_sphere_points(n_points=n_points)

  let mut total_sasa = 0.0
  let mut i = 0
  while i < n {
    let atom_i = atoms[i]
    let radius_i = atom_i.vdw_radius + probe_radius
    let r_sq = radius_i * radius_i

    let mut n_accessible = 0
    for sp in sphere_points {
      let px = atom_i.x + sp.x * radius_i
      let py = atom_i.y + sp.y * radius_i
      let pz = atom_i.z + sp.z * radius_i

      let mut buried = false
      let mut j = 0
      while j < n {
        if j != i {
          let atom_j = atoms[j]
          let dx = px - atom_j.x
          let dy = py - atom_j.y
          let dz = pz - atom_j.z
          let dist_sq = dx * dx + dy * dy + dz * dz
          let radius_j = atom_j.vdw_radius + probe_radius
          if dist_sq < radius_j * radius_j {
            buried = true
            break
          }
        }
        j = j + 1
      }

      if !buried {
        n_accessible = n_accessible + 1
      }
    }

    let sphere_area = 4.0 * 3.14159265358979323846 * r_sq
    let atom_sasa = (n_accessible.to_double() / n_points.to_double()) * sphere_area
    total_sasa = total_sasa + atom_sasa

    i = i + 1
  }

  total_sasa
}

///|
/// Normalize packing density values relative to a threshold.
/// Values above threshold are set to 1.0, values below are
/// scaled linearly.
pub fn normalize_packing_density(
  density : PackingResult,
  threshold? : Double = 1.0,
) -> PackingResult {
  let norm_val = if threshold > 0.0 {
    if density.density > threshold {
      1.0
    } else {
      density.density / threshold
    }
  } else {
    density.density
  }

  PackingResult::new(
    residue_name=density.residue_name,
    residue_seq=density.residue_seq,
    chain_id=density.chain_id,
    density=norm_val,
    n_contacting_atoms=density.n_contacting_atoms,
    n_shell_points=density.n_shell_points,
  )
}

///|
/// Identify residues with low packing density.
/// Returns an array of residue identifier strings (chain+seq+name)
/// for residues whose density falls below the cutoff value.
pub fn identify_low_packing(
  analysis : PackingAnalysisResult,
  cutoff? : Double = 0.5,
) -> Array[String] {
  let low_residues : Array[String] = Array::new()
  let results = analysis.results

  for r in results {
    if r.density < cutoff {
      low_residues.push(r.get_residue_id())
    }
  }

  low_residues
}

///|
/// Calculate the radius of a sphere with given volume.
pub fn sphere_radius_from_volume(volume : Double) -> Double {
  if volume <= 0.0 {
    return 0.0
  }
  @math.pow((3.0 * volume) / (4.0 * 3.14159265358979323846), 1.0 / 3.0)
}

///|
/// Calculate the volume of a sphere with given radius.
pub fn sphere_volume(radius : Double) -> Double {
  (4.0 / 3.0) * 3.14159265358979323846 * radius * radius * radius
}

///|
/// Calculate packing efficiency as the ratio of
/// actual sphere volume to total volume.
pub fn calculate_packing_efficiency(
  atoms : Array[PackingAtom],
  container_radius? : Double = 100.0,
) -> Double {
  let n = atoms.length()
  if n == 0 {
    return 0.0
  }

  let total_volume = sphere_volume(container_radius)

  let mut atom_volume = 0.0
  for atom in atoms {
    atom_volume = atom_volume + sphere_volume(atom.vdw_radius)
  }

  if total_volume <= 0.0 {
    return 0.0
  }

  atom_volume / total_volume
}

///|
/// Create a demo set of PackingAtom arrays for testing.
pub fn create_demo_packing_atoms() -> Array[PackingAtom] {
  let atoms : Array[PackingAtom] = [
    PackingAtom::new(
      x=0.0,
      y=0.0,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
      chain_id='A',
      atom_name="N",
      element="N",
      vdw_radius=1.55,
    ),
    PackingAtom::new(
      x=1.47,
      y=0.0,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
      chain_id='A',
      atom_name="CA",
      element="C",
      vdw_radius=1.7,
    ),
    PackingAtom::new(
      x=2.01,
      y=1.42,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
      chain_id='A',
      atom_name="C",
      element="C",
      vdw_radius=1.7,
    ),
    PackingAtom::new(
      x=1.23,
      y=2.35,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
      chain_id='A',
      atom_name="O",
      element="O",
      vdw_radius=1.52,
    ),
    PackingAtom::new(
      x=3.33,
      y=1.53,
      z=0.0,
      residue_name="GLY",
      residue_seq=2,
      chain_id='A',
      atom_name="N",
      element="N",
      vdw_radius=1.55,
    ),
    PackingAtom::new(
      x=3.93,
      y=2.87,
      z=0.0,
      residue_name="GLY",
      residue_seq=2,
      chain_id='A',
      atom_name="CA",
      element="C",
      vdw_radius=1.7,
    ),
    PackingAtom::new(
      x=5.43,
      y=2.87,
      z=0.0,
      residue_name="GLY",
      residue_seq=2,
      chain_id='A',
      atom_name="C",
      element="C",
      vdw_radius=1.7,
    ),
    PackingAtom::new(
      x=6.02,
      y=3.92,
      z=0.0,
      residue_name="GLY",
      residue_seq=2,
      chain_id='A',
      atom_name="O",
      element="O",
      vdw_radius=1.52,
    ),
    PackingAtom::new(
      x=6.13,
      y=1.74,
      z=0.0,
      residue_name="VAL",
      residue_seq=3,
      chain_id='A',
      atom_name="N",
      element="N",
      vdw_radius=1.55,
    ),
    PackingAtom::new(
      x=7.54,
      y=1.42,
      z=0.0,
      residue_name="VAL",
      residue_seq=3,
      chain_id='A',
      atom_name="CA",
      element="C",
      vdw_radius=1.7,
    ),
  ]
  atoms
}

///|
/// Get Van der Waals radius for common atom names.
pub fn get_packing_vdw_radius(atom_name : String) -> Double {
  let name = atom_name.to_upper()
  if name == "H" || (contains(name, "H") && starts_with(name, 0, "H")) {
    return 1.2
  }
  if name == "C" || name == "CA" || name == "CB" || name == "CG" ||
     name == "CD" || name == "CE" || name == "CZ" || name == "CH2" || name == "CH3" {
    return 1.7
  }
  if name == "N" || name == "NZ" || name == "ND" || name == "NE" || name == "NH1" || name == "NH2" {
    return 1.55
  }
  if name == "O" || name == "OH" || name == "OG" || name == "OD1" || name == "OD2" ||
     name == "OE1" || name == "OE2" || name == "OT1" || name == "OT2" {
    return 1.52
  }
  if name == "S" || name == "SG" || name == "SD" {
    return 1.8
  }
  if name == "P" || name == "PA" || name == "PB" {
    return 1.8
  }
  1.7
}