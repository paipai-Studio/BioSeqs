/// biopython_sasa.mbt — Solvent Accessible Surface Area (SASA)
/// Inspired by Biopython Bio.PDB.SASA: Shrake-Rupley algorithm

///|
/// Atom representation for SASA calculation.
pub struct BiopythonAtom {
  // Struct for SASA calculation

  name : String
  element : String
  x : Double
  y : Double
  z : Double
  radius : Double // van der Waals radius
}

///|
/// Create a new atom.
pub fn BiopythonAtom::new(
  name~ : String,
  element~ : String,
  x~ : Double,
  y~ : Double,
  z~ : Double,
  radius~ : Double = 1.7,
) -> BiopythonAtom {
  { name, element, x, y, z, radius, }
}

///|
/// Calculate SASA using the Shrake-Rupley algorithm.
/// Places N test points on a sphere around each atom, counts points
/// not buried inside other atoms.
pub fn biopython_sasa_shrake_rupley(
  atoms : Array[BiopythonAtom],
  n_points~ : Int32 = 960, // number of test points per atom
  probe_radius~ : Double = 1.4, // water probe radius
) -> Array[Double] {
  let n_atoms = atoms.length()
  let sasa_values : Array[Double] = Array::make(n_atoms, 0.0)
  if n_atoms == 0 {
    return sasa_values
  }
  // Generate Fibonacci sphere points
  let sphere_points = biopython_fibonacci_sphere(n_points)
  for i = 0; i < n_atoms; i = i + 1 {
    let atom = atoms[i]
    let r = atom.radius + probe_radius
    let mut accessible = 0
    for p = 0; p < n_points; p = p + 1 {
      // Test point on sphere around atom i
      let tx = atom.x + r * sphere_points[p].0
      let ty = atom.y + r * sphere_points[p].1
      let tz = atom.z + r * sphere_points[p].2
      let mut buried = false
      for j = 0; j < n_atoms; j = j + 1 {
        if i == j {
          continue
        }
        let other = atoms[j]
        let dx = tx - other.x
        let dy = ty - other.y
        let dz = tz - other.z
        let dist_sq = dx * dx + dy * dy + dz * dz
        let cutoff = other.radius + probe_radius
        if dist_sq < cutoff * cutoff {
          buried = true
          break
        }
      }
      if !buried {
        accessible = accessible + 1
      }
    }
    // SASA = (accessible / n_points) * 4 * pi * r^2
    let area = 4.0 * 3.14159265358979323846 * r * r
    sasa_values[i] = accessible.to_double() / n_points.to_double() * area
  }
  sasa_values
}

///|
/// Generate points on a sphere using Fibonacci spiral.
fn biopython_fibonacci_sphere(n : Int32) -> Array[(Double, Double, Double)] {
  let points : Array[(Double, Double, Double)] = Array::make(n, (0.0, 0.0, 0.0))
  let golden = (1.0 + @math.pow(5.0, 0.5)) / 2.0
  for i = 0; i < n; i = i + 1 {
    let t = (i.to_double() + 0.5) / n.to_double()
    let phi = @math.acos(1.0 - 2.0 * t) // polar angle
    let theta = 2.0 * 3.14159265358979323846 * i.to_double() / golden // azimuthal
    let x = biopython_local_sin(phi) * biopython_local_cos(theta)
    let y = biopython_local_sin(phi) * biopython_local_sin(theta)
    let z = biopython_local_cos(phi)
    points[i] = (x, y, z)
  }
  points
}

///|
/// Get van der Waals radius for an element.
pub fn biopython_vdw_radius(element : String) -> Double {
  match element {
    "H" => 1.20
    "C" => 1.70
    "N" => 1.55
    "O" => 1.52
    "F" => 1.47
    "P" => 1.80
    "S" => 1.80
    "CL" | "Cl" => 1.75
    "BR" | "Br" => 1.85
    "I" => 1.98
    "FE" | "Fe" => 2.0
    "ZN" | "Zn" => 1.39
    "CA" | "Ca" => 2.0
    "MG" | "Mg" => 1.72
    "NA" | "Na" => 2.27
    "K" => 2.75
    _ => 1.70
  }
}

///|
/// Total SASA of all atoms.
pub fn biopython_total_sasa(sasa_values : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < sasa_values.length(); i = i + 1 {
    total = total + sasa_values[i]
  }
  total
}

///|
/// Relative SASA: ratio of observed to maximum possible SASA.
pub fn biopython_relative_sasa(observed : Double, max_sasa : Double) -> Double {
  if max_sasa == 0.0 {
    0.0
  } else {
    observed / max_sasa
  }
}

///|
/// Classify residue burial based on relative SASA.
/// Returns: "buried" (< 0.2), "intermediate" (0.2-0.5), "exposed" (> 0.5).
pub fn biopython_burial_class(relative_sasa : Double) -> String {
  if relative_sasa < 0.2 {
    "buried"
  } else if relative_sasa < 0.5 {
    "intermediate"
  } else {
    "exposed"
  }
}

///|
fn biopython_local_sin(x : Double) -> Double {
  @math.sin(x)
}

///|
fn biopython_local_cos(x : Double) -> Double {
  @math.cos(x)
}
