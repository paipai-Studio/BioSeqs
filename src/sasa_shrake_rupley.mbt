///|
/// Shrake-Rupley SASA with Fibonacci sphere points and RSA classification.
///
/// Complements the existing `sasa.mbt` with an alternative implementation
/// using Fibonacci sphere point distribution and relative SASA analysis.

///|
/// Atom radius lookup (van der Waals radii in Angstroms).
fn sr_atom_radius(element : String) -> Double {
  match element {
    "H" => 1.2
    "C" => 1.7
    "N" => 1.55
    "O" => 1.52
    "S" => 1.8
    "P" => 1.8
    "F" => 1.47
    "CL" => 1.75
    "BR" => 1.85
    "I" => 1.98
    "FE" => 1.94
    "ZN" => 1.39
    "CU" => 1.4
    "MG" => 1.73
    "CA" => 2.31
    "NA" => 2.27
    "K" => 2.75
    "SE" => 1.9
    _ => 1.5
  }
}

///|
/// Generate Fibonacci sphere points (evenly distributed on unit sphere).
fn fibonacci_sphere(n : Int32) -> Array[Vector3] {
  let points : Array[Vector3] = []
  let golden_ratio = (1.0 + 2.23606797749979) / 2.0
  let angle_increment = 2.0 * 3.141592653589793 / golden_ratio
  for i in 0..<n {
    let t = (i.to_double() + 0.5) / n.to_double()
    let inclination = @math.acos(1.0 - 2.0 * t)
    let azimuth = angle_increment * i.to_double()
    let x = @math.sin(inclination) * @math.cos(azimuth)
    let y = @math.sin(inclination) * @math.sin(azimuth)
    let z = @math.cos(inclination)
    points.push(Vector3::new(x, y, z))
  }
  points
}

///|
/// Calculate SASA using the Shrake-Rupley algorithm with Fibonacci points.
pub fn shrake_rupley_sasa(
  coords : Array[Vector3],
  elements : Array[String],
  n_points~ : Int32 = 960,
  probe_radius~ : Double = 1.4,
) -> Array[Double] {
  let n = coords.length()
  let sasa_values = Array::make(n, 0.0)
  let sphere_points = fibonacci_sphere(n_points)
  let radii : Array[Double] = []
  for i in 0..<n {
    radii.push(sr_atom_radius(elements[i]) + probe_radius)
  }
  for i in 0..<n {
    let r = radii[i]
    let center = coords[i]
    let mut accessible = 0
    for sp in sphere_points {
      let test_point = Vector3::new(
        center.x + r * sp.x,
        center.y + r * sp.y,
        center.z + r * sp.z,
      )
      let mut is_accessible = true
      for j in 0..<n {
        if j != i {
          let dist = test_point.distance(coords[j])
          if dist < radii[j] {
            is_accessible = false
            break
          }
        }
      }
      if is_accessible {
        accessible = accessible + 1
      }
    }
    let area = 4.0 * 3.141592653589793 * r * r
    sasa_values[i] = area * accessible.to_double() / n_points.to_double()
  }
  sasa_values
}

///|
/// Calculate total SASA.
pub fn total_sasa(sasa_values : Array[Double]) -> Double {
  let mut total = 0.0
  for s in sasa_values {
    total = total + s
  }
  total
}

///|
/// Reference maximum SASA values (Tien et al. 2013).
fn reference_max_sasa(resname : String) -> Double {
  match resname {
    "ALA" => 129.0
    "ARG" => 274.0
    "ASN" => 195.0
    "ASP" => 193.0
    "CYS" => 167.0
    "GLN" => 224.0
    "GLU" => 223.0
    "GLY" => 104.0
    "HIS" => 224.0
    "ILE" => 197.0
    "LEU" => 201.0
    "LYS" => 236.0
    "MET" => 224.0
    "PHE" => 240.0
    "PRO" => 159.0
    "SER" => 155.0
    "THR" => 172.0
    "TRP" => 285.0
    "TYR" => 263.0
    "VAL" => 174.0
    _ => 200.0
  }
}

///|
/// Calculate relative SASA (RSA) for each residue.
pub fn calculate_rsa(residue_sasa : Array[(String, Double)]) -> Array[Double] {
  let result : Array[Double] = []
  for pair in residue_sasa {
    let max_sasa = reference_max_sasa(pair.0)
    if max_sasa > 0.0 {
      result.push(pair.1 / max_sasa)
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Residue exposure classification.
pub enum SasaExposure {
  Exposed
  PartiallyBuried
  Buried
}

///|
/// Classify a single residue's exposure.
pub fn classify_sasa_exposure(rsa : Double) -> SasaExposure {
  if rsa > 0.25 {
    Exposed
  } else if rsa > 0.05 {
    PartiallyBuried
  } else {
    Buried
  }
}
