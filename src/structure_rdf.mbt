///|
/// Radial distribution function g(r) — histogram of pairwise distances
/// normalized by the ideal-gas expectation.
/// Ported from biotite.structure.rdf + MDAnalysis.analysis.rdf.
///
/// g(r) = N(r, r+dr) / (4π·r²·dr·ρ·N)
/// where ρ = N/V (number density), N = number of reference atoms.

///|
/// RDF result: bin centers and g(r) values.
pub struct RDFResult {
  r : Array[Double]     // bin centers (Å)
  g_r : Array[Double]   // g(r) values
  n_bins : Int
  dr : Double
}

///|
/// Compute the radial distribution function for a set of atoms.
///
/// Parameters:
///   `atoms` — atom array (typically one element type from a trajectory frame)
///   `r_max` — maximum distance to histogram (Å), default 12.0
///   `n_bins` — number of histogram bins, default 100
///   `box_size` — optional periodic box dimensions [lx, ly, lz] for PBC
pub fn compute_rdf(
  atoms : Array[Atom],
  r_max~ : Double = 12.0,
  n_bins~ : Int = 100,
  box_size~ : Array[Double]? = None,
) -> RDFResult {
  let n = atoms.length()
  let dr = r_max / n_bins.to_double()
  // Histogram counts
  let counts : Array[Int] = Array::make(n_bins, 0)
  // Count pairs in each bin
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut dist = atoms[i].coord.distance(atoms[j].coord)
      // Apply minimum image convention if PBC box is provided
      match box_size {
        Some(box) => {
          if box.length() >= 3 {
            let dx = atoms[i].coord.x - atoms[j].coord.x
            let dy = atoms[i].coord.y - atoms[j].coord.y
            let dz = atoms[i].coord.z - atoms[j].coord.z
            // Minimum image
            let mx = dx - box[0] * round(dx / box[0])
            let my = dy - box[1] * round(dy / box[1])
            let mz = dz - box[2] * round(dz / box[2])
            dist = (mx * mx + my * my + mz * mz).sqrt()
          }
        }
        None => ()
      }
      if dist < r_max {
        let bin_idx = (dist / dr).to_int()
        if bin_idx >= 0 && bin_idx < n_bins {
          counts[bin_idx] = counts[bin_idx] + 2 // count both (i,j) and (j,i)
        }
      }
    }
  }
  // Normalize to get g(r)
  // Volume estimate: use cubic box or bounding box
  let volume = match box_size {
    Some(box) => {
      if box.length() >= 3 {
        box[0] * box[1] * box[2]
      } else {
        estimate_volume(atoms)
      }
    }
    None => estimate_volume(atoms)
  }
  let density = n.to_double() / volume
  let r_arr : Array[Double] = Array(capacity=n_bins)
  let g_arr : Array[Double] = Array(capacity=n_bins)
  for b = 0; b < n_bins; b = b + 1 {
    let r = (b.to_double() + 0.5) * dr
    r_arr.push(r)
    // Shell volume: 4π·r²·dr
    let shell_vol = 4.0 * 3.141592653589793 * r * r * dr
    // Expected count in ideal gas: ρ·shell_vol·N
    let expected = density * shell_vol * n.to_double()
    if expected > 0.0 {
      g_arr.push(counts[b].to_double() / expected)
    } else {
      g_arr.push(0.0)
    }
  }
  RDFResult::{ r: r_arr, g_r: g_arr, n_bins, dr }
}

///|
/// Round to nearest integer (for minimum image convention).
fn round(x : Double) -> Double {
  if x >= 0.0 {
    (x + 0.5).to_int().to_double()
  } else {
    -((-x + 0.5).to_int().to_double())
  }
}

///|
/// Estimate volume from atom bounding box (with 1Å padding).
fn estimate_volume(atoms : Array[Atom]) -> Double {
  if atoms.length() == 0 {
    return 1.0
  }
  let mut min_x = atoms[0].coord.x
  let mut max_x = atoms[0].coord.x
  let mut min_y = atoms[0].coord.y
  let mut max_y = atoms[0].coord.y
  let mut min_z = atoms[0].coord.z
  let mut max_z = atoms[0].coord.z
  for a in atoms {
    if a.coord.x < min_x {
      min_x = a.coord.x
    }
    if a.coord.x > max_x {
      max_x = a.coord.x
    }
    if a.coord.y < min_y {
      min_y = a.coord.y
    }
    if a.coord.y > max_y {
      max_y = a.coord.y
    }
    if a.coord.z < min_z {
      min_z = a.coord.z
    }
    if a.coord.z > max_z {
      max_z = a.coord.z
    }
  }
  (max_x - min_x + 2.0) * (max_y - min_y + 2.0) * (max_z - min_z + 2.0)
}
