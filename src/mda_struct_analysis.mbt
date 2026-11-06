/// mda_struct_analysis.mbt — ANM, simplified DSSP, Voronoi, bond geometry, protein-membrane
/// Items 11-15

///|
pub struct MdaAnmResult {
  n_modes : Int32
  eigenvalues : Array[Double]
  hessian : Array[Array[Double]]
  mode_fluctuations : Array[Double]
}

///|
pub fn mda_anm(
  coords : Array[(Double, Double, Double)],
  cutoff~ : Double = 12.0,
  gamma~ : Double = 1.0,
) -> MdaAnmResult {
  let n = coords.length()
  let dim = n * 3
  let hessian : Array[Array[Double]] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let row : Array[Double] = Array(capacity=dim)
    for j = 0; j < dim; j = j + 1 {
      row.push(0.0)
    }
    hessian.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let d2 = dx * dx + dy * dy + dz * dz
      let d = @math.pow(d2, 0.5)
      if d > cutoff || d < 0.0000000001 {
        continue
      }
      let g = -gamma / d2
      let hxx = g * dx * dx
      let hxy = g * dx * dy
      let hxz = g * dx * dz
      let hyy = g * dy * dy
      let hyz = g * dy * dz
      let hzz = g * dz * dz
      let i3 = i * 3
      let j3 = j * 3
      mda_anm_fill(hessian, i3, j3, hxx, hxy, hxz, hyy, hyz, hzz)
    }
  }
  let eigenvalues : Array[Double] = Array(capacity=dim)
  let mode_fluctuations : Array[Double] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let ev = hessian[i][i]
    eigenvalues.push(ev)
    if mda_abs(ev) > 0.0000000001 {
      mode_fluctuations.push(1.0 / ev)
    } else {
      mode_fluctuations.push(0.0)
    }
  }
  MdaAnmResult::{ n_modes: dim, eigenvalues, hessian, mode_fluctuations, }
}

///|
fn mda_anm_fill(
  h : Array[Array[Double]],
  i3 : Int32,
  j3 : Int32,
  hxx : Double,
  hxy : Double,
  hxz : Double,
  hyy : Double,
  hyz : Double,
  hzz : Double,
) -> Unit {
  h[i3][j3] = hxx
  h[i3][j3 + 1] = hxy
  h[i3][j3 + 2] = hxz
  h[i3 + 1][j3] = hxy
  h[i3 + 1][j3 + 1] = hyy
  h[i3 + 1][j3 + 2] = hyz
  h[i3 + 2][j3] = hxz
  h[i3 + 2][j3 + 1] = hyz
  h[i3 + 2][j3 + 2] = hzz
  h[j3][i3] = hxx
  h[j3][i3 + 1] = hxy
  h[j3][i3 + 2] = hxz
  h[j3 + 1][i3] = hxy
  h[j3 + 1][i3 + 1] = hyy
  h[j3 + 1][i3 + 2] = hyz
  h[j3 + 2][i3] = hxz
  h[j3 + 2][i3 + 1] = hyz
  h[j3 + 2][i3 + 2] = hzz
  h[i3][i3] = h[i3][i3] - hxx
  h[i3][i3 + 1] = h[i3][i3 + 1] - hxy
  h[i3][i3 + 2] = h[i3][i3 + 2] - hxz
  h[i3 + 1][i3] = h[i3 + 1][i3] - hxy
  h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] - hyy
  h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] - hyz
  h[i3 + 2][i3] = h[i3 + 2][i3] - hxz
  h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] - hyz
  h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] - hzz
  h[j3][j3] = h[j3][j3] - hxx
  h[j3][j3 + 1] = h[j3][j3 + 1] - hxy
  h[j3][j3 + 2] = h[j3][j3 + 2] - hxz
  h[j3 + 1][j3] = h[j3 + 1][j3] - hxy
  h[j3 + 1][j3 + 1] = h[j3 + 1][j3 + 1] - hyy
  h[j3 + 1][j3 + 2] = h[j3 + 1][j3 + 2] - hyz
  h[j3 + 2][j3] = h[j3 + 2][j3] - hxz
  h[j3 + 2][j3 + 1] = h[j3 + 2][j3 + 1] - hyz
  h[j3 + 2][j3 + 2] = h[j3 + 2][j3 + 2] - hzz
}

///|
/// DSSP-like secondary structure assignment (simplified).
pub fn mda_dssp_assign(
  ca_coords : Array[(Double, Double, Double)],
  _res_names : Array[String],
) -> Array[String] {
  let n = ca_coords.length()
  let ss : Array[String] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    ss.push("C")
  }
  for i = 0; i < n - 3; i = i + 1 {
    let d = mda_sa_dist(ca_coords[i], ca_coords[i + 3])
    if d > 4.5 && d < 6.0 {
      ss[i] = "H"
      ss[i + 1] = "H"
      ss[i + 2] = "H"
      ss[i + 3] = "H"
    }
  }
  for i = 0; i < n - 2; i = i + 1 {
    if ss[i] == "H" {
      continue
    }
    let d = mda_sa_dist(ca_coords[i], ca_coords[i + 2])
    if d > 6.5 && d < 8.0 {
      if ss[i] == "C" {
        ss[i] = "E"
      }
      if ss[i + 2] == "C" {
        ss[i + 2] = "E"
      }
    }
  }
  for i = 0; i < n - 2; i = i + 1 {
    if ss[i] == "C" {
      let d = mda_sa_dist(ca_coords[i], ca_coords[i + 2])
      if d < 5.5 {
        ss[i] = "T"
      }
    }
  }
  ss
}

///|
/// Voronoi cell volume estimation (simplified).
pub fn mda_voronoi_volumes(
  coords : Array[(Double, Double, Double)],
  cutoff~ : Double = 5.0,
) -> Array[Double] {
  let n = coords.length()
  let volumes : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let dists : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let d = mda_sa_dist(coords[i], coords[j])
      if d <= cutoff {
        dists.push(d)
      }
    }
    dists.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
    if dists.length() >= 3 {
      volumes.push(dists[0] * dists[1] * dists[2] / 6.0)
    } else if dists.length() > 0 {
      volumes.push(@math.pow(dists[0], 3.0) / 6.0)
    } else {
      volumes.push(0.0)
    }
  }
  volumes
}

///|
pub struct MdaBondGeomStats {
  mean : Double
  std : Double
  min_ : Double
  max_ : Double
  n_samples : Int32
}

///|
pub fn mda_bond_length_stats(
  frames : Array[Array[(Double, Double, Double)]],
  bond_pairs : Array[(Int32, Int32)],
) -> MdaBondGeomStats {
  let values : Array[Double] = Array(capacity=0)
  for f = 0; f < frames.length(); f = f + 1 {
    for b = 0; b < bond_pairs.length(); b = b + 1 {
      let (a1, a2) = bond_pairs[b]
      if a1 < frames[f].length() && a2 < frames[f].length() {
        values.push(mda_sa_dist(frames[f][a1], frames[f][a2]))
      }
    }
  }
  mda_sa_compute_stats(values)
}

///|
pub fn mda_bond_angle_stats(
  frames : Array[Array[(Double, Double, Double)]],
  angle_triples : Array[(Int32, Int32, Int32)],
) -> MdaBondGeomStats {
  let values : Array[Double] = Array(capacity=0)
  for f = 0; f < frames.length(); f = f + 1 {
    for t = 0; t < angle_triples.length(); t = t + 1 {
      let (a1, a2, a3) = angle_triples[t]
      if a1 < frames[f].length() &&
        a2 < frames[f].length() &&
        a3 < frames[f].length() {
        values.push(mda_sa_angle(frames[f][a1], frames[f][a2], frames[f][a3]))
      }
    }
  }
  mda_sa_compute_stats(values)
}

///|
pub struct MdaProteinMembraneResult {
  insertion_depth : Double
  tilt_angle : Double
  n_transmembrane : Int32
  transmembrane_residues : Array[Int32]
}

///|
pub fn mda_protein_membrane(
  ca_coords : Array[(Double, Double, Double)],
  res_ids : Array[Int32],
  membrane_center~ : Double = 0.0,
  membrane_thickness~ : Double = 30.0,
) -> MdaProteinMembraneResult {
  let n = ca_coords.length()
  if n == 0 {
    return MdaProteinMembraneResult::{
      insertion_depth: 0.0,
      tilt_angle: 0.0,
      n_transmembrane: 0,
      transmembrane_residues: Array(capacity=0),
    }
  }
  let mut avg_z = 0.0
  for i = 0; i < n; i = i + 1 {
    avg_z = avg_z + ca_coords[i].2
  }
  avg_z = avg_z / n.to_double()
  let insertion_depth = avg_z - membrane_center
  let half_thick = membrane_thickness / 2.0
  let tm_residues : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let z = ca_coords[i].2 - membrane_center
    if mda_abs(z) <= half_thick {
      tm_residues.push(res_ids[i])
    }
  }
  let mut x_min = 999999.0
  let mut x_max = -999999.0
  let mut z_min = 999999.0
  let mut z_max = -999999.0
  for i = 0; i < n; i = i + 1 {
    if ca_coords[i].0 < x_min {
      x_min = ca_coords[i].0
    }
    if ca_coords[i].0 > x_max {
      x_max = ca_coords[i].0
    }
    if ca_coords[i].2 < z_min {
      z_min = ca_coords[i].2
    }
    if ca_coords[i].2 > z_max {
      z_max = ca_coords[i].2
    }
  }
  let x_range = x_max - x_min
  let z_range = z_max - z_min
  let tilt_angle = if z_range > 0.0000000001 {
    @math.atan(x_range / z_range) * 180.0 / 3.141592653589793
  } else {
    90.0
  }
  MdaProteinMembraneResult::{
    insertion_depth,
    tilt_angle,
    n_transmembrane: tm_residues.length(),
    transmembrane_residues: tm_residues,
  }
}

///|
fn mda_sa_dist(
  a : (Double, Double, Double),
  b : (Double, Double, Double),
) -> Double {
  let dx = a.0 - b.0
  let dy = a.1 - b.1
  let dz = a.2 - b.2
  @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
}

///|
fn mda_sa_angle(
  a : (Double, Double, Double),
  b : (Double, Double, Double),
  c : (Double, Double, Double),
) -> Double {
  let v1x = a.0 - b.0
  let v1y = a.1 - b.1
  let v1z = a.2 - b.2
  let v2x = c.0 - b.0
  let v2y = c.1 - b.1
  let v2z = c.2 - b.2
  let dot = v1x * v2x + v1y * v2y + v1z * v2z
  let l1 = @math.pow(v1x * v1x + v1y * v1y + v1z * v1z, 0.5)
  let l2 = @math.pow(v2x * v2x + v2y * v2y + v2z * v2z, 0.5)
  if l1 < 0.0000000001 || l2 < 0.0000000001 {
    return 0.0
  }
  let cos_angle = dot / (l1 * l2)
  let clamped = if cos_angle > 1.0 {
    1.0
  } else if cos_angle < -1.0 {
    -1.0
  } else {
    cos_angle
  }
  @math.acos(clamped) * 180.0 / 3.141592653589793
}

///|
fn mda_sa_compute_stats(values : Array[Double]) -> MdaBondGeomStats {
  let n = values.length()
  if n == 0 {
    return MdaBondGeomStats::{
      mean: 0.0,
      std: 0.0,
      min_: 0.0,
      max_: 0.0,
      n_samples: 0,
    }
  }
  let mut sum = 0.0
  let mut min_v = values[0]
  let mut max_v = values[0]
  for i = 0; i < n; i = i + 1 {
    sum = sum + values[i]
    if values[i] < min_v {
      min_v = values[i]
    }
    if values[i] > max_v {
      max_v = values[i]
    }
  }
  let mean = sum / n.to_double()
  let mut variance = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = values[i] - mean
    variance = variance + diff * diff
  }
  variance = variance / n.to_double()
  MdaBondGeomStats::{
    mean,
    std: @math.pow(variance, 0.5),
    min_: min_v,
    max_: max_v,
    n_samples: n,
  }
}

///|
pub fn mda_struct_analysis_dummy() -> Unit {
  ignore(MdaAnmResult::{
    n_modes: 0,
    eigenvalues: Array(capacity=0),
    hessian: Array(capacity=0),
    mode_fluctuations: Array(capacity=0),
  })
  ignore(MdaBondGeomStats::{
    mean: 0.0,
    std: 0.0,
    min_: 0.0,
    max_: 0.0,
    n_samples: 0,
  })
  ignore(MdaProteinMembraneResult::{
    insertion_depth: 0.0,
    tilt_angle: 0.0,
    n_transmembrane: 0,
    transmembrane_residues: Array(capacity=0),
  })
}
