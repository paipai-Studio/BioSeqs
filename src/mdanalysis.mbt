///|
/// MDAnalysis-style trajectory analysis on top of the minimal MD engine and
/// the Trajectory type: distance arrays, RMSD/RMSF, radial distribution
/// function and mean-squared displacement.

///|
fn mda_sq_dist(a : Array[Double], b : Array[Double]) -> Double {
  let dx = a[0] - b[0]
  let dy = a[1] - b[1]
  let dz = a[2] - b[2]
  dx * dx + dy * dy + dz * dz
}

///|
fn mda_euclid(a : Array[Double], b : Array[Double]) -> Double {
  mda_sq_dist(a, b).sqrt()
}

///|
fn mda_pts_to_qs(points : Array[Array[Double]]) -> Array[QsPoint] {
  let out : Array[QsPoint] = Array(capacity=0)
  for p in points {
    out.push(QsPoint::new(p[0], p[1], p[2]))
  }
  out
}

///|
/// All-pairs Euclidean distance between two coordinate sets
/// (`lib.distances.distance_array`): returns an N1 × N2 matrix.
pub fn mdanalysis_distance_array(
  coords1 : Array[Array[Double]],
  coords2 : Array[Array[Double]],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array::make(coords1.length(), [])
  for i = 0; i < coords1.length(); i = i + 1 {
    let row = Array::make(coords2.length(), 0.0)
    for j = 0; j < coords2.length(); j = j + 1 {
      row[j] = mda_euclid(coords1[i], coords2[j])
    }
    out[i] = row
  }
  out
}

///|
/// RMSD between two coordinate sets after optimal (quaternion) superposition
/// (`analysis.rms.rmsd`): `b` is rotated/translated onto `a` first.
pub fn mdanalysis_rmsd(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Double {
  if a.length() == 0 || a.length() != b.length() {
    return 0.0
  }
  let qa = mda_pts_to_qs(a)
  let qb = mda_pts_to_qs(b)
  let r = qs_superimpose(qa, qb)
  let b_sup = qs_apply_transform(qb, r.rotation, r.translation)
  let n = a.length()
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let bx = b_sup[i].x
    let by = b_sup[i].y
    let bz = b_sup[i].z
    let dx = a[i][0] - bx
    let dy = a[i][1] - by
    let dz = a[i][2] - bz
    sum = sum + dx * dx + dy * dy + dz * dz
  }
  (sum / n.to_double()).sqrt()
}

///|
/// Root mean square fluctuation per atom over a trajectory
/// (`analysis.rms.RMSF`): the per-atom deviation from its mean position.
pub fn mdanalysis_rmsf(traj : Array[Array[Array[Double]]]) -> Array[Double] {
  let n_frames = traj.length()
  let n_atoms = if n_frames > 0 { traj[0].length() } else { 0 }
  let mean : Array[Array[Double]] = Array::make(n_atoms, [])
  for a = 0; a < n_atoms; a = a + 1 {
    mean[a] = [0.0, 0.0, 0.0]
  }
  for t = 0; t < n_frames; t = t + 1 {
    for a = 0; a < n_atoms; a = a + 1 {
      mean[a][0] = mean[a][0] + traj[t][a][0]
      mean[a][1] = mean[a][1] + traj[t][a][1]
      mean[a][2] = mean[a][2] + traj[t][a][2]
    }
  }
  for a = 0; a < n_atoms; a = a + 1 {
    mean[a][0] = mean[a][0] / n_frames.to_double()
    mean[a][1] = mean[a][1] / n_frames.to_double()
    mean[a][2] = mean[a][2] / n_frames.to_double()
  }
  let acc = Array::make(n_atoms, 0.0)
  for t = 0; t < n_frames; t = t + 1 {
    for a = 0; a < n_atoms; a = a + 1 {
      acc[a] = acc[a] + mda_sq_dist(traj[t][a], mean[a])
    }
  }
  let out : Array[Double] = Array(capacity=0)
  for a = 0; a < n_atoms; a = a + 1 {
    out.push((acc[a] / n_frames.to_double()).sqrt())
  }
  out
}

///|
/// Radial distribution function g(r) between two coordinate sets
/// (`analysis.rdf.InterRDF`), normalized against a uniform background of
/// density `n2` in a `range_max` sphere. Returns `(bin_centers, g)`.
pub fn mdanalysis_rdf(
  coords1 : Array[Array[Double]],
  coords2 : Array[Array[Double]],
  nbins : Int,
  range_max : Double,
) -> (Array[Double], Array[Double]) {
  let n1 = coords1.length()
  let n2 = coords2.length()
  let dr = range_max / nbins.to_double()
  let hist = Array::make(nbins, 0.0)
  for i = 0; i < n1; i = i + 1 {
    for j = 0; j < n2; j = j + 1 {
      let d = mda_euclid(coords1[i], coords2[j])
      if d < range_max {
        let b = (d / dr).to_int()
        if b < nbins {
          hist[b] = hist[b] + 1.0
        }
      }
    }
  }
  let v_total = 4.1887902047863905 * range_max * range_max * range_max // 4π/3
  let centers = Array::make(nbins, 0.0)
  let g = Array::make(nbins, 0.0)
  for b = 0; b < nbins; b = b + 1 {
    let r_in = b.to_double() * dr
    let r_out = (b + 1).to_double() * dr
    let shell = 4.1887902047863905 * (r_out * r_out * r_out - r_in * r_in * r_in)
    centers[b] = (r_in + r_out) / 2.0
    if shell > 0.0 && n1 > 0 && n2 > 0 {
      g[b] = hist[b] * v_total / (n1.to_double() * n2.to_double() * shell)
    }
  }
  (centers, g)
}

///|
/// Mean-squared displacement over a trajectory (`analysis.msd.EinsteinMSD`),
/// via the windowed estimate: for lag `t`, MSD(t) = mean over time origins and
/// particles of `|x(τ+t) − x(τ)|²`. Returns one value per frame (index 0 = 0).
pub fn mdanalysis_msd(traj : Array[Array[Array[Double]]]) -> Array[Double] {
  let n_frames = traj.length()
  let n_atoms = if n_frames > 0 { traj[0].length() } else { 0 }
  let out = Array::make(n_frames, 0.0)
  for lag = 1; lag < n_frames; lag = lag + 1 {
    let mut sum = 0.0
    let mut count = 0
    let t_max = n_frames - lag
    for t = 0; t < t_max; t = t + 1 {
      let mut s = 0.0
      for a = 0; a < n_atoms; a = a + 1 {
        s = s + mda_sq_dist(traj[t + lag][a], traj[t][a])
      }
      sum = sum + s / n_atoms.to_double()
      count = count + 1
    }
    out[lag] = sum / count.to_double()
  }
  out
}
