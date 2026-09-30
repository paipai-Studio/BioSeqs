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

///|
/// Hydrogen-bond detection (`analysis.hydrogenbonds.HydrogenBondAnalysis`):
/// for each donor–hydrogen pair and acceptor, a hydrogen bond exists when the
/// donor–acceptor distance ≤ `d_a_cutoff` (default 3.0 Å) and the D-H···A
/// angle at the donor ≥ `angle_cutoff_deg` (default 150°).
/// Returns `(donor_idx, hydrogen_idx, acceptor_idx, d_a_dist, dha_angle_deg)`.
pub fn mdanalysis_hydrogen_bonds(
  coords : Array[Array[Double]],
  donor_hydrogen_pairs : Array[(Int, Int)],
  acceptors : Array[Int],
  d_a_cutoff : Double,
  angle_cutoff_deg : Double,
) -> Array[(Int, Int, Int, Double, Double)] {
  let result : Array[(Int, Int, Int, Double, Double)] = Array(capacity=0)
  for pair in donor_hydrogen_pairs {
    let (d_idx, h_idx) = pair
    let d = coords[d_idx]
    let h = coords[h_idx]
    for a_idx in acceptors {
      if a_idx != d_idx {
        let a = coords[a_idx]
        let da = mda_euclid(d, a)
        if da <= d_a_cutoff {
          let dh_x = d[0] - h[0]
          let dh_y = d[1] - h[1]
          let dh_z = d[2] - h[2]
          let ah_x = a[0] - h[0]
          let ah_y = a[1] - h[1]
          let ah_z = a[2] - h[2]
          let dot = dh_x * ah_x + dh_y * ah_y + dh_z * ah_z
          let mag_dh = (dh_x * dh_x + dh_y * dh_y + dh_z * dh_z).sqrt()
          let mag_ah = (ah_x * ah_x + ah_y * ah_y + ah_z * ah_z).sqrt()
          if mag_dh > 0.0 && mag_ah > 0.0 {
            let cos_t = dot / (mag_dh * mag_ah)
            let cos_c = if cos_t > 1.0 {
              1.0
            } else if cos_t < -1.0 {
              -1.0
            } else {
              cos_t
            }
            let angle_deg = @math.acos(cos_c) * 57.29577951308232
            if angle_deg >= angle_cutoff_deg {
              result.push((d_idx, h_idx, a_idx, da, angle_deg))
            }
          }
        }
      }
    }
  }
  result
}

///|
/// Protein backbone φ/ψ dihedrals (`analysis.dihedrals.Ramachandran`):
/// φ(i) = C'(i-1)–N–CA–C, ψ(i) = N–CA–C–N'(i+1), in degrees.
/// Input is one (N, CA, C) coordinate triple per residue; the first residue's
/// φ and the last residue's ψ are `None` (no flanking atom).
pub fn mdanalysis_phi_psi(
  backbone : Array[(Array[Double], Array[Double], Array[Double])],
) -> Array[(Double?, Double?)] {
  let n = backbone.length()
  let result : Array[(Double?, Double?)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let (n_i, ca_i, c_i) = backbone[i]
    let phi : Double? = if i > 0 {
      let (_, _, c_prev) = backbone[i - 1]
      Some(chem_dihedral_angle(c_prev, n_i, ca_i, c_i))
    } else {
      None
    }
    let psi : Double? = if i < n - 1 {
      let (n_next, _, _) = backbone[i + 1]
      Some(chem_dihedral_angle(n_i, ca_i, c_i, n_next))
    } else {
      None
    }
    result.push((phi, psi))
  }
  result
}

// ============================================================================
// Contact Analysis (MDAnalysis.analysis.contacts)
// ============================================================================

///|
/// Native contact analysis between two coordinate sets.
/// Returns the list of index pairs `(i, j)` where `coords_a[i]` and
/// `coords_b[j]` are within `cutoff` Å of each other.
/// Mirrors `MDAnalysis.analysis.contacts.Contacts` initial-contact detection.
pub fn mdanalysis_contacts(
  coords_a : Array[Array[Double]],
  coords_b : Array[Array[Double]],
  cutoff : Double,
) -> Array[(Int, Int)] {
  let cutoff_sq = cutoff * cutoff
  let pairs : Array[(Int, Int)] = Array(capacity=0)
  for i = 0; i < coords_a.length(); i = i + 1 {
    for j = 0; j < coords_b.length(); j = j + 1 {
      if mda_sq_dist(coords_a[i], coords_b[j]) <= cutoff_sq {
        pairs.push((i, j))
      }
    }
  }
  pairs
}

///|
/// Fraction of native contacts Q(t) for a trajectory.
/// `reference` defines the initial contacts (pairs within `cutoff`).
/// For each frame, returns the fraction of reference contacts still present.
/// `hard_cut = true` uses the same cutoff for every frame (hard-cut metric);
/// `hard_cut = false` uses a soft switching function (1 / (1 + (r/r0)^6)).
pub fn mdanalysis_contact_fraction(
  reference_a : Array[Array[Double]],
  reference_b : Array[Array[Double]],
  frames_a : Array[Array[Array[Double]]],
  frames_b : Array[Array[Array[Double]]],
  cutoff : Double,
  hard_cut : Bool,
) -> Array[Double] {
  // Determine initial contacts from reference
  let native = mdanalysis_contacts(reference_a, reference_b, cutoff)
  let n_native = native.length()
  if n_native == 0 {
    return Array::make(frames_a.length(), 0.0)
  }
  let cutoff_sq = cutoff * cutoff
  let result : Array[Double] = Array(capacity=frames_a.length())
  for f = 0; f < frames_a.length(); f = f + 1 {
    let fa = frames_a[f]
    let fb = frames_b[f]
    let mut q = 0.0
    for k = 0; k < n_native; k = k + 1 {
      let (i, j) = native[k]
      let r_sq = mda_sq_dist(fa[i], fb[j])
      if hard_cut {
        if r_sq <= cutoff_sq {
          q = q + 1.0
        }
      } else {
        // soft: r0 from reference
        let r0_sq = mda_sq_dist(reference_a[i], reference_b[j])
        let ratio = r0_sq / r_sq
        q = q + (1.0 / (1.0 + ratio * ratio * ratio))
      }
    }
    result.push(q / n_native.to_double())
  }
  result
}

// ============================================================================
// Ramachandran Plot (MDAnalysis.analysis.dihedrals.Ramachandran)
// ============================================================================

///|
/// Ramachandran region classification for a single (phi, psi) pair.
/// Returns:
/// - `"core"`      — favoured region (right-handed α-helix / β-sheet)
/// - `"allowed"`   — marginally allowed region
/// - `"outlier"`   — disallowed region
/// Classification uses a simplified rectangular approximation of the
/// Lovell et al. (2003) Ramachandran contours.
pub fn mdanalysis_ramachandran_region(phi : Double, psi : Double) -> String {
  // Right-handed α-helix:  φ ≈ -60°, ψ ≈ -45°
  // β-sheet:               φ ≈ -120°, ψ ≈ 120°
  // Left-handed α-helix:   φ ≈ 60°, ψ ≈ 45°
  // Core regions (generous rectangles)
  
  let abs_psi = psi.abs()
  // αR core: φ ∈ [-90, -30], ψ ∈ [-70, -10]
  if phi >= -90.0 && phi <= -30.0 && psi >= -70.0 && psi <= -10.0 {
    return "core"
  }
  // β core: φ ∈ [-150, -90], ψ ∈ [90, 180] or ψ ∈ [-180, -150]
  if phi >= -150.0 && phi <= -90.0 && (psi >= 90.0 || psi <= -150.0) {
    return "core"
  }
  // αL core: φ ∈ [30, 90], ψ ∈ [10, 70]
  if phi >= 30.0 && phi <= 90.0 && psi >= 10.0 && psi <= 70.0 {
    return "core"
  }
  // Allowed regions (wider rectangles)
  // αR allowed
  if phi >= -120.0 && phi <= 0.0 && psi >= -100.0 && psi <= 30.0 {
    return "allowed"
  }
  // β allowed
  if phi >= -180.0 && phi <= -60.0 && abs_psi >= 80.0 {
    return "allowed"
  }
  // αL allowed
  if phi >= 0.0 && phi <= 120.0 && psi >= -30.0 && psi <= 100.0 {
    return "allowed"
  }
  "outlier"
}

///|
/// Multi-frame Ramachandran analysis.
/// `backbone_frames` is an array of frames, each frame being the backbone
/// `[(N, CA, C), ...]` for each residue (same format as `mdanalysis_phi_psi`).
/// Returns `(angles, regions)` where:
/// - `angles[f]`  = array of `(phi, psi)` for frame f (same as `mdanalysis_phi_psi`)
/// - `regions[f]` = array of `"core"` / `"allowed"` / `"outlier"` for each residue
pub fn mdanalysis_ramachandran(
  backbone_frames : Array[Array[(Array[Double], Array[Double], Array[Double])]],
) -> (
    Array[Array[(Double?, Double?)]],
    Array[Array[String?]],
  ) {
  let all_angles : Array[Array[(Double?, Double?)]] = Array(capacity=
    backbone_frames.length()
  )
  let all_regions : Array[Array[String?]] = Array(capacity=
    backbone_frames.length()
  )
  for frame in backbone_frames {
    let angles = mdanalysis_phi_psi(frame)
    let regions : Array[String?] = Array(capacity=angles.length())
    for pair in angles {
      match pair {
        (Some(phi), Some(psi)) => regions.push(
            Some(mdanalysis_ramachandran_region(phi, psi)),
          )
        _ => regions.push(None)
      }
    }
    all_angles.push(angles)
    all_regions.push(regions)
  }
  (all_angles, all_regions)
}

// ============================================================================
// Linear Density (MDAnalysis.analysis.lineardensity.LinearDensity)
// ============================================================================

///|
/// 1-D linear density histogram along a single axis.
/// `coords` are atomic positions; `axis` is 0 (x), 1 (y), or 2 (z).
/// `nbins` bins are spread over `[xmin, xmax]`.
/// If `weights` is provided (e.g. masses or charges), each atom contributes
/// its weight to the bin instead of 1.  Returns the bin counts/weights.
pub fn mdanalysis_linear_density_1d(
  coords : Array[Array[Double]],
  axis : Int,
  nbins : Int,
  xmin : Double,
  xmax : Double,
  weights : Array[Double]?,
) -> Array[Double] {
  let hist : Array[Double] = Array::make(nbins, 0.0)
  let bin_width = (xmax - xmin) / nbins.to_double()
  if bin_width <= 0.0 {
    return hist
  }
  let has_weights = weights.is_some()
  let w = weights.unwrap_or([])
  for i = 0; i < coords.length(); i = i + 1 {
    let val = coords[i][axis]
    if val < xmin || val >= xmax {
      continue
    }
    let bin = ((val - xmin) / bin_width).to_int()
    let b = if bin >= nbins { nbins - 1 } else { bin }
    let wt = if has_weights { w[i] } else { 1.0 }
    hist[b] = hist[b] + wt
  }
  hist
}

///|
/// Linear density along all three axes (x, y, z).
/// Returns `(hist_x, hist_y, hist_z)`, each an array of `nbins` values.
/// `box` is `[lx, ly, lz]`; bins span `[0, l_dim]` per axis.
/// If `weights` is provided, it is used for all three histograms.
pub fn mdanalysis_linear_density(
  coords : Array[Array[Double]],
  nbins : Int,
  box : Array[Double],
  weights : Array[Double]?,
) -> (Array[Double], Array[Double], Array[Double]) {
  let hx = mdanalysis_linear_density_1d(coords, 0, nbins, 0.0, box[0], weights)
  let hy = mdanalysis_linear_density_1d(coords, 1, nbins, 0.0, box[1], weights)
  let hz = mdanalysis_linear_density_1d(coords, 2, nbins, 0.0, box[2], weights)
  (hx, hy, hz)
}

// ============================================================================
// Persistence Length (MDAnalysis.analysis.polymer.PersistenceLength)
// ============================================================================

///|
/// Bond-vector autocorrelation for a single polymer chain.
/// `positions` are ordered atomic positions along the chain.
/// Returns `C(n)` for n = 0 .. N-2, where
///   C(n) = ⟨ cos θ_{i,i+n} ⟩ = ⟨ â_i · â_{i+n} ⟩
/// averaged over all valid i.  Also returns the average bond length `lb`.
pub fn mdanalysis_bond_autocorrelation(
  positions : Array[Array[Double]],
) -> (Array[Double], Double) {
  let n = positions.length()
  if n < 2 {
    return ([], 0.0)
  }
  // Compute unit bond vectors
  let n_bonds = n - 1
  let bonds : Array[Array[Double]] = Array(capacity=n_bonds)
  let mut total_bl = 0.0
  for i = 0; i < n_bonds; i = i + 1 {
    let dx = positions[i + 1][0] - positions[i][0]
    let dy = positions[i + 1][1] - positions[i][1]
    let dz = positions[i + 1][2] - positions[i][2]
    let bl = (dx * dx + dy * dy + dz * dz).sqrt()
    total_bl = total_bl + bl
    if bl > 0.0 {
      bonds.push([dx / bl, dy / bl, dz / bl])
    } else {
      bonds.push([0.0, 0.0, 0.0])
    }
  }
  let lb = total_bl / n_bonds.to_double()
  // C(n) for n = 0 .. n_bonds-1
  let autocorr : Array[Double] = Array(capacity=n_bonds)
  for sep = 0; sep < n_bonds; sep = sep + 1 {
    let mut sum = 0.0
    let mut count = 0
    for i = 0; i + sep < n_bonds; i = i + 1 {
      let a = bonds[i]
      let b = bonds[i + sep]
      sum = sum + (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
      count = count + 1
    }
    if count > 0 {
      autocorr.push(sum / count.to_double())
    } else {
      autocorr.push(0.0)
    }
  }
  (autocorr, lb)
}

///|
/// Persistence length from bond autocorrelation.
/// Fits `C(n) ≈ exp(-n * lb / lp)` by linear regression on `ln(C(n))` vs `n`.
/// Returns `(lp, lb, autocorr, fit)` where `fit` is the fitted exponential
/// curve evaluated at each n.
pub fn mdanalysis_persistence_length(
  positions : Array[Array[Double]],
) -> (Double, Double, Array[Double], Array[Double]) {
  let (autocorr, lb) = mdanalysis_bond_autocorrelation(positions)
  let n = autocorr.length()
  if n < 2 || lb <= 0.0 {
    return (0.0, lb, autocorr, [])
  }
  // Linear regression: ln(C(n)) = -n * lb / lp
  // Only use points where C(n) > 0
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xx = 0.0
  let mut sum_xy = 0.0
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    let c = autocorr[i]
    if c > 0.0 {
      let x = i.to_double()
      let y = @math.ln(c)
      sum_x = sum_x + x
      sum_y = sum_y + y
      sum_xx = sum_xx + x * x
      sum_xy = sum_xy + x * y
      count = count + 1
    }
  }
  let lp = if count > 1 {
    let m = count.to_double()
    let denom = m * sum_xx - sum_x * sum_x
    if denom.abs() > 1.0e-12 {
      let slope = (m * sum_xy - sum_x * sum_y) / denom
      // slope = -lb / lp  →  lp = -lb / slope
      if slope < 0.0 {
        -lb / slope
      } else {
        0.0
      }
    } else {
      0.0
    }
  } else {
    0.0
  }
  // Fitted curve
  let fit : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let x = i.to_double() * lb
    fit.push(@math.exp(-x / lp))
  }
  (lp, lb, autocorr, fit)
}

// ============================================================================
// PCA (MDAnalysis.analysis.pca.PCA)
// ============================================================================

///|
/// Principal Component Analysis of a trajectory.
/// `frames` is an array of frames, each frame being N atoms × 3D coordinates.
/// Builds the N×N covariance matrix of atomic fluctuations (each atom's
/// position minus its time-average), then performs eigendecomposition.
///
/// Returns `(eigenvalues, eigenvectors, variance_fraction)` where:
/// - `eigenvalues`     — sorted descending (variance per component)
/// - `eigenvectors[i]` — i-th principal component (N-dimensional unit vector)
/// - `variance_fraction` — cumulative fraction of variance explained
///
/// For very large N this is O(N³); use on CA atoms or other reduced sets.
pub fn mdanalysis_pca(
  frames : Array[Array[Array[Double]]],
) -> (Array[Double], Array[Array[Double]], Array[Double]) {
  let n_frames = frames.length()
  if n_frames < 2 {
    return ([], [], [])
  }
  let n_atoms = frames[0].length()
  if n_atoms == 0 {
    return ([], [], [])
  }
  // Compute mean position per atom
  let mean : Array[Array[Double]] = Array::make(n_atoms, [0.0, 0.0, 0.0])
  for f in frames {
    for i = 0; i < n_atoms; i = i + 1 {
      mean[i][0] = mean[i][0] + f[i][0]
      mean[i][1] = mean[i][1] + f[i][1]
      mean[i][2] = mean[i][2] + f[i][2]
    }
  }
  for i = 0; i < n_atoms; i = i + 1 {
    mean[i][0] = mean[i][0] / n_frames.to_double()
    mean[i][1] = mean[i][1] / n_frames.to_double()
    mean[i][2] = mean[i][2] / n_frames.to_double()
  }
  // Build N×N covariance matrix: C[i][j] = <Δr_i · Δr_j>
  let cov : Array[Array[Double]] = Array::makei(n_atoms, fn(_) {
    Array::make(n_atoms, 0.0)
  })
  for f in frames {
    for i = 0; i < n_atoms; i = i + 1 {
      let dxi = f[i][0] - mean[i][0]
      let dyi = f[i][1] - mean[i][1]
      let dzi = f[i][2] - mean[i][2]
      for j = 0; j < n_atoms; j = j + 1 {
        let dxj = f[j][0] - mean[j][0]
        let dyj = f[j][1] - mean[j][1]
        let dzj = f[j][2] - mean[j][2]
        cov[i][j] = cov[i][j] + dxi * dxj + dyi * dyj + dzi * dzj
      }
    }
  }
  let denom = (n_frames - 1).to_double()
  for i = 0; i < n_atoms; i = i + 1 {
    for j = 0; j < n_atoms; j = j + 1 {
      cov[i][j] = cov[i][j] / denom
    }
  }
  // Eigendecomposition (symmetric → Jacobi)
  let (evals, evects) = eigh_symmetric(cov, n_atoms)
  // Sort descending
  let indices : Array[Int] = []
  for i = 0; i < n_atoms; i = i + 1 {
    indices.push(i)
  }
  indices.sort_by(fn(a, b) {
    if evals[b] > evals[a] {
      1
    } else if evals[b] < evals[a] {
      -1
    } else {
      0
    }
  })
  let sorted_evals : Array[Double] = Array(capacity=n_atoms)
  let sorted_evects : Array[Array[Double]] = Array(capacity=n_atoms)
  for idx in indices {
    sorted_evals.push(evals[idx])
    sorted_evects.push(evects[idx])
  }
  // Cumulative variance fraction
  let total = {
    let mut s = 0.0
    for v in sorted_evals {
      if v > 0.0 {
        s = s + v
      }
    }
    s
  }
  let var_frac : Array[Double] = Array(capacity=n_atoms)
  let mut cumul = 0.0
  for v in sorted_evals {
    if v > 0.0 {
      cumul = cumul + v
    }
    var_frac.push(if total > 0.0 { cumul / total } else { 0.0 })
  }
  (sorted_evals, sorted_evects, var_frac)
}

// ============================================================================
// GNM (MDAnalysis.analysis.gnm.GNMAnalysis)
// ============================================================================

///|
/// Build the Kirchhoff (connectivity) matrix for a Gaussian Network Model.
/// `positions` are CA atom coordinates; atoms within `cutoff` Å are connected.
/// The Kirchhoff matrix Γ has:
/// - Γ[i][j] = -1  if i≠j and atoms i,j are within cutoff
/// - Γ[i][i] = number of neighbors of atom i
/// - Γ[i][j] = 0   otherwise
pub fn mdanalysis_gnm_kirchhoff(
  positions : Array[Array[Double]],
  cutoff : Double,
) -> Array[Array[Double]] {
  let n = positions.length()
  let matrix : Array[Array[Double]] = Array::makei(n, fn(_) { Array::make(n, 0.0) })
  let cutoff_sq = cutoff * cutoff
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if mda_sq_dist(positions[i], positions[j]) < cutoff_sq {
        matrix[i][j] = -1.0
        matrix[j][i] = -1.0
        matrix[i][i] = matrix[i][i] + 1.0
        matrix[j][j] = matrix[j][j] + 1.0
      }
    }
  }
  matrix
}

///|
/// Gaussian Network Model analysis.
/// Builds the Kirchhoff matrix from `positions` with `cutoff`, then
/// performs eigendecomposition.  The smallest eigenvalue is ~0 (rigid-body
/// translation); the remaining modes describe internal fluctuations.
///
/// Returns `(eigenvalues, eigenvectors)` sorted ascending.
/// B-factors can be estimated as B_i ∝ Σ_k (1/λ_k) * v_k[i]²  (k > 0).
pub fn mdanalysis_gnm(
  positions : Array[Array[Double]],
  cutoff : Double,
) -> (Array[Double], Array[Array[Double]]) {
  let n = positions.length()
  if n < 2 {
    return ([], [])
  }
  let matrix = mdanalysis_gnm_kirchhoff(positions, cutoff)
  let (evals, evects) = eigh_symmetric(matrix, n)
  // Sort ascending (smallest mode first — mode 0 ≈ 0)
  let indices : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    indices.push(i)
  }
  indices.sort_by(fn(a, b) {
    if evals[a] < evals[b] {
      -1
    } else if evals[a] > evals[b] {
      1
    } else {
      0
    }
  })
  let sorted_evals : Array[Double] = Array(capacity=n)
  let sorted_evects : Array[Array[Double]] = Array(capacity=n)
  for idx in indices {
    sorted_evals.push(evals[idx])
    sorted_evects.push(evects[idx])
  }
  (sorted_evals, sorted_evects)
}

///|
/// Estimate B-factors from GNM eigenvalues/eigenvectors.
/// B_i ∝ Σ_k (1/λ_k) * v_k[i]²  (sum over non-zero modes k > 0).
/// Returns an array of N B-factor values (arbitrary units).
pub fn mdanalysis_gnm_bfactors(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
) -> Array[Double] {
  let n_modes = eigenvalues.length()
  if n_modes == 0 {
    return []
  }
  let n_atoms = eigenvectors[0].length()
  let bfactors : Array[Double] = Array::make(n_atoms, 0.0)
  for k = 1; k < n_modes; k = k + 1 {
    let lambda = eigenvalues[k]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    let inv_lambda = 1.0 / lambda
    let vk = eigenvectors[k]
    for i = 0; i < n_atoms; i = i + 1 {
      bfactors[i] = bfactors[i] + inv_lambda * vk[i] * vk[i]
    }
  }
  bfactors
}

// ============================================================================
// Helix Analysis (MDAnalysis.analysis.helix_analysis.helix_analysis)
// ============================================================================

///|
/// Helix analysis result.
/// All `local_*` arrays have length N-3 (sliding window of 4 atoms),
/// except `local_bends` which has length N-6.
pub struct MdanalysisHelixResult {
  /// Twist angle (degrees) per local window
  local_twists : Array[Double]
  /// Residues per turn per local window
  local_nres_per_turn : Array[Double]
  /// Helix axis unit vector per local window
  local_axes : Array[Array[Double]]
  /// Bending angle between axes 3 windows apart (degrees)
  local_bends : Array[Double]
  /// Rise per residue per local window (Å)
  local_heights : Array[Double]
  /// Radius of local cylinder per window (Å)
  local_radii : Array[Double]
}

///|
/// Calculate helix properties from atomic coordinates.
/// Mirrors `MDAnalysis.analysis.helix_analysis.helix_analysis`.
/// Each property is calculated from a sliding window of 4 atoms.
/// `positions` should be ordered backbone atoms (e.g. CA) of a helix.
pub fn mdanalysis_helix_analysis(
  positions : Array[Array[Double]],
) -> MdanalysisHelixResult {
  let n = positions.length()
  let empty = MdanalysisHelixResult::{
    local_twists: [],
    local_nres_per_turn: [],
    local_axes: [],
    local_bends: [],
    local_heights: [],
    local_radii: [],
  }
  if n < 4 {
    return empty
  }
  // vectors: positions[1:] - positions[:-1]  → (N-1, 3)
  let vectors : Array[Array[Double]] = Array(capacity=n - 1)
  for i = 0; i < n - 1; i = i + 1 {
    vectors.push([
      positions[i + 1][0] - positions[i][0],
      positions[i + 1][1] - positions[i][1],
      positions[i + 1][2] - positions[i][2],
    ])
  }
  // bisectors: vectors[:-1] - vectors[1:]  → (N-2, 3)
  let bisectors : Array[Array[Double]] = Array(capacity=n - 2)
  for i = 0; i < n - 2; i = i + 1 {
    bisectors.push([
      vectors[i][0] - vectors[i + 1][0],
      vectors[i][1] - vectors[i + 1][1],
      vectors[i][2] - vectors[i + 1][2],
    ])
  }
  // bimags: norm of bisectors → (N-2,)
  let bimags : Array[Double] = Array(capacity=n - 2)
  for b in bisectors {
    bimags.push((b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt())
  }
  let n_windows = n - 3 // number of 4-atom windows
  // cos_theta between adjacent bisectors → (N-3,)
  let cos_theta : Array[Double] = Array(capacity=n_windows)
  let adjacent_mag : Array[Double] = Array(capacity=n_windows)
  for i = 0; i < n_windows; i = i + 1 {
    let bi = bisectors[i]
    let bi1 = bisectors[i + 1]
    let dot = bi[0] * bi1[0] + bi[1] * bi1[1] + bi[2] * bi1[2]
    let mag = bimags[i] * bimags[i + 1]
    adjacent_mag.push(mag)
    let ct = if mag > 0.0 { dot / mag } else { 1.0 }
    // clip to [-1, 1]
    cos_theta.push(if ct > 1.0 { 1.0 } else if ct < -1.0 { -1.0 } else { ct })
  }
  // twists (radians) → local_twists (degrees)
  let local_twists : Array[Double] = Array(capacity=n_windows)
  let local_nres_per_turn : Array[Double] = Array(capacity=n_windows)
  let twists_rad : Array[Double] = Array(capacity=n_windows)
  for i = 0; i < n_windows; i = i + 1 {
    let twist = @math.acos(cos_theta[i])
    twists_rad.push(twist)
    local_twists.push(twist * 180.0 / @math.pi)
    local_nres_per_turn.push(if twist > 0.0 { 2.0 * @math.pi / twist } else { 0.0 })
  }
  // local_axes: cross product of adjacent bisectors, normalized → (N-3, 3)
  let local_axes : Array[Array[Double]] = Array(capacity=n_windows)
  for i = 0; i < n_windows; i = i + 1 {
    let bi = bisectors[i]
    let bi1 = bisectors[i + 1]
    let cx = bi[1] * bi1[2] - bi[2] * bi1[1]
    let cy = bi[2] * bi1[0] - bi[0] * bi1[2]
    let cz = bi[0] * bi1[1] - bi[1] * bi1[0]
    let mag = (cx * cx + cy * cy + cz * cz).sqrt()
    if mag > 0.0 {
      local_axes.push([cx / mag, cy / mag, cz / mag])
    } else {
      local_axes.push([0.0, 0.0, 0.0])
    }
  }
  // local_bends: angles between axes 3 windows apart → (N-6,)
  let local_bends : Array[Double] = Array(capacity=0)
  for i = 0; i + 3 < n_windows; i = i + 1 {
    let a = local_axes[i]
    let b = local_axes[i + 3]
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    let clipped = if dot > 1.0 { 1.0 } else if dot < -1.0 { -1.0 } else { dot }
    local_bends.push(@math.acos(clipped) * 180.0 / @math.pi)
  }
  // radii: (adjacent_mag^0.5) / (2 * (1 - cos_theta)) → (N-3,)
  let local_radii : Array[Double] = Array(capacity=n_windows)
  for i = 0; i < n_windows; i = i + 1 {
    let am = adjacent_mag[i]
    let ct = cos_theta[i]
    if (1.0 - ct).abs() > 1.0e-10 {
      local_radii.push(am.sqrt() / (2.0 * (1.0 - ct)))
    } else {
      local_radii.push(am.sqrt() / 2.0)
    }
  }
  // heights: |dot(vectors[1:-1], local_axes)| → (N-3,)
  let local_heights : Array[Double] = Array(capacity=n_windows)
  for i = 0; i < n_windows; i = i + 1 {
    let v = vectors[i + 1]
    let ax = local_axes[i]
    let h = v[0] * ax[0] + v[1] * ax[1] + v[2] * ax[2]
    local_heights.push(h.abs())
  }
  MdanalysisHelixResult::{
    local_twists,
    local_nres_per_turn,
    local_axes,
    local_bends,
    local_heights,
    local_radii,
  }
}

// ============================================================================
// Nucleic Acid Backbone Torsions (MDAnalysis.analysis.nuclinfo)
// ============================================================================

///|
/// Compute a dihedral angle modulo 360 (matching MDAnalysis convention).
fn mda_dihedral_mod360(
  a1 : Array[Double],
  a2 : Array[Double],
  a3 : Array[Double],
  a4 : Array[Double],
) -> Double {
  let d = chem_dihedral_angle(a1, a2, a3, a4)
  let r = d % 360.0
  if r < 0.0 { r + 360.0 } else { r }
}

///|
/// Nucleic acid backbone torsion alpha: O3'(i-1) - P(i) - O5'(i) - C5'(i)
pub fn mdanalysis_tors_alpha(
  o3p_prev : Array[Double],
  p : Array[Double],
  o5p : Array[Double],
  c5p : Array[Double],
) -> Double {
  mda_dihedral_mod360(o3p_prev, p, o5p, c5p)
}

///|
/// Nucleic acid backbone torsion beta: P(i) - O5'(i) - C5'(i) - C4'(i)
pub fn mdanalysis_tors_beta(
  p : Array[Double],
  o5p : Array[Double],
  c5p : Array[Double],
  c4p : Array[Double],
) -> Double {
  mda_dihedral_mod360(p, o5p, c5p, c4p)
}

///|
/// Nucleic acid backbone torsion gamma: O5'(i) - C5'(i) - C4'(i) - C3'(i)
pub fn mdanalysis_tors_gamma(
  o5p : Array[Double],
  c5p : Array[Double],
  c4p : Array[Double],
  c3p : Array[Double],
) -> Double {
  mda_dihedral_mod360(o5p, c5p, c4p, c3p)
}

///|
/// Nucleic acid backbone torsion delta: C5'(i) - C4'(i) - C3'(i) - O3'(i)
pub fn mdanalysis_tors_delta(
  c5p : Array[Double],
  c4p : Array[Double],
  c3p : Array[Double],
  o3p : Array[Double],
) -> Double {
  mda_dihedral_mod360(c5p, c4p, c3p, o3p)
}

///|
/// Nucleic acid backbone torsion epsilon: C4'(i) - C3'(i) - O3'(i) - P(i+1)
pub fn mdanalysis_tors_epsilon(
  c4p : Array[Double],
  c3p : Array[Double],
  o3p : Array[Double],
  p_next : Array[Double],
) -> Double {
  mda_dihedral_mod360(c4p, c3p, o3p, p_next)
}

///|
/// Nucleic acid backbone torsion zeta: C3'(i) - O3'(i) - P(i+1) - O5'(i+1)
pub fn mdanalysis_tors_zeta(
  c3p : Array[Double],
  o3p : Array[Double],
  p_next : Array[Double],
  o5p_next : Array[Double],
) -> Double {
  mda_dihedral_mod360(c3p, o3p, p_next, o5p_next)
}

///|
/// Nucleic acid glycosidic torsion chi (purine): O4'(i) - C1'(i) - N9(i) - C4(i)
pub fn mdanalysis_tors_chi(
  o4p : Array[Double],
  c1p : Array[Double],
  n9 : Array[Double],
  c4 : Array[Double],
) -> Double {
  mda_dihedral_mod360(o4p, c1p, n9, c4)
}

///|
/// Compute all 7 backbone torsions at once.
/// `atoms` maps atom name to 3D coords for residue i,
/// `atoms_prev` for residue i-1 (needs O3'),
/// `atoms_next` for residue i+1 (needs P, O5').
pub fn mdanalysis_tors_all(
  atoms_prev : Map[String, Array[Double]],
  atoms : Map[String, Array[Double]],
  atoms_next : Map[String, Array[Double]],
) -> (Double, Double, Double, Double, Double, Double, Double) {
  let alpha = mda_dihedral_mod360(atoms_prev["O3'"], atoms["P"], atoms["O5'"], atoms["C5'"])
  let beta = mda_dihedral_mod360(atoms["P"], atoms["O5'"], atoms["C5'"], atoms["C4'"])
  let gamma = mda_dihedral_mod360(atoms["O5'"], atoms["C5'"], atoms["C4'"], atoms["C3'"])
  let delta = mda_dihedral_mod360(atoms["C5'"], atoms["C4'"], atoms["C3'"], atoms["O3'"])
  let eps = mda_dihedral_mod360(atoms["C4'"], atoms["C3'"], atoms["O3'"], atoms_next["P"])
  let zeta = mda_dihedral_mod360(atoms["C3'"], atoms["O3'"], atoms_next["P"], atoms_next["O5'"])
  let chi = mda_dihedral_mod360(atoms["O4'"], atoms["C1'"], atoms["N9"], atoms["C4"])
  (alpha, beta, gamma, delta, eps, zeta, chi)
}

// ============================================================================
// Sugar Pucker Phase (MDAnalysis.analysis.nuclinfo.phase_cp / phase_as)
// ============================================================================

///|
/// Sugar pucker phase angle using the CP (Cremer-Pople) method.
/// `sugar_atoms` = [O4', C1', C2', C3', C4'] — the 5 ribose ring atoms.
/// Returns phase angle in degrees [0, 360).
pub fn mdanalysis_phase_cp(
  sugar_atoms : Array[Array[Double]],
) -> Double {
  if sugar_atoms.length() < 5 {
    return 0.0
  }
  let r1 = sugar_atoms[0]
  let r2 = sugar_atoms[1]
  let r3 = sugar_atoms[2]
  let r4 = sugar_atoms[3]
  let r5 = sugar_atoms[4]
  // Centroid
  let r0 = [
    (r1[0] + r2[0] + r3[0] + r4[0] + r5[0]) / 5.0,
    (r1[1] + r2[1] + r3[1] + r4[1] + r5[1]) / 5.0,
    (r1[2] + r2[2] + r3[2] + r4[2] + r5[2]) / 5.0,
  ]
  // Centered coordinates
  let pts : Array[Array[Double]] = []
  for r in [r1, r2, r3, r4, r5] {
    pts.push([r[0] - r0[0], r[1] - r0[1], r[2] - r0[2]])
  }
  // R1 = sum r_j * sin(2pi*j/5), R2 = sum r_j * cos(2pi*j/5)
  let mut r1x = 0.0
  let mut r1y = 0.0
  let mut r1z = 0.0
  let mut r2x = 0.0
  let mut r2y = 0.0
  let mut r2z = 0.0
  for j = 0; j < 5; j = j + 1 {
    let angle = 2.0 * @math.pi * j.to_double() / 5.0
    let s = @math.sin(angle)
    let c = @math.cos(angle)
    r1x = r1x + pts[j][0] * s
    r1y = r1y + pts[j][1] * s
    r1z = r1z + pts[j][2] * s
    r2x = r2x + pts[j][0] * c
    r2y = r2y + pts[j][1] * c
    r2z = r2z + pts[j][2] * c
  }
  // n = (R1 x R2) normalized
  let nx = r1y * r2z - r1z * r2y
  let ny = r1z * r2x - r1x * r2z
  let nz = r1x * r2y - r1y * r2x
  let nmag = (nx * nx + ny * ny + nz * nz).sqrt()
  if nmag < 1.0e-10 {
    return 0.0
  }
  let n = [nx / nmag, ny / nmag, nz / nmag]
  // Project onto n
  let ds : Array[Double] = []
  for p in pts {
    ds.push(p[0] * n[0] + p[1] * n[1] + p[2] * n[2])
  }
  // D = -sqrt(2/5) * sum d_j * sin(4pi*j/5)
  // C =  sqrt(2/5) * sum d_j * cos(4pi*j/5)
  let mut d_val = 0.0
  let mut c_val = 0.0
  for j = 0; j < 5; j = j + 1 {
    let angle = 4.0 * @math.pi * j.to_double() / 5.0
    d_val = d_val + ds[j] * @math.sin(angle)
    c_val = c_val + ds[j] * @math.cos(angle)
  }
  let sqrt_2_5 = (2.0 / 5.0).sqrt()
  d_val = d_val * -1.0 * sqrt_2_5
  c_val = c_val * sqrt_2_5
  let phase = (@math.atan2(d_val, c_val) + @math.pi / 2.0) * 180.0 / @math.pi
  let result = phase % 360.0
  if result < 0.0 { result + 360.0 } else { result }
}

///|
/// Sugar pucker phase angle using the AS (Altona-Sundaralingam) method.
/// `sugar_atoms` = [O4', C1', C2', C3', C4'] — the 5 ribose ring atoms.
/// Computes 5 pseudo-torsion angles from consecutive quadruples,
/// then Fourier-transforms them. Returns phase angle in degrees [0, 360).
pub fn mdanalysis_phase_as(
  sugar_atoms : Array[Array[Double]],
) -> Double {
  if sugar_atoms.length() < 5 {
    return 0.0
  }
  let o4p = sugar_atoms[0]
  let c1p = sugar_atoms[1]
  let c2p = sugar_atoms[2]
  let c3p = sugar_atoms[3]
  let c4p = sugar_atoms[4]
  // 5 pseudo-torsions nu0..nu4
  let nu0 = chem_dihedral_angle(c4p, o4p, c1p, c2p)
  let nu1 = chem_dihedral_angle(o4p, c1p, c2p, c3p)
  let nu2 = chem_dihedral_angle(c1p, c2p, c3p, c4p)
  let nu3 = chem_dihedral_angle(c2p, c3p, c4p, o4p)
  let nu4 = chem_dihedral_angle(c3p, c4p, o4p, c1p)
  let nus = [nu0, nu1, nu2, nu3, nu4]
  // B = -2/5 * sum nu_j * sin(4pi*j/5)
  // A =  2/5 * sum nu_j * cos(4pi*j/5)
  let mut b_val = 0.0
  let mut a_val = 0.0
  for j = 0; j < 5; j = j + 1 {
    let angle = 4.0 * @math.pi * j.to_double() / 5.0
    b_val = b_val + nus[j] * @math.sin(angle)
    a_val = a_val + nus[j] * @math.cos(angle)
  }
  b_val = b_val * -2.0 / 5.0
  a_val = a_val * 2.0 / 5.0
  let phase = @math.atan2(b_val, a_val) * 180.0 / @math.pi
  let result = phase % 360.0
  if result < 0.0 { result + 360.0 } else { result }
}

// ============================================================================
// Leaflet Finder (MDAnalysis.analysis.leaflet.LeafletFinder)
// ============================================================================

///|
/// Identify lipid bilayer leaflets using a distance-based connectivity graph.
/// `head_positions` are the 3D coordinates of lipid head-group atoms.
/// Atoms within `cutoff` Angstroms of each other are considered connected;
/// connected components define the leaflets.
/// Returns an array of leaflet indices (0-based) per atom.
pub fn mdanalysis_leaflet_finder(
  head_positions : Array[Array[Double]],
  cutoff : Double,
) -> Array[Int] {
  let n = head_positions.length()
  if n == 0 {
    return []
  }
  // Union-Find
  let parent : Array[Int] = Array::makei(n, fn(i) { i })
  let find = fn(parent : Array[Int], x : Int) -> Int {
    let mut root = x
    while parent[root] != root {
      root = parent[root]
    }
    let mut curr = x
    while parent[curr] != root {
      let next = parent[curr]
      parent[curr] = root
      curr = next
    }
    root
  }
  let cutoff_sq = cutoff * cutoff
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if mda_sq_dist(head_positions[i], head_positions[j]) < cutoff_sq {
        let ri = find(parent, i)
        let rj = find(parent, j)
        if ri != rj {
          parent[ri] = rj
        }
      }
    }
  }
  // Assign leaflet indices
  let labels : Array[Int] = Array::make(n, -1)
  let mut n_leaflets = 0
  for i = 0; i < n; i = i + 1 {
    let root = find(parent, i)
    if labels[root] == -1 {
      labels[root] = n_leaflets
      n_leaflets = n_leaflets + 1
    }
    labels[i] = labels[root]
  }
  labels
}

// ============================================================================
// BAT: Bond-Angle-Torsion Internal Coordinates
// (MDAnalysis.analysis.bat.BAT)
// ============================================================================

///|
/// BAT root-atom parameters from 3 root atoms.
/// Returns (phi, theta, omega, r01, r12, a012) where:
/// - phi, theta: polar angles of the bond vector p0->p1
/// - omega: rotation angle about the bond axis
/// - r01: distance p0-p1
/// - r12: distance p1-p2
/// - a012: angle p0-p1-p2 (radians)
pub fn mdanalysis_bat_root(
  p0 : Array[Double],
  p1 : Array[Double],
  p2 : Array[Double],
) -> (Double, Double, Double, Double, Double, Double) {
  let v01 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]
  let v21 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]]
  let r01 = (v01[0] * v01[0] + v01[1] * v01[1] + v01[2] * v01[2]).sqrt()
  let r12 = (v21[0] * v21[0] + v21[1] * v21[1] + v21[2] * v21[2]).sqrt()
  let dot = v01[0] * v21[0] + v01[1] * v21[1] + v01[2] * v21[2]
  let cos_a = if r01 > 0.0 && r12 > 0.0 {
    let c = dot / (r01 * r12)
    if c > 1.0 { 1.0 } else if c < -1.0 { -1.0 } else { c }
  } else {
    1.0
  }
  let a012 = @math.acos(cos_a)
  let phi = @math.atan2(v01[1], v01[0])
  let theta = @math.acos(if r01 > 0.0 { v01[2] / r01 } else { 1.0 })
  // Rotation to z axis: Rz . (p2 - p1)
  let cp = @math.cos(phi)
  let sp = @math.sin(phi)
  let ct = @math.cos(theta)
  let st = @math.sin(theta)
  let dp = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]
  let pos2_0 = cp * ct * dp[0] + ct * sp * dp[1] - st * dp[2]
  let pos2_1 = -sp * dp[0] + cp * dp[1]
  let omega = @math.atan2(pos2_1, pos2_0)
  (phi, theta, omega, r01, r12, a012)
}

///|
/// BAT internal coordinates for a chain of atoms.
/// Given 4 arrays of consecutive atom positions (ag1, ag2, ag3, ag4),
/// computes bonds (ag1-ag2), angles (ag1-ag2-ag3), and torsions (ag1-ag2-ag3-ag4).
/// Returns (bonds, angles, torsions) where each array has length = len(ag1).
pub fn mdanalysis_bat_internal(
  ag1 : Array[Array[Double]],
  ag2 : Array[Array[Double]],
  ag3 : Array[Array[Double]],
  ag4 : Array[Array[Double]],
) -> (Array[Double], Array[Double], Array[Double]) {
  let n = ag1.length()
  let bonds : Array[Double] = Array(capacity=n)
  let angles : Array[Double] = Array(capacity=n)
  let torsions : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    // Bond: dist(ag1[i], ag2[i])
    let dx = ag2[i][0] - ag1[i][0]
    let dy = ag2[i][1] - ag1[i][1]
    let dz = ag2[i][2] - ag1[i][2]
    bonds.push((dx * dx + dy * dy + dz * dz).sqrt())
    // Angle: ag1[i]-ag2[i]-ag3[i]
    let v1 = [ag1[i][0] - ag2[i][0], ag1[i][1] - ag2[i][1], ag1[i][2] - ag2[i][2]]
    let v2 = [ag3[i][0] - ag2[i][0], ag3[i][1] - ag2[i][1], ag3[i][2] - ag2[i][2]]
    let m1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt()
    let m2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt()
    let cos_ang = if m1 > 0.0 && m2 > 0.0 {
      let c = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m1 * m2)
      if c > 1.0 { 1.0 } else if c < -1.0 { -1.0 } else { c }
    } else {
      1.0
    }
    angles.push(@math.acos(cos_ang))
    // Torsion: ag1[i]-ag2[i]-ag3[i]-ag4[i]
    torsions.push(chem_dihedral_angle(ag1[i], ag2[i], ag3[i], ag4[i]) * @math.pi / 180.0)
  }
  (bonds, angles, torsions)
}
