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
