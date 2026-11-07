/// rdkit_3d_align.mbt — 3D alignment & shape descriptors (RDKit gap fill, A3)
///
/// GetBestRMS-style optimal-alignment RMSD with symmetry permutation
/// search, a grid-based shape alignment (rdShapeAlign analogue with
/// Gaussian atom densities), mass-weighted principal moments of inertia
/// with NPR ratios, a planarity descriptor (PBF proxy) and USR moments.

///|
/// Mass-weighted principal moments of inertia (smallest, middle, largest)
/// computed from real atomic masses.
pub fn rdkit_pmi_mass_weighted(
  conf : RdkitConformer,
) -> (Double, Double, Double) {
  let n = conf.coords.length()
  if n == 0 {
    return (0.0, 0.0, 0.0)
  }
  // center of mass
  let mut total_mass = 0.0
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for i = 0; i < n; i = i + 1 {
    let m = rdkit_atomic_weight(conf.atoms[i])
    let (x, y, z) = conf.coords[i]
    total_mass = total_mass + m
    cx = cx + m * x
    cy = cy + m * y
    cz = cz + m * z
  }
  if total_mass <= 0.0 {
    return (0.0, 0.0, 0.0)
  }
  cx = cx / total_mass
  cy = cy / total_mass
  cz = cz / total_mass
  // inertia tensor (symmetric 3x3)
  let mut ixx = 0.0
  let mut iyy = 0.0
  let mut izz = 0.0
  let mut ixy = 0.0
  let mut ixz = 0.0
  let mut iyz = 0.0
  for i = 0; i < n; i = i + 1 {
    let m = rdkit_atomic_weight(conf.atoms[i])
    let (x, y, z) = conf.coords[i]
    let dx = x - cx
    let dy = y - cy
    let dz = z - cz
    ixx = ixx + m * (dy * dy + dz * dz)
    iyy = iyy + m * (dx * dx + dz * dz)
    izz = izz + m * (dx * dx + dy * dy)
    ixy = ixy - m * dx * dy
    ixz = ixz - m * dx * dz
    iyz = iyz - m * dy * dz
  }
  let tensor : Array[Array[Double]] = [
    [ixx, ixy, ixz],
    [ixy, iyy, iyz],
    [ixz, iyz, izz],
  ]
  let eigs = rdkit_desc_jacobi_eigenvalues(tensor)
  // sort ascending
  for i = 0; i < 3; i = i + 1 {
    for j = i + 1; j < 3; j = j + 1 {
      if eigs[j] < eigs[i] {
        let tmp = eigs[i]
        eigs[i] = eigs[j]
        eigs[j] = tmp
      }
    }
  }
  (eigs[0], eigs[1], eigs[2])
}

///|
/// Normalized principal moment ratios: npr1 = I1/I3, npr2 = I2/I3,
/// npr3 = 1 (returned for completeness).
pub fn rdkit_npr_ratios(conf : RdkitConformer) -> (Double, Double, Double) {
  let (i1, i2, i3) = rdkit_pmi_mass_weighted(conf)
  if i3 <= 0.0000000001 {
    return (0.0, 0.0, 0.0)
  }
  (i1 / i3, i2 / i3, 1.0)
}

///|
/// Planarity descriptor (PBF proxy): 3 * lambda_min / trace of the
/// (unweighted) gyration tensor eigenvalues. ~0 for planar point sets,
/// ~1/3 for a sphere-like distribution.
pub fn rdkit_pbf(conf : RdkitConformer) -> Double {
  let n = conf.coords.length()
  if n < 3 {
    return 0.0
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for i = 0; i < n; i = i + 1 {
    let (x, y, z) = conf.coords[i]
    cx = cx + x
    cy = cy + y
    cz = cz + z
  }
  cx = cx / n.to_double()
  cy = cy / n.to_double()
  cz = cz / n.to_double()
  let mut sxx = 0.0
  let mut syy = 0.0
  let mut szz = 0.0
  let mut sxy = 0.0
  let mut sxz = 0.0
  let mut syz = 0.0
  for i = 0; i < n; i = i + 1 {
    let (x, y, z) = conf.coords[i]
    let dx = x - cx
    let dy = y - cy
    let dz = z - cz
    sxx = sxx + dx * dx
    syy = syy + dy * dy
    szz = szz + dz * dz
    sxy = sxy + dx * dy
    sxz = sxz + dx * dz
    syz = syz + dy * dz
  }
  let nf = n.to_double()
  let tensor : Array[Array[Double]] = [
    [sxx / nf, sxy / nf, sxz / nf],
    [sxy / nf, syy / nf, syz / nf],
    [sxz / nf, syz / nf, szz / nf],
  ]
  let eigs = rdkit_desc_jacobi_eigenvalues(tensor)
  let trace = eigs[0] + eigs[1] + eigs[2]
  if trace <= 0.0000000001 {
    return 0.0
  }
  let lambda_min = if eigs[0] < eigs[1] {
    if eigs[0] < eigs[2] { eigs[0] } else { eigs[2] }
  } else if eigs[1] < eigs[2] {
    eigs[1]
  } else {
    eigs[2]
  }
  3.0 * lambda_min / trace
}

///|
/// USR (Ultrafast Shape Recognition) moments: 4 reference points
/// (centroid, farthest from centroid, farthest from that, farthest from
/// the midpoint pair), each described by (mean, variance, skewness) of
/// the distance distribution -> 12 dimensions.
pub fn rdkit_usr_descriptors(conf : RdkitConformer) -> Array[Double] {
  let n = conf.coords.length()
  let out : Array[Double] = Array::make(12, 0.0)
  if n == 0 {
    return out
  }
  // centroid
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for i = 0; i < n; i = i + 1 {
    let (x, y, z) = conf.coords[i]
    cx = cx + x
    cy = cy + y
    cz = cz + z
  }
  cx = cx / n.to_double()
  cy = cy / n.to_double()
  cz = cz / n.to_double()
  let centroid = (cx, cy, cz)
  let ctd = usr_farthest(conf, centroid)
  let cst = usr_farthest(conf, ctd)
  let midx = (ctd.0 + cst.0) / 2.0
  let midy = (ctd.1 + cst.1) / 2.0
  let midz = (ctd.2 + cst.2) / 2.0
  let fct = usr_farthest(conf, (midx, midy, midz))
  let refs = [centroid, ctd, cst, fct]
  for r = 0; r < 4; r = r + 1 {
    // distance distribution to reference point r
    let dists : Array[Double] = Array(capacity=0)
    let (rx, ry, rz) = refs[r]
    for i = 0; i < n; i = i + 1 {
      let (x, y, z) = conf.coords[i]
      let dx = x - rx
      let dy = y - ry
      let dz = z - rz
      dists.push(rdkit_sqrt(dx * dx + dy * dy + dz * dz))
    }
    let mut mean = 0.0
    for d in dists {
      mean = mean + d
    }
    mean = mean / n.to_double()
    let mut m2 = 0.0
    let mut m3 = 0.0
    for d in dists {
      let dd = d - mean
      m2 = m2 + dd * dd
      m3 = m3 + dd * dd * dd
    }
    m2 = m2 / n.to_double()
    m3 = m3 / n.to_double()
    out[r * 3] = mean
    out[r * 3 + 1] = m2
    out[r * 3 + 2] = m3
  }
  out
}

///|
fn usr_farthest(
  conf : RdkitConformer,
  point : (Double, Double, Double),
) -> (Double, Double, Double) {
  let n = conf.coords.length()
  let mut best = 0
  let mut best_d = -1.0
  for i = 0; i < n; i = i + 1 {
    let (x, y, z) = conf.coords[i]
    let dx = x - point.0
    let dy = y - point.1
    let dz = z - point.2
    let d = dx * dx + dy * dy + dz * dz
    if d > best_d {
      best_d = d
      best = i
    }
  }
  conf.coords[best]
}

///|
/// Euclidean distance between two USR descriptor vectors (the USR shape
/// similarity is 1 / (1 + distance)).
pub fn rdkit_usr_distance(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let d = a[i] - b[i]
    sum = sum + d * d
  }
  rdkit_sqrt(sum)
}

///|
/// Optimal-alignment RMSD with symmetry-aware atom permutation search:
/// identical-element atom groups are permuted (capped at 5040 total
/// permutations) and the rigid-body superposition RMSD is evaluated for
/// each mapping; the minimum is returned together with the transformed
/// probe coordinates. Returns None when atom counts or element
/// compositions differ.
pub fn rdkit_get_best_rms(
  probe : RdkitConformer,
  reference : RdkitConformer,
) -> Option[(Double, RdkitConformer)] {
  let n = probe.coords.length()
  if n != reference.coords.length() || n == 0 {
    return None
  }
  // element-wise quick check: both conformers come from the same molecule
  // so the symbol sequences should line up index by index
  let symbols : Array[String] = Array(capacity=0)
  for s in reference.atoms {
    symbols.push(s)
  }
  for i = 0; i < n; i = i + 1 {
    let sym = if i < probe.atoms.length() { probe.atoms[i] } else { "C" }
    if sym != symbols[i] {
      return None
    }
  }
  // rebuild multiset check separately (index identity is too strict for
  // symmetric molecules)
  let probe_counts : Map[String, Int32] = Map([], capacity=0)
  for s in probe.atoms {
    probe_counts[s] = match probe_counts.get(s) {
      Some(c) => c + 1
      None => 1
    }
  }
  let ref_counts : Map[String, Int32] = Map([], capacity=0)
  for s in reference.atoms {
    ref_counts[s] = match ref_counts.get(s) {
      Some(c) => c + 1
      None => 1
    }
  }
  for s, c in probe_counts {
    match ref_counts.get(s) {
      Some(c2) => if c != c2 {
        return None
      }
      None => return None
    }
  }
  // permutation groups per element (positions of that element in probe)
  let by_symbol : Map[String, Array[Int32]] = Map([], capacity=0)
  for i = 0; i < n; i = i + 1 {
    let sym = probe.atoms[i]
    match by_symbol.get(sym) {
      Some(arr) => arr.push(i)
      None => by_symbol[sym] = [i]
    }
  }
  // count total permutations
  let mut total_perms = 1
  for s, idxs in by_symbol {
    let k = idxs.length()
    let mut fact = 1
    for i = 2; i <= k; i = i + 1 {
      fact = fact * i
    }
    total_perms = total_perms * fact
    let _ = s
  }
  let mut best_rms = 1.0e300
  let mut best_coords : Array[(Double, Double, Double)] = probe.coords
  if total_perms <= 5040 {
    let result = rms_permute(probe, reference, by_symbol, n, 1.0e300)
    best_rms = result.0
    best_coords = result.1
  } else {
    // identity mapping fallback
    best_rms = aligned_rms(probe.coords, reference.coords)
    best_coords = transform_coords(probe.coords, reference.coords)
  }
  Some((
    best_rms,
    RdkitConformer::{ atoms: probe.atoms, coords: best_coords, energy: probe.energy, },
  ))
}

///|
fn rms_permute(
  probe : RdkitConformer,
  reference : RdkitConformer,
  by_symbol : Map[String, Array[Int32]],
  n : Int32,
  best_so_far : Double,
) -> (Double, Array[(Double, Double, Double)]) {
  // gather symbols deterministically
  let symbols : Array[String] = Array(capacity=0)
  for s, _idxs in by_symbol {
    symbols.push(s)
  }
  let mut best_rms = best_so_far
  let mut best_coords = probe.coords
  // generate cartesian product of permutations for each symbol group
  let per_group : Array[Array[Array[Int32]]] = Array(capacity=0)
  for s in symbols {
    per_group.push(permute_array(by_symbol[s]))
  }
  // iterate over cartesian product via odometer
  let counters : Array[Int32] = Array::make(symbols.length(), 0)
  while true {
    // build full mapping: probe position -> probe index to use
    let mapping : Array[Int32] = Array::make(n, 0)
    for gi = 0; gi < symbols.length(); gi = gi + 1 {
      let perm = per_group[gi][counters[gi]]
      let idxs = by_symbol[symbols[gi]]
      for k = 0; k < idxs.length(); k = k + 1 {
        mapping[idxs[k]] = perm[k]
      }
    }
    // apply mapping: reordered probe coordinates
    let mapped : Array[(Double, Double, Double)] = Array(capacity=0)
    for i = 0; i < n; i = i + 1 {
      mapped.push(probe.coords[mapping[i]])
    }
    let rms = aligned_rms(mapped, reference.coords)
    if rms < best_rms {
      best_rms = rms
      best_coords = transform_coords(mapped, reference.coords)
    }
    // odometer increment
    let mut carry = symbols.length() - 1
    let mut overflow = true
    while carry >= 0 {
      counters[carry] = counters[carry] + 1
      if counters[carry] < per_group[carry].length() {
        overflow = false
        break
      }
      counters[carry] = 0
      carry = carry - 1
    }
    if overflow {
      break
    }
  }
  (best_rms, best_coords)
}

///|
/// Generate all permutations of the input array (n! entries).
fn permute_array(items : Array[Int32]) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = Array(capacity=0)
  permute_rec(items, [], out)
  out
}

///|
fn permute_rec(
  remaining : Array[Int32],
  current : Array[Int32],
  out : Array[Array[Int32]],
) -> Unit {
  if remaining.length() == 0 {
    out.push(current.copy())
    return
  }
  for i = 0; i < remaining.length(); i = i + 1 {
    let rest : Array[Int32] = Array(capacity=0)
    for j = 0; j < remaining.length(); j = j + 1 {
      if j != i {
        rest.push(remaining[j])
      }
    }
    let next = current.copy()
    next.push(remaining[i])
    permute_rec(rest, next, out)
  }
}

///|
/// Kabsch superposition (SVD-style, numerically stable): the rotation is
/// rebuilt from the eigendecomposition of M^T M where M = sum q' p'^T.
/// Small singular directions are completed with cross products instead
/// of inverting near-zero eigenvalues, so coplanar point sets do not
/// blow up. Returns (rotation 3x3, translation, rmsd) mapping probe
/// onto reference.
fn kabsch_fit(
  probe : Array[(Double, Double, Double)],
  reference : Array[(Double, Double, Double)],
) -> (Array[Array[Double]], (Double, Double, Double), Double) {
  let n = probe.length()
  let identity : Array[Array[Double]] = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
  ]
  if n == 0 || n != reference.length() {
    return (identity, (0.0, 0.0, 0.0), 0.0)
  }
  // centroids
  let mut cpx = 0.0
  let mut cpy = 0.0
  let mut cpz = 0.0
  let mut cqx = 0.0
  let mut cqy = 0.0
  let mut cqz = 0.0
  for i = 0; i < n; i = i + 1 {
    let (px, py, pz) = probe[i]
    let (qx, qy, qz) = reference[i]
    cpx = cpx + px
    cpy = cpy + py
    cpz = cpz + pz
    cqx = cqx + qx
    cqy = cqy + qy
    cqz = cqz + qz
  }
  let nf = n.to_double()
  cpx = cpx / nf
  cpy = cpy / nf
  cpz = cpz / nf
  cqx = cqx / nf
  cqy = cqy / nf
  cqz = cqz / nf
  // cross-covariance M[i][j] = sum_k Q'_k[i] * P'_k[j]
  let m : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for i = 0; i < n; i = i + 1 {
    let (px, py, pz) = probe[i]
    let (qx, qy, qz) = reference[i]
    let pd = [px - cpx, py - cpy, pz - cpz]
    let qd = [qx - cqx, qy - cqy, qz - cqz]
    for a = 0; a < 3; a = a + 1 {
      for b = 0; b < 3; b = b + 1 {
        let _ = m[a].set(b, m[a][b] + qd[a] * pd[b])
      }
    }
  }
  // N = M^T M, eigendecomposed with eigenvectors
  let nmat : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for a = 0; a < 3; a = a + 1 {
    for b = 0; b < 3; b = b + 1 {
      let mut s = 0.0
      for l = 0; l < 3; l = l + 1 {
        s = s + m[l][a] * m[l][b]
      }
      let _ = nmat[a].set(b, s)
    }
  }
  let (vals, vecs) = rdkit_desc_jacobi_full(nmat)
  // sort eigenpairs descending by eigenvalue
  for i = 0; i < 3; i = i + 1 {
    for j = i + 1; j < 3; j = j + 1 {
      if vals[j] > vals[i] {
        let tv = vals[i]
        vals[i] = vals[j]
        vals[j] = tv
        for k = 0; k < 3; k = k + 1 {
          let tk = vecs[k][i]
          vecs[k][i] = vecs[k][j]
          vecs[k][j] = tk
        }
      }
    }
  }
  let tol = 0.00000001 * rdkit_sqrt(rdkit_desc_max0(vals[0]))
  // singular values
  let s1 = rdkit_sqrt(rdkit_desc_max0(vals[0]))
  let s2 = rdkit_sqrt(rdkit_desc_max0(vals[1]))
  // build U columns: u_k = M v_k / sigma_k when sigma is large enough
  let mut u1 = [0.0, 0.0, 0.0]
  let mut u2 = [0.0, 0.0, 0.0]
  if s1 > tol {
    for a = 0; a < 3; a = a + 1 {
      u1[a] = (m[a][0] * vecs[0][0] + m[a][1] * vecs[1][0] + m[a][2] * vecs[2][0]) / s1
    }
  } else {
    u1 = [1.0, 0.0, 0.0]
  }
  if s2 > tol {
    for a = 0; a < 3; a = a + 1 {
      u2[a] = (m[a][0] * vecs[0][1] + m[a][1] * vecs[1][1] + m[a][2] * vecs[2][1]) / s2
    }
  } else {
    // orthogonal complement of u1
    u2 = ff_orthogonal_complement(u1)
  }
  // re-orthonormalize u2 against u1 (numerical hygiene)
  let mut dot12 = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2]
  let u1n = rdkit_sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2])
  if u1n > 0.0 {
    for a = 0; a < 3; a = a + 1 {
      let _ = u1.set(a, u1[a] / u1n)
    }
    dot12 = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2]
    for a = 0; a < 3; a = a + 1 {
      let _ = u2.set(a, u2[a] - dot12 * u1[a])
    }
  }
  let u2n = rdkit_sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2])
  if u2n > 0.0 {
    for a = 0; a < 3; a = a + 1 {
      let _ = u2.set(a, u2[a] / u2n)
    }
  }
  // u3 = u1 x u2 completes the right-handed frame
  let u3 = [
    u1[1] * u2[2] - u1[2] * u2[1],
    u1[2] * u2[0] - u1[0] * u2[2],
    u1[0] * u2[1] - u1[1] * u2[0],
  ]
  // V columns are vecs[:, k]; assemble R = U * D * V^T with
  // D = diag(1, 1, sign(det(U V^T))) to force det(R) = +1
  let uv00 = u1[0] * vecs[0][0] + u2[0] * vecs[0][1] + u3[0] * vecs[0][2]
  let uv01 = u1[1] * vecs[0][0] + u2[1] * vecs[0][1] + u3[1] * vecs[0][2]
  let uv02 = u1[2] * vecs[0][0] + u2[2] * vecs[0][1] + u3[2] * vecs[0][2]
  let uv10 = u1[0] * vecs[1][0] + u2[0] * vecs[1][1] + u3[0] * vecs[1][2]
  let uv11 = u1[1] * vecs[1][0] + u2[1] * vecs[1][1] + u3[1] * vecs[1][2]
  let uv12 = u1[2] * vecs[1][0] + u2[2] * vecs[1][1] + u3[2] * vecs[1][2]
  let uv20 = u1[0] * vecs[2][0] + u2[0] * vecs[2][1] + u3[0] * vecs[2][2]
  let uv21 = u1[1] * vecs[2][0] + u2[1] * vecs[2][1] + u3[1] * vecs[2][2]
  let uv22 = u1[2] * vecs[2][0] + u2[2] * vecs[2][1] + u3[2] * vecs[2][2]
  // det of U * V^T
  let det_uvt = uv00 * (uv11 * uv22 - uv12 * uv21) -
    uv01 * (uv10 * uv22 - uv12 * uv20) +
    uv02 * (uv10 * uv21 - uv11 * uv20)
  let d3 = if det_uvt < 0.0 { -1.0 } else { 1.0 }
  // R = U * diag(1,1,d3) * V^T
  // (U D)[a][k] = U[a][k] * (k==2 ? d3 : 1)
  let rot : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for a = 0; a < 3; a = a + 1 {
    let ua = [u1[a], u2[a], u3[a]]
    for b = 0; b < 3; b = b + 1 {
      let mut s = 0.0
      for k = 0; k < 3; k = k + 1 {
        let uk = if k == 2 { ua[k] * d3 } else { ua[k] }
        s = s + uk * vecs[b][k]
      }
      let _ = rot[a].set(b, s)
    }
  }
  // rmsd
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let (px, py, pz) = probe[i]
    let (qx, qy, qz) = reference[i]
    let tx = rot[0][0] * (px - cpx) + rot[0][1] * (py - cpy) + rot[0][2] * (pz - cpz)
    let ty = rot[1][0] * (px - cpx) + rot[1][1] * (py - cpy) + rot[1][2] * (pz - cpz)
    let tz = rot[2][0] * (px - cpx) + rot[2][1] * (py - cpy) + rot[2][2] * (pz - cpz)
    let dx = cqx + tx - qx
    let dy = cqy + ty - qy
    let dz = cqz + tz - qz
    sum = sum + dx * dx + dy * dy + dz * dz
  }
  (rot, (cqx, cqy, cqz), rdkit_sqrt(sum / nf))
}

///|
fn rdkit_desc_max0(x : Double) -> Double {
  if x < 0.0 { 0.0 } else { x }
}

///|
/// Arbitrary unit vector orthogonal to v.
fn ff_orthogonal_complement(v : Array[Double]) -> Array[Double] {
  // pick the coordinate axis least aligned with v
  let ax = rdkit_desc_abs(v[0])
  let ay = rdkit_desc_abs(v[1])
  let az = rdkit_desc_abs(v[2])
  let other = if ax <= ay && ax <= az {
    [1.0, 0.0, 0.0]
  } else if ay <= az {
    [0.0, 1.0, 0.0]
  } else {
    [0.0, 0.0, 1.0]
  }
  // cross product v x other
  let c = [
    v[1] * other[2] - v[2] * other[1],
    v[2] * other[0] - v[0] * other[2],
    v[0] * other[1] - v[1] * other[0],
  ]
  let n = rdkit_sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
  if n < 0.000000000001 {
    [0.0, 1.0, 0.0]
  } else {
    [c[0] / n, c[1] / n, c[2] / n]
  }
}

///|
/// Kabsch RMSD between equal-length coordinate sets.
fn aligned_rms(
  probe : Array[(Double, Double, Double)],
  reference : Array[(Double, Double, Double)],
) -> Double {
  let (_, _, rmsd) = kabsch_fit(probe, reference)
  rmsd
}

///|
/// Rigid-body transform of probe onto reference; returns transformed
/// probe coordinates.
fn transform_coords(
  probe : Array[(Double, Double, Double)],
  reference : Array[(Double, Double, Double)],
) -> Array[(Double, Double, Double)] {
  let (rot, _t, _rmsd) = kabsch_fit(probe, reference)
  // recompute both centroids here: the transform is R*(p - cp) + cq
  let n = probe.length()
  let mut cpx = 0.0
  let mut cpy = 0.0
  let mut cpz = 0.0
  let mut cqx = 0.0
  let mut cqy = 0.0
  let mut cqz = 0.0
  for i = 0; i < n; i = i + 1 {
    let (px, py, pz) = probe[i]
    let (qx, qy, qz) = reference[i]
    cpx = cpx + px
    cpy = cpy + py
    cpz = cpz + pz
    cqx = cqx + qx
    cqy = cqy + qy
    cqz = cqz + qz
  }
  let nf = n.to_double()
  cpx = cpx / nf
  cpy = cpy / nf
  cpz = cpz / nf
  cqx = cqx / nf
  cqy = cqy / nf
  cqz = cqz / nf
  let out : Array[(Double, Double, Double)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let (px, py, pz) = probe[i]
    let dx = px - cpx
    let dy = py - cpy
    let dz = pz - cpz
    out.push((
      rot[0][0] * dx + rot[0][1] * dy + rot[0][2] * dz + cqx,
      rot[1][0] * dx + rot[1][1] * dy + rot[1][2] * dz + cqy,
      rot[2][0] * dx + rot[2][1] * dy + rot[2][2] * dz + cqz,
    ))
  }
  out
}

///|
/// Grid-based shape alignment (rdShapeAlign analogue): the probe is
/// optimally superposed onto the reference, then both are sampled on a
/// 0.5 A grid with Gaussian atom densities. Returns the shape Tanimoto
/// and the volume ratio (probe occupied / reference occupied).
pub struct RdkitShapeAlignResult {
  shape_tanimoto : Double
  protrude : Double
  grid_points : Int32
}

///|
pub fn rdkit_shape_align(
  probe : RdkitConformer,
  reference : RdkitConformer,
) -> RdkitShapeAlignResult {
  // superpose probe onto reference first
  let transformed = transform_coords(probe.coords, reference.coords)
  // bounding box around the reference (+/- 3 A padding)
  let mut minx = 1.0e300
  let mut miny = 1.0e300
  let mut minz = 1.0e300
  let mut maxx = -1.0e300
  let mut maxy = -1.0e300
  let mut maxz = -1.0e300
  for (x, y, z) in reference.coords {
    if x < minx { minx = x }
    if y < miny { miny = y }
    if z < minz { minz = z }
    if x > maxx { maxx = x }
    if y > maxy { maxy = y }
    if z > maxz { maxz = z }
  }
  minx = minx - 3.0
  miny = miny - 3.0
  minz = minz - 3.0
  maxx = maxx + 3.0
  maxy = maxy + 3.0
  maxz = maxz + 3.0
  let spacing = 0.5
  let nx = ((maxx - minx) / spacing).ceil().to_int().max(1)
  let ny = ((maxy - miny) / spacing).ceil().to_int().max(1)
  let nz = ((maxz - minz) / spacing).ceil().to_int().max(1)
  // cap the grid for safety
  let nxc = if nx > 40 { 40 } else { nx }
  let nyc = if ny > 40 { 40 } else { ny }
  let nzc = if nz > 40 { 40 } else { nz }
  let alpha = 1.78 // gaussian steepness
  let mut sum_p = 0.0
  let mut sum_r = 0.0
  let mut sum_min = 0.0
  let mut occ_p = 0
  let mut occ_r = 0
  for gx = 0; gx < nxc; gx = gx + 1 {
    for gy = 0; gy < nyc; gy = gy + 1 {
      for gz = 0; gz < nzc; gz = gz + 1 {
        let px = minx + gx.to_double() * spacing
        let py = miny + gy.to_double() * spacing
        let pz = minz + gz.to_double() * spacing
        // gaussian density of each molecule
        let mut dp = 0.0
        for i = 0; i < probe.atoms.length() && i < transformed.length(); i = i + 1 {
          let (x, y, z) = transformed[i]
          let dx = x - px
          let dy = y - py
          let dz = z - pz
          dp = dp + @math.exp(-alpha * (dx * dx + dy * dy + dz * dz))
        }
        let mut dr = 0.0
        for (x, y, z) in reference.coords {
          let dx = x - px
          let dy = y - py
          let dz = z - pz
          dr = dr + @math.exp(-alpha * (dx * dx + dy * dy + dz * dz))
        }
        sum_p = sum_p + dp
        sum_r = sum_r + dr
        if dp < dr {
          sum_min = sum_min + dp
        } else {
          sum_min = sum_min + dr
        }
        if dp > 0.35 {
          occ_p = occ_p + 1
        }
        if dr > 0.35 {
          occ_r = occ_r + 1
        }
      }
    }
  }
  let denom = sum_p + sum_r - sum_min
  let tanimoto = if denom > 0.0000000001 { sum_min / denom } else { 0.0 }
  let protrude = if occ_r > 0 {
    occ_p.to_double() / occ_r.to_double()
  } else {
    0.0
  }
  RdkitShapeAlignResult::{
    shape_tanimoto: tanimoto,
    protrude,
    grid_points: nxc * nyc * nzc,
  }
}
