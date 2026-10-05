///|
/// Helix geometry: best-fit line/plane, helix axis, helix kink analysis.
///
/// Used for analyzing alpha-helices in protein structures, detecting
/// bends/kinks (e.g. at proline positions), and characterizing
/// transmembrane helix geometry.

///|
/// Helix kink analysis result.
pub(all) struct HelixKink {
  kink_position : Int       // residue index where kink occurs
  bend_angle : Double       // angle between pre/post axes (degrees)
  wobble_angle : Double     // rotation around perpendicular axis (degrees)
  face_shift : Double       // rotation around helix axis (degrees)
  pre_kink_axis : Vector3  // helix axis direction before kink
  post_kink_axis : Vector3 // helix axis direction after kink
}

// ─── 3x3 symmetric eigendecomposition (Jacobi) ─────────────

///|
fn hg_sqrt(x : Double) -> Double {
  if x < 0.0 { 0.0 } else { x.sqrt() }
}

///|
/// Compute eigenvalues/eigenvectors of 3x3 symmetric matrix.
/// Returns (eigenvalues, eigenvectors) sorted by eigenvalue descending.
fn jacobi_eigen3x3(
  m : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let a = Array(capacity=0)
  let mut i = 0
  while i < 3 {
    let row = Array::make(3, 0.0)
    let mut j = 0
    while j < 3 { row[j] = m[i][j]; j = j + 1 }
    a.push(row)
    i = i + 1
  }
  let v = Array(capacity=0)
  i = 0
  while i < 3 {
    let row = Array::make(3, 0.0)
    row[i] = 1.0
    v.push(row)
    i = i + 1
  }
  let mut sweep = 0
  while sweep < 100 {
    let mut max_off = 0.0
    let mut mi = 0
    let mut mj = 1
    let mut p = 0
    while p < 3 {
      let mut q = 0
      while q < 3 {
        if p != q && a[p][q].abs() > max_off {
          max_off = a[p][q].abs()
          mi = p
          mj = q
        }
        q = q + 1
      }
      p = p + 1
    }
    if max_off < 0.0000000001 { break }
    let ai = a[mi][mi]
    let aj = a[mj][mj]
    let aij = a[mi][mj]
    let theta = 0.5 * @math.atan2(2.0 * aij, aj - ai)
    let c = @math.cos(theta)
    let s = @math.sin(theta)
    // Rotate a
    let mut k = 0
    while k < 3 {
      let ak_mi = a[k][mi]
      let ak_mj = a[k][mj]
      a[k][mi] = c * ak_mi - s * ak_mj
      a[k][mj] = s * ak_mi + c * ak_mj
      a[mi][k] = a[k][mi]
      a[mj][k] = a[k][mj]
      k = k + 1
    }
    a[mi][mi] = c * c * ai - 2.0 * s * c * aij + s * s * aj
    a[mj][mj] = s * s * ai + 2.0 * s * c * aij + c * c * aj
    a[mi][mj] = 0.0
    a[mj][mi] = 0.0
    // Rotate v
    k = 0
    while k < 3 {
      let vk_mi = v[k][mi]
      let vk_mj = v[k][mj]
      v[k][mi] = c * vk_mi - s * vk_mj
      v[k][mj] = s * vk_mi + c * vk_mj
      k = k + 1
    }
    sweep = sweep + 1
  }
  // Collect and sort by eigenvalue descending
  let eigs = [a[0][0], a[1][1], a[2][2]]
  let indices = [0, 1, 2]
  let sorted = Array(capacity=0)
  for idx in indices {
    sorted.push((eigs[idx], idx))
  }
  sorted.sort_by(fn(a, b) { if a.0 > b.0 { -1 } else if a.0 < b.0 { 1 } else { 0 } })
  let evals = Array(capacity=3)
  let evecs = Array(capacity=3)
  for (eval, idx) in sorted {
    evals.push(eval)
    let vec = Array(capacity=3)
    let mut k = 0
    while k < 3 { vec.push(v[k][idx]); k = k + 1 }
    evecs.push(vec)
  }
  (evals, evecs)
}

///|
/// Compute centroid of a set of points.
fn centroid(coords : Array[Vector3]) -> Vector3 {
  let n = coords.length()
  if n == 0 { return Vector3::new(0.0, 0.0, 0.0) }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  for v in coords {
    sx = sx + v.x
    sy = sy + v.y
    sz = sz + v.z
  }
  Vector3::new(sx / n.to_double(), sy / n.to_double(), sz / n.to_double())
}

///|
/// Compute 3x3 covariance matrix of points around their centroid.
fn helix_covariance_matrix(coords : Array[Vector3], center : Vector3) -> Array[Array[Double]] {
  let n = coords.length()
  let mut cxx = 0.0
  let mut cyy = 0.0
  let mut czz = 0.0
  let mut cxy = 0.0
  let mut cxz = 0.0
  let mut cyz = 0.0
  for v in coords {
    let dx = v.x - center.x
    let dy = v.y - center.y
    let dz = v.z - center.z
    cxx = cxx + dx * dx
    cyy = cyy + dy * dy
    czz = czz + dz * dz
    cxy = cxy + dx * dy
    cxz = cxz + dx * dz
    cyz = cyz + dy * dz
  }
  let nf = n.to_double()
  [
    [cxx / nf, cxy / nf, cxz / nf],
    [cxy / nf, cyy / nf, cyz / nf],
    [cxz / nf, cyz / nf, czz / nf],
  ]
}

///|
/// Best-fit line through a set of 3D points.
/// Returns (origin=centroid, direction=principal axis).
/// The principal axis is the eigenvector with the largest eigenvalue.
pub fn calc_best_fit_line(
  coords : Array[Vector3],
) -> (Vector3, Vector3) {
  let n = coords.length()
  if n == 0 {
    return (Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0))
  }
  if n == 1 {
    return (coords[0], Vector3::new(0.0, 0.0, 1.0))
  }
  let center = centroid(coords)
  let cov = helix_covariance_matrix(coords, center)
  let (evals, evecs) = jacobi_eigen3x3(cov)
  // Largest eigenvalue is first (sorted descending)
  let dir = Vector3::new(evecs[0][0], evecs[0][1], evecs[0][2])
  // Normalize
  let mag = hg_sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z)
  let direction = if mag > 0.0 {
    Vector3::new(dir.x / mag, dir.y / mag, dir.z / mag)
  } else {
    Vector3::new(0.0, 0.0, 1.0)
  }
  (center, direction)
}

///|
/// Best-fit plane through a set of 3D points.
/// Returns (origin=centroid, normal=plane normal).
/// The normal is the eigenvector with the smallest eigenvalue.
pub fn calc_best_fit_plane(
  coords : Array[Vector3],
) -> (Vector3, Vector3) {
  let n = coords.length()
  if n == 0 {
    return (Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0))
  }
  if n == 1 {
    return (coords[0], Vector3::new(0.0, 0.0, 1.0))
  }
  let center = centroid(coords)
  let cov = helix_covariance_matrix(coords, center)
  let (evals, evecs) = jacobi_eigen3x3(cov)
  // Smallest eigenvalue is last (sorted descending)
  let normal_raw = Vector3::new(evecs[2][0], evecs[2][1], evecs[2][2])
  let mag = hg_sqrt(normal_raw.x * normal_raw.x + normal_raw.y * normal_raw.y + normal_raw.z * normal_raw.z)
  let normal = if mag > 0.0 {
    Vector3::new(normal_raw.x / mag, normal_raw.y / mag, normal_raw.z / mag)
  } else {
    Vector3::new(0.0, 0.0, 1.0)
  }
  (center, normal)
}

///|
/// Compute helix axis from CA coordinates.
/// Wrapper around calc_best_fit_line.
pub fn calc_helix_axis(
  ca_coords : Array[Vector3],
) -> (Vector3, Vector3) {
  calc_best_fit_line(ca_coords)
}

// ─── Helix kink analysis ───────────────────────────────────

///|
fn vec_sub(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

///|
fn vec_cross(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

///|
fn helix_vec_dot(a : Vector3, b : Vector3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
fn vec_normalize(v : Vector3) -> Vector3 {
  let mag = hg_sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
  if mag > 0.0 {
    Vector3::new(v.x / mag, v.y / mag, v.z / mag)
  } else {
    Vector3::new(0.0, 0.0, 0.0)
  }
}

///|
/// Compute helix kink at a given position.
///
/// Splits the helix into pre-kink and post-kink segments,
/// fits helix axes to each, and computes:
/// - bend_angle: angle between the two axes
/// - wobble_angle: rotation component perpendicular to the average axis
/// - face_shift: rotation component around the average axis
pub fn calc_helix_kink(
  ca_coords : Array[Vector3],
  kink_position : Int,
) -> HelixKink? {
  let n = ca_coords.length()
  // Need at least 3 residues on each side
  if kink_position < 6 || kink_position > n - 6 { return None }
  let pre = Array(capacity=kink_position)
  let mut i = 0
  while i < kink_position {
    pre.push(ca_coords[i])
    i = i + 1
  }
  let post = Array(capacity=n - kink_position)
  i = kink_position
  while i < n {
    post.push(ca_coords[i])
    i = i + 1
  }
  let (_, pre_dir) = calc_best_fit_line(pre)
  let (_, post_dir) = calc_best_fit_line(post)
  // Bend angle
  let dot = helix_vec_dot(pre_dir, post_dir)
  let clamped = if dot > 1.0 { 1.0 } else if dot < -1.0 { -1.0 } else { dot }
  let bend_angle = @math.acos(clamped) * 180.0 / 3.14159265358979323846
  // Average axis
  let avg_axis = vec_normalize(Vector3::new(
    pre_dir.x + post_dir.x,
    pre_dir.y + post_dir.y,
    pre_dir.z + post_dir.z,
  ))
  // Wobble: component of post_dir perpendicular to avg_axis
  let perp = vec_sub(post_dir, Vector3::new(
    avg_axis.x * helix_vec_dot(post_dir, avg_axis),
    avg_axis.y * helix_vec_dot(post_dir, avg_axis),
    avg_axis.z * helix_vec_dot(post_dir, avg_axis),
  ))
  let perp_mag = hg_sqrt(perp.x * perp.x + perp.y * perp.y + perp.z * perp.z)
  let wobble_angle = @math.atan2(perp_mag, helix_vec_dot(post_dir, avg_axis).abs()) * 180.0 / 3.14159265358979323846
  // Face-shift: rotation around the helix axis
  // Computed as the angle between the perpendicular components
  // of pre_dir and post_dir projected onto the plane normal to avg_axis
  let pre_perp = vec_sub(pre_dir, Vector3::new(
    avg_axis.x * helix_vec_dot(pre_dir, avg_axis),
    avg_axis.y * helix_vec_dot(pre_dir, avg_axis),
    avg_axis.z * helix_vec_dot(pre_dir, avg_axis),
  ))
  let post_perp = perp
  let pre_perp_n = vec_normalize(pre_perp)
  let post_perp_n = vec_normalize(post_perp)
  let fs_dot = helix_vec_dot(pre_perp_n, post_perp_n)
  let fs_clamped = if fs_dot > 1.0 { 1.0 } else if fs_dot < -1.0 { -1.0 } else { fs_dot }
  let fs_cross = helix_vec_dot(vec_cross(pre_perp_n, post_perp_n), avg_axis)
  let face_shift = @math.atan2(fs_cross, fs_clamped) * 180.0 / 3.14159265358979323846
  Some(HelixKink::{
    kink_position,
    bend_angle,
    wobble_angle,
    face_shift,
    pre_kink_axis: pre_dir,
    post_kink_axis: post_dir,
  })
}

///|
/// Scan a helix for kinks by sliding a window.
///
/// Computes the local bend angle at each position and returns
/// positions where the bend exceeds a threshold.
pub fn analyze_helix_kinks(
  ca_coords : Array[Vector3],
  window_size? : Int = 4,
) -> Array[HelixKink] {
  let n = ca_coords.length()
  let kinks = Array(capacity=0)
  let mut i = window_size
  while i < n - window_size {
    match calc_helix_kink(ca_coords, i) {
      Some(kink) =>
        if kink.bend_angle > 20.0 {
          kinks.push(kink)
        }
      None => ()
    }
    i = i + 1
  }
  kinks
}
