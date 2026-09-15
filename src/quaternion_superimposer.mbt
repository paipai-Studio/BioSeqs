///|
/// quaternion_superimposer - Protein structure superposition using Quaternion method.
/// Inspired by Biopython Bio.PDB.Superimposer (Quaternion method).
/// Alternative to SVD superimposer; uses quaternion-based rotation optimization.

///|
/// 3D point.
pub struct QsPoint {
  x : Double
  y : Double
  z : Double
}

///|
/// Quaternion: (w, x, y, z) = (scalar, vector)
pub struct Quaternion {
  w : Double
  x : Double
  y : Double
  z : Double
}

///|
/// Superposition result: rotation matrix, translation, RMSD.
pub struct QsSuperimposeResult {
  rotation : Array[Array[Double]]
  translation : Array[Double]
  rmsd : Double
}

///|
/// Construct a point.
pub fn QsPoint::new(x : Double, y : Double, z : Double) -> QsPoint {
  QsPoint::{ x, y, z }
}

///|
/// Vector subtraction.
pub fn qs_point_sub(a : QsPoint, b : QsPoint) -> QsPoint {
  QsPoint::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

///|
/// Vector addition.
pub fn qs_point_add(a : QsPoint, b : QsPoint) -> QsPoint {
  QsPoint::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

///|
/// Dot product.
pub fn qs_point_dot(a : QsPoint, b : QsPoint) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
/// Cross product.
pub fn qs_point_cross(a : QsPoint, b : QsPoint) -> QsPoint {
  QsPoint::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  )
}

///|
/// Vector norm.
pub fn qs_point_norm(a : QsPoint) -> Double {
  (a.x * a.x + a.y * a.y + a.z * a.z).sqrt()
}

///|
/// Construct a quaternion.
pub fn Quaternion::new(w : Double, x : Double, y : Double, z : Double) -> Quaternion {
  Quaternion::{ w, x, y, z }
}

///|
/// Identity quaternion.
pub fn quaternion_identity() -> Quaternion {
  Quaternion::new(1.0, 0.0, 0.0, 0.0)
}

///|
/// Quaternion multiplication q1 * q2.
pub fn quaternion_multiply(q1 : Quaternion, q2 : Quaternion) -> Quaternion {
  Quaternion::new(
    q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
    q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
    q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
    q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
  )
}

///|
/// Quaternion conjugate.
pub fn quaternion_conjugate(q : Quaternion) -> Quaternion {
  Quaternion::new(q.w, -q.x, -q.y, -q.z)
}

///|
/// Quaternion norm.
pub fn quaternion_norm(q : Quaternion) -> Double {
  (q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z).sqrt()
}

///|
/// Normalize quaternion.
pub fn quaternion_normalize(q : Quaternion) -> Quaternion {
  let n = quaternion_norm(q)
  if n < 1.0e-12 {
    return q
  }
  Quaternion::new(q.w / n, q.x / n, q.y / n, q.z / n)
}

///|
/// Convert quaternion to 3x3 rotation matrix.
pub fn quaternion_to_matrix(q : Quaternion) -> Array[Array[Double]] {
  let qn = quaternion_normalize(q)
  let w = qn.w
  let x = qn.x
  let y = qn.y
  let z = qn.z

  let xx = x * x
  let yy = y * y
  let zz = z * z
  let xy = x * y
  let xz = x * z
  let yz = y * z
  let wx = w * x
  let wy = w * y
  let wz = w * z

  let mat : Array[Array[Double]] = []
  mat.push([1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy)])
  mat.push([2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)])
  mat.push([2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)])
  mat
}

///|
/// Compute centroid of a set of points.
pub fn qs_centroid(points : Array[QsPoint]) -> QsPoint {
  let n = points.length()
  if n == 0 {
    return QsPoint::new(0.0, 0.0, 0.0)
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  let mut i = 0
  while i < n {
    cx = cx + points[i].x
    cy = cy + points[i].y
    cz = cz + points[i].z
    i = i + 1
  }
  QsPoint::new(cx / n.to_double(), cy / n.to_double(), cz / n.to_double())
}

///|
/// Translate points by a vector.
pub fn qs_translate(points : Array[QsPoint], v : QsPoint) -> Array[QsPoint] {
  let result : Array[QsPoint] = []
  let mut i = 0
  while i < points.length() {
    result.push(qs_point_add(points[i], v))
    i = i + 1
  }
  result
}

///|
/// Compute the 4x4 N matrix for quaternion-based superposition.
/// Given point sets x (moving) and y (target), define
/// S_ij = sum_k x_k[i] * y_k[j].
/// The eigenvector (w, x, y, z) of N associated with its largest eigenvalue
/// yields the rotation matrix that maps x onto y (Coutsias et al. 2004).
fn qs_compute_n_matrix(
  x : Array[QsPoint],
  y : Array[QsPoint]
) -> Array[Array[Double]] {
  let n = x.length()
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut sxz = 0.0
  let mut syx = 0.0
  let mut syy = 0.0
  let mut syz = 0.0
  let mut szx = 0.0
  let mut szy = 0.0
  let mut szz = 0.0

  let mut i = 0
  while i < n {
    sxx = sxx + x[i].x * y[i].x
    sxy = sxy + x[i].x * y[i].y
    sxz = sxz + x[i].x * y[i].z
    syx = syx + x[i].y * y[i].x
    syy = syy + x[i].y * y[i].y
    syz = syz + x[i].y * y[i].z
    szx = szx + x[i].z * y[i].x
    szy = szy + x[i].z * y[i].y
    szz = szz + x[i].z * y[i].z
    i = i + 1
  }

  let n_mat : Array[Array[Double]] = []
  n_mat.push([sxx + syy + szz, syz - szy, szx - sxz, sxy - syx])
  n_mat.push([syz - szy, sxx - syy - szz, sxy + syx, szx + sxz])
  n_mat.push([szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy])
  n_mat.push([sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz])
  n_mat
}

///|
/// Jacobi eigenvalue decomposition of a 4x4 symmetric matrix.
/// Returns the eigenvalues (diagonal) and the matrix whose columns are the
/// corresponding eigenvectors. More robust than the power method, which fails
/// when the extreme eigenvalues have equal magnitude.
fn qs_jacobi_eigen(
  mat : Array[Array[Double]],
  max_sweeps~ : Int = 100
) -> (Array[Double], Array[Array[Double]]) {
  let n = 4
  let a : Array[Array[Double]] = []
  let v : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row_a : Array[Double] = []
    let row_v : Array[Double] = []
    let mut j = 0
    while j < n {
      row_a.push(mat[i][j])
      row_v.push(if i == j { 1.0 } else { 0.0 })
      j = j + 1
    }
    a.push(row_a)
    v.push(row_v)
    i = i + 1
  }

  let mut sweep = 0
  while sweep < max_sweeps {
    // Sum of squares of the off-diagonal entries.
    let mut off = 0.0
    let mut p = 0
    while p < n {
      let mut q = p + 1
      while q < n {
        off = off + a[p][q] * a[p][q]
        q = q + 1
      }
      p = p + 1
    }
    if off < 1.0e-20 {
      break
    }

    p = 0
    while p < n {
      let mut q = p + 1
      while q < n {
        if a[p][q].abs() > 1.0e-15 {
          let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
          let t = if theta >= 0.0 {
            1.0 / (theta + (theta * theta + 1.0).sqrt())
          } else {
            -1.0 / (-theta + (theta * theta + 1.0).sqrt())
          }
          let c = 1.0 / (t * t + 1.0).sqrt()
          let s = t * c
          let a_pp = a[p][p]
          let a_qq = a[q][q]
          let a_pq = a[p][q]
          a[p][p] = a_pp - t * a_pq
          a[q][q] = a_qq + t * a_pq
          a[p][q] = 0.0
          a[q][p] = 0.0
          let mut k = 0
          while k < n {
            if k != p && k != q {
              let a_kp = a[k][p]
              let a_kq = a[k][q]
              a[k][p] = c * a_kp - s * a_kq
              a[p][k] = a[k][p]
              a[k][q] = s * a_kp + c * a_kq
              a[q][k] = a[k][q]
            }
            k = k + 1
          }
          k = 0
          while k < n {
            let v_kp = v[k][p]
            let v_kq = v[k][q]
            v[k][p] = c * v_kp - s * v_kq
            v[k][q] = s * v_kp + c * v_kq
            k = k + 1
          }
        }
        q = q + 1
      }
      p = p + 1
    }
    sweep = sweep + 1
  }

  let eig : Array[Double] = []
  i = 0
  while i < n {
    eig.push(a[i][i])
    i = i + 1
  }
  (eig, v)
}

///|
/// Return the eigenvector of a 4x4 symmetric matrix associated with its
/// largest eigenvalue.
fn qs_largest_eigenvector(mat : Array[Array[Double]]) -> Array[Double] {
  let (eig, vecs) = qs_jacobi_eigen(mat)
  let mut best = 0
  let mut i = 1
  while i < 4 {
    if eig[i] > eig[best] {
      best = i
    }
    i = i + 1
  }
  [vecs[0][best], vecs[1][best], vecs[2][best], vecs[3][best]]
}

///|
/// Superimpose two sets of points using the Quaternion method.
/// Returns rotation matrix, translation, and RMSD.
pub fn qs_superimpose(
  a : Array[QsPoint],
  b : Array[QsPoint]
) -> QsSuperimposeResult {
  let n = a.length()
  if n == 0 {
    return QsSuperimposeResult::{
      rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
      translation: [0.0, 0.0, 0.0],
      rmsd: 0.0
    }
  }

  // Compute centroids
  let centroid_a = qs_centroid(a)
  let centroid_b = qs_centroid(b)

  // Center points
  let a_centered : Array[QsPoint] = []
  let b_centered : Array[QsPoint] = []
  let mut i = 0
  while i < n {
    a_centered.push(QsPoint::new(a[i].x - centroid_a.x, a[i].y - centroid_a.y, a[i].z - centroid_a.z))
    b_centered.push(QsPoint::new(b[i].x - centroid_b.x, b[i].y - centroid_b.y, b[i].z - centroid_b.z))
    i = i + 1
  }

  // N matrix for mapping b (moving) onto a (target)
  let n_mat = qs_compute_n_matrix(b_centered, a_centered)

  // Quaternion = eigenvector of N with the largest eigenvalue
  let q_vec = qs_largest_eigenvector(n_mat)

  // Convert quaternion to rotation matrix
  let q = Quaternion::new(q_vec[0], q_vec[1], q_vec[2], q_vec[3])
  let rot = quaternion_to_matrix(q)

  // Compute RMSD
  // RMSD = sqrt( sum_i |a_i - (R·b_i + t)|^2 / n )
  // where t = centroid_a - R·centroid_b
  let mut sum_sq = 0.0
  let mut j = 0
  while j < n {
    // Rotate b_centered[j] by rot
    let x_rot = rot[0][0] * b_centered[j].x + rot[0][1] * b_centered[j].y + rot[0][2] * b_centered[j].z
    let y_rot = rot[1][0] * b_centered[j].x + rot[1][1] * b_centered[j].y + rot[1][2] * b_centered[j].z
    let z_rot = rot[2][0] * b_centered[j].x + rot[2][1] * b_centered[j].y + rot[2][2] * b_centered[j].z

    let dx = a_centered[j].x - x_rot
    let dy = a_centered[j].y - y_rot
    let dz = a_centered[j].z - z_rot
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    j = j + 1
  }
  let rmsd = (sum_sq / n.to_double()).sqrt()

  // Translation: t = centroid_a - R·centroid_b
  let tx = centroid_a.x - (rot[0][0] * centroid_b.x + rot[0][1] * centroid_b.y + rot[0][2] * centroid_b.z)
  let ty = centroid_a.y - (rot[1][0] * centroid_b.x + rot[1][1] * centroid_b.y + rot[1][2] * centroid_b.z)
  let tz = centroid_a.z - (rot[2][0] * centroid_b.x + rot[2][1] * centroid_b.y + rot[2][2] * centroid_b.z)

  QsSuperimposeResult::{
    rotation: rot,
    translation: [tx, ty, tz],
    rmsd
  }
}

///|
/// Apply transformation to a set of points.
pub fn qs_apply_transform(
  points : Array[QsPoint],
  rotation : Array[Array[Double]],
  translation : Array[Double]
) -> Array[QsPoint] {
  let result : Array[QsPoint] = []
  let mut i = 0
  while i < points.length() {
    let x_rot = rotation[0][0] * points[i].x + rotation[0][1] * points[i].y + rotation[0][2] * points[i].z
    let y_rot = rotation[1][0] * points[i].x + rotation[1][1] * points[i].y + rotation[1][2] * points[i].z
    let z_rot = rotation[2][0] * points[i].x + rotation[2][1] * points[i].y + rotation[2][2] * points[i].z
    result.push(QsPoint::new(x_rot + translation[0], y_rot + translation[1], z_rot + translation[2]))
    i = i + 1
  }
  result
}

///|
/// Compute RMSD between two sets of points (already aligned).
pub fn qs_rmsd(a : Array[QsPoint], b : Array[QsPoint]) -> Double {
  let n = a.length()
  if n == 0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let dx = a[i].x - b[i].x
    let dy = a[i].y - b[i].y
    let dz = a[i].z - b[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Example: two simple structures for testing.
pub fn qs_example_structures() -> (Array[QsPoint], Array[QsPoint]) {
  let a = [
    QsPoint::new(1.0, 0.0, 0.0),
    QsPoint::new(0.0, 1.0, 0.0),
    QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b = [
    QsPoint::new(0.0, 1.0, 0.0),
    QsPoint::new(-1.0, 0.0, 0.0),
    QsPoint::new(0.0, 0.0, 1.0),
  ]
  (a, b)
}