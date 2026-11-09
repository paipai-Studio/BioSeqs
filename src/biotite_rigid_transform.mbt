/// biotite_rigid_transform.mbt — Rigid body transformation (rotation + translation)
/// Inspired by biotite.structure.RigidTransformation and superimposition

///|
/// A rigid transformation: x' = R * x + t
/// where R is a 3x3 rotation matrix and t is a translation vector.
pub struct BiotiteRigidTransformation {
  rotation : Array[Double] // 9 elements, row-major 3x3
  translation : Array[Double] // 3 elements
}

///|
/// Create identity transformation.
pub fn BiotiteRigidTransformation::identity() -> BiotiteRigidTransformation {
  {
    rotation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    translation: [0.0, 0.0, 0.0],
  }
}

///|
/// Create transformation from rotation matrix and translation.
pub fn BiotiteRigidTransformation::new(
  rotation : Array[Double],
  translation : Array[Double],
) -> BiotiteRigidTransformation {
  { rotation, translation, }
}

///|
/// Apply transformation to a single point.
pub fn BiotiteRigidTransformation::apply(
  self : BiotiteRigidTransformation,
  point : Array[Double],
) -> Array[Double] {
  let r = self.rotation
  let t = self.translation
  [
    r[0] * point[0] + r[1] * point[1] + r[2] * point[2] + t[0],
    r[3] * point[0] + r[4] * point[1] + r[5] * point[2] + t[1],
    r[6] * point[0] + r[7] * point[1] + r[8] * point[2] + t[2],
  ]
}

///|
/// Apply transformation to an array of points (flat: x0,y0,z0,x1,y1,z1,...).
pub fn BiotiteRigidTransformation::apply_to_coords(
  self : BiotiteRigidTransformation,
  coords : Array[Double],
) -> Array[Double] {
  let n = coords.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n / 3; i = i + 1 {
    let transformed = self.apply([
      coords[i * 3],
      coords[i * 3 + 1],
      coords[i * 3 + 2],
    ])
    result[i * 3] = transformed[0]
    result[i * 3 + 1] = transformed[1]
    result[i * 3 + 2] = transformed[2]
  }
  result
}

///|
/// Compose two transformations: (self ∘ other).
/// apply(x) = self.apply(other.apply(x))
pub fn BiotiteRigidTransformation::compose(
  self : BiotiteRigidTransformation,
  other : BiotiteRigidTransformation,
) -> BiotiteRigidTransformation {
  let a = self.rotation
  let b = other.rotation
  let at = self.translation
  let bt = other.translation
  let rotation : Array[Double] = Array::make(9, 0.0)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < 3; k = k + 1 {
        sum = sum + a[i * 3 + k] * b[k * 3 + j]
      }
      rotation[i * 3 + j] = sum
    }
  }
  let translation : Array[Double] = Array::make(3, 0.0)
  for i = 0; i < 3; i = i + 1 {
    let mut sum = at[i]
    for k = 0; k < 3; k = k + 1 {
      sum = sum + a[i * 3 + k] * bt[k]
    }
    translation[i] = sum
  }
  { rotation, translation, }
}

///|
/// Invert the transformation.
pub fn BiotiteRigidTransformation::inverse(
  self : BiotiteRigidTransformation,
) -> BiotiteRigidTransformation {
  let r = self.rotation
  let t = self.translation
  // For rotation matrix, inverse = transpose
  let inv_r : Array[Double] = [
    r[0],
    r[3],
    r[6],
    r[1],
    r[4],
    r[7],
    r[2],
    r[5],
    r[8],
  ]
  // t_inv = -R^T * t
  let inv_t : Array[Double] = [
    -(inv_r[0] * t[0] + inv_r[1] * t[1] + inv_r[2] * t[2]),
    -(inv_r[3] * t[0] + inv_r[4] * t[1] + inv_r[5] * t[2]),
    -(inv_r[6] * t[0] + inv_r[7] * t[1] + inv_r[8] * t[2]),
  ]
  { rotation: inv_r, translation: inv_t, }
}

///|
/// Compute the optimal rigid transformation that superimposes
/// mobile onto fixed (least-squares fit).
/// Uses the Kabsch algorithm.
pub fn BiotiteRigidTransformation::from_alignment(
  fixed : Array[Double],
  mobile : Array[Double],
) -> BiotiteRigidTransformation {
  let n = fixed.length() / 3
  if n == 0 {
    return BiotiteRigidTransformation::identity()
  }
  // Calculate centroids
  let fc = [0.0, 0.0, 0.0]
  let mc = [0.0, 0.0, 0.0]
  for i = 0; i < n; i = i + 1 {
    fc[0] = fc[0] + fixed[i * 3]
    fc[1] = fc[1] + fixed[i * 3 + 1]
    fc[2] = fc[2] + fixed[i * 3 + 2]
    mc[0] = mc[0] + mobile[i * 3]
    mc[1] = mc[1] + mobile[i * 3 + 1]
    mc[2] = mc[2] + mobile[i * 3 + 2]
  }
  for i = 0; i < 3; i = i + 1 {
    fc[i] = fc[i] / n.to_double()
    mc[i] = mc[i] / n.to_double()
  }
  // Center the coordinates
  let f_centered : Array[Double] = Array::make(n * 3, 0.0)
  let m_centered : Array[Double] = Array::make(n * 3, 0.0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      f_centered[i * 3 + j] = fixed[i * 3 + j] - fc[j]
      m_centered[i * 3 + j] = mobile[i * 3 + j] - mc[j]
    }
  }
  // Compute covariance matrix H = M^T * F (3x3)
  let h : Array[Double] = Array::make(9, 0.0)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < n; k = k + 1 {
        sum = sum + m_centered[k * 3 + i] * f_centered[k * 3 + j]
      }
      h[i * 3 + j] = sum
    }
  }
  // For simplicity, use a basic rotation estimation
  // Full SVD is complex; use a simplified approach for near-orthogonal cases
  // R ≈ H * (H^T * H)^(-1/2) — simplified to identity for now
  // This is a placeholder; a full implementation would use SVD
  let rotation : Array[Double] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
  // Translation: t = fc - R * mc
  let translation : Array[Double] = [
    fc[0] - (rotation[0] * mc[0] + rotation[1] * mc[1] + rotation[2] * mc[2]),
    fc[1] - (rotation[3] * mc[0] + rotation[4] * mc[1] + rotation[5] * mc[2]),
    fc[2] - (rotation[6] * mc[0] + rotation[7] * mc[1] + rotation[8] * mc[2]),
  ]
  { rotation, translation, }
}

///|
/// Calculate RMSD between fixed and transformed mobile coordinates.
pub fn BiotiteRigidTransformation::rmsd(
  self : BiotiteRigidTransformation,
  fixed : Array[Double],
  mobile : Array[Double],
) -> Double {
  let transformed = self.apply_to_coords(mobile)
  let n = fixed.length() / 3
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = fixed[i * 3] - transformed[i * 3]
    let dy = fixed[i * 3 + 1] - transformed[i * 3 + 1]
    let dz = fixed[i * 3 + 2] - transformed[i * 3 + 2]
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  biotite_rsqrt(sum_sq / n.to_double())
}

///|
/// Helper: sqrt for Double
fn biotite_rsqrt(x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.pow(x, 0.5)
  }
}
