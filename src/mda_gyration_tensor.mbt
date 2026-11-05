/// Radius of gyration tensor and shape descriptors.
///
/// MDAnalysis equivalent: `analysis.RadiusOfGyration` (tensor form)
///
/// Computes the full gyration tensor, its eigenvalues (principal moments),
/// and shape descriptors: asphericity, acylindricity, relative shape anisotropy.

///|
/// Gyration tensor result
pub struct MdaGyrationTensor {
  tensor : Array[Array[Double]] // 3x3 tensor
  eigenvalues : Array[Double] // 3 principal moments (sorted descending)
  rg : Double // scalar radius of gyration
  asphericity : Double
  acylindricity : Double
  relative_shape_anisotropy : Double
}

///|
pub fn MdaGyrationTensor::tensor(
  self : MdaGyrationTensor,
) -> Array[Array[Double]] {
  self.tensor
}

///|
pub fn MdaGyrationTensor::eigenvalues(
  self : MdaGyrationTensor,
) -> Array[Double] {
  self.eigenvalues
}

///|
pub fn MdaGyrationTensor::rg(self : MdaGyrationTensor) -> Double {
  self.rg
}

///|
pub fn MdaGyrationTensor::asphericity(self : MdaGyrationTensor) -> Double {
  self.asphericity
}

///|
pub fn MdaGyrationTensor::acylindricity(self : MdaGyrationTensor) -> Double {
  self.acylindricity
}

///|
/// Compute the gyration tensor from coordinates.
///
/// coords: atom coordinates [atom][xyz]
/// masses: optional atom masses (uniform if empty)
pub fn mda_gyration_tensor(
  coords : Array[Array[Double]],
  masses~ : Array[Double] = [],
) -> MdaGyrationTensor {
  let n = coords.length()
  if n == 0 {
    return {
      tensor: [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
      eigenvalues: [0.0, 0.0, 0.0],
      rg: 0.0,
      asphericity: 0.0,
      acylindricity: 0.0,
      relative_shape_anisotropy: 0.0,
    }
  }
  // Compute center of mass
  let mut total_mass = 0.0
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for i = 0; i < n; i = i + 1 {
    let m = if masses.length() > 0 { masses[i] } else { 1.0 }
    total_mass = total_mass + m
    cx = cx + m * coords[i][0]
    cy = cy + m * coords[i][1]
    cz = cz + m * coords[i][2]
  }
  cx = cx / total_mass
  cy = cy / total_mass
  cz = cz / total_mass
  // Compute gyration tensor S = (1/M) * sum(m_i * (r_i - r_cm)(r_i - r_cm)^T)
  let mut sxx = 0.0
  let mut syy = 0.0
  let mut szz = 0.0
  let mut sxy = 0.0
  let mut sxz = 0.0
  let mut syz = 0.0
  for i = 0; i < n; i = i + 1 {
    let m = if masses.length() > 0 { masses[i] } else { 1.0 }
    let dx = coords[i][0] - cx
    let dy = coords[i][1] - cy
    let dz = coords[i][2] - cz
    sxx = sxx + m * dx * dx
    syy = syy + m * dy * dy
    szz = szz + m * dz * dz
    sxy = sxy + m * dx * dy
    sxz = sxz + m * dx * dz
    syz = syz + m * dy * dz
  }
  sxx = sxx / total_mass
  syy = syy / total_mass
  szz = szz / total_mass
  sxy = sxy / total_mass
  sxz = sxz / total_mass
  syz = syz / total_mass
  let tensor = [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
  // Eigenvalues of symmetric 3x3 (characteristic polynomial)
  let p1 = sxx * syy + sxx * szz + syy * szz - sxy * sxy - sxz * sxz - syz * syz
  let p2 = sxx * syy * szz +
    2.0 * sxy * sxz * syz -
    sxx * syz * syz -
    syy * sxz * sxz -
    szz * sxy * sxy
  let q = (sxx + syy + szz) / 3.0
  let p2_val = p1 / 3.0 - q * q
  let p3_val = p2 - q * p1 + 2.0 * q * q * q
  // If degenerate, all eigenvalues equal
  if p2_val >= 0.0 {
    return {
      tensor,
      eigenvalues: [q, q, q],
      rg: (3.0 * q).sqrt(),
      asphericity: 0.0,
      acylindricity: 0.0,
      relative_shape_anisotropy: 0.0,
    }
  }
  let r = (-p2_val * p2_val * p2_val).sqrt()
  let phi = -p3_val / (2.0 * r)
  let phi_clamped = if phi > 1.0 {
    1.0
  } else if phi < -1.0 {
    -1.0
  } else {
    phi
  }
  let theta = @math.acos(phi_clamped) / 3.0
  let e1 = q + 2.0 * (-p2_val).sqrt() * @math.cos(theta)
  let e2 = q +
    2.0 * (-p2_val).sqrt() * @math.cos(theta + 2.0 * 3.141592653589793 / 3.0)
  let e3 = q +
    2.0 * (-p2_val).sqrt() * @math.cos(theta + 4.0 * 3.141592653589793 / 3.0)
  // Sort descending
  let mut l1 = e1
  let mut l2 = e2
  let mut l3 = e3
  if l1 < l2 {
    let tmp = l1
    l1 = l2
    l2 = tmp
  }
  if l1 < l3 {
    let tmp = l1
    l1 = l3
    l3 = tmp
  }
  if l2 < l3 {
    let tmp = l2
    l2 = l3
    l3 = tmp
  }
  let rg = (l1 + l2 + l3).sqrt()
  // Shape descriptors
  // Asphericity: b = l1 - (l2 + l3)/2
  let asph = l1 - 0.5 * (l2 + l3)
  // Acylindricity: c = l2 - l3
  let acyl = l2 - l3
  // Relative shape anisotropy: kappa^2 = (3/2) * (b^2 + (3/4)c^2) / (l1+l2+l3)^2
  let trace = l1 + l2 + l3
  let rsa = if trace > 0.0 {
    1.5 * (asph * asph + 0.75 * acyl * acyl) / (trace * trace)
  } else {
    0.0
  }
  {
    tensor,
    eigenvalues: [l1, l2, l3],
    rg,
    asphericity: asph,
    acylindricity: acyl,
    relative_shape_anisotropy: rsa,
  }
}

///|
/// Compute shape descriptors from coordinates directly.
pub fn mda_shape_descriptors(
  coords : Array[Array[Double]],
) -> (Double, Double, Double) {
  let gt = mda_gyration_tensor(coords)
  (gt.asphericity, gt.acylindricity, gt.relative_shape_anisotropy)
}
