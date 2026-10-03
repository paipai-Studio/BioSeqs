///|
/// ccANM — Crystal Contact ANM.
///
/// In a crystal, neighboring asymmetric units (symmetry mates) exert
/// additional forces on the protein.  ccANM extends standard ANM by
/// adding Hessian contributions from these crystal contacts, yielding
/// more realistic dynamics for crystal structures.
///
/// The Hessian remains 3N×3N (only asymmetric-unit atoms are nodes),
/// but each contact between atom i and a symmetry mate of atom k
/// contributes:
///   H_ii += (γ/d²) * n n^T     (diagonal, positive)
///   H_ik -= (γ/d²) * n n^T     (off-diagonal, negative)
///
/// Reference: Phillips, Izhaky, Phillips (2006) "Optimization of the
///            ANM for large-scale fluctuations."
///            ProDy: Bakan, Meireles, Bahar (2014) Bioinformatics 30:1782.

///|
/// Symmetry operation in fractional coordinates: x' = R·x + t.
pub(all) struct SymmetryOp {
  rot : Array[Array[Double]]  // 3×3 rotation in fractional coords
  trans : Array[Double]       // 3D translation in fractional coords
}

///|
/// Crystal unit cell + space-group symmetry operations.
pub(all) struct CrystalInfo {
  a : Double           // unit cell length (Å)
  b : Double
  c : Double
  alpha : Double       // unit cell angle (degrees)
  beta : Double
  gamma : Double
  sym_ops : Array[SymmetryOp]
}

///|
pub struct CcANM {
  n : Int
  cutoff : Double
  gamma : Double
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build fractional-to-Cartesian transformation matrix from unit cell.
/// Returns (M, M_inv) where x_cart = M * x_frac.
fn build_orthogonal_matrix(
  a : Double,
  b : Double,
  c : Double,
  alpha_deg : Double,
  beta_deg : Double,
  gamma_deg : Double,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let pi = 3.14159265358979323846
  let alpha = alpha_deg * pi / 180.0
  let beta = beta_deg * pi / 180.0
  let gamma = gamma_deg * pi / 180.0
  let ca = @math.cos(alpha)
  let cb = @math.cos(beta)
  let cg = @math.cos(gamma)
  let sg = @math.sin(gamma)
  // M (fractional → Cartesian)
  let m11 = a
  let m12 = b * cg
  let m13 = c * cb
  let m22 = b * sg
  let m23 = c * (ca - cb * cg) / sg
  let m33 = c * ((1.0 - cb * cb - (ca - cb * cg) * (ca - cb * cg) / (sg * sg)).sqrt())
  let m : Array[Array[Double]] = [
    [m11, m12, m13],
    [0.0, m22, m23],
    [0.0, 0.0, m33],
  ]
  // M_inv (Cartesian → fractional) — inverse of upper-triangular
  let inv11 = 1.0 / m11
  let inv12 = -m12 / (m11 * m22)
  let inv13 = -(m13 * m22 - m12 * m23) / (m11 * m22 * m33)
  let inv22 = 1.0 / m22
  let inv23 = -m23 / (m22 * m33)
  let inv33 = 1.0 / m33
  let m_inv : Array[Array[Double]] = [
    [inv11, inv12, inv13],
    [0.0, inv22, inv23],
    [0.0, 0.0, inv33],
  ]
  (m, m_inv)
}

///|
/// Apply symmetry operation in Cartesian coordinates.
/// x_cart' = M * (R * M_inv * x_cart + t)
fn apply_symmetry(
  x : Double,
  y : Double,
  z : Double,
  rot : Array[Array[Double]],
  trans : Array[Double],
  m : Array[Array[Double]],
  m_inv : Array[Array[Double]],
) -> (Double, Double, Double) {
  // Convert to fractional
  let fx = m_inv[0][0] * x + m_inv[0][1] * y + m_inv[0][2] * z
  let fy = m_inv[1][0] * x + m_inv[1][1] * y + m_inv[1][2] * z
  let fz = m_inv[2][0] * x + m_inv[2][1] * y + m_inv[2][2] * z
  // Apply rotation in fractional
  let rx = rot[0][0] * fx + rot[0][1] * fy + rot[0][2] * fz + trans[0]
  let ry = rot[1][0] * fx + rot[1][1] * fy + rot[1][2] * fz + trans[1]
  let rz = rot[2][0] * fx + rot[2][1] * fy + rot[2][2] * fz + trans[2]
  // Convert back to Cartesian
  let cx = m[0][0] * rx + m[0][1] * ry + m[0][2] * rz
  let cy = m[1][0] * rx + m[1][1] * ry + m[1][2] * rz
  let cz = m[2][0] * rx + m[2][1] * ry + m[2][2] * rz
  (cx, cy, cz)
}

///|
/// Build a ccANM from a protein chain and crystal information.
///
/// * `cutoff` — interaction range in Å (default 15.0)
/// * `gamma`  — spring constant (default 1.0)
/// * `ncell`  — unit-cell translation range to search (default 2)
pub fn CcANM::new(
  chain : Chain,
  crystal : CrystalInfo,
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
  ncell? : Int = 2,
) -> CcANM? {
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  if n < 3 {
    return None
  }
  let dim = 3 * n
  let h : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(dim, 0.0) })

  // Step 1: Standard intra-chain ANM Hessian
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dx = ca_coords[j].x - ca_coords[i].x
      let dy = ca_coords[j].y - ca_coords[i].y
      let dz = ca_coords[j].z - ca_coords[i].z
      let d2 = dx * dx + dy * dy + dz * dz
      let d = d2.sqrt()
      if d <= cutoff && d > 0.001 {
        let g = gamma / d2
        add_hessian_contact(h, i, j, n, g, dx, dy, dz)
      }
      j = j + 1
    }
    i = i + 1
  }

  // Step 2: Crystal contact contributions
  let (m, m_inv) = build_orthogonal_matrix(
    crystal.a, crystal.b, crystal.c,
    crystal.alpha, crystal.beta, crystal.gamma,
  )
  for sym_idx = 0; sym_idx < crystal.sym_ops.length(); sym_idx = sym_idx + 1 {
    let op = crystal.sym_ops[sym_idx]
    // For each unit cell translation (na, nb, nc)
    for na = -ncell; na <= ncell; na = na + 1 {
      for nb = -ncell; nb <= ncell; nb = nb + 1 {
        for nc = -ncell; nc <= ncell; nc = nc + 1 {
          // Skip identity (na=0, nb=0, nc=0, sym=identity) handled by step 1
          // For non-identity symmetry ops with zero translation, we skip the identity check
          let is_identity = (op.rot[0][0] == 1.0 && op.rot[1][1] == 1.0 && op.rot[2][2] == 1.0 &&
            op.trans[0] == 0.0 && op.trans[1] == 0.0 && op.trans[2] == 0.0 &&
            na == 0 && nb == 0 && nc == 0)
          if is_identity {
            continue
          }
          // For each original atom k, compute its symmetry mate position
          // Then check contacts with all original atoms i
          for k = 0; k < n; k = k + 1 {
            let (sx, sy, sz) = apply_symmetry(
              ca_coords[k].x, ca_coords[k].y, ca_coords[k].z,
              op.rot, op.trans, m, m_inv,
            )
            // Add unit cell translation
            let na_d = na.to_double()
            let nb_d = nb.to_double()
            let nc_d = nc.to_double()
            let sx2 = sx + na_d * m[0][0] + nb_d * m[0][1] + nc_d * m[0][2]
            let sy2 = sy + na_d * m[1][0] + nb_d * m[1][1] + nc_d * m[1][2]
            let sz2 = sz + na_d * m[2][0] + nb_d * m[2][1] + nc_d * m[2][2]
            // Check contacts with all original atoms
            for ci = 0; ci < n; ci = ci + 1 {
              let dx = ca_coords[ci].x - sx2
              let dy = ca_coords[ci].y - sy2
              let dz = ca_coords[ci].z - sz2
              let d2 = dx * dx + dy * dy + dz * dz
              let d = d2.sqrt()
              if d <= cutoff && d > 0.001 {
                let g = gamma / d2
                // H_ii += g * n n^T (diagonal)
                // H_ik -= g * n n^T (off-diagonal)
                add_hessian_contact(h, ci, k, n, g, dx, dy, dz)
              }
            }
          }
        }
      }
    }
  }

  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(CcANM::{ n, cutoff, gamma, hessian: h, eigenvalues, eigenvectors })
}

///|
/// Add a Hessian contact between atoms i and j.
/// H_ij -= g * [dx*dx ...], H_ii +=, H_jj +=
fn add_hessian_contact(
  h : Array[Array[Double]],
  i : Int,
  j : Int,
  _n : Int,
  g : Double,
  dx : Double,
  dy : Double,
  dz : Double,
) -> Unit {
  let hxx = g * dx * dx
  let hxy = g * dx * dy
  let hxz = g * dx * dz
  let hyy = g * dy * dy
  let hyz = g * dy * dz
  let hzz = g * dz * dz
  let i3 = 3 * i
  let j3 = 3 * j
  // Off-diagonal H_ij (negative)
  h[i3][j3] = h[i3][j3] - hxx
  h[i3][j3 + 1] = h[i3][j3 + 1] - hxy
  h[i3][j3 + 2] = h[i3][j3 + 2] - hxz
  h[i3 + 1][j3] = h[i3 + 1][j3] - hxy
  h[i3 + 1][j3 + 1] = h[i3 + 1][j3 + 1] - hyy
  h[i3 + 1][j3 + 2] = h[i3 + 1][j3 + 2] - hyz
  h[i3 + 2][j3] = h[i3 + 2][j3] - hxz
  h[i3 + 2][j3 + 1] = h[i3 + 2][j3 + 1] - hyz
  h[i3 + 2][j3 + 2] = h[i3 + 2][j3 + 2] - hzz
  // Off-diagonal H_ji (negative, symmetric)
  h[j3][i3] = h[j3][i3] - hxx
  h[j3][i3 + 1] = h[j3][i3 + 1] - hxy
  h[j3][i3 + 2] = h[j3][i3 + 2] - hxz
  h[j3 + 1][i3] = h[j3 + 1][i3] - hxy
  h[j3 + 1][i3 + 1] = h[j3 + 1][i3 + 1] - hyy
  h[j3 + 1][i3 + 2] = h[j3 + 1][i3 + 2] - hyz
  h[j3 + 2][i3] = h[j3 + 2][i3] - hxz
  h[j3 + 2][i3 + 1] = h[j3 + 2][i3 + 1] - hyz
  h[j3 + 2][i3 + 2] = h[j3 + 2][i3 + 2] - hzz
  // Diagonal H_ii (positive)
  h[i3][i3] = h[i3][i3] + hxx
  h[i3][i3 + 1] = h[i3][i3 + 1] + hxy
  h[i3][i3 + 2] = h[i3][i3 + 2] + hxz
  h[i3 + 1][i3] = h[i3 + 1][i3] + hxy
  h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + hyy
  h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] + hyz
  h[i3 + 2][i3] = h[i3 + 2][i3] + hxz
  h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] + hyz
  h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + hzz
  // Diagonal H_jj (positive)
  h[j3][j3] = h[j3][j3] + hxx
  h[j3][j3 + 1] = h[j3][j3 + 1] + hxy
  h[j3][j3 + 2] = h[j3][j3 + 2] + hxz
  h[j3 + 1][j3] = h[j3 + 1][j3] + hxy
  h[j3 + 1][j3 + 1] = h[j3 + 1][j3 + 1] + hyy
  h[j3 + 1][j3 + 2] = h[j3 + 1][j3 + 2] + hyz
  h[j3 + 2][j3] = h[j3 + 2][j3] + hxz
  h[j3 + 2][j3 + 1] = h[j3 + 2][j3 + 1] + hyz
  h[j3 + 2][j3 + 2] = h[j3 + 2][j3 + 2] + hzz
}

///|
/// Get the Hessian matrix (3N × 3N).
pub fn CcANM::hessian_matrix(self : CcANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn CcANM::get_eigenvalues(self : CcANM) -> Array[Double] {
  let dim = 3 * self.n
  let vals = Array::make(dim, 0.0)
  for i = 0; i < dim; i = i + 1 {
    vals[i] = self.eigenvalues[i]
  }
  for i = 1; i < dim; i = i + 1 {
    let key = vals[i]
    let mut j = i - 1
    while j >= 0 && vals[j] > key {
      vals[j + 1] = vals[j]
      j = j - 1
    }
    vals[j + 1] = key
  }
  vals
}

///|
/// Get eigenvector matrix.
pub fn CcANM::get_eigenvectors(self : CcANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors: B_i = sum_k (1/λ_k) * |v_k[3i:3i+3]|²
pub fn CcANM::predict_bfactors(self : CcANM) -> Array[Double] {
  let n = self.n
  let dim = 3 * n
  let bfactors = Array::make(n, 0.0)
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  for k = 6; k < dim; k = k + 1 {
    let lambda = self.eigenvalues[indices[k]]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      let vx = self.eigenvectors[3 * i][indices[k]]
      let vy = self.eigenvectors[3 * i + 1][indices[k]]
      let vz = self.eigenvectors[3 * i + 2][indices[k]]
      bfactors[i] = bfactors[i] + (vx * vx + vy * vy + vz * vz) / lambda
    }
  }
  for i = 0; i < n; i = i + 1 {
    if bfactors[i] < 0.0 {
      bfactors[i] = -bfactors[i]
    }
  }
  bfactors
}

///|
/// Mean-square fluctuations.
pub fn CcANM::fluctuations(self : CcANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes (skip 6 zero modes).
pub fn CcANM::slow_modes(
  self : CcANM,
  k : Int,
) -> Array[(Int, Double, Array[Double])] {
  let n = self.n
  let dim = 3 * n
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  let result : Array[(Int, Double, Array[Double])] = Array(capacity=0)
  let count = if k < dim - 6 { k } else { dim - 6 }
  for m = 6; m < 6 + count; m = m + 1 {
    let idx = indices[m]
    let eigval = self.eigenvalues[idx]
    let eigvec = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      eigvec[i] = self.eigenvectors[i][idx]
    }
    result.push((idx, eigval, eigvec))
  }
  result
}

///|
/// Cross-correlation matrix between residues.
pub fn CcANM::cross_correlation(self : CcANM) -> Array[Array[Double]] {
  let n = self.n
  let dim = 3 * n
  let corr : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  for k = 6; k < dim; k = k + 1 {
    let lambda = self.eigenvalues[indices[k]]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let vix = self.eigenvectors[3 * i][indices[k]]
        let viy = self.eigenvectors[3 * i + 1][indices[k]]
        let viz = self.eigenvectors[3 * i + 2][indices[k]]
        let vjx = self.eigenvectors[3 * j][indices[k]]
        let vjy = self.eigenvectors[3 * j + 1][indices[k]]
        let vjz = self.eigenvectors[3 * j + 2][indices[k]]
        corr[i][j] = corr[i][j] + (vix * vjx + viy * vjy + viz * vjz) / lambda
      }
    }
  }
  let mut max_val = 1.0e-10
  for i = 0; i < n; i = i + 1 {
    if corr[i][i].abs() > max_val {
      max_val = corr[i][i].abs()
    }
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      corr[i][j] = corr[i][j] / max_val
    }
  }
  corr
}

///|
/// Number of nodes.
pub fn CcANM::num_nodes(self : CcANM) -> Int {
  self.n
}

// ===== Common space group symmetry operations =====

///|
/// Identity operation (always present).
pub fn sym_identity() -> SymmetryOp {
  SymmetryOp::{
    rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    trans: [0.0, 0.0, 0.0],
  }
}

///|
/// P1 space group — only identity.
pub fn spacegroup_p1() -> Array[SymmetryOp] {
  [sym_identity()]
}

///|
/// P21 space group — identity + 2-fold screw along b.
pub fn spacegroup_p21() -> Array[SymmetryOp] {
  [
    sym_identity(),
    SymmetryOp::{
      rot: [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
      trans: [0.0, 0.5, 0.0],
    },
  ]
}

///|
/// P212121 space group — 4 operations (common for protein crystals).
pub fn spacegroup_p212121() -> Array[SymmetryOp] {
  [
    sym_identity(),
    SymmetryOp::{
      rot: [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
      trans: [0.5, 0.0, 0.5],
    },
    SymmetryOp::{
      rot: [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
      trans: [0.0, 0.5, 0.5],
    },
    SymmetryOp::{
      rot: [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
      trans: [0.5, 0.5, 0.0],
    },
  ]
}
