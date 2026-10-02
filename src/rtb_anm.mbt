///|
/// Rotational-Translational Block (RTB) ANM.
/// Reduces the 3N x 3N Hessian to a 6N_b x 6N_b effective Hessian
/// by treating blocks of residues as rigid bodies (6 DOF each).
/// Based on ProDy's RTB method for large proteins.

///|
/// RTB-ANM result.
pub struct RTBANM {
  n_residues : Int // total residues
  n_blocks : Int // number of blocks
  block_sizes : Array[Int] // size of each block
  block_assignments : Array[Int] // block index for each residue
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]] // in reduced space (6N_b)
  full_eigenvectors : Array[Array[Double]] // projected back to 3N
  hessian_reduced : Array[Array[Double]] // 6N_b x 6N_b
}

///|
/// Create RTB-ANM from CA coordinates.
/// Partitions residues into blocks of `block_size` residues each.
/// Each block contributes 6 degrees of freedom (3 translation + 3 rotation).
pub fn RTBANM::new(
  ca_coords : Array[Vector3],
  block_size? : Int = 10,
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
) -> RTBANM? {
  let n = ca_coords.length()
  if n < 3 {
    return None
  }
  // Partition into blocks
  let block_assignments : Array[Int] = Array::make(n, 0)
  let block_sizes : Array[Int] = Array(capacity=0)
  let n_blocks = (n + block_size - 1) / block_size
  let mut current_block = 0
  let mut count_in_block = 0
  for i = 0; i < n; i = i + 1 {
    block_assignments[i] = current_block
    count_in_block = count_in_block + 1
    if count_in_block >= block_size {
      block_sizes.push(count_in_block)
      current_block = current_block + 1
      count_in_block = 0
    }
  }
  if count_in_block > 0 {
    block_sizes.push(count_in_block)
  }
  let nb = block_sizes.length()
  let dim_reduced = 6 * nb
  let dim_full = 3 * n
  // Build full Hessian (3N x 3N)
  let h : Array[Array[Double]] = Array::makei(dim_full, { _ =>
    Array::make(dim_full, 0.0)
  })
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
        let hxx = g * dx * dx
        let hxy = g * dx * dy
        let hxz = g * dx * dz
        let hyy = g * dy * dy
        let hyz = g * dy * dz
        let hzz = g * dz * dz
        let i3 = 3 * i
        let j3 = 3 * j
        h[i3][j3] = h[i3][j3] - hxx
        h[i3][j3 + 1] = h[i3][j3 + 1] - hxy
        h[i3][j3 + 2] = h[i3][j3 + 2] - hxz
        h[i3 + 1][j3] = h[i3 + 1][j3] - hxy
        h[i3 + 1][j3 + 1] = h[i3 + 1][j3 + 1] - hyy
        h[i3 + 1][j3 + 2] = h[i3 + 1][j3 + 2] - hyz
        h[i3 + 2][j3] = h[i3 + 2][j3] - hxz
        h[i3 + 2][j3 + 1] = h[i3 + 2][j3 + 1] - hyz
        h[i3 + 2][j3 + 2] = h[i3 + 2][j3 + 2] - hzz
        h[j3][i3] = h[j3][i3] - hxx
        h[j3][i3 + 1] = h[j3][i3 + 1] - hxy
        h[j3][i3 + 2] = h[j3][i3 + 2] - hxz
        h[j3 + 1][i3] = h[j3 + 1][i3] - hxy
        h[j3 + 1][i3 + 1] = h[j3 + 1][i3 + 1] - hyy
        h[j3 + 1][i3 + 2] = h[j3 + 1][i3 + 2] - hyz
        h[j3 + 2][i3] = h[j3 + 2][i3] - hxz
        h[j3 + 2][i3 + 1] = h[j3 + 2][i3 + 1] - hyz
        h[j3 + 2][i3 + 2] = h[j3 + 2][i3 + 2] - hzz
        h[i3][i3] = h[i3][i3] + hxx
        h[i3][i3 + 1] = h[i3][i3 + 1] + hxy
        h[i3][i3 + 2] = h[i3][i3 + 2] + hxz
        h[i3 + 1][i3] = h[i3 + 1][i3] + hxy
        h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + hyy
        h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] + hyz
        h[i3 + 2][i3] = h[i3 + 2][i3] + hxz
        h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] + hyz
        h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + hzz
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
      j = j + 1
    }
    i = i + 1
  }
  // Build projection matrix P (3N x 6N_b)
  // For each block, 6 columns: [tx, ty, tz, rx, ry, rz]
  let p = Array::makei(dim_full, { _ =>
    Array::make(dim_reduced, 0.0)
  })
  for b = 0; b < nb; b = b + 1 {
    // Collect residues in this block
    let block_residues : Array[Int] = Array(capacity=0)
    for r = 0; r < n; r = r + 1 {
      if block_assignments[r] == b {
        block_residues.push(r)
      }
    }
    let sz = block_residues.length()
    if sz == 0 {
      continue
    }
    // Compute centroid
    let mut cx = 0.0
    let mut cy = 0.0
    let mut cz = 0.0
    for r in block_residues {
      cx = cx + ca_coords[r].x
      cy = cy + ca_coords[r].y
      cz = cz + ca_coords[r].z
    }
    cx = cx / sz.to_double()
    cy = cy / sz.to_double()
    cz = cz / sz.to_double()
    // Build 6 basis vectors for this block
    let col_base = 6 * b
    let inv_sqrt_sz = 1.0 / sz.to_double().sqrt()
    for r in block_residues {
      let r3 = 3 * r
      // Translation: unit vectors / sqrt(sz)
      p[r3][col_base] = inv_sqrt_sz
      p[r3 + 1][col_base + 1] = inv_sqrt_sz
      p[r3 + 2][col_base + 2] = inv_sqrt_sz
      // Rotation: (r - centroid) x e_k / norm
      let dx = ca_coords[r].x - cx
      let dy = ca_coords[r].y - cy
      let dz = ca_coords[r].z - cz
      // rx: rotation about x-axis -> (0, -dz, dy)
      p[r3 + 1][col_base + 3] = -dz * inv_sqrt_sz
      p[r3 + 2][col_base + 3] = dy * inv_sqrt_sz
      // ry: rotation about y-axis -> (dz, 0, -dx)
      p[r3][col_base + 4] = dz * inv_sqrt_sz
      p[r3 + 2][col_base + 4] = -dx * inv_sqrt_sz
      // rz: rotation about z-axis -> (-dy, dx, 0)
      p[r3][col_base + 5] = -dy * inv_sqrt_sz
      p[r3 + 1][col_base + 5] = dx * inv_sqrt_sz
    }
  }
  // Reduced Hessian: H_red = P^T * H * P
  // First compute HP = H * P (3N x 6N_b)
  let hp = Array::makei(dim_full, { _ =>
    Array::make(dim_reduced, 0.0)
  })
  for row = 0; row < dim_full; row = row + 1 {
    for col = 0; col < dim_reduced; col = col + 1 {
      let mut s = 0.0
      for kk = 0; kk < dim_full; kk = kk + 1 {
        s = s + h[row][kk] * p[kk][col]
      }
      hp[row][col] = s
    }
  }
  // H_red = P^T * HP (6N_b x 6N_b)
  let h_red = Array::makei(dim_reduced, { _ =>
    Array::make(dim_reduced, 0.0)
  })
  for i2 = 0; i2 < dim_reduced; i2 = i2 + 1 {
    for j2 = 0; j2 < dim_reduced; j2 = j2 + 1 {
      let mut s = 0.0
      for kk = 0; kk < dim_full; kk = kk + 1 {
        s = s + p[kk][i2] * hp[kk][j2]
      }
      h_red[i2][j2] = s
    }
  }
  // Eigen-decomposition of reduced Hessian
  let (eigenvalues, eigenvectors) = linalg_jacobi(h_red, dim_reduced, 200)
  // Sort by eigenvalue (ascending) using insertion sort
  let indices = Array::makei(dim_reduced, { i => i })
  for sort_i = 1; sort_i < dim_reduced; sort_i = sort_i + 1 {
    let key = indices[sort_i]
    let key_val = eigenvalues[key]
    let mut sort_j = sort_i - 1
    while sort_j >= 0 && eigenvalues[indices[sort_j]] > key_val {
      indices[sort_j + 1] = indices[sort_j]
      sort_j = sort_j - 1
    }
    indices[sort_j + 1] = key
  }
  let sorted_evals = Array::make(dim_reduced, 0.0)
  let sorted_evecs = Array::makei(dim_reduced, { _ =>
    Array::make(dim_reduced, 0.0)
  })
  for idx = 0; idx < dim_reduced; idx = idx + 1 {
    let src = indices[idx]
    sorted_evals[idx] = eigenvalues[src]
    for j = 0; j < dim_reduced; j = j + 1 {
      sorted_evecs[idx][j] = eigenvectors[src][j]
    }
  }
  // Project eigenvectors back to full space: V_full = P * V_red
  let full_evecs = Array::makei(dim_reduced, { _ =>
    Array::make(dim_full, 0.0)
  })
  for mode = 0; mode < dim_reduced; mode = mode + 1 {
    for j = 0; j < dim_full; j = j + 1 {
      let mut s = 0.0
      for kk = 0; kk < dim_reduced; kk = kk + 1 {
        s = s + p[j][kk] * sorted_evecs[mode][kk]
      }
      full_evecs[mode][j] = s
    }
  }
  Some(RTBANM::{
    n_residues: n,
    n_blocks: nb,
    block_sizes,
    block_assignments,
    eigenvalues: sorted_evals,
    eigenvectors: sorted_evecs,
    full_eigenvectors: full_evecs,
    hessian_reduced: h_red,
  })
}

///|
/// Get number of blocks.
pub fn RTBANM::num_blocks(self : RTBANM) -> Int {
  self.n_blocks
}

///|
/// Get number of residues.
pub fn RTBANM::num_residues(self : RTBANM) -> Int {
  self.n_residues
}

///|
/// Get eigenvalues.
pub fn RTBANM::get_eigenvalues(self : RTBANM) -> Array[Double] {
  self.eigenvalues
}

///|
/// Get eigenvectors (reduced space).
pub fn RTBANM::get_eigenvectors(self : RTBANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Get full-space eigenvectors (3N).
pub fn RTBANM::get_full_eigenvectors(self : RTBANM) -> Array[Array[Double]] {
  self.full_eigenvectors
}

///|
/// Predict B-factors from RTB modes.
/// B_i = sum_k (1/lambda_k) * |v_k(3i)|^2 + |v_k(3i+1)|^2 + |v_k(3i+2)|^2
pub fn RTBANM::predict_bfactors(self : RTBANM) -> Array[Double] {
  let n = self.n_residues
  let bf = Array::make(n, 0.0)
  let n_modes = self.eigenvalues.length()
  for k = 6; k < n_modes; k = k + 1 { // skip 6 zero modes
    let lam = self.eigenvalues[k]
    if lam < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      let i3 = 3 * i
      let v = self.full_eigenvectors[k]
      bf[i] = bf[i] + (v[i3] * v[i3] + v[i3 + 1] * v[i3 + 1] + v[i3 + 2] * v[i3 + 2]) / lam
    }
  }
  bf
}

///|
/// Get block sizes.
pub fn RTBANM::get_block_sizes(self : RTBANM) -> Array[Int] {
  self.block_sizes
}
