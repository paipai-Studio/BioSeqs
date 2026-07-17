///|
/// SVD Superimposer module - Protein structure superposition using SVD.
/// Inspired by Biopython Bio.PDB.SVDSuperimposer.
///
/// Provides functions for aligning protein structures and calculating RMSD.

///|
/// Represents 3D coordinates of an atom.
pub struct AtomCoordinate {
  x : Double
  y : Double
  z : Double
}

///|
/// Create a new AtomCoordinate.
pub fn AtomCoordinate::new(
  x : Double,
  y : Double,
  z : Double,
) -> AtomCoordinate {
  AtomCoordinate::{ x, y, z }
}

///|
/// Convert to array.
pub fn AtomCoordinate::to_array(self : AtomCoordinate) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Calculate distance between two atoms.
pub fn AtomCoordinate::distance(
  self : AtomCoordinate,
  other : AtomCoordinate,
) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// SVD Superimposer struct.
pub struct SVDSuperimposer {
  rot : Array[Array[Double]]
  tran : Array[Double]
  rmsd : Double
}

///|
/// Create a new SVDSuperimposer.
pub fn SVDSuperimposer::new() -> SVDSuperimposer {
  SVDSuperimposer::{
    rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    tran: [0.0, 0.0, 0.0],
    rmsd: 0.0,
  }
}

///|
/// Center coordinates to origin.
fn center_coords(
  coords : Array[AtomCoordinate],
) -> (Array[AtomCoordinate], Array[Double]) {
  let n = coords.length()
  if n == 0 {
    return (coords, [0.0, 0.0, 0.0])
  }

  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_z = 0.0
  let mut i = 0
  while i < n {
    sum_x = sum_x + coords[i].x
    sum_y = sum_y + coords[i].y
    sum_z = sum_z + coords[i].z
    i = i + 1
  }

  let cx = sum_x / n.to_double()
  let cy = sum_y / n.to_double()
  let cz = sum_z / n.to_double()

  let centered : Array[AtomCoordinate] = Array::new()
  i = 0
  while i < n {
    centered.push(
      AtomCoordinate::new(coords[i].x - cx, coords[i].y - cy, coords[i].z - cz),
    )
    i = i + 1
  }

  (centered, [cx, cy, cz])
}

///|
/// 3x3 matrix multiplication.
fn matmul_3x3(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < 3 {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < 3 {
      let mut sum = 0.0
      let mut k = 0
      while k < 3 {
        sum = sum + a[i][k] * b[k][j]
        k = k + 1
      }
      row.push(sum)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Multiply 3x3 matrix with 3D vector.
fn matmul_vec_3x3(m : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < 3 {
    let mut sum = 0.0
    let mut j = 0
    while j < 3 {
      sum = sum + m[i][j] * v[j]
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

///|
/// Transpose 3x3 matrix.
fn mat_transpose(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < 3 {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < 3 {
      row.push(m[j][i])
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Compute cross-correlation matrix.
fn compute_cross_correlation(
  coords1 : Array[AtomCoordinate],
  coords2 : Array[AtomCoordinate],
) -> Array[Array[Double]] {
  let c_0 = Array::make(3, 0.0)
  let c_1 = Array::make(3, 0.0)
  let c_2 = Array::make(3, 0.0)

  let n = coords1.length()

  let mut i = 0
  while i < n {
    let a = coords1[i]
    let b = coords2[i]
    c_0[0] = c_0[0] + a.x * b.x
    c_0[1] = c_0[1] + a.x * b.y
    c_0[2] = c_0[2] + a.x * b.z
    c_1[0] = c_1[0] + a.y * b.x
    c_1[1] = c_1[1] + a.y * b.y
    c_1[2] = c_1[2] + a.y * b.z
    c_2[0] = c_2[0] + a.z * b.x
    c_2[1] = c_2[1] + a.z * b.y
    c_2[2] = c_2[2] + a.z * b.z
    i = i + 1
  }

  [c_0, c_1, c_2]
}

///|
/// Compute determinant of 3x3 matrix.
fn determinant_3x3(m : Array[Array[Double]]) -> Double {
  let a = m[0][0]
  let b = m[0][1]
  let c = m[0][2]
  let d = m[1][0]
  let e = m[1][1]
  let f = m[1][2]
  let g = m[2][0]
  let h = m[2][1]
  let i = m[2][2]
  a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

///|
/// Compute eigenvalues and eigenvectors of 3x3 symmetric matrix using Jacobi method.
fn jacobi_eigenvalue_3x3(
  m : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let mut a = Array::new()
  let mut i = 0
  while i < 3 {
    let row = Array::make(3, 0.0)
    let mut j = 0
    while j < 3 {
      row[j] = m[i][j]
      j = j + 1
    }
    a.push(row)
    i = i + 1
  }

  let mut v = Array::new()
  i = 0
  while i < 3 {
    let row = Array::make(3, 0.0)
    row[i] = 1.0
    v.push(row)
    i = i + 1
  }

  let mut sweep = 0
  while sweep < 100 {
    let mut max_off_diag = 0.0
    let mut max_i = 0
    let mut max_j = 1

    let mut p = 0
    while p < 3 {
      let mut q = 0
      while q < 3 {
        if p != q && a[p][q].abs() > max_off_diag {
          max_off_diag = a[p][q].abs()
          max_i = p
          max_j = q
        }
        q = q + 1
      }
      p = p + 1
    }

    if max_off_diag < 0.0000000001 {
      break
    }

    let ai = a[max_i][max_i]
    let aj = a[max_j][max_j]
    let aij = a[max_i][max_j]

    let theta = 0.5 * @math.atan((aj - ai) / (2.0 * aij))
    let c = @math.cos(theta)
    let s = @math.sin(theta)

    let new_a = Array::new()
    let mut r = 0
    while r < 3 {
      let row = Array::make(3, 0.0)
      let mut t = 0
      while t < 3 {
        row[t] = a[r][t]
        t = t + 1
      }
      new_a.push(row)
      r = r + 1
    }

    new_a[max_i][max_i] = ai * c * c + aj * s * s + 2.0 * aij * c * s
    new_a[max_j][max_j] = ai * s * s + aj * c * c - 2.0 * aij * c * s
    new_a[max_i][max_j] = (aj - ai) * c * s + aij * (c * c - s * s)
    new_a[max_j][max_i] = new_a[max_i][max_j]

    let mut k = 0
    while k < 3 {
      if k != max_i && k != max_j {
        let aik = a[max_i][k]
        let ajk = a[max_j][k]
        new_a[max_i][k] = aik * c + ajk * s
        new_a[k][max_i] = new_a[max_i][k]
        new_a[max_j][k] = -aik * s + ajk * c
        new_a[k][max_j] = new_a[max_j][k]
      }
      k = k + 1
    }

    a = new_a

    let new_v = Array::new()
    r = 0
    while r < 3 {
      let row = Array::make(3, 0.0)
      let mut t = 0
      while t < 3 {
        row[t] = v[r][t]
        t = t + 1
      }
      new_v.push(row)
      r = r + 1
    }

    k = 0
    while k < 3 {
      let vik = v[k][max_i]
      let vjk = v[k][max_j]
      new_v[k][max_i] = vik * c + vjk * s
      new_v[k][max_j] = -vik * s + vjk * c
      k = k + 1
    }

    v = new_v
    sweep = sweep + 1
  }

  let eigenvalues = [a[0][0], a[1][1], a[2][2]]
  (eigenvalues, v)
}

///|
/// Compute SVD of 3x3 matrix.
fn svd_3x3(
  m : Array[Array[Double]],
) -> (Array[Array[Double]], Array[Double], Array[Array[Double]]) {
  let m_t = mat_transpose(m)
  let cov = matmul_3x3(m_t, m)

  let (eigenvalues, v) = jacobi_eigenvalue_3x3(cov)

  let sorted_indices = Array::make(3, 0)
  sorted_indices[0] = 0
  sorted_indices[1] = 1
  sorted_indices[2] = 2

  let mut i = 0
  while i < 2 {
    let mut j = i + 1
    while j < 3 {
      if eigenvalues[sorted_indices[j]] > eigenvalues[sorted_indices[i]] {
        let tmp = sorted_indices[i]
        sorted_indices[i] = sorted_indices[j]
        sorted_indices[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  let s : Array[Double] = Array::new()
  let v_sorted : Array[Array[Double]] = Array::new()
  i = 0
  while i < 3 {
    let idx = sorted_indices[i]
    let val = if eigenvalues[idx] > 0.0 { eigenvalues[idx].sqrt() } else { 0.0 }
    s.push(val)

    let col = Array::make(3, 0.0)
    let mut j = 0
    while j < 3 {
      col[j] = v[j][idx]
      j = j + 1
    }
    v_sorted.push(col)
    i = i + 1
  }

  let u : Array[Array[Double]] = Array::new()
  i = 0
  while i < 3 {
    if s[i] > 0.0000000001 {
      let col = matmul_vec_3x3(m, v_sorted[i])
      let col_normalized = Array::make(3, 0.0)
      let mut j = 0
      while j < 3 {
        col_normalized[j] = col[j] / s[i]
        j = j + 1
      }
      u.push(col_normalized)
    } else {
      let col = Array::make(3, 0.0)
      col[0] = 1.0
      u.push(col)
    }
    i = i + 1
  }

  let v_t = mat_transpose(v_sorted)
  (u, s, v_t)
}

///|
/// Set atoms for superimposition.
pub fn SVDSuperimposer::set(
  atoms1 : Array[AtomCoordinate],
  atoms2 : Array[AtomCoordinate],
) -> SVDSuperimposer {
  let n = atoms1.length()
  if n != atoms2.length() || n == 0 {
    return SVDSuperimposer::new()
  }

  let (centered1, center1) = center_coords(atoms1)
  let (centered2, center2) = center_coords(atoms2)

  let c = compute_cross_correlation(centered1, centered2)

  let (u, _, v_t) = svd_3x3(c)

  let mut rot = matmul_3x3(u, v_t)

  let det = determinant_3x3(rot)
  if det < 0.0 {
    let u_adj = Array::new()
    let mut i = 0
    while i < 3 {
      let row = Array::make(3, 0.0)
      let mut j = 0
      while j < 3 {
        row[j] = u[i][j]
        j = j + 1
      }
      u_adj.push(row)
      i = i + 1
    }
    let u_adj_row2 = Array::make(3, 0.0)
    u_adj_row2[0] = -u_adj[2][0]
    u_adj_row2[1] = -u_adj[2][1]
    u_adj_row2[2] = -u_adj[2][2]
    let u_adj_fixed = Array::new()
    u_adj_fixed.push(u_adj[0])
    u_adj_fixed.push(u_adj[1])
    u_adj_fixed.push(u_adj_row2)
    rot = matmul_3x3(u_adj_fixed, v_t)
  }

  let rot_center1 = matmul_vec_3x3(rot, center1)
  let tran : Array[Double] = Array::new()
  let mut k = 0
  while k < 3 {
    tran.push(center2[k] - rot_center1[k])
    k = k + 1
  }

  let mut rmsd_sum = 0.0
  let mut m = 0
  while m < n {
    let transformed = matmul_vec_3x3(rot, centered1[m].to_array())
    let dx = centered2[m].x - transformed[0]
    let dy = centered2[m].y - transformed[1]
    let dz = centered2[m].z - transformed[2]
    rmsd_sum = rmsd_sum + dx * dx + dy * dy + dz * dz
    m = m + 1
  }
  let rmsd = (rmsd_sum / n.to_double()).sqrt()

  SVDSuperimposer::{ rot, tran, rmsd }
}

///|
/// Get rotation matrix.
pub fn SVDSuperimposer::get_rotation(
  self : SVDSuperimposer,
) -> Array[Array[Double]] {
  self.rot
}

///|
/// Get translation vector.
pub fn SVDSuperimposer::get_translation(
  self : SVDSuperimposer,
) -> Array[Double] {
  self.tran
}

///|
/// Get RMSD.
pub fn SVDSuperimposer::get_rmsd(self : SVDSuperimposer) -> Double {
  self.rmsd
}

///|
/// Apply transformation to coordinates.
pub fn SVDSuperimposer::apply(
  self : SVDSuperimposer,
  coords : Array[AtomCoordinate],
) -> Array[AtomCoordinate] {
  let result : Array[AtomCoordinate] = Array::new()
  let mut i = 0
  while i < coords.length() {
    let transformed = matmul_vec_3x3(self.rot, coords[i].to_array())
    result.push(
      AtomCoordinate::new(
        transformed[0] + self.tran[0],
        transformed[1] + self.tran[1],
        transformed[2] + self.tran[2],
      ),
    )
    i = i + 1
  }
  result
}

///|
/// Superimpose two sets of atoms.
pub fn superimpose(
  atoms1 : Array[AtomCoordinate],
  atoms2 : Array[AtomCoordinate],
) -> (Array[AtomCoordinate], Double) {
  let superimposer = SVDSuperimposer::set(atoms1, atoms2)
  let transformed = superimposer.apply(atoms1)
  (transformed, superimposer.rmsd)
}

///|
/// Calculate RMSD between two sets of coordinates.
pub fn calculate_rmsd(
  coords1 : Array[AtomCoordinate],
  coords2 : Array[AtomCoordinate],
) -> Double {
  let n = coords1.length()
  if n != coords2.length() || n == 0 {
    return 0.0
  }

  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let dx = coords1[i].x - coords2[i].x
    let dy = coords1[i].y - coords2[i].y
    let dz = coords1[i].z - coords2[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }

  (sum_sq / n.to_double()).sqrt()
}
