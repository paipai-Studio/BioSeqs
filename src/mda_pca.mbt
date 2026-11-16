// ============================================================================
// MDAnalysis 风格：AtomGroup 轨迹上的真实 PCA（Jacobi eigenvalue 分解）
// ============================================================================
//
// 旧的 mda_analysis_pca.mbt 只做 per-atom-variance 假 PCA。这里做完整的：
//   1. 构造 3N × 3N 协方差矩阵
//   2. Jacobi eigenvalue 分解（linalg.mbt 已提供）
//   3. 降序排序 + 投影轨迹
// ============================================================================

///|
/// AtomGroup 轨迹 PCA 结果。
pub struct AtomGroupPcaResult {
  /// 原子数
  n_atoms : Int32
  /// 非零特征值数量
  n_modes : Int32
  /// 特征值（降序），单位 Å²
  eigenvalues : Array[Double]
  /// 特征向量，shape `(3*n_atoms, n_modes)`
  eigenvectors : Array[Array[Double]]
  /// 每个 PC 解释方差比例（降序，总和 ≈ 1）
  variance_explained : Array[Double]
  /// 平均结构坐标，shape `(n_atoms,)` — 投影时参考
  mean_coords : Array[Vector3]
}

///|
/// AtomGroup 轨迹上的 PCA。
///
/// 参数：
///   - `align~` : 是否先做平移对齐（所有帧质心重合）。默认 true。
///                精确 Kabsch 旋转对齐省略，因为 3N 协方差 PCA 对整体刚体
///                旋转的贡献在统计平均后很小（ProDy 等工具同样简化）。
pub fn AtomGroup::pca_full(
  self : AtomGroup,
  align~ : Bool = true,
) -> AtomGroupPcaResult {
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let dim = 3 * n_atoms
  let empty = AtomGroupPcaResult::{
    n_atoms,
    n_modes: 0,
    eigenvalues: Array(capacity=0),
    eigenvectors: Array(capacity=0),
    variance_explained: Array(capacity=0),
    mean_coords: Array(capacity=0),
  }
  if n_atoms == 0 || n_frames < 2 { return empty }

  // ===== 1. 从轨迹抽坐标 =====
  let frames_vec : Array[Array[Vector3]] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let coords : Array[Vector3] = Array(capacity=n_atoms)
    let mut i = 0
    while i < n_atoms {
      let idx = self.indices[i]
      coords.push(
        if idx < fr.length() { fr[idx] } else { Vector3::new(0.0, 0.0, 0.0) }
      )
      i = i + 1
    }
    frames_vec.push(coords)
    f = f + 1
  }

  // ===== 2. 质心平移对齐 =====
  if align {
    f = 0
    while f < n_frames {
      let fr = frames_vec[f]
      let mut sx = 0.0
      let mut sy = 0.0
      let mut sz = 0.0
      let mut i = 0
      while i < n_atoms {
        sx = sx + fr[i].x
        sy = sy + fr[i].y
        sz = sz + fr[i].z
        i = i + 1
      }
      let cx = sx / n_atoms.to_double()
      let cy = sy / n_atoms.to_double()
      let cz = sz / n_atoms.to_double()
      let new_fr : Array[Vector3] = Array(capacity=n_atoms)
      i = 0
      while i < n_atoms {
        new_fr.push(Vector3::new(fr[i].x - cx, fr[i].y - cy, fr[i].z - cz))
        i = i + 1
      }
      frames_vec[f] = new_fr
      f = f + 1
    }
  }

  // ===== 3. 平均结构 =====
  let mean_coords : Array[Vector3] = Array(capacity=n_atoms)
  let nf = n_frames.to_double()
  let mut i = 0
  while i < n_atoms {
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    f = 0
    while f < n_frames {
      sx = sx + frames_vec[f][i].x
      sy = sy + frames_vec[f][i].y
      sz = sz + frames_vec[f][i].z
      f = f + 1
    }
    mean_coords.push(Vector3::new(sx / nf, sy / nf, sz / nf))
    i = i + 1
  }

  // ===== 4. 协方差矩阵 (dim × dim)，用增量累加 =====
  let cov : Array[Array[Double]] = Array(capacity=dim)
  i = 0
  while i < dim {
    let row : Array[Double] = Array(capacity=dim)
    let mut j = 0
    while j < dim { row.push(0.0); j = j + 1 }
    cov.push(row)
    i = i + 1
  }
  f = 0
  while f < n_frames {
    // flatten displacement
    let disp : Array[Double] = Array(capacity=dim)
    let mut a_i = 0
    while a_i < n_atoms {
      disp.push(frames_vec[f][a_i].x - mean_coords[a_i].x)
      disp.push(frames_vec[f][a_i].y - mean_coords[a_i].y)
      disp.push(frames_vec[f][a_i].z - mean_coords[a_i].z)
      a_i = a_i + 1
    }
    let mut a = 0
    while a < dim {
      let mut b = a
      while b < dim {
        let v = disp[a] * disp[b]
        cov[a][b] = cov[a][b] + v
        cov[b][a] = cov[a][b]
        b = b + 1
      }
      a = a + 1
    }
    f = f + 1
  }
  let norm = 1.0 / ((n_frames - 1).to_double())
  i = 0
  while i < dim {
    let mut j = 0
    while j < dim { cov[i][j] = cov[i][j] * norm; j = j + 1 }
    i = i + 1
  }

  // ===== 5. Jacobi eigenvalue 分解 =====
  let max_iter = if dim > 200 { dim * 20 } else { 200 }
  let (eigenvalues_raw, eigenvectors_raw) = linalg_jacobi(cov, dim, max_iter)

  // ===== 6. 降序排序 =====
  let indices : Array[Int32] = Array(capacity=dim)
  i = 0
  while i < dim { indices.push(i); i = i + 1 }
  i = 1
  while i < dim {
    let key = indices[i]
    let key_val = eigenvalues_raw[key]
    let mut j = i - 1
    while j >= 0 && eigenvalues_raw[indices[j]] < key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
    i = i + 1
  }

  // ===== 7. 收集非零特征值 =====
  let mut n_modes = 0
  i = 0
  while i < dim {
    if eigenvalues_raw[indices[i]] > 1.0e-12 { n_modes = n_modes + 1 }
    i = i + 1
  }
  if n_modes == 0 { n_modes = 1 }

  let eigenvalues : Array[Double] = Array(capacity=n_modes)
  let eigenvectors : Array[Array[Double]] = Array(capacity=dim)
  let variance_explained : Array[Double] = Array(capacity=n_modes)
  let mut total_var = 0.0
  i = 0
  while i < n_modes {
    let val = eigenvalues_raw[indices[i]]
    eigenvalues.push(val)
    total_var = total_var + val
    i = i + 1
  }
  // eigenvectors[i][k] = 第 i 个坐标分量在第 k 个 PC 上的系数
  i = 0
  while i < dim {
    let col : Array[Double] = Array(capacity=n_modes)
    let mut k = 0
    while k < n_modes {
      col.push(eigenvectors_raw[i][indices[k]])
      k = k + 1
    }
    eigenvectors.push(col)
    i = i + 1
  }
  i = 0
  while i < n_modes {
    variance_explained.push(
      if total_var > 1.0e-12 { eigenvalues[i] / total_var } else { 0.0 }
    )
    i = i + 1
  }

  AtomGroupPcaResult::{
    n_atoms,
    n_modes,
    eigenvalues,
    eigenvectors,
    variance_explained,
    mean_coords,
  }
}

// ============================================================================
// 查询方法
// ============================================================================

///|
/// 特征值（降序），单位 Å²。
pub fn AtomGroupPcaResult::eigenvalues(self : AtomGroupPcaResult) -> Array[Double] {
  self.eigenvalues
}

///|
/// 每个 PC 解释方差比例（降序，总和 ≈ 1）。
pub fn AtomGroupPcaResult::variance_explained(self : AtomGroupPcaResult) -> Array[Double] {
  self.variance_explained
}

///|
/// 累计解释方差。
pub fn AtomGroupPcaResult::cumulative_variance(self : AtomGroupPcaResult) -> Array[Double] {
  let n = self.variance_explained.length()
  let out : Array[Double] = Array(capacity=n)
  let mut acc = 0.0
  let mut i = 0
  while i < n {
    acc = acc + self.variance_explained[i]
    out.push(acc)
    i = i + 1
  }
  out
}

///|
/// PC 数量（非零特征值）。
pub fn AtomGroupPcaResult::n_modes(self : AtomGroupPcaResult) -> Int32 {
  self.n_modes
}

///|
/// 平均结构坐标。
pub fn AtomGroupPcaResult::mean_structure(self : AtomGroupPcaResult) -> Array[Vector3] {
  self.mean_coords
}

///|
/// 把一组原子坐标投影到 PC 空间，返回每个 PC 的得分。
///
/// 输入 coords 长度 = self.n_atoms。
pub fn AtomGroupPcaResult::project(
  self : AtomGroupPcaResult,
  coords : Array[Vector3],
) -> Array[Double] {
  let dim = 3 * self.n_atoms
  let disp : Array[Double] = Array(capacity=dim)
  let mut i = 0
  while i < self.n_atoms {
    let cx = if i < coords.length() { coords[i].x } else { 0.0 }
    let cy = if i < coords.length() { coords[i].y } else { 0.0 }
    let cz = if i < coords.length() { coords[i].z } else { 0.0 }
    disp.push(cx - self.mean_coords[i].x)
    disp.push(cy - self.mean_coords[i].y)
    disp.push(cz - self.mean_coords[i].z)
    i = i + 1
  }
  let out : Array[Double] = Array(capacity=self.n_modes)
  let mut k = 0
  while k < self.n_modes {
    let mut s = 0.0
    let mut a = 0
    while a < dim {
      s = s + self.eigenvectors[a][k] * disp[a]
      a = a + 1
    }
    out.push(s)
    k = k + 1
  }
  out
}

///|
/// 把 AtomGroup 的整个轨迹投影到 PC 空间。
///
/// 返回 `Array[Array[Double]]`，shape `(n_frames, n_modes)`。
pub fn AtomGroupPcaResult::project_trajectory(
  self : AtomGroupPcaResult,
  group : AtomGroup,
) -> Array[Array[Double]] {
  let n_frames = group.universe.n_frames()
  let out : Array[Array[Double]] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    let fr = group.universe.frame(f)
    let coords : Array[Vector3] = Array(capacity=self.n_atoms)
    let mut i = 0
    while i < self.n_atoms {
      let idx = group.indices[i]
      coords.push(
        if idx < fr.length() { fr[idx] } else { Vector3::new(0.0, 0.0, 0.0) }
      )
      i = i + 1
    }
    out.push(self.project(coords))
    f = f + 1
  }
  out
}
