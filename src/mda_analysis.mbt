// ============================================================================
// MDAnalysis 风格分析方法层 — 基础分析
// ============================================================================
//
// 这里给 AtomGroup 加方法，复用 BioSeqs 已有的底层工具：
//   - RMSD / RMSF / mean_structure    → 自己实现（Vector3 输入，复用 PBC helper）
//   - distance_matrix                → O(n²) 两两距离，支持 PBC
//   - radius_of_gyration_trajectory → 逐帧回转半径
//   - MSD                             → Einstein relation ⟨|r_i(t) - r_i(0)|²⟩
//   - SASA（简化）                    → 球近似
// ============================================================================

///|
/// 原子-原子距离（非 PBC）。
fn mda_analysis_dist2(a : Vector3, b : Vector3) -> Double {
  let dx = a.x - b.x
  let dy = a.y - b.y
  let dz = a.z - b.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// 最小镜像距离（orthogonal PBC box）。
fn mda_analysis_min_dist(a : Vector3, b : Vector3, box : PBCBox?) -> Double {
  let dx0 = a.x - b.x
  let dy0 = a.y - b.y
  let dz0 = a.z - b.z
  match box {
    Some(bx) => {
      let ax = bx.v1.x
      let by = bx.v2.y
      let cz = bx.v3.z
      let dx = if ax > 0.0 {
        mda_analysis_min_image(dx0, ax)
      } else { dx0 }
      let dy = if by > 0.0 {
        mda_analysis_min_image(dy0, by)
      } else { dy0 }
      let dz = if cz > 0.0 {
        mda_analysis_min_image(dz0, cz)
      } else { dz0 }
      (dx * dx + dy * dy + dz * dz).sqrt()
    }
    None => mda_analysis_dist2(a, b)
  }
}

fn mda_analysis_min_image(d : Double, l : Double) -> Double {
  if l <= 0.0 { return d }
  let k = (d / l + 0.5).floor()
  d - k * l
}

// ============================================================================
// AtomGroup → 轨迹级分析
// ============================================================================

///|
/// 两帧间 RMSD（无对齐，直接坐标差）。
pub fn AtomGroup::rmsd_between(
  self : AtomGroup,
  frame_a : Int32,
  frame_b : Int32,
) -> Double {
  let n = self.indices.length()
  let fr_a = self.universe.frame(frame_a)
  let fr_b = self.universe.frame(frame_b)
  if fr_a.length() == 0 || fr_b.length() == 0 || n == 0 { return 0.0 }
  let mut ssd = 0.0
  let mut i = 0
  while i < n {
    let idx = self.indices[i]
    if idx < fr_a.length() && idx < fr_b.length() {
      let dx = fr_a[idx].x - fr_b[idx].x
      let dy = fr_a[idx].y - fr_b[idx].y
      let dz = fr_a[idx].z - fr_b[idx].z
      ssd = ssd + dx * dx + dy * dy + dz * dz
    }
    i = i + 1
  }
  (ssd / n.to_double()).sqrt()
}

///|
/// 每帧相对于第 0 帧的 RMSD 轨迹。
pub fn AtomGroup::rmsd_trajectory(self : AtomGroup) -> Array[Double] {
  let n_frames = self.universe.n_frames()
  let result : Array[Double] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    result.push(self.rmsd_between(0, f))
    f = f + 1
  }
  result
}

///|
/// 平均结构（所有帧的原子坐标均值）。
///
/// 返回 `Array[Vector3]`，长度等于组内原子数。
pub fn AtomGroup::mean_structure(self : AtomGroup) -> Array[Vector3] {
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let out : Array[Vector3] = Array(capacity=n_atoms)
  if n_atoms == 0 || n_frames == 0 { return out }
  let nf = n_frames.to_double()
  let mut i = 0
  while i < n_atoms {
    let idx = self.indices[i]
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    let mut f = 0
    while f < n_frames {
      let fr = self.universe.frame(f)
      if idx < fr.length() {
        sx = sx + fr[idx].x
        sy = sy + fr[idx].y
        sz = sz + fr[idx].z
      }
      f = f + 1
    }
    out.push(Vector3::new(sx / nf, sy / nf, sz / nf))
    i = i + 1
  }
  out
}

///|
/// 每帧相对于 mean structure 的 RMSD 轨迹。
pub fn AtomGroup::rmsd_to_mean(self : AtomGroup) -> Array[Double] {
  let mean = self.mean_structure()
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let result : Array[Double] = Array(capacity=n_frames)
  if n_atoms == 0 || n_frames == 0 { return result }
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let mut ssd = 0.0
    let mut i = 0
    while i < n_atoms {
      let idx = self.indices[i]
      if idx < fr.length() {
        let dx = fr[idx].x - mean[i].x
        let dy = fr[idx].y - mean[i].y
        let dz = fr[idx].z - mean[i].z
        ssd = ssd + dx * dx + dy * dy + dz * dz
      }
      i = i + 1
    }
    result.push((ssd / n_atoms.to_double()).sqrt())
    f = f + 1
  }
  result
}

///|
/// 每原子的 RMSF（原子位置相对均值的波动，所有帧的 mean-square-fluctuation）。
pub fn AtomGroup::rmsf(self : AtomGroup) -> Array[Double] {
  let mean = self.mean_structure()
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let out : Array[Double] = Array(capacity=n_atoms)
  if n_atoms == 0 || n_frames == 0 { return out }
  let mut i = 0
  while i < n_atoms {
    let idx = self.indices[i]
    let mut msf = 0.0
    let mut f = 0
    while f < n_frames {
      let fr = self.universe.frame(f)
      if idx < fr.length() {
        let dx = fr[idx].x - mean[i].x
        let dy = fr[idx].y - mean[i].y
        let dz = fr[idx].z - mean[i].z
        msf = msf + dx * dx + dy * dy + dz * dz
      }
      f = f + 1
    }
    out.push((msf / n_frames.to_double()).sqrt())
    i = i + 1
  }
  out
}

// ============================================================================
// 距离矩阵 & 接触
// ============================================================================

///|
/// 当前帧下原子两两距离矩阵（n×n，下三角镜像）。
///
/// 可选 PBC。`cutoff` > 0 时返回 `Option` 矩阵（None = 超过 cutoff）。
pub fn AtomGroup::distance_matrix(
  self : AtomGroup,
  frame~ : Int32 = 0,
  pbc~ : Bool = false,
) -> Array[Array[Double]] {
  let n = self.indices.length()
  let fr = self.universe.frame(frame)
  let box = if pbc { Some(self.universe.box) } else { None }
  let out : Array[Array[Double]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array(capacity=n)
    let idx_i = self.indices[i]
    let ci = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let mut j = 0
    while j < n {
      let idx_j = self.indices[j]
      let cj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
      if pbc {
        row.push(mda_analysis_min_dist(ci, cj, box))
      } else {
        row.push(mda_analysis_dist2(ci, cj))
      }
      j = j + 1
    }
    out.push(row)
    i = i + 1
  }
  out
}

///|
/// 某帧内的接触原子对（cutoff 距离内）。
///
/// 返回 `Array[(Int32, Int32)]`，每个元素是 (idx_i, idx_j) 且 i < j。
/// 可选 PBC。
pub fn AtomGroup::contacts(
  self : AtomGroup,
  cutoff : Double,
  frame~ : Int32 = 0,
  pbc~ : Bool = false,
) -> Array[(Int32, Int32)] {
  let n = self.indices.length()
  let fr = self.universe.frame(frame)
  let box = if pbc { Some(self.universe.box) } else { None }
  let out : Array[(Int32, Int32)] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let idx_i = self.indices[i]
    let ci = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let mut j = i + 1
    while j < n {
      let idx_j = self.indices[j]
      let cj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
      let d = if pbc {
        mda_analysis_min_dist(ci, cj, box)
      } else {
        mda_analysis_dist2(ci, cj)
      }
      if d <= cutoff {
        out.push((idx_i, idx_j))
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

///|
/// 两 AtomGroup 之间的接触原子对（A 中某个原子与 B 中某个原子距离 < cutoff）。
///
/// 可选 PBC。
pub fn AtomGroup::contacts_with(
  self : AtomGroup,
  other : AtomGroup,
  cutoff : Double,
  frame~ : Int32 = 0,
  pbc~ : Bool = false,
) -> Array[(Int32, Int32)] {
  let n_a = self.indices.length()
  let n_b = other.indices.length()
  let fr = self.universe.frame(frame)
  let box = if pbc { Some(self.universe.box) } else { None }
  let out : Array[(Int32, Int32)] = Array(capacity=n_a)
  let mut i = 0
  while i < n_a {
    let idx_i = self.indices[i]
    let ci = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let mut j = 0
    while j < n_b {
      let idx_j = other.indices[j]
      let cj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
      let d = if pbc {
        mda_analysis_min_dist(ci, cj, box)
      } else {
        mda_analysis_dist2(ci, cj)
      }
      if d <= cutoff {
        out.push((idx_i, idx_j))
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

// ============================================================================
// 轨迹级的回转半径 + MSD
// ============================================================================

///|
/// 每帧回转半径轨迹（加权 mass，默认）。
pub fn AtomGroup::radius_of_gyration_trajectory(self : AtomGroup) -> Array[Double] {
  let n_frames = self.universe.n_frames()
  let result : Array[Double] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    result.push(self.radius_of_gyration(frame=f))
    f = f + 1
  }
  result
}

///|
/// 均方位移（MSD）相对于第 0 帧。
///
/// 对每组原子计算 ⟨|r_i(t) - r_i(0)|²⟩（时间平均），结果长度 = n_frames。
pub fn AtomGroup::msd(self : AtomGroup) -> Array[Double] {
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let out : Array[Double] = Array(capacity=n_frames)
  if n_atoms == 0 || n_frames < 2 {
    let mut f = 0
    while f < n_frames { out.push(0.0); f = f + 1 }
    return out
  }
  // 第 0 帧参考坐标
  let fr0 = self.universe.frame(0)
  // 预存每个原子的参考坐标
  let ref_coords : Array[Vector3] = Array(capacity=n_atoms)
  let mut i = 0
  while i < n_atoms {
    let idx = self.indices[i]
    ref_coords.push(
      if idx < fr0.length() { fr0[idx] } else { Vector3::new(0.0, 0.0, 0.0) }
    )
    i = i + 1
  }
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let mut total = 0.0
    let mut i = 0
    while i < n_atoms {
      let idx = self.indices[i]
      if idx < fr.length() {
        let dx = fr[idx].x - ref_coords[i].x
        let dy = fr[idx].y - ref_coords[i].y
        let dz = fr[idx].z - ref_coords[i].z
        total = total + dx * dx + dy * dy + dz * dz
      }
      i = i + 1
    }
    out.push(total / n_atoms.to_double())
    f = f + 1
  }
  out
}

// ============================================================================
// 简单 SASA（基于 probe radius 的原子-原子接触近似）
// ============================================================================

///|
/// 某帧 SASA（简化版：radius = vdw_radius + probe）。
///
/// 用 Shrake-Rupley 思路的近似：原子被邻居原子遮挡的部分不暴露。
/// 不做真正的球面网格化，只做 pairwise van-der-Waals overlap 计算。
/// 结果是每个原子的近似暴露表面积（Å²）。
pub fn AtomGroup::sasa_per_atom(
  self : AtomGroup,
  frame~ : Int32 = 0,
  probe~ : Double = 1.4,
) -> Array[Double] {
  let n = self.indices.length()
  let fr = self.universe.frame(frame)
  if n == 0 { return Array(capacity=0) }
  // vdW radii（近似，常见元素）
  let vdw : Array[Double] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let idx = self.indices[i]
    let radius = if idx < self.universe.atoms.length() {
      let el = self.universe.atoms[idx].element.to_upper()
      mda_vdw_radius(el)
    } else { 1.7 }
    vdw.push(radius + probe)
    i = i + 1
  }
  // 对每个原子，计算它的邻居遮挡比例
  let out : Array[Double] = Array(capacity=n)
  i = 0
  while i < n {
    let idx_i = self.indices[i]
    let ci = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let ri = vdw[i]
    let area = 4.0 * 3.141592653589793 * ri * ri
    // 检查每个邻居原子，若 overlap 则减去遮挡的球面帽面积
    let mut exposed = 1.0
    let mut j = 0
    while j < n {
      if i != j {
        let idx_j = self.indices[j]
        let cj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
        let d = mda_analysis_dist2(ci, cj)
        let rj = vdw[j]
        // 两球重叠的球面帽面积（仅近似扣除遮挡的一半）
        if d < ri + rj && d > 0.001 {
          let cos_theta = (ri * ri + d * d - rj * rj) / (2.0 * ri * d)
          if cos_theta > -1.0 && cos_theta < 1.0 {
            let hat_area = 2.0 * 3.141592653589793 * ri * ri * (1.0 - cos_theta)
            // 近似：每个邻居最多遮挡一半（避免超扣）
            let frac = if hat_area / area < 0.5 { hat_area / area } else { 0.5 }
            exposed = exposed - frac
          }
        }
      }
      j = j + 1
    }
    out.push((if exposed > 0.0 { exposed } else { 0.0 }) * area)
    i = i + 1
  }
  out
}

///|
/// 某帧总 SASA。
pub fn AtomGroup::sasa(self : AtomGroup, frame~ : Int32 = 0, probe~ : Double = 1.4) -> Double {
  let per_atom = self.sasa_per_atom(frame=frame, probe=probe)
  let mut total = 0.0
  let mut i = 0
  while i < per_atom.length() {
    total = total + per_atom[i]
    i = i + 1
  }
  total
}

// ============================================================================
// vdW radii 表（Bondi 1964）
// ============================================================================

fn mda_vdw_radius(element : String) -> Double {
  match element {
    "H"  => 1.2
    "C"  => 1.7
    "N"  => 1.55
    "O"  => 1.52
    "F"  => 1.47
    "P"  => 1.8
    "S"  => 1.8
    "CL" => 1.75
    "BR" => 1.85
    "I"  => 1.98
    "FE" => 1.94
    "ZN" => 1.39
    "CA" => 2.31
    "NA" => 2.27
    "MG" => 1.73
    // 默认碳
    _    => 1.7
  }
}
