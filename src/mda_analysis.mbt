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
// A5: Distance / Angle / Dihedral 时间序列（AtomGroup 接口）
// ============================================================================

///|
/// 两个指定原子（self.indices[i] 和 self.indices[j]）的距离时间序列。
///
/// 返回 Array[Double]，长度 = n_frames。
/// i, j 是 AtomGroup 内的**局部索引**（0-based，不是 Universe 全局索引）。
pub fn AtomGroup::distance_trajectory(
  self : AtomGroup,
  i : Int32,
  j : Int32,
) -> Array[Double] {
  let n_frames = self.universe.n_frames()
  let idx_i = self.indices[i]
  let idx_j = self.indices[j]
  let out : Array[Double] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let di = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let dj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
    out.push(mda_analysis_dist2(di, dj))
    f = f + 1
  }
  out
}

///|
/// 三点夹角（j 为顶点）时间序列，单位度。
///
/// self.indices[i] — self.indices[j] — self.indices[k] 的夹角。
pub fn AtomGroup::angle_trajectory(
  self : AtomGroup,
  i : Int32,
  j : Int32,
  k : Int32,
) -> Array[Double] {
  let n_frames = self.universe.n_frames()
  let out : Array[Double] = Array(capacity=n_frames)
  let idx_i = self.indices[i]
  let idx_j = self.indices[j]
  let idx_k = self.indices[k]
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let pi = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
    let pj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
    let pk = if idx_k < fr.length() { fr[idx_k] } else { Vector3::new(0.0, 0.0, 0.0) }
    let v1 = Vector3::new(pi.x - pj.x, pi.y - pj.y, pi.z - pj.z)
    let v2 = Vector3::new(pk.x - pj.x, pk.y - pj.y, pk.z - pj.z)
    let dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
    let l1 = (v1.x * v1.x + v1.y * v1.y + v1.z * v1.z).sqrt()
    let l2 = (v2.x * v2.x + v2.y * v2.y + v2.z * v2.z).sqrt()
    let cos = if l1 > 1.0e-10 && l2 > 1.0e-10 {
      let c = dot / (l1 * l2)
      if c > 1.0 { 1.0 } else if c < -1.0 { -1.0 } else { c }
    } else { 0.0 }
    out.push(@math.acos(cos) * 180.0 / 3.141592653589793)
    f = f + 1
  }
  out
}

///|
/// 四原子二面角时间序列（度）。
///
/// 输入是 self.indices 中的**局部索引**。复用 pdb_vectors.mbt 的
/// `vector_dihedral_angle_degrees`（已 public）。
pub fn AtomGroup::dihedral_trajectory(
  self : AtomGroup,
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
) -> Array[Double] {
  let n_frames = self.universe.n_frames()
  let out : Array[Double] = Array(capacity=n_frames)
  let ia = self.indices[a]
  let ib = self.indices[b]
  let ic = self.indices[c]
  let id = self.indices[d]
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let pa = if ia < fr.length() { fr[ia] } else { Vector3::new(0.0, 0.0, 0.0) }
    let pb = if ib < fr.length() { fr[ib] } else { Vector3::new(0.0, 0.0, 0.0) }
    let pc = if ic < fr.length() { fr[ic] } else { Vector3::new(0.0, 0.0, 0.0) }
    let pd = if id < fr.length() { fr[id] } else { Vector3::new(0.0, 0.0, 0.0) }
    out.push(vector_dihedral_angle_degrees(pa, pb, pc, pd))
    f = f + 1
  }
  out
}

// ============================================================================
// Backbone phi/psi 时间序列（按残基自动识别）
// ============================================================================

///|
/// 内部：按 resid + chain 建立残基边界索引。
fn mda_analysis_residue_boundaries(
  atoms : Array[UnivAtom],
  indices : Array[Int32],
) -> (Array[Int32], Array[Int32]) {
  // starts[i] = AtomGroup 局部起始原子索引
  // resids_out[i] = 残基 ID
  let starts : Array[Int32] = Array(capacity=indices.length())
  let resids_out : Array[Int32] = Array(capacity=indices.length())
  let mut n_res = 0
  let mut cur_resid = -1
  let mut cur_chain = ""
  let mut cur_resname = ""
  let mut k = 0
  while k < indices.length() {
    let idx = indices[k]
    let a = if idx < atoms.length() { atoms[idx] } else { UnivAtom::new() }
    let is_new = (a.resid != cur_resid) || (a.chain != cur_chain) || (a.resname != cur_resname)
    if is_new {
      starts.push(k)
      resids_out.push(a.resid)
      cur_resid = a.resid
      cur_chain = a.chain
      cur_resname = a.resname
      n_res = n_res + 1
    }
    k = k + 1
  }
  (starts, resids_out)
}

///|
/// 在 AtomGroup 的一个残基边界内按名字找原子（返回 Universe 全局索引）。
fn mda_analysis_find_atom(
  atoms : Array[UnivAtom],
  indices : Array[Int32],
  res_start : Int32,   // AtomGroup 局部起始
  res_end : Int32,     // AtomGroup 局部结束（exclusive）
  name : String,
) -> Int32? {
  let mut k = res_start
  while k < res_end && k < indices.length() {
    let idx = indices[k]
    if idx < atoms.length() && atoms[idx].name.to_upper() == name.to_upper() {
      return Some(idx)
    }
    k = k + 1
  }
  None
}

///|
/// Backbone φ (phi) 时间序列 — C(i-1)-N(i)-CA(i)-C(i)。
///
/// 返回 Array[Array[Double?]]，shape `(n_res - 1, n_frames)`；
/// 首残基没有 φ；跨 chain 断裂点也为 None。
pub fn AtomGroup::backbone_phi_trajectory(self : AtomGroup) -> Array[Array[Double?]] {
  let atoms = self.universe.atoms
  let n_frames = self.universe.n_frames()
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let n_res = starts.length()
  let out : Array[Array[Double?]] = Array(capacity=n_res)
  let mut ri = 0
  while ri < n_res {
    let row : Array[Double?] = Array(capacity=n_frames)
    // 残基边界
    let cur_start = starts[ri]
    let cur_end = if ri + 1 < n_res { starts[ri + 1] } else { self.indices.length() }
    // 找 N(i), CA(i), C(i)
    let n_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "N")
    let ca_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "CA")
    let c_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "C")
    // 前一个残基的 C
    let c_prev = if ri > 0 {
      let prev_s = starts[ri - 1]
      let prev_e = cur_start
      mda_analysis_find_atom(atoms, self.indices, prev_s, prev_e, "C")
    } else { None }
    let mut f = 0
    while f < n_frames {
      match (c_prev, n_i, ca_i, c_i) {
        (Some(ic), Some(in1), Some(ica), Some(ic2)) => {
          let fr = self.universe.frame(f)
          let pc = if ic < fr.length() { fr[ic] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pn = if in1 < fr.length() { fr[in1] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pca = if ica < fr.length() { fr[ica] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pci = if ic2 < fr.length() { fr[ic2] } else { Vector3::new(0.0, 0.0, 0.0) }
          row.push(Some(vector_dihedral_angle_degrees(pc, pn, pca, pci)))
        }
        _ => row.push(None)
      }
      f = f + 1
    }
    out.push(row)
    ri = ri + 1
  }
  out
}

///|
/// Backbone ψ (psi) 时间序列 — N(i)-CA(i)-C(i)-N(i+1)。
pub fn AtomGroup::backbone_psi_trajectory(self : AtomGroup) -> Array[Array[Double?]] {
  let atoms = self.universe.atoms
  let n_frames = self.universe.n_frames()
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let n_res = starts.length()
  let out : Array[Array[Double?]] = Array(capacity=n_res)
  let mut ri = 0
  while ri < n_res {
    let row : Array[Double?] = Array(capacity=n_frames)
    let cur_start = starts[ri]
    let cur_end = if ri + 1 < n_res { starts[ri + 1] } else { self.indices.length() }
    let n_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "N")
    let ca_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "CA")
    let c_i = mda_analysis_find_atom(atoms, self.indices, cur_start, cur_end, "C")
    let n_next = if ri + 1 < n_res {
      let next_s = starts[ri + 1]
      let next_e = if ri + 2 < n_res { starts[ri + 2] } else { self.indices.length() }
      mda_analysis_find_atom(atoms, self.indices, next_s, next_e, "N")
    } else { None }
    let mut f = 0
    while f < n_frames {
      match (n_i, ca_i, c_i, n_next) {
        (Some(in1), Some(ica), Some(ic), Some(in2)) => {
          let fr = self.universe.frame(f)
          let pn = if in1 < fr.length() { fr[in1] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pca = if ica < fr.length() { fr[ica] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pci = if ic < fr.length() { fr[ic] } else { Vector3::new(0.0, 0.0, 0.0) }
          let pnn = if in2 < fr.length() { fr[in2] } else { Vector3::new(0.0, 0.0, 0.0) }
          row.push(Some(vector_dihedral_angle_degrees(pn, pca, pci, pnn)))
        }
        _ => row.push(None)
      }
      f = f + 1
    }
    out.push(row)
    ri = ri + 1
  }
  out
}

// ============================================================================
// A6: Radial Distribution Function（Vector3 版，支持多帧平均）
// ============================================================================

///|
/// RDF 结果（与 structure_rdf.mbt 的 RDFResult 字段一致，类型简化）。
pub struct MdaRdfResult {
  r : Array[Double]
  g_r : Array[Double]
  n_bins : Int32
  dr : Double
}

///|
/// 两群原子的径向分布函数 g(r) — 单帧。
///
/// 可选 PBC（orthogonal box + minimum image convention）。
/// rmax 默认 12 Å，nbins 默认 100。
pub fn AtomGroup::rdf(
  self : AtomGroup,
  other : AtomGroup,
  rmax~ : Double = 12.0,
  nbins~ : Int32 = 100,
  frame~ : Int32 = 0,
  pbc~ : Bool = false,
) -> MdaRdfResult {
  let fr_a = self.universe.frame(frame)
  let fr_b = other.universe.frame(frame)
  let na = self.indices.length()
  let nb = other.indices.length()
  let dr = rmax / nbins.to_double()
  let box = if pbc { Some(self.universe.box) } else { None }
  let counts : Array[Int32] = Array(capacity=nbins)
  let mut b = 0
  while b < nbins { counts.push(0); b = b + 1 }
  let mut i = 0
  while i < na {
    let idx_a = self.indices[i]
    let ca = if idx_a < fr_a.length() { fr_a[idx_a] } else { Vector3::new(0.0, 0.0, 0.0) }
    let mut j = 0
    while j < nb {
      let idx_b = other.indices[j]
      let cb = if idx_b < fr_b.length() { fr_b[idx_b] } else { Vector3::new(0.0, 0.0, 0.0) }
      let d = if pbc {
        mda_analysis_min_dist(ca, cb, box)
      } else {
        mda_analysis_dist2(ca, cb)
      }
      if d > 0.0 && d < rmax {
        let bin = (d / dr).to_int()
        if bin >= 0 && bin < nbins {
          counts[bin] = counts[bin] + 1
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  // 归一化：g(r) = counts / (na * nb * shell_volume * rho_bulk)
  let box_vol = if pbc && self.universe.box.v1.x > 0.0 && self.universe.box.v2.y > 0.0 && self.universe.box.v3.z > 0.0 {
    self.universe.box.v1.x * self.universe.box.v2.y * self.universe.box.v3.z
  } else {
    let cx = (na + nb).to_double().sqrt() * 20.0
    cx * cx * cx
  }
  let rho = nb.to_double() / box_vol
  let g_r : Array[Double] = Array(capacity=nbins)
  let r : Array[Double] = Array(capacity=nbins)
  let mut k = 0
  while k < nbins {
    let r_lo = k.to_double() * dr
    let r_hi = r_lo + dr
    let shell_vol = (4.0 / 3.0) * 3.141592653589793 * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo)
    let ideal = na.to_double() * rho * shell_vol
    g_r.push(if ideal > 1.0e-10 { counts[k].to_double() / ideal } else { 0.0 })
    r.push((r_lo + r_hi) / 2.0)
    k = k + 1
  }
  MdaRdfResult::{ r, g_r, n_bins: nbins, dr }
}

///|
/// 整个轨迹上的 RDF 平均（降低单帧噪声）。
///
/// 逐帧累加 counts，最后一起归一化。比每帧算完平均更精确。
pub fn AtomGroup::rdf_averaged(
  self : AtomGroup,
  other : AtomGroup,
  rmax~ : Double = 12.0,
  nbins~ : Int32 = 100,
  pbc~ : Bool = false,
) -> MdaRdfResult {
  let n_frames = self.universe.n_frames()
  let na = self.indices.length()
  let nb = other.indices.length()
  let dr = rmax / nbins.to_double()
  let box = if pbc { Some(self.universe.box) } else { None }
  let total_counts : Array[Double] = Array(capacity=nbins)
  let mut b = 0
  while b < nbins { total_counts.push(0.0); b = b + 1 }
  let mut f = 0
  while f < n_frames {
    let fr_a = self.universe.frame(f)
    let fr_b = other.universe.frame(f)
    let mut i = 0
    while i < na {
      let idx_a = self.indices[i]
      let ca = if idx_a < fr_a.length() { fr_a[idx_a] } else { Vector3::new(0.0, 0.0, 0.0) }
      let mut j = 0
      while j < nb {
        let idx_b = other.indices[j]
        let cb = if idx_b < fr_b.length() { fr_b[idx_b] } else { Vector3::new(0.0, 0.0, 0.0) }
        let d = if pbc {
          mda_analysis_min_dist(ca, cb, box)
        } else {
          mda_analysis_dist2(ca, cb)
        }
        if d > 0.0 && d < rmax {
          let bin = (d / dr).to_int()
          if bin >= 0 && bin < nbins {
            total_counts[bin] = total_counts[bin] + 1.0
          }
        }
        j = j + 1
      }
      i = i + 1
    }
    f = f + 1
  }
  let box_vol = if pbc && self.universe.box.v1.x > 0.0 && self.universe.box.v2.y > 0.0 && self.universe.box.v3.z > 0.0 {
    self.universe.box.v1.x * self.universe.box.v2.y * self.universe.box.v3.z
  } else {
    let cx = (na + nb).to_double().sqrt() * 20.0
    cx * cx * cx
  }
  let rho = nb.to_double() / box_vol
  let g_r : Array[Double] = Array(capacity=nbins)
  let r : Array[Double] = Array(capacity=nbins)
  let nf = n_frames.to_double()
  let mut k = 0
  while k < nbins {
    let r_lo = k.to_double() * dr
    let r_hi = r_lo + dr
    let shell_vol = (4.0 / 3.0) * 3.141592653589793 * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo)
    let ideal = na.to_double() * rho * shell_vol * nf
    g_r.push(if ideal > 1.0e-10 { total_counts[k] / ideal } else { 0.0 })
    r.push((r_lo + r_hi) / 2.0)
    k = k + 1
  }
  MdaRdfResult::{ r, g_r, n_bins: nbins, dr }
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
