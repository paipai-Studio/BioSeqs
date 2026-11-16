// ============================================================================
// C1: Ramachandran — 轨迹版 phi/psi 散点 + 多帧叠加密度 + 分类
// ============================================================================

pub struct RamachandranResult {
  /// 残基 × 帧，phi 值（度，None 表示无 backbone 原子）
  phi_frames : Array[Array[Double?]]
  /// 残基 × 帧，psi 值（度）
  psi_frames : Array[Array[Double?]]
  /// 残基名 + resid 标签
  res_labels : Array[String]
  /// 帧 × 残基的区域分类（0=coil/general, 1=alpha, 2=beta, 3=polyproline, -1=disallowed）
  region_frames : Array[Array[Int32]]
}

fn mda_rama_classify(phi : Double, psi : Double) -> Int32 {
  // 经典 Ramachandran 区域（简化版）
  // α-helix        φ∈[-90,-30]  ψ∈[-70,-10]
  if phi >= -90.0 && phi <= -30.0 && psi >= -70.0 && psi <= -10.0 { return 1 }
  // β-sheet        φ∈[-160,-100] ψ∈[90,180] 或 ψ∈[-180,-120]
  if phi >= -160.0 && phi <= -100.0 && ((psi >= 90.0 && psi <= 180.0) || (psi >= -180.0 && psi <= -120.0)) { return 2 }
  // left-handed α  φ∈[30,90]  ψ∈[-30,70]
  if phi >= 30.0 && phi <= 90.0 && psi >= -30.0 && psi <= 70.0 { return 1 }
  // polyproline II φ∈[-90,-40]  ψ∈[120,180]
  if phi >= -90.0 && phi <= -40.0 && psi >= 120.0 && psi <= 180.0 { return 3 }
  return 0
}

fn mda_rama_region_name(r : Int32) -> String {
  if r == 1 { "alpha" } else if r == 2 { "beta" } else if r == 3 { "ppII" } else if r < 0 { "disallowed" } else { "coil" }
}

///|
/// AtomGroup 全轨迹 Ramachandran 分析。
/// 自动按残基边界计算每残基每帧的 φ/ψ 二面角。
pub fn AtomGroup::ramachandran_trajectory(self : AtomGroup) -> RamachandranResult {
  let atoms = self.universe.atoms
  let nf = self.universe.n_frames()
  let (starts, resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let nr = starts.length()
  let labels : Array[String] = Array(capacity=nr)
  let mut i = 0
  while i < nr {
    let rid = if i < resids.length() { resids[i] } else { i + 1 }
    labels.push(rid.to_string())
    i = i + 1
  }
  let phi_frames : Array[Array[Double?]] = Array(capacity=nr)
  let psi_frames : Array[Array[Double?]] = Array(capacity=nr)
  let region_frames : Array[Array[Int32]] = Array(capacity=nr)
  let mut ri = 0
  while ri < nr {
    let cs = starts[ri]
    let ce = if ri + 1 < nr { starts[ri + 1] } else { self.indices.length() }
    let na = mda_analysis_find_atom(atoms, self.indices, cs, ce, "N")
    let caa = mda_analysis_find_atom(atoms, self.indices, cs, ce, "CA")
    let ca2 = mda_analysis_find_atom(atoms, self.indices, cs, ce, "C")
    let c_prev = if ri > 0 { mda_analysis_find_atom(atoms, self.indices, starts[ri - 1], cs, "C") } else { None }
    let n_next = if ri + 1 < nr {
      mda_analysis_find_atom(atoms, self.indices, starts[ri + 1],
        if ri + 2 < nr { starts[ri + 2] } else { self.indices.length() }, "N")
    } else { None }
    let phi_row : Array[Double?] = Array(capacity=nf)
    let psi_row : Array[Double?] = Array(capacity=nf)
    let reg_row : Array[Int32] = Array(capacity=nf)
    let mut f = 0
    while f < nf {
      let fr = self.universe.frame(f)
      // φ = C(i-1)-N(i)-CA(i)-C(i)
      let ph : Double? = match (c_prev, na, caa, ca2) {
        (Some(cp), Some(nd), Some(cad), Some(cc)) => {
          if cp < fr.length() && nd < fr.length() && cad < fr.length() && cc < fr.length() {
            Some(vector_dihedral_angle_degrees(fr[cp], fr[nd], fr[cad], fr[cc]))
          } else { None }
        }
        _ => None
      }
      // ψ = N(i)-CA(i)-C(i)-N(i+1)
      let ps : Double? = match (na, caa, ca2, n_next) {
        (Some(nd), Some(cad), Some(cc), Some(nn)) => {
          if nd < fr.length() && cad < fr.length() && cc < fr.length() && nn < fr.length() {
            Some(vector_dihedral_angle_degrees(fr[nd], fr[cad], fr[cc], fr[nn]))
          } else { None }
        }
        _ => None
      }
      phi_row.push(ph); psi_row.push(ps)
      match (ph, ps) {
        (Some(pp), Some(ss)) => reg_row.push(mda_rama_classify(pp, ss))
        _ => reg_row.push(-1)
      }
      f = f + 1
    }
    phi_frames.push(phi_row); psi_frames.push(psi_row); region_frames.push(reg_row)
    ri = ri + 1
  }
  RamachandranResult::{ phi_frames, psi_frames, res_labels: labels, region_frames }
}

///|
/// Ramachandran 2D 散点密度（phi × psi 网格，用于可视化）。
///
/// n_bins 默认 36（10° 分辨率）。返回 (n_bins × n_bins) 概率密度矩阵。
pub fn RamachandranResult::density_grid(self : RamachandranResult, n_bins~ : Int32 = 36) -> Array[Array[Double]] {
  let total = self.phi_frames.length() * (if self.phi_frames.length() > 0 { self.phi_frames[0].length() } else { 0 })
  if total == 0 { return Array(capacity=n_bins) }
  let hist : Array[Array[Int32]] = Array(capacity=n_bins)
  let mut i = 0
  while i < n_bins {
    let row : Array[Int32] = Array(capacity=n_bins)
    let mut j = 0
    while j < n_bins { row.push(0); j = j + 1 }
    hist.push(row); i = i + 1
  }
  let bin = 360.0 / n_bins.to_double()
  let mut n_total = 0
  let n_res = self.phi_frames.length()
  let nf = if n_res > 0 { self.phi_frames[0].length() } else { 0 }
  let mut ri = 0
  while ri < n_res {
    let mut f = 0
    while f < nf {
      match (self.phi_frames[ri][f], self.psi_frames[ri][f]) {
        (Some(pp), Some(ss)) => {
          let px = (pp + 180.0).to_int() / bin.to_int()
          let py = (ss + 180.0).to_int() / bin.to_int()
          let bx = if px < 0 { 0 } else if px >= n_bins { n_bins - 1 } else { px }
          let by = if py < 0 { 0 } else if py >= n_bins { n_bins - 1 } else { py }
          hist[by][bx] = hist[by][bx] + 1
          n_total = n_total + 1
        }
        _ => { let _ = () }
      }
      f = f + 1
    }
    ri = ri + 1
  }
  let out : Array[Array[Double]] = Array(capacity=n_bins)
  i = 0
  while i < n_bins {
    let row : Array[Double] = Array(capacity=n_bins)
    let mut j = 0
    while j < n_bins {
      row.push(if n_total > 0 { hist[i][j].to_double() / n_total.to_double() } else { 0.0 })
      j = j + 1
    }
    out.push(row); i = i + 1
  }
  out
}

///|
/// 区域分布（整个轨迹 alpha/beta/ppII/coil 的残基-帧对比例）。
pub fn RamachandranResult::region_fractions(self : RamachandranResult) -> Array[Double] {
  let counts : Array[Int32] = [0, 0, 0, 0, 0] // alpha, beta, ppII, coil, disallowed
  let mut total = 0
  let n_res = self.region_frames.length()
  let nf = if n_res > 0 { self.region_frames[0].length() } else { 0 }
  let mut ri = 0
  while ri < n_res {
    let mut f = 0
    while f < nf {
      let r = self.region_frames[ri][f]
      if r >= 0 {
        counts[r] = counts[r] + 1
        total = total + 1
      } else {
        counts[4] = counts[4] + 1; total = total + 1
      }
      f = f + 1
    }
    ri = ri + 1
  }
  let out : Array[Double] = Array(capacity=5)
  let t = total.to_double()
  let mut i = 0
  while i < 5 { out.push(if t > 0.0 { counts[i].to_double() / t } else { 0.0 }); i = i + 1 }
  out
}
