// ============================================================================
// B4: Helix Analysis — α / 3₁₀ / π 螺旋检测
// ============================================================================
//
// 基于 backbone φ/ψ 二面角：
//   α-helix:   φ ∈ [-150°, -30°],  ψ ∈ [-90°, 30°]     (≥ 3 残基)
//   3₁₀-helix: φ ∈ [-150°, -40°],  ψ ∈ [-80°, -10°]    (≥ 3 残基)
//   π-helix:   φ ∈ [-150°, -30°],  ψ ∈ [ 30°, 100°]    (≥ 5 残基)
// ============================================================================

///|
/// 检测指定帧的螺旋。
/// 返回 Array[(start_res_idx, length, helix_type)]。
pub fn AtomGroup::helix_detection(
  self : AtomGroup,
  frame~ : Int32 = 0,
) -> Array[(Int32, Int32, String)] {
  let atoms = self.universe.atoms
  let fr = self.universe.frame(frame)
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let n_res = starts.length()
  if n_res < 3 {
    return Array(capacity=0)
  }

  // 每残基 φ/ψ
  let phi : Array[Double?] = Array(capacity=n_res)
  let psi : Array[Double?] = Array(capacity=n_res)
  let mut ri = 0
  while ri < n_res {
    let cs = starts[ri]
    let ce = if ri + 1 < n_res { starts[ri + 1] } else { self.indices.length() }
    let na = mda_analysis_find_atom(atoms, self.indices, cs, ce, "N")
    let caa = mda_analysis_find_atom(atoms, self.indices, cs, ce, "CA")
    let ca2 = mda_analysis_find_atom(atoms, self.indices, cs, ce, "C")
    // φ
    match
      (
        if ri > 0 {
          mda_analysis_find_atom(atoms, self.indices, starts[ri - 1], cs, "C")
        } else {
          None
        },
        na,
        caa,
        ca2,
      ) {
      (Some(cp), Some(nd), Some(cad), Some(cc)) =>
        if cp < fr.length() &&
          nd < fr.length() &&
          cad < fr.length() &&
          cc < fr.length() {
          phi.push(
            Some(vector_dihedral_angle_degrees(fr[cp], fr[nd], fr[cad], fr[cc])),
          )
        } else {
          phi.push(None)
        }
      _ => phi.push(None)
    }
    // ψ
    match
      (
        na,
        caa,
        ca2,
        if ri + 1 < n_res {
          mda_analysis_find_atom(
            atoms,
            self.indices,
            starts[ri + 1],
            if ri + 2 < n_res {
              starts[ri + 2]
            } else {
              self.indices.length()
            },
            "N",
          )
        } else {
          None
        },
      ) {
      (Some(nd), Some(cad), Some(cc), Some(nn)) =>
        if nd < fr.length() &&
          cad < fr.length() &&
          cc < fr.length() &&
          nn < fr.length() {
          psi.push(
            Some(vector_dihedral_angle_degrees(fr[nd], fr[cad], fr[cc], fr[nn])),
          )
        } else {
          psi.push(None)
        }
      _ => psi.push(None)
    }
    ri = ri + 1
  }

  // 每残基螺旋类型（0=无, 1=3₁₀, 2=α, 3=π）
  let helix_flag : Array[Int32] = Array(capacity=n_res)
  ri = 0
  while ri < n_res {
    match (phi[ri], psi[ri]) {
      (Some(ph), Some(ps)) => {
        let t : Int32 = if ph >= -150.0 &&
          ph <= -30.0 &&
          ps >= -90.0 &&
          ps <= 30.0 {
          2
        } else if ph >= -150.0 && ph <= -40.0 && ps >= -80.0 && ps <= -10.0 {
          1
        } else if ph >= -150.0 && ph <= -30.0 && ps >= 30.0 && ps <= 100.0 {
          3
        } else {
          0
        }
        helix_flag.push(t)
      }
      _ => helix_flag.push(0)
    }
    ri = ri + 1
  }

  // 合并连续段
  let out : Array[(Int32, Int32, String)] = Array(capacity=5)
  ri = 0
  while ri < n_res {
    let t = helix_flag[ri]
    if t == 0 {
      ri = ri + 1
      continue
    }
    let start = ri
    let mut len = 1
    ri = ri + 1
    while ri < n_res && helix_flag[ri] == t {
      len = len + 1
      ri = ri + 1
    }
    let name = if t == 2 { "alpha" } else if t == 1 { "3_10" } else { "pi" }
    let ok = if t == 3 { len >= 5 } else { len >= 3 }
    if ok {
      out.push((start, len, name))
    }
  }
  out
}

///|
/// 螺旋比例时间序列（螺旋残基数 / 总残基数）。
pub fn AtomGroup::helix_fraction_trajectory(self : AtomGroup) -> Array[Double] {
  let nf = self.universe.n_frames()
  let atoms = self.universe.atoms
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let nr = starts.length()
  let out : Array[Double] = Array(capacity=nf)
  if nr < 3 {
    let mut f = 0
    while f < nf {
      out.push(0.0)
      f = f + 1
    }
    return out
  }
  let mut f = 0
  while f < nf {
    let hs = self.helix_detection(frame=f)
    let mut total = 0
    for h in hs {
      total = total + h.1
    }
    out.push(total.to_double() / nr.to_double())
    f = f + 1
  }
  out
}
