// ============================================================================
// B1: Salt Bridge Analysis — 正-负残基间 4.5 Å 接触
// ============================================================================

///|
fn mda_is_pos_res(name : String) -> Bool {
  let u = name.to_upper()
  u == "LYS" || u == "ARG" || u == "HIS"
}

///|
fn mda_is_neg_res(name : String) -> Bool {
  let u = name.to_upper()
  u == "ASP" || u == "GLU"
}

///|
/// 指定帧内检测盐桥（正-负残基对，最近原子距离 ≤ cutoff Å）。
/// 返回 Array[(res_i, res_j, min_distance)]。
pub fn AtomGroup::salt_bridges(
  self : AtomGroup,
  frame~ : Int32 = 0,
  cutoff~ : Double = 4.5,
) -> Array[(Int32, Int32, Double)] {
  let atoms = self.universe.atoms
  let fr = self.universe.frame(frame)
  let n = self.indices.length()
  let cutoff2 = cutoff * cutoff
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let n_res = starts.length()
  let out : Array[(Int32, Int32, Double)] = Array(capacity=10)
  let mut i = 0
  while i < n_res {
    let cs = starts[i]
    let ce = if i + 1 < n_res { starts[i + 1] } else { n }
    let rn = if cs < n && self.indices[cs] < atoms.length() {
      atoms[self.indices[cs]].resname.to_upper()
    } else {
      ""
    }
    let pos_i = mda_is_pos_res(rn)
    let neg_i = mda_is_neg_res(rn)
    if !pos_i && !neg_i {
      i = i + 1
      continue
    }
    let mut j = i + 1
    while j < n_res {
      let ns = starts[j]
      let ne = if j + 1 < n_res { starts[j + 1] } else { n }
      let rn2 = if ns < n && self.indices[ns] < atoms.length() {
        atoms[self.indices[ns]].resname.to_upper()
      } else {
        ""
      }
      let pos_j = mda_is_pos_res(rn2)
      let neg_j = mda_is_neg_res(rn2)
      if !((pos_i && neg_j) || (neg_i && pos_j)) {
        j = j + 1
        continue
      }
      let mut min_d2 = 1.7976931348623157e308
      let mut a = cs
      while a < ce {
        let ia = self.indices[a]
        let ca = if ia < fr.length() {
          fr[ia]
        } else {
          Vector3::new(0.0, 0.0, 0.0)
        }
        let mut b = ns
        while b < ne {
          let ib = self.indices[b]
          let cb = if ib < fr.length() {
            fr[ib]
          } else {
            Vector3::new(0.0, 0.0, 0.0)
          }
          let dx = ca.x - cb.x
          let dy = ca.y - cb.y
          let dz = ca.z - cb.z
          let d2 = dx * dx + dy * dy + dz * dz
          if d2 < min_d2 {
            min_d2 = d2
          }
          b = b + 1
        }
        a = a + 1
      }
      if min_d2 <= cutoff2 {
        out.push((i, j, min_d2.sqrt()))
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

///|
/// 盐桥数量时间序列。
pub fn AtomGroup::salt_bridge_trajectory(
  self : AtomGroup,
  cutoff~ : Double = 4.5,
) -> Array[Int32] {
  let nf = self.universe.n_frames()
  let out : Array[Int32] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    out.push(self.salt_bridges(frame=f, cutoff~).length())
    f = f + 1
  }
  out
}
