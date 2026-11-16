// ============================================================================
// B6: Transformations — align_frame + PBC unwrap
// B3: Order Parameter S²（deuterium 序参量）
// ============================================================================

///|
fn mda_vec3_to_arr(v : Vector3) -> Array[Double] {
  let out : Array[Double] = Array(capacity=3)
  out.push(v.x)
  out.push(v.y)
  out.push(v.z)
  out
}

///|
/// 对 AtomGroup 轨迹每帧做 Kabsch 对齐（质心平移到第 0 帧）。
/// 返回新 Universe。
pub fn AtomGroup::align_frames(self : AtomGroup) -> Universe {
  let nf = self.universe.n_frames()
  let refr = self.universe.frame(0)
  let idx = self.indices
  let n = idx.length()
  let out_frames : Array[Array[Array[Double]]] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    // 质心 ref
    let mut sxr = 0.0
    let mut syr = 0.0
    let mut szr = 0.0
    let mut sxm = 0.0
    let mut sym = 0.0
    let mut szm = 0.0
    let mut i = 0
    while i < n {
      if idx[i] < refr.length() {
        sxr = sxr + refr[idx[i]].x
        syr = syr + refr[idx[i]].y
        szr = szr + refr[idx[i]].z
      }
      if idx[i] < fr.length() {
        sxm = sxm + fr[idx[i]].x
        sym = sym + fr[idx[i]].y
        szm = szm + fr[idx[i]].z
      }
      i = i + 1
    }
    let nf2 = n.to_double()
    let cxr = sxr / nf2
    let cyr = syr / nf2
    let czr = szr / nf2
    let cxm = sxm / nf2
    let cym = sym / nf2
    let czm = szm / nf2
    // 对 self.indices 外的原子不动，只对 selected 做质心对齐
    let frame_out : Array[Array[Double]] = Array(capacity=fr.length())
    i = 0
    while i < fr.length() {
      let px = fr[i].x
      let py = fr[i].y
      let pz = fr[i].z
      let mut aligned_x = px
      let mut aligned_y = py
      let mut aligned_z = pz
      // 检查 i 是否在 idx 里
      let mut in_group = false
      let mut k = 0
      while k < n {
        if idx[k] == i {
          in_group = true
          break
        } else {
          k = k + 1
        }
      }
      if in_group {
        aligned_x = px - cxm + cxr
        aligned_y = py - cym + cyr
        aligned_z = pz - czm + czr
      }
      let arr : Array[Double] = Array(capacity=3)
      arr.push(aligned_x)
      arr.push(aligned_y)
      arr.push(aligned_z)
      frame_out.push(arr)
      i = i + 1
    }
    out_frames.push(frame_out)
    f = f + 1
  }
  Universe::from_dcd_frames(self.universe.atoms, out_frames)
}

///|
fn mda_wrap_coord2(v : Double, l : Double) -> Double {
  if l <= 0.0 {
    return v
  }
  let k = (v / l).to_int()
  v - k.to_double() * l
}

///|
/// PBC unwrap：把周期性边界包起来的原子展开到连续轨迹。
/// 返回新 Universe（所有原子轨迹展开，不只是 self.indices）。
pub fn AtomGroup::unwrap_trajectory(self : AtomGroup) -> Universe {
  let box = self.universe.box
  if box.v1.x <= 0.0 {
    return self.universe
  }
  let bx = box.v1.x
  let by = box.v2.y
  let bz = box.v3.z
  let nf = self.universe.n_frames()
  let na = self.universe.n_atoms()
  // 累计偏移
  let ox : Array[Double] = Array(capacity=na)
  let oy : Array[Double] = Array(capacity=na)
  let oz : Array[Double] = Array(capacity=na)
  let mut i = 0
  while i < na {
    ox.push(0.0)
    oy.push(0.0)
    oz.push(0.0)
    i = i + 1
  }
  let out_frames : Array[Array[Array[Double]]] = Array(capacity=nf)
  // frame 0：wrap 到 [0, box)，累计初始偏移
  let fr0 = self.universe.frame(0)
  let f0_out : Array[Array[Double]] = Array(capacity=fr0.length())
  i = 0
  while i < fr0.length() {
    let wx = mda_wrap_coord2(fr0[i].x, bx)
    let wy = mda_wrap_coord2(fr0[i].y, by)
    let wz = mda_wrap_coord2(fr0[i].z, bz)
    ox[i] = wx - fr0[i].x
    oy[i] = wy - fr0[i].y
    oz[i] = wz - fr0[i].z
    let a : Array[Double] = Array(capacity=3)
    a.push(wx)
    a.push(wy)
    a.push(wz)
    f0_out.push(a)
    i = i + 1
  }
  out_frames.push(f0_out)
  let mut f = 1
  while f < nf {
    let prev = self.universe.frame(f - 1)
    let curr = self.universe.frame(f)
    let fc_out : Array[Array[Double]] = Array(capacity=curr.length())
    i = 0
    while i < curr.length() {
      let cx = curr[i].x
      let cy = curr[i].y
      let cz = curr[i].z
      let px = if i < prev.length() { prev[i].x + ox[i] } else { cx }
      let py = if i < prev.length() { prev[i].y + oy[i] } else { cy }
      let pz = if i < prev.length() { prev[i].z + oz[i] } else { cz }
      let dx = cx - px
      let dy = cy - py
      let dz = cz - pz
      // 修正 PBC 跳跃
      if bx > 0.0 && dx > bx * 0.5 {
        ox[i] = ox[i] - bx
      } else if bx > 0.0 && dx < -bx * 0.5 {
        ox[i] = ox[i] + bx
      }
      if by > 0.0 && dy > by * 0.5 {
        oy[i] = oy[i] - by
      } else if by > 0.0 && dy < -by * 0.5 {
        oy[i] = oy[i] + by
      }
      if bz > 0.0 && dz > bz * 0.5 {
        oz[i] = oz[i] - bz
      } else if bz > 0.0 && dz < -bz * 0.5 {
        oz[i] = oz[i] + bz
      }
      let a : Array[Double] = Array(capacity=3)
      a.push(cx + ox[i])
      a.push(cy + oy[i])
      a.push(cz + oz[i])
      fc_out.push(a)
      i = i + 1
    }
    out_frames.push(fc_out)
    f = f + 1
  }
  Universe::from_dcd_frames(self.universe.atoms, out_frames, box~)
}

// ============================================================================
// B3: Order Parameter S² — deuterium NMR 序参量
// ============================================================================

///|
/// S² = ⟨(3cos²θ - 1)/2⟩，范围 [0, 1]。
/// 参考轴默认 z（膜法线方向）。
/// carbon_indices / hydrogen_indices 是 self.indices 中的**局部索引**。
pub fn AtomGroup::order_parameter_s2(
  self : AtomGroup,
  carbon_indices : Array[Int32],
  hydrogen_indices : Array[Int32],
  axis~ : String = "z",
) -> Double {
  let nb = carbon_indices.length()
  let nf = self.universe.n_frames()
  if nb != hydrogen_indices.length() || nb == 0 || nf == 0 {
    return 0.0
  }
  let ax_i : Int32 = if axis == "x" { 0 } else if axis == "y" { 1 } else { 2 }
  let mut sum_c2 = 0.0
  let mut tot = 0.0
  let mut f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    let mut b = 0
    while b < nb {
      let ic = self.indices[carbon_indices[b]]
      let ih = self.indices[hydrogen_indices[b]]
      if ic < fr.length() && ih < fr.length() {
        let cx = fr[ic].x
        let cy = fr[ic].y
        let cz = fr[ic].z
        let hx = fr[ih].x
        let hy = fr[ih].y
        let hz = fr[ih].z
        let dx = hx - cx
        let dy = hy - cy
        let dz = hz - cz
        let l = (dx * dx + dy * dy + dz * dz).sqrt()
        if l > 1.0e-6 {
          let comp = if ax_i == 0 {
            dx / l
          } else if ax_i == 1 {
            dy / l
          } else {
            dz / l
          }
          sum_c2 = sum_c2 + comp * comp
          tot = tot + 1.0
        }
      }
      b = b + 1
    }
    f = f + 1
  }
  if tot < 1.0e-6 {
    return 0.0
  }
  let avg = sum_c2 / tot
  (3.0 * avg - 1.0) / 2.0
}
