// ============================================================================
// MDAnalysis 风格：氢键检测（几何准则）
// ============================================================================
//
// 准则（Baker & Hubbard 1984, biotite 同款）：
//   D-A   距离 ≤ max_distance（默认 3.5 Å）
//   H-A   距离 ≤ max_h_distance（默认 2.5 Å）
//   D-H-A 角度 ≥ min_angle（默认 120°）
//
// donor_type / acceptor_type 白名单（标准蛋白）：
//   donors:   N, NE, NH, NZ, OG, OH, NE2, ND2, SD
//   acceptors: O, OD, OG, OE, OD1, OD2, OE1, OE2
// ============================================================================

///|
/// 检测标准蛋白残基中可能作为 HBond donor 的原子名。
pub fn mda_hbond_is_donor(name : String) -> Bool {
  match name.to_upper() {
    "N" | "NE" | "NH" | "NZ" | "OG" | "OH" | "NE2" | "ND2" | "SD" => true
    _ => false
  }
}

///|
/// 检测标准蛋白残基中可能作为 HBond acceptor 的原子名。
pub fn mda_hbond_is_acceptor(name : String) -> Bool {
  match name.to_upper() {
    "O" | "OD" | "OG" | "OE" | "OD1" | "OD2" | "OE1" | "OE2" => true
    _ => false
  }
}

///|
/// 检测 HBond 几何：返回 Some(dist, angle) 或 None。
///
/// 需要 H 原子坐标（若不显式提供，用 D-H-A 近似角度）。
fn mda_hbond_geom(
  donor : Vector3,
  acceptor : Vector3,
  hydrogen : Vector3?,
  max_da : Double,
  max_ha : Double,
  min_angle : Double,
) -> (Double, Double)? {
  let dx = donor.x - acceptor.x
  let dy = donor.y - acceptor.y
  let dz = donor.z - acceptor.z
  let da = (dx * dx + dy * dy + dz * dz).sqrt()
  if da > max_da { return None }
  match hydrogen {
    Some(h) => {
      let hdx = h.x - acceptor.x
      let hdy = h.y - acceptor.y
      let hdz = h.z - acceptor.z
      let ha = (hdx * hdx + hdy * hdy + hdz * hdz).sqrt()
      if ha > max_ha { return None }
      // D-H-A 角度（在 H 处）
      let hh_dx = donor.x - h.x
      let hh_dy = donor.y - h.y
      let hh_dz = donor.z - h.z
      let ha_dx = acceptor.x - h.x
      let ha_dy = acceptor.y - h.y
      let ha_dz = acceptor.z - h.z
      let len_hh = (hh_dx * hh_dx + hh_dy * hh_dy + hh_dz * hh_dz).sqrt()
      let len_ha = (ha_dx * ha_dx + ha_dy * ha_dy + ha_dz * ha_dz).sqrt()
      if len_hh < 1.0e-6 || len_ha < 1.0e-6 { return None }
      let cos_a = (hh_dx * ha_dx + hh_dy * ha_dy + hh_dz * ha_dz) / (len_hh * len_ha)
      let cos_clamped = if cos_a > 1.0 { 1.0 } else if cos_a < -1.0 { -1.0 } else { cos_a }
      let angle_deg = @math.acos(cos_clamped) * 180.0 / 3.141592653589793
      if angle_deg < min_angle { return None }
      Some((da, angle_deg))
    }
    None => {
      // 无 H：用 D-A 距离 + 近似角度（默认 180°，宽松匹配）
      Some((da, 180.0))
    }
  }
}

///|
/// 检测一帧内的所有 HBond。
///
/// 返回的 HBond 包含 donor/acceptor 原子索引 + 距离 + 角度。
/// 参数完全匹配 MDAnalysis 的默认值。
pub fn AtomGroup::hbonds(
  self : AtomGroup,
  frame~ : Int32 = 0,
  max_distance~ : Double = 3.5,
  max_h_distance~ : Double = 2.5,
  min_angle~ : Double = 120.0,
  consider_h~ : Bool = true,
) -> Array[(Int32, Int32, Double, Double)] {
  let atoms = self.universe.atoms
  let fr = self.universe.frame(frame)
  let n = atoms.length()
  let donor_idx : Array[Int32] = Array(capacity=n)
  let acceptor_idx : Array[Int32] = Array(capacity=n)
  let hydrogen_idx : Array[Int32] = Array(capacity=n)
  let mut d_count = 0
  let mut a_count = 0
  let mut h_count = 0
  let mut i = 0
  while i < n {
    let name = atoms[i].name.to_upper()
    let el = atoms[i].element.to_upper()
    if mda_hbond_is_donor(name) || (el == "N" || el == "O") {
      donor_idx.push(i); d_count = d_count + 1
    }
    if mda_hbond_is_acceptor(name) || el == "O" {
      acceptor_idx.push(i); a_count = a_count + 1
    }
    if consider_h && (el == "H") {
      hydrogen_idx.push(i); h_count = h_count + 1
    }
    i = i + 1
  }
  let out : Array[(Int32, Int32, Double, Double)] = Array(capacity=10)
  let mut di = 0
  while di < donor_idx.length() {
    let d = donor_idx[di]
    let cd = if d < fr.length() { fr[d] } else { Vector3::new(0.0, 0.0, 0.0) }
    let mut ai = 0
    while ai < acceptor_idx.length() {
      let a = acceptor_idx[ai]
      if d == a { ai = ai + 1; continue }
      let ca = if a < fr.length() { fr[a] } else { Vector3::new(0.0, 0.0, 0.0) }
      let geom = if consider_h {
        // 找最近的 H 原子（不严格属于同一个 donor；简化匹配）
        let mut best_h : Vector3? = None
        let mut best_h_dist = 1.2
        let mut hi = 0
        while hi < hydrogen_idx.length() {
          let h = hydrogen_idx[hi]
          let ch = if h < fr.length() { fr[h] } else { Vector3::new(0.0, 0.0, 0.0) }
          let dx = ch.x - cd.x
          let dy = ch.y - cd.y
          let dz = ch.z - cd.z
          let dh = (dx * dx + dy * dy + dz * dz).sqrt()
          if dh < best_h_dist && dh > 0.5 {
            best_h_dist = dh
            best_h = Some(ch)
          }
          hi = hi + 1
        }
        mda_hbond_geom(cd, ca, best_h, max_distance, max_h_distance, min_angle)
      } else {
        mda_hbond_geom(cd, ca, None, max_distance, max_h_distance, min_angle)
      }
      match geom {
        Some((da, ang)) => out.push((d, a, da, ang))
        None => ()
      }
      ai = ai + 1
    }
    di = di + 1
  }
  out
}

///|
/// 统计整个轨迹的 HBond（每帧计数）。
pub fn AtomGroup::hbond_trajectory(
  self : AtomGroup,
  max_distance~ : Double = 3.5,
  max_h_distance~ : Double = 2.5,
  min_angle~ : Double = 120.0,
  consider_h~ : Bool = true,
) -> Array[Int32] {
  let n_frames = self.universe.n_frames()
  let out : Array[Int32] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    let hb = self.hbonds(
      frame=f,
      max_distance=max_distance,
      max_h_distance=max_h_distance,
      min_angle=min_angle,
      consider_h=consider_h,
    )
    out.push(hb.length())
    f = f + 1
  }
  out
}
