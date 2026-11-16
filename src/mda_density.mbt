// ============================================================================
// B5: Density Map — 3D 空间概率密度网格
// ============================================================================

///|
pub struct MdaDensityMap {
  origin : Array[Double]
  shape : (Int32, Int32, Int32)
  spacing : Array[Double]
  values : Array[Double]
  total_samples : Int32
}

///|
fn mda_wrap_coord(v : Double, l : Double) -> Double {
  if l <= 0.0 {
    return v
  }
  let k = (v / l).to_int()
  v - k.to_double() * l
}

///|
/// AtomGroup 所有帧的 3D 网格概率密度。
pub fn AtomGroup::density_map(
  self : AtomGroup,
  cell_size~ : Double = 1.0,
  pad~ : Double = 2.0,
  pbc~ : Bool = false,
) -> MdaDensityMap {
  let nf = self.universe.n_frames()
  let na = self.indices.length()
  let bx = self.universe.box.v1.x
  let by = self.universe.box.v2.y
  let bz = self.universe.box.v3.z
  let mut min_x = 1.7976931348623157e308
  let mut min_y = 1.7976931348623157e308
  let mut min_z = 1.7976931348623157e308
  let mut max_x = -1.7976931348623157e308
  let mut max_y = -1.7976931348623157e308
  let mut max_z = -1.7976931348623157e308
  let mut f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    let mut i = 0
    while i < na {
      let idx = self.indices[i]
      if idx < fr.length() {
        let c = fr[idx]
        if c.x < min_x {
          min_x = c.x
        }
        if c.x > max_x {
          max_x = c.x
        }
        if c.y < min_y {
          min_y = c.y
        }
        if c.y > max_y {
          max_y = c.y
        }
        if c.z < min_z {
          min_z = c.z
        }
        if c.z > max_z {
          max_z = c.z
        }
      }
      i = i + 1
    }
    f = f + 1
  }
  if pbc && bx > 0.0 {
    min_x = 0.0
    max_x = bx
    min_y = 0.0
    max_y = by
    min_z = 0.0
    max_z = bz
  } else {
    min_x = min_x - pad
    max_x = max_x + pad
    min_y = min_y - pad
    max_y = max_y + pad
    min_z = min_z - pad
    max_z = max_z + pad
  }
  let nx = ((max_x - min_x) / cell_size).to_int() + 1
  let ny = ((max_y - min_y) / cell_size).to_int() + 1
  let nz = ((max_z - min_z) / cell_size).to_int() + 1
  let total = nx * ny * nz
  let counts : Array[Int32] = Array(capacity=total)
  let mut k = 0
  while k < total {
    counts.push(0)
    k = k + 1
  }
  let mut tot_samp = 0
  f = 0
  while f < nf {
    let fr = self.universe.frame(f)
    let mut i = 0
    while i < na {
      let idx = self.indices[i]
      if idx < fr.length() {
        let mut c = fr[idx]
        if pbc {
          c = Vector3::new(
            mda_wrap_coord(c.x, max_x),
            mda_wrap_coord(c.y, max_y),
            mda_wrap_coord(c.z, max_z),
          )
        }
        let bx2 = ((c.x - min_x) / cell_size).to_int()
        let by2 = ((c.y - min_y) / cell_size).to_int()
        let bz2 = ((c.z - min_z) / cell_size).to_int()
        if bx2 >= 0 && bx2 < nx && by2 >= 0 && by2 < ny && bz2 >= 0 && bz2 < nz {
          let flat = bz2 * nx * ny + by2 * nx + bx2
          counts[flat] = counts[flat] + 1
          tot_samp = tot_samp + 1
        }
      }
      i = i + 1
    }
    f = f + 1
  }
  let values : Array[Double] = Array(capacity=total)
  let ts = tot_samp.to_double()
  k = 0
  while k < total {
    values.push(if ts > 0.0 { counts[k].to_double() / ts } else { 0.0 })
    k = k + 1
  }
  MdaDensityMap::{
    origin: [min_x, min_y, min_z],
    shape: (nx, ny, nz),
    spacing: [cell_size, cell_size, cell_size],
    values,
    total_samples: tot_samp,
  }
}

///|
/// 密度最大值（最可能位置）。
pub fn MdaDensityMap::max_probability(self : MdaDensityMap) -> Double {
  let mut m = 0.0
  let mut i = 0
  while i < self.values.length() {
    if self.values[i] > m {
      m = self.values[i]
    }
    i = i + 1
  }
  m
}
