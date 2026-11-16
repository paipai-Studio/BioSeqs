// ============================================================================
// C5: Spatial Hash Grid — 快速 O(1) 邻域搜索（Array-based，无 Map）
// ============================================================================
//
// 对 N 个原子按 cell_size 均匀分桶，邻居查询只扫描相邻 27 个桶。
// 比 O(N²) 蛮力快 10-100×，用于 HBond/Contacts/RDF 加速。
// ============================================================================

///|
pub struct SpatialGrid {
  cell_size : Double
  /// 所有原子坐标
  coords : Array[Vector3]
  /// 所有原子原始全局索引
  global_indices : Array[Int32]
  n_atoms : Int32
  /// cell 边界（为避免 Map 依赖，用 sorted (key, bucket) pair）
  bucket_keys : Array[String]
  bucket_values : Array[Array[Int32]]
}

///|
fn mda_sg_ix(v : Double, cs : Double) -> Int32 {
  if cs <= 0.0 {
    0
  } else {
    (v / cs).to_int()
  }
}

///|
fn mda_sg_key(ix : Int32, iy : Int32, iz : Int32) -> String {
  ix.to_string() + "," + iy.to_string() + "," + iz.to_string()
}

///|
fn mda_sg_find_idx(keys : Array[String], target : String) -> Int32 {
  let mut i = 0
  while i < keys.length() {
    if keys[i] == target {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// 从 AtomGroup 构建 SpatialGrid。
pub fn AtomGroup::build_spatial_grid(
  self : AtomGroup,
  cell_size~ : Double = 4.5,
  frame~ : Int32 = 0,
) -> SpatialGrid {
  let idx = self.indices
  let n = idx.length()
  let fr = self.universe.frame(frame)
  let coords : Array[Vector3] = Array(capacity=n)
  let gi : Array[Int32] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let id = idx[i]
    coords.push(
      if id < fr.length() {
        fr[id]
      } else {
        Vector3::new(0.0, 0.0, 0.0)
      },
    )
    gi.push(id)
    i = i + 1
  }
  let bucket_keys : Array[String] = Array(capacity=32)
  let bucket_values : Array[Array[Int32]] = Array(capacity=32)
  i = 0
  while i < n {
    let c = coords[i]
    let key = mda_sg_key(
      mda_sg_ix(c.x, cell_size),
      mda_sg_ix(c.y, cell_size),
      mda_sg_ix(c.z, cell_size),
    )
    let bi = mda_sg_find_idx(bucket_keys, key)
    if bi >= 0 {
      bucket_values[bi].push(i)
    } else {
      let new_arr : Array[Int32] = Array(capacity=4)
      new_arr.push(i)
      bucket_keys.push(key)
      bucket_values.push(new_arr)
    }
    i = i + 1
  }
  SpatialGrid::{
    cell_size,
    coords,
    global_indices: gi,
    n_atoms: n,
    bucket_keys,
    bucket_values,
  }
}

///|
/// 查找所有与 local_idx 距离 ≤ cutoff 的原子（返回局部索引）。
pub fn SpatialGrid::neighbors(
  self : SpatialGrid,
  local_idx : Int32,
  cutoff~ : Double = 4.5,
) -> Array[Int32] {
  if local_idx < 0 || local_idx >= self.n_atoms {
    return Array(capacity=0)
  }
  let c = self.coords[local_idx]
  let cs = self.cell_size
  let ix = mda_sg_ix(c.x, cs)
  let iy = mda_sg_ix(c.y, cs)
  let iz = mda_sg_ix(c.z, cs)
  let cutoff2 = cutoff * cutoff
  let out : Array[Int32] = Array(capacity=16)
  let mut dx = -1
  while dx <= 1 {
    let mut dy = -1
    while dy <= 1 {
      let mut dz = -1
      while dz <= 1 {
        let key = mda_sg_key(ix + dx, iy + dy, iz + dz)
        let bi = mda_sg_find_idx(self.bucket_keys, key)
        if bi >= 0 {
          let bucket = self.bucket_values[bi]
          let mut k = 0
          while k < bucket.length() {
            let other = bucket[k]
            if other == local_idx {
              k = k + 1
              continue
            }
            let oc = self.coords[other]
            let ddx = c.x - oc.x
            let ddy = c.y - oc.y
            let ddz = c.z - oc.z
            let d2 = ddx * ddx + ddy * ddy + ddz * ddz
            if d2 <= cutoff2 {
              out.push(other)
            }
            k = k + 1
          }
        }
        dz = dz + 1
      }
      dy = dy + 1
    }
    dx = dx + 1
  }
  out
}

///|
/// 返回 self 内所有原子对距离 ≤ cutoff 的边（无向唯一对）。
pub fn SpatialGrid::all_pairs_within(
  self : SpatialGrid,
  cutoff~ : Double = 4.5,
) -> Array[(Int32, Int32)] {
  let cutoff2 = cutoff * cutoff
  let out : Array[(Int32, Int32)] = Array(capacity=64)
  let cs = self.cell_size
  let mut i = 0
  while i < self.n_atoms {
    let c = self.coords[i]
    let ix = mda_sg_ix(c.x, cs)
    let iy = mda_sg_ix(c.y, cs)
    let iz = mda_sg_ix(c.z, cs)
    let mut dx = 0
    while dx <= 1 {
      let mut dy = if dx == 0 { -1 } else { -1 }
      while dy <= 1 {
        let mut dz = if dx == 0 && dy == 0 { 1 } else { -1 }
        while dz <= 1 {
          let key = mda_sg_key(ix + dx, iy + dy, iz + dz)
          let bi = mda_sg_find_idx(self.bucket_keys, key)
          if bi >= 0 {
            let bucket = self.bucket_values[bi]
            let mut k = 0
            while k < bucket.length() {
              let other = bucket[k]
              if other > i {
                let oc = self.coords[other]
                let ddx = c.x - oc.x
                let ddy = c.y - oc.y
                let ddz = c.z - oc.z
                let d2 = ddx * ddx + ddy * ddy + ddz * ddz
                if d2 <= cutoff2 {
                  out.push((i, other))
                }
              }
              k = k + 1
            }
          }
          dz = dz + 1
        }
        dy = dy + 1
      }
      dx = dx + 1
    }
    i = i + 1
  }
  out
}

///|
pub fn SpatialGrid::global_index(
  self : SpatialGrid,
  local_idx : Int32,
) -> Int32 {
  if local_idx >= 0 && local_idx < self.global_indices.length() {
    self.global_indices[local_idx]
  } else {
    -1
  }
}
