///|
/// docking.mbt — Vina-lite 风格蛋白质-配体刚体对接
/// (本地确定性实现).
///
/// 诚实声明(与 AutoDock Vina 的差距):
/// - 构象搜索: 口袋内平移格点穷举 + 均匀欧拉角旋转网格
///   (刚体, 无可旋转键扭转自由度; Vina 为遗传算法 + 柔性键);
/// - 打分: Lennard-Jones 12-6 + Coulomb (元素参数表组合规则)
///   + 硬冲突惩罚; Vina 为经验项(gauss/repulsion/hydrophobic/hbond);
/// - 口袋检测: 参考中心 + 半径, 或格点原子密度峰值(简化凹陷).
///
/// 功能:
/// - docking_find_pocket: 指定中心邻域原子收集, 或无中心时的
///   密度峰值口袋(3A 格点, 8A 邻域计数最大者);
/// - docking_search: 平移 x 旋转网格构象枚举, 每构象 LJ+Coulomb
///   + 冲突惩罚打分, 返回能量升序 top-k;
/// - 姿势输出: 配体变换后坐标 + PDBQT 风格文本 + JSON.
///
/// 数据不出本机.

// ---------------------------------------------------------------------------
// 原子模型与参数
// ---------------------------------------------------------------------------

///|
/// 对接原子(受体或配体统一模型).
pub(all) struct DockAtom {
  x : Double
  y : Double
  z : Double
  /// 元素符号("C"/"N"/"O"/"S"/"H"/...).
  element : String
  /// 部分电荷(e; 默认 0 时仅 LJ 项).
  charge : Double
}

///|
/// 元素 LJ 参数(sigma A, epsilon kcal/mol).
fn dock_lj_params(element : String) -> (Double, Double) {
  match element {
    "C" => (3.5, 0.1)
    "N" => (3.3, 0.15)
    "O" => (3.2, 0.15)
    "S" => (3.9, 0.2)
    "P" => (4.0, 0.2)
    "F" => (3.0, 0.1)
    "CL" => (4.1, 0.2)
    "BR" => (4.2, 0.25)
    _ => (3.0, 0.05)
  }
}

// ---------------------------------------------------------------------------
// 口袋检测
// ---------------------------------------------------------------------------

///|
/// 检测到的口袋.
pub(all) struct DockPocket {
  /// 口袋中心.
  cx : Double
  cy : Double
  cz : Double
  /// 口袋邻域原子数.
  n_atoms : Int32
  /// 口内受体原子(供打分).
  atoms : Array[DockAtom]
}

///|
/// 口袋检测.
///
/// center 为 Some 时收集中心 radius 邻域的受体原子(参考配体/已知
/// 活性位点模式); center 为 None 时做格点密度峰值检测: 3A 格点上
/// 8A 邻域原子计数最大者为口袋中心.
pub fn docking_find_pocket(
  receptor : Array[DockAtom],
  center~ : (Double, Double, Double)? = None,
  radius~ : Double = 8.0,
) -> DockPocket? {
  if receptor.length() == 0 {
    return None
  }
  match center {
    Some((cx, cy, cz)) => {
      let atoms = dock_atoms_within(receptor, cx, cy, cz, radius)
      if atoms.length() == 0 {
        return None
      }
      Some(DockPocket::{ cx, cy, cz, n_atoms: atoms.length(), atoms, })
    }
    None => {
      // 包围盒
      let mut min_x = 1.0e30
      let mut min_y = 1.0e30
      let mut min_z = 1.0e30
      let mut max_x = -1.0e30
      let mut max_y = -1.0e30
      let mut max_z = -1.0e30
      for a in receptor {
        if a.x < min_x {
          min_x = a.x
        }
        if a.y < min_y {
          min_y = a.y
        }
        if a.z < min_z {
          min_z = a.z
        }
        if a.x > max_x {
          max_x = a.x
        }
        if a.y > max_y {
          max_y = a.y
        }
        if a.z > max_z {
          max_z = a.z
        }
      }
      // 3A 格点密度峰值
      let grid = 3.0
      let mut best_cx = 0.0
      let mut best_cy = 0.0
      let mut best_cz = 0.0
      let mut best_count = -1
      let mut gx = min_x
      while gx <= max_x {
        let mut gy = min_y
        while gy <= max_y {
          let mut gz = min_z
          while gz <= max_z {
            let cnt = dock_count_within(receptor, gx, gy, gz, 8.0)
            if cnt > best_count {
              best_count = cnt
              best_cx = gx
              best_cy = gy
              best_cz = gz
            }
            gz = gz + grid
          }
          gy = gy + grid
        }
        gx = gx + grid
      }
      if best_count <= 0 {
        return None
      }
      let atoms = dock_atoms_within(receptor, best_cx, best_cy, best_cz, radius)
      Some(DockPocket::{
        cx: best_cx,
        cy: best_cy,
        cz: best_cz,
        n_atoms: atoms.length(),
        atoms,
      })
    }
  }
}

///|
fn dock_atoms_within(
  atoms : Array[DockAtom],
  cx : Double,
  cy : Double,
  cz : Double,
  radius : Double,
) -> Array[DockAtom] {
  let out : Array[DockAtom] = Array(capacity=0)
  for a in atoms {
    let dx = a.x - cx
    let dy = a.y - cy
    let dz = a.z - cz
    if dx * dx + dy * dy + dz * dz <= radius * radius {
      out.push(a)
    }
  }
  out
}

///|
fn dock_count_within(
  atoms : Array[DockAtom],
  cx : Double,
  cy : Double,
  cz : Double,
  radius : Double,
) -> Int32 {
  let mut n = 0
  for a in atoms {
    let dx = a.x - cx
    let dy = a.y - cy
    let dz = a.z - cz
    if dx * dx + dy * dy + dz * dz <= radius * radius {
      n = n + 1
    }
  }
  n
}

// ---------------------------------------------------------------------------
// 打分
// ---------------------------------------------------------------------------

///|
/// 一个构象的受体-配体相互作用能(kcal/mol, 越低越稳).
/// LJ 12-6(组合规则) + Coulomb(k_e = 332.0637) + 硬冲突惩罚
/// (r < 0.65 * sigma_ij 时 +100, 模拟 Vina repulsion).
pub fn docking_score(
  receptor : Array[DockAtom],
  ligand : Array[DockAtom],
) -> Double {
  let ke = 332.0637
  let mut e = 0.0
  for la in ligand {
    let (ls, le) = dock_lj_params(la.element)
    for ra in receptor {
      let dx = la.x - ra.x
      let dy = la.y - ra.y
      let dz = la.z - ra.z
      let r2 = dx * dx + dy * dy + dz * dz
      if r2 == 0.0 {
        continue
      }
      let r = r2.sqrt()
      let (rs, re) = dock_lj_params(ra.element)
      let sigma = (ls + rs) / 2.0
      let eps = (le * re).sqrt()
      // 硬冲突: 原子间距过近
      if r < 0.65 * sigma {
        e = e + 100.0
        continue
      }
      let sr = sigma / r
      let sr6 = sr * sr * sr * sr * sr * sr
      let sr12 = sr6 * sr6
      e = e + 4.0 * eps * (sr12 - sr6)
      // 截断长程(12A)加速
      if r < 12.0 && (la.charge != 0.0 || ra.charge != 0.0) {
        e = e + ke * la.charge * ra.charge / r
      }
    }
  }
  e
}

// ---------------------------------------------------------------------------
// 构象搜索
// ---------------------------------------------------------------------------

///|
/// 搜索参数.
pub(all) struct DockParams {
  /// 平移搜索半径(A, 绕口袋中心).
  search_radius : Double
  /// 平移格点步长(A).
  translate_step : Double
  /// 每轴旋转采样数(欧拉角 theta/phi/psi 各 N 份).
  rotations : Int32
  /// 返回前 k 个姿势.
  top_k : Int32
}

///|
/// 默认参数.
pub fn dock_params() -> DockParams {
  { search_radius: 4.0, translate_step: 2.0, rotations: 6, top_k: 10, }
}

///|
/// 一个对接姿势.
pub(all) struct DockPose {
  /// 配体平移(相对初始坐标).
  tx : Double
  ty : Double
  tz : Double
  /// 欧拉角(弧度).
  theta : Double
  phi : Double
  psi : Double
  /// 相互作用能(kcal/mol).
  energy : Double
  /// 变换后的配体原子坐标.
  ligand : Array[DockAtom]
}

///|
/// 绕轴旋转点 (绕中心 cx,cy,cz, 用欧拉角 theta/phi/psi 的
/// Z-Y-Z 旋转矩阵).
fn dock_rotate_atom(
  a : DockAtom,
  cx : Double,
  cy : Double,
  cz : Double,
  theta : Double,
  phi : Double,
  psi : Double,
) -> DockAtom {
  // 平移到原点
  let x = a.x - cx
  let y = a.y - cy
  let z = a.z - cz
  // Z-Y-Z 欧拉旋转
  let (c1, s1) = (@math.cos(theta), @math.sin(theta))
  let (c2, s2) = (@math.cos(phi), @math.sin(phi))
  let (c3, s3) = (@math.cos(psi), @math.sin(psi))
  // R = Rz(psi) * Ry(phi) * Rz(theta)
  let r00 = c3 * c2 * c1 - s3 * s1
  let r01 = -c3 * c2 * s1 - s3 * c1
  let r02 = c3 * s2
  let r10 = s3 * c2 * c1 + c3 * s1
  let r11 = -s3 * c2 * s1 + c3 * c1
  let r12 = s3 * s2
  let r20 = -s2 * c1
  let r21 = s2 * s1
  let r22 = c2
  let nx = r00 * x + r01 * y + r02 * z
  let ny = r10 * x + r11 * y + r12 * z
  let nz = r20 * x + r21 * y + r22 * z
  // 平移回中心
  DockAtom::{
    x: nx + cx,
    y: ny + cy,
    z: nz + cz,
    element: a.element,
    charge: a.charge,
  }
}

///|
/// 平移配体原子(绕质心先居中再放至目标中心).
fn dock_place_ligand(
  ligand : Array[DockAtom],
  from_cx : Double,
  from_cy : Double,
  from_cz : Double,
  to_cx : Double,
  to_cy : Double,
  to_cz : Double,
) -> Array[DockAtom] {
  let dx = to_cx - from_cx
  let dy = to_cy - from_cy
  let dz = to_cz - from_cz
  let out : Array[DockAtom] = Array(capacity=ligand.length())
  for a in ligand {
    out.push(DockAtom::{
      x: a.x + dx,
      y: a.y + dy,
      z: a.z + dz,
      element: a.element,
      charge: a.charge,
    })
  }
  out
}

///|
/// 配体质心.
pub fn docking_ligand_center(
  ligand : Array[DockAtom],
) -> (Double, Double, Double) {
  if ligand.length() == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  for a in ligand {
    sx = sx + a.x
    sy = sy + a.y
    sz = sz + a.z
  }
  let n = ligand.length().to_double()
  (sx / n, sy / n, sz / n)
}

///|
/// 刚体对接搜索: 平移格点 x 欧拉角旋转网格穷举,
/// 返回能量升序前 k 个姿势.
pub fn docking_search(
  pocket : DockPocket,
  ligand : Array[DockAtom],
  params~ : DockParams = dock_params(),
) -> Array[DockPose] {
  if ligand.length() == 0 || pocket.atoms.length() == 0 {
    return []
  }
  let (lcx, lcy, lcz) = docking_ligand_center(ligand)
  // 先把配体质心放到口袋中心
  let centered = dock_place_ligand(
    ligand,
    lcx,
    lcy,
    lcz,
    pocket.cx,
    pocket.cy,
    pocket.cz,
  )
  let poses : Array[DockPose] = Array(capacity=0)
  let rot = params.rotations
  let d2pi = 6.283185307179586
  // 平移格点(绕口袋中心, 步长 grid, 半径 search_radius)
  let grid = params.translate_step
  let r = params.search_radius
  let n_steps = (r / grid).to_int()
  let mut ix = -n_steps
  while ix <= n_steps {
    let mut iy = -n_steps
    while iy <= n_steps {
      let mut iz = -n_steps
      while iz <= n_steps {
        let tx = ix.to_double() * grid
        let ty = iy.to_double() * grid
        let tz = iz.to_double() * grid
        // 平移后的质心
        let ncx = pocket.cx + tx
        let ncy = pocket.cy + ty
        let ncz = pocket.cz + tz
        let moved = dock_place_ligand(
          centered,
          pocket.cx,
          pocket.cy,
          pocket.cz,
          ncx,
          ncy,
          ncz,
        )
        // 旋转网格(绕新质心)
        let mut it = 0
        while it < rot {
          let theta = d2pi * it.to_double() / rot.to_double()
          let mut ip = 0
          while ip < rot {
            let phi = d2pi * ip.to_double() / rot.to_double()
            let mut ipsi = 0
            while ipsi < rot {
              let psi = d2pi * ipsi.to_double() / rot.to_double()
              let rotated : Array[DockAtom] = Array(capacity=moved.length())
              for a in moved {
                rotated.push(
                  dock_rotate_atom(a, ncx, ncy, ncz, theta, phi, psi),
                )
              }
              let energy = docking_score(pocket.atoms, rotated)
              poses.push(DockPose::{
                tx,
                ty,
                tz,
                theta,
                phi,
                psi,
                energy,
                ligand: rotated,
              })
              ipsi = ipsi + 1
            }
            ip = ip + 1
          }
          it = it + 1
        }
        iz = iz + 1
      }
      iy = iy + 1
    }
    ix = ix + 1
  }
  // 能量升序排序(确定性 tie-break 保持枚举序)
  poses.sort_by(fn(a, b) {
    if a.energy < b.energy {
      -1
    } else if a.energy > b.energy {
      1
    } else {
      0
    }
  })
  let k = if params.top_k < poses.length() {
    params.top_k
  } else {
    poses.length()
  }
  let out : Array[DockPose] = Array(capacity=k)
  for i = 0; i < k; i = i + 1 {
    out.push(poses[i])
  }
  out
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 姿势配体坐标 -> PDBQT 风格文本(ATOMIC 记录, 简化为
/// 元素+电荷标注).
pub fn docking_pose_pdbqt(pose : DockPose) -> String {
  let sb = StringBuilder()
  sb.write_string("REMARK vina-lite pose energy=\{pose.energy}\n")
  sb.write_string("REMARK translate=\{pose.tx} \{pose.ty} \{pose.tz}\n")
  sb.write_string("ROOT\n")
  for a in pose.ligand {
    sb.write_string(
      "ATOM\t\t\{a.element}\t\{a.x}\t\{a.y}\t\{a.z}\t\{a.charge}\n",
    )
  }
  sb.write_string("ENDROOT\nTORSDOF 0\n")
  sb.to_string()
}

///|
/// 姿势列表 -> JSON 值.
pub fn docking_value(poses : Array[DockPose]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=poses.length())
  for i = 0; i < poses.length(); i = i + 1 {
    let p = poses[i]
    let lig : Array[JsonValue] = Array(capacity=p.ligand.length())
    for a in p.ligand {
      lig.push(
        JsonArr([
          JsonNum(a.x),
          JsonNum(a.y),
          JsonNum(a.z),
          JsonStr(a.element),
          JsonNum(a.charge),
        ]),
      )
    }
    arr.push(
      JsonObj(
        Map([
          ("rank", JsonNum((i + 1).to_double())),
          ("energy", JsonNum(p.energy)),
          ("tx", JsonNum(p.tx)),
          ("ty", JsonNum(p.ty)),
          ("tz", JsonNum(p.tz)),
          ("theta", JsonNum(p.theta)),
          ("phi", JsonNum(p.phi)),
          ("psi", JsonNum(p.psi)),
          ("ligand", JsonArr(lig)),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 姿势列表 -> JSON 文本.
pub fn docking_json(poses : Array[DockPose]) -> String {
  json_to_string(docking_value(poses))
}

///|
/// 摘要行.
pub fn docking_summary_line(poses : Array[DockPose]) -> String {
  if poses.length() == 0 {
    return "docking: no poses"
  }
  "docking: poses=\{poses.length()}, best_energy=\{poses[0].energy}, worst_energy=\{poses[poses.length() - 1].energy}"
}
