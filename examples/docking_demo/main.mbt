///|
/// docking_demo — Vina-lite 刚体对接演示:
/// 原子墙受体口袋 + 4C 小配体 -> 口袋检测 -> 平移 x 旋转网格
/// 构象搜索 -> LJ+Coulomb 打分 -> 姿势排名(PDBQT/JSON).

///|
fn dkd_receptor() -> Array[@src.DockAtom] {
  let atoms : Array[@src.DockAtom] = Array(capacity=0)
  // x=10 处 7x7 原子墙(口袋朝向 x-)
  for i = -3; i <= 3; i = i + 1 {
    for j = -3; j <= 3; j = j + 1 {
      atoms.push(@src.DockAtom::{
        x: 10.0,
        y: i.to_double(),
        z: j.to_double(),
        element: "C",
        charge: 0.0,
      })
    }
  }
  // 带正电的极性位点(墙中央)
  atoms.push(@src.DockAtom::{
    x: 10.0,
    y: 0.0,
    z: 0.0,
    element: "N",
    charge: 0.4,
  })
  atoms
}

///|
/// 6 原子配体(苯环简化 + O 受体), 初始在原点.
fn dkd_ligand() -> Array[@src.DockAtom] {
  [
    @src.DockAtom::{ x: -0.7, y: -0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: 0.7, y: -0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: 0.7, y: 0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: -0.7, y: 0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: 1.9, y: 0.0, z: 0.0, element: "O", charge: -0.5, },
    @src.DockAtom::{ x: -1.9, y: 0.0, z: 0.0, element: "H", charge: 0.2, },
  ]
}

///|
fn main {
  println(
    "1. Vina-lite 刚体对接演示(平移格点 x 欧拉角网格, LJ+Coulomb 打分)",
  )
  let receptor = dkd_receptor()
  println(
    "2. 受体: \{receptor.length()} 原子(7x7 原子墙 + 极性位点)",
  )
  // 3. 口袋检测(参考位点模式: 墙前 LJ 势阱区中心)
  let pocket = @src.docking_find_pocket(
    receptor,
    center=Some((5.0, 0.0, 0.0)),
    radius=6.0,
  ).unwrap()
  println(
    "3. 口袋检测(参考位点): 中心=(\{pocket.cx}, \{pocket.cy}, \{pocket.cz}), 原子数=\{pocket.n_atoms}",
  )
  // 4. 构象搜索
  let params : @src.DockParams = {
    search_radius: 2.0,
    translate_step: 1.0,
    rotations: 4,
    top_k: 5,
  }
  let poses = @src.docking_search(pocket, dkd_ligand(), params~)
  println(
    "4. 构象搜索(半径\{params.search_radius}A/步\{params.translate_step}A/每轴\{params.rotations}旋转): 前 5 姿势",
  )
  for i = 0; i < poses.length(); i = i + 1 {
    let p = poses[i]
    let (lcx, lcy, lcz) = @src.docking_ligand_center(p.ligand)
    println(
      "  #\{i + 1}: energy=\{p.energy}, 配体质心=(\{lcx}, \{lcy}, \{lcz})",
    )
  }
  // 5. 导出
  let pdbqt = @src.docking_pose_pdbqt(poses[0])
  let json = @src.docking_json(poses)
  println(
    "5. 最优姿势 PDBQT \{pdbqt.length()} 字节, JSON \{json.length()} 字节",
  )
  println(@src.docking_summary_line(poses))
  println(
    "演示完成: 受体口袋 -> 刚体构象枚举 -> 能量排名 -> 姿势输出闭环",
  )
}
