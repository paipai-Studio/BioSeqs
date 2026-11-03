// ===========================================================================
// docking_test.mbt — Vina-lite 刚体对接测试.
// ===========================================================================

///|
/// 构造受体"口袋墙": x=10 平面附近的原子墙(口袋中心 (8,0,0)).
fn dk_receptor() -> Array[@src.DockAtom] {
  let atoms : Array[@src.DockAtom] = Array(capacity=0)
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
  // 远处 decoy 原子簇(不构成主口袋)
  atoms.push(@src.DockAtom::{
    x: -30.0,
    y: 0.0,
    z: 0.0,
    element: "C",
    charge: 0.0,
  })
  atoms
}

///|
/// 4 原子小配体(苯环简化: 4C 方形, 边长 1.4A), 初始在原点.
fn dk_ligand() -> Array[@src.DockAtom] {
  [
    @src.DockAtom::{ x: -0.7, y: -0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: 0.7, y: -0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: 0.7, y: 0.7, z: 0.0, element: "C", charge: 0.0, },
    @src.DockAtom::{ x: -0.7, y: 0.7, z: 0.0, element: "C", charge: 0.0, },
  ]
}

///|
test "docking: 口袋检测(指定中心)" {
  let pocket = @src.docking_find_pocket(
    dk_receptor(),
    center=Some((8.0, 0.0, 0.0)),
    radius=6.0,
  )
  assert_true(pocket is Some(_))
  let p = pocket.unwrap()
  assert_eq(p.cx, 8.0)
  assert_eq(p.cy, 0.0)
  // 7x7 = 49 墙原子在 6A 半径内
  inspect(p.n_atoms, content="49")
  inspect(p.atoms.length(), content="49")
}

///|
test "docking: 口口袋检测(密度峰值)" {
  // 无中心: 3A 格点密度峰值(8A 邻域计数)应落在原子墙附近
  let pocket = @src.docking_find_pocket(dk_receptor(), radius=6.0)
  assert_true(pocket is Some(_))
  let p = pocket.unwrap()
  // 墙(x=10)前方同计数格点先到先得(8A 密度球覆盖墙前 4-10A)
  assert_true(p.cx >= 3.0)
  assert_true(p.cx <= 13.0)
  assert_true(p.n_atoms > 0)
}

///|
test "docking: 打分函数物理性" {
  let wall = dk_receptor()
  // 配体远离墙: LJ 吸引弱 -> 能量接近 0
  let far = [
    @src.DockAtom::{ x: 0.0, y: 0.0, z: 0.0, element: "C", charge: 0.0, },
  ]
  let e_far = @src.docking_score(wall, far)
  assert_true(e_far < 0.0)
  assert_true(e_far > -2.0)
  // 配体贴墙(~4.4A, LJ 势阱底): 能量更低(更强吸引)
  let near = [
    @src.DockAtom::{ x: 5.6, y: 0.0, z: 0.0, element: "C", charge: 0.0, },
  ]
  let e_near = @src.docking_score(wall, near)
  assert_true(e_near < e_far)
  // 硬冲突(穿墙 0.5A): 大正罚
  let clash = [
    @src.DockAtom::{ x: 10.5, y: 0.0, z: 0.0, element: "C", charge: 0.0, },
  ]
  let e_clash = @src.docking_score(wall, clash)
  assert_true(e_clash > 0.0)
  // 库仑: 异性电荷增强吸引
  let neg = [
    @src.DockAtom::{ x: 5.6, y: 0.0, z: 0.0, element: "O", charge: -0.5, },
  ]
  let pos_wall : Array[@src.DockAtom] = Array(capacity=0)
  for a in wall {
    if a.charge == 0.0 && a.y == 0.0 && a.z == 0.0 {
      pos_wall.push(@src.DockAtom::{
        x: a.x,
        y: a.y,
        z: a.z,
        element: a.element,
        charge: 0.3,
      })
    }
  }
  let neutral = [
    @src.DockAtom::{ x: 5.6, y: 0.0, z: 0.0, element: "O", charge: 0.0, },
  ]
  let e_coul = @src.docking_score(pos_wall, neg)
  let e_neutral = @src.docking_score(pos_wall, neutral)
  assert_true(e_coul < e_neutral)
}

///|
test "docking: 构象搜索与排名" {
  let pocket = @src.docking_find_pocket(
    dk_receptor(),
    center=Some((8.0, 0.0, 0.0)),
    radius=6.0,
  ).unwrap()
  let params : @src.DockParams = {
    search_radius: 2.0,
    translate_step: 1.0,
    rotations: 3,
    top_k: 5,
  }
  let poses = @src.docking_search(pocket, dk_ligand(), params~)
  inspect(poses.length(), content="5")
  // 能量升序
  let mut sorted_ok = true
  for i = 1; i < poses.length(); i = i + 1 {
    if poses[i].energy < poses[i - 1].energy {
      sorted_ok = false
      break
    }
  }
  assert_eq(sorted_ok, true)
  // 最优姿势能量有限(非硬冲突主导)
  assert_true(poses[0].energy < 0.0)
  // 每姿势配体 4 原子
  inspect(poses[0].ligand.length(), content="4")
}

///|
test "docking: 配体质心放置" {
  let (cx, cy, cz) = @src.docking_ligand_center(dk_ligand())
  assert_eq(cx, 0.0)
  assert_eq(cy, 0.0)
  assert_eq(cz, 0.0)
  // 空配体
  let (ex, ey, ez) = @src.docking_ligand_center([])
  assert_eq(ex, 0.0)
  ignore(ey)
  ignore(ez)
}

///|
test "docking: PDBQT 与 JSON 导出" {
  let pocket = @src.docking_find_pocket(
    dk_receptor(),
    center=Some((8.0, 0.0, 0.0)),
    radius=6.0,
  ).unwrap()
  let params : @src.DockParams = {
    search_radius: 1.0,
    translate_step: 1.0,
    rotations: 2,
    top_k: 3,
  }
  let poses = @src.docking_search(pocket, dk_ligand(), params~)
  let pdbqt = @src.docking_pose_pdbqt(poses[0])
  assert_eq(pdbqt.contains("REMARK vina-lite pose energy="), true)
  assert_eq(pdbqt.contains("ROOT"), true)
  assert_eq(pdbqt.contains("TORSDOF 0"), true)
  assert_eq(pdbqt.contains("ATOM"), true)
  let json = @src.docking_json(poses)
  assert_eq(json.contains("\"rank\""), true)
  assert_eq(json.contains("\"energy\""), true)
  assert_eq(json.contains("\"ligand\""), true)
  let line = @src.docking_summary_line(poses)
  assert_eq(line.contains("poses=3"), true)
  assert_eq(line.contains("best_energy="), true)
}

///|
test "docking: 空输入与无姿势" {
  assert_true(@src.docking_find_pocket([]) is None)
  let pocket = @src.docking_find_pocket(
    dk_receptor(),
    center=Some((8.0, 0.0, 0.0)),
  ).unwrap()
  let poses = @src.docking_search(pocket, [], params=@src.dock_params())
  inspect(poses.length(), content="0")
  inspect(@src.docking_summary_line(poses), content="docking: no poses")
  // 远离口袋的中心 -> 空口袋 None
  let far = @src.docking_find_pocket(
    dk_receptor(),
    center=Some((100.0, 100.0, 100.0)),
    radius=2.0,
  )
  assert_true(far is None)
}
