///|
/// Tests for the extended MDAnalysis-style selection language.
test "mda_select_universe_basic_keywords" {
  let atoms = test_sel_atoms()
  // name / resname
  assert_eq(@src.mda_select_universe("name CA", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("resname HOH", atoms).length(), 1)
  // segid / type / element
  assert_eq(@src.mda_select_universe("segid B", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("type CT1", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("element O", atoms).length(), 2)
  // resid 单值与区间
  assert_eq(@src.mda_select_universe("resid 2", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("resid 1:2", atoms).length(), 4)
  // 组合
  assert_eq(@src.mda_select_universe("protein and name CA", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("not water", atoms).length(), 4)
  assert_eq(@src.mda_select_universe("all", atoms).length(), 5)
  assert_eq(@src.mda_select_universe("none", atoms).length(), 0)
}

///|
test "mda_select_universe_bynum_byres" {
  let atoms = test_sel_atoms()
  // bynum：1-based 编号，区间与多值
  let r = @src.mda_select_universe("bynum 1:2", atoms)
  assert_eq(r.length(), 2)
  assert_eq(r[0], 0)
  assert_eq(r[1], 1)
  let r2 = @src.mda_select_universe("bynum 1 3 5", atoms)
  assert_eq(r2.length(), 3)
  assert_eq(r2[2], 4)
  // byres：CA 所在残基扩展为整个残基（resid 1 的 3 个原子）
  let r3 = @src.mda_select_universe("byres name CA", atoms)
  assert_eq(r3.length(), 3)
  assert_eq(r3[0], 0)
  assert_eq(r3[2], 2)
  // (byres name O) and (not protein)：name O 命中 resid 2 的 O 与水的 O，
  // 前者是蛋白质被排除，只剩水（索引 4）
  let r4 = @src.mda_select_universe("byres name O and not protein", atoms)
  assert_eq(r4.length(), 1)
  assert_eq(r4[0], 4)
}

///|
test "mda_select_universe_around_point" {
  let atoms = test_sel_atoms()
  let pos = test_sel_positions()
  // CA 在 (1.5,0,0)：原子 0 距 1.5、原子 1 距 1.0 均在 2.0 内
  let r = @src.mda_select_universe(
    "around 2.0 name CA",
    atoms,
    positions=Some(pos),
  )
  assert_eq(r.length(), 2)
  assert_eq(r[0], 0)
  assert_eq(r[1], 1)
  // PBC 下原子 3 (9,0,0) 镜像距离 2.5 > 2.0 不命中；原子 4 (9.5,0,0)
  // 镜像距离恰为 2.0 命中
  let r_pbc = @src.mda_select_universe(
    "around 2.0 name CA",
    atoms,
    positions=Some(pos),
    box=Some([10.0, 10.0, 10.0]),
  )
  assert_eq(r_pbc.length(), 3)
  // 放宽到 2.5 后 PBC 下原子 3 也命中（镜像距离 2.5），共 4 个
  let r_pbc2 = @src.mda_select_universe(
    "around 2.5 name CA",
    atoms,
    positions=Some(pos),
    box=Some([10.0, 10.0, 10.0]),
  )
  assert_eq(r_pbc2.length(), 4)
  // point：距 (0,0,0) 2.0 以内的原子 0/1/2（含 CA 距 1.5）
  let p = @src.mda_select_universe(
    "point 0.0 0.0 0.0 2.0",
    atoms,
    positions=Some(pos),
  )
  assert_eq(p.length(), 3)
  assert_eq(p[0], 0)
  assert_eq(p[1], 1)
  assert_eq(p[2], 2)
}

///|
test "mda_select_universe_prop" {
  let atoms = test_sel_atoms()
  // prop mass：H 1.008, C 12.011, O 15.999
  assert_eq(@src.mda_select_universe("prop mass < 10.0", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("prop mass >= 12.0", atoms).length(), 4)
  // prop charge：H 0.2719 与 CA 0.0213 均为正
  assert_eq(@src.mda_select_universe("prop charge > 0.0", atoms).length(), 2)
  assert_eq(@src.mda_select_universe("prop charge == 0.0", atoms).length(), 1)
  assert_eq(@src.mda_select_universe("prop charge <= -0.4", atoms).length(), 2)
}

///|
test "mda_select_universe_precedence_parens" {
  let atoms = test_sel_atoms()
  // 优先级：not > and > or
  // protein and (not name H) = {0, 2, 3}
  let r1 = @src.mda_select_universe("protein and not name H", atoms)
  assert_eq(r1.length(), 3)
  assert_eq(r1[0], 0)
  assert_eq(r1[1], 2)
  assert_eq(r1[2], 3)
  // 括号分组：(name N or name O) and not water = {0, 3}
  let r2 = @src.mda_select_universe("(name N or name O) and not water", atoms)
  assert_eq(r2.length(), 2)
  assert_eq(r2[0], 0)
  assert_eq(r2[1], 3)
  // not 作用于括号表达式：not (protein and backbone) = {1, 4}
  let r3 = @src.mda_select_universe("not (protein and backbone)", atoms)
  assert_eq(r3.length(), 2)
  assert_eq(r3[0], 1)
  assert_eq(r3[1], 4)
  // or 链：water or type CT1 = {2, 4}
  let r4 = @src.mda_select_universe("water or type CT1", atoms)
  assert_eq(r4.length(), 2)
  assert_eq(r4[0], 2)
  assert_eq(r4[1], 4)
}

///|
fn test_sel_atoms() -> Array[@src.UnivAtom] {
  // 0: N (ALA, resid 1, chain A)
  // 1: H (ALA, resid 1, chain A)
  // 2: CA (ALA, resid 1, chain A, type CT1)
  // 3: O (ALA, resid 2, chain A)
  // 4: O (HOH, resid 3, chain B, segid B)
  let base = @src.UnivAtom::{
    name: "N",
    resname: "ALA",
    resid: 1,
    chain: "A",
    segid: "A",
    element: "N",
    atom_type: "NH1",
    charge: -0.4157,
    mass: 14.007,
    is_protein: true,
    is_water: false,
    is_backbone: true,
  }
  let h = @src.UnivAtom::{
    ..base,
    name: "H",
    element: "H",
    atom_type: "H",
    charge: 0.2719,
    mass: 1.008,
    is_backbone: false,
  }
  let ca = @src.UnivAtom::{
    ..base,
    name: "CA",
    element: "C",
    atom_type: "CT1",
    charge: 0.0213,
    mass: 12.011,
    is_backbone: true,
  }
  let o2 = @src.UnivAtom::{
    ..base,
    name: "O",
    resid: 2,
    element: "O",
    charge: 0.0,
    mass: 15.999,
  }
  let water = @src.UnivAtom::{
    name: "O",
    resname: "HOH",
    resid: 3,
    chain: "B",
    segid: "B",
    element: "O",
    atom_type: "OH2",
    charge: -0.834,
    mass: 15.999,
    is_protein: false,
    is_water: true,
    is_backbone: false,
  }
  [base, h, ca, o2, water]
}

///|
fn test_sel_positions() -> Array[Array[Double]] {
  [
    [0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0],
    [1.5, 0.0, 0.0],
    [9.0, 0.0, 0.0],
    [9.5, 0.0, 0.0],
  ]
}
