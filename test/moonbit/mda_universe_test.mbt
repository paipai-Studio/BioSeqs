///|
/// Tests for the MDAnalysis-style Universe / AtomGroup object model.
test "universe_from_pdb_text_multi_model" {
  // 两帧三原子：N / CA / O 沿 x 排布，第二帧 z 偏移 1.0
  let pdb = "MODEL        1\n" +
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n" +
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000\n" +
    "ATOM      3  O   ALA A   1       3.000   0.000   0.000\n" +
    "ENDMDL\n" +
    "MODEL        2\n" +
    "ATOM      1  N   ALA A   1       0.000   0.000   1.000\n" +
    "ATOM      2  CA  ALA A   1       1.500   0.000   1.000\n" +
    "ATOM      3  O   ALA A   1       3.000   0.000   1.000\n" +
    "ENDMDL\n" +
    "END\n"
  let u = @src.universe_from_pdb_text(pdb)
  assert_eq(u.n_atoms(), 3)
  assert_eq(u.n_frames(), 2)
  // 元素从原子名推断：CA -> C, O -> O, N -> N
  assert_eq(u.atom_info()[0].element, "N")
  assert_eq(u.atom_info()[1].element, "C")
  assert_eq(u.atom_info()[2].element, "O")
  // 蛋白 / 骨架分类
  assert_true(u.atom_info()[0].is_protein)
  assert_true(u.atom_info()[0].is_backbone)
  assert_true(u.atom_info()[1].is_backbone)
  assert_false(u.atom_info()[0].is_water)
  // 帧切换
  assert_eq(u.frame_index(), 0)
  let p0 = u.positions()
  assert_true(@moonbit.assert_num_close(p0[1][0], 1.5, 1.0e-9))
  u.set_frame(1)
  assert_eq(u.frame_index(), 1)
  let p1 = u.positions()
  assert_true(@moonbit.assert_num_close(p1[0][2], 1.0, 1.0e-9))
  // set_frame 越界钳制
  u.set_frame(99)
  assert_eq(u.frame_index(), 1)
  u.set_frame(-5)
  assert_eq(u.frame_index(), 0)
}

///|
test "universe_from_psf_dcd_roundtrip" {
  // PSF：3 原子（ALA 链 A + HOH 链 B），1 个键
  let psf = "PSF EXT\n\n" +
    "       1 !NTITLE\n REMARKS test\n\n" +
    "       3 !NATOM\n" +
    "   1 A     1  ALA  N    NH1   -0.415700       14.0070           0\n" +
    "   2 A     1  ALA  CA   CT1    0.021300       12.0110           0\n" +
    "   3 B     2  HOH  O    OH2   -0.834000       15.9990           0\n\n" +
    "       1 !NBOND: bonds\n   1     2\n\n"
  // DCD：2 帧 3 原子
  let f1 = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [10.0, 0.0, 0.0]]
  let f2 = [[0.0, 0.0, 2.0], [1.5, 0.0, 2.0], [10.0, 0.0, 2.0]]
  let dcd = @src.write_dcd(@src.Trajectory::new([f1, f2]))
  let u = @src.universe_from_psf_dcd(psf, dcd)
  assert_eq(u.n_atoms(), 3)
  assert_eq(u.n_frames(), 2)
  // PSF 元数据
  assert_eq(u.atom_info()[0].segid, "A")
  assert_eq(u.atom_info()[2].atom_type, "OH2")
  assert_true(
    @moonbit.assert_num_close(u.atom_info()[0].charge, -0.4157, 1.0e-9),
  )
  assert_true(u.atom_info()[0].is_backbone)
  assert_true(u.atom_info()[2].is_water)
  assert_false(u.atom_info()[2].is_protein)
  // DCD 坐标
  assert_true(@moonbit.assert_num_close(u.positions()[1][0], 1.5, 1.0e-5))
  u.set_frame(1)
  assert_true(@moonbit.assert_num_close(u.positions()[0][2], 2.0, 1.0e-5))
}

///|
test "universe_from_pqr_text" {
  let pqr = "REMARK test\n" +
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  -0.415700  14.0100\n" +
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000   0.021300  12.0110\n"
  let u = @src.universe_from_pqr_text(pqr)
  assert_eq(u.n_atoms(), 2)
  assert_eq(u.n_frames(), 1)
  assert_true(
    @moonbit.assert_num_close(u.atom_info()[0].charge, -0.4157, 1.0e-9),
  )
  assert_true(@moonbit.assert_num_close(u.positions()[1][0], 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(u.atom_info()[1].mass, 12.011, 1.0e-9))
}

///|
test "atomgroup_geometry_and_mass" {
  // 等质量 3 原子沿 x：(0,0,0),(2,0,0),(4,0,0)
  let a1 = @src.UnivAtom::{
    name: "C",
    resname: "ALA",
    resid: 1,
    chain: "A",
    segid: "A",
    element: "C",
    atom_type: "",
    charge: 0.0,
    mass: 12.011,
    is_protein: true,
    is_water: false,
    is_backbone: false,
  }
  let a2 = a1
  let a3 = a1
  let traj = @src.Trajectory::new([
    [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
  ])
  let u = @src.Universe::new([a1, a2, a3], traj)
  let g = u.all()
  assert_eq(g.n_atoms(), 3)
  let cog = g.center_of_geometry()
  assert_true(@moonbit.assert_num_close(cog[0], 2.0, 1.0e-9))
  let com = g.center_of_mass()
  assert_true(@moonbit.assert_num_close(com[0], 2.0, 1.0e-9))
  // Rg = sqrt((4+0+4)/3)
  assert_true(
    @moonbit.assert_num_close(g.radius_of_gyration(), 1.632993, 1.0e-5),
  )
  // 总电荷 / 总质量
  assert_true(@moonbit.assert_num_close(g.total_mass(), 36.033, 1.0e-9))
  // names / resids
  assert_eq(g.names().length(), 3)
  assert_eq(g.resids()[0], 1)
}

///|
test "atomgroup_subselect_and_select_within" {
  let pdb = "MODEL        1\n" +
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n" +
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000\n" +
    "ATOM      3  O   ALA A   1       3.000   0.000   0.000\n" +
    "ATOM      4  O   HOH B   2       9.000   0.000   0.000\n" +
    "ENDMDL\n"
  let u = @src.universe_from_pdb_text(pdb)
  // 基础选择
  let cas = u.select_atoms("name CA")
  assert_eq(cas.n_atoms(), 1)
  assert_eq(cas.indices()[0], 1)
  // 复合选择
  let bb = u.select_atoms("protein and backbone")
  assert_eq(bb.n_atoms(), 3)
  let nowater = u.select_atoms("not water")
  assert_eq(nowater.n_atoms(), 3)
  // 子组内再选择
  let sub = bb.select("name O")
  assert_eq(sub.n_atoms(), 1)
  assert_eq(sub.indices()[0], 2)
  // positions_at 指定帧
  assert_eq(cas.positions().length(), 1)
  assert_true(@moonbit.assert_num_close(cas.positions()[0][0], 1.5, 1.0e-9))
  // 排除水后 positions 仍是全集坐标
  assert_true(@moonbit.assert_num_close(nowater.positions()[0][0], 0.0, 1.0e-9))
}

///|
test "mda_infer_element_and_mass" {
  assert_eq(@src.mda_infer_element("CA"), "C")
  assert_eq(@src.mda_infer_element("N"), "N")
  assert_eq(@src.mda_infer_element("1HG2"), "H")
  assert_eq(@src.mda_infer_element("CL1"), "Cl")
  assert_eq(@src.mda_infer_element("FE"), "Fe")
  assert_true(
    @moonbit.assert_num_close(@src.mda_element_mass("H"), 1.008, 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(@src.mda_element_mass("FE"), 55.845, 1.0e-9),
  )
  // 未知元素回退碳
  assert_true(
    @moonbit.assert_num_close(@src.mda_element_mass("Xx"), 12.011, 1.0e-9),
  )
  // 水判定
  assert_true(@src.mda_is_water_resname("HOH"))
  assert_true(@src.mda_is_water_resname("tip3"))
  assert_false(@src.mda_is_water_resname("ALA"))
}
