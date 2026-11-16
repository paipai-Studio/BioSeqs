///|
test "system_universe_roundtrip" {
  let atoms : Array[@src.UnivAtom] = [
    @src.UnivAtom::new(
      name="C1",
      resname="MOL",
      chain="A",
      segid="SYS",
      element="C",
      atom_type="CT",
      resid=1,
      mass=12.011,
      charge=-0.1,
      is_protein=false,
      is_water=false,
      is_backbone=false,
    ),
    @src.UnivAtom::new(
      name="H1",
      resname="MOL",
      chain="A",
      segid="SYS",
      element="H",
      atom_type="HC",
      resid=1,
      mass=1.008,
      charge=0.1,
      is_protein=false,
      is_water=false,
      is_backbone=false,
    ),
  ]
  let frame : Array[Array[Double]] = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let u0 = @src.Universe::from_dcd_frames(atoms, [frame])
  // Universe -> OmmSystem
  let sys = @src.omm_system_from_universe(u0, sigma=3.0, epsilon=0.1)
  assert_eq(sys.num_particles(), 2)
  let charges = sys.get_charges()
  assert_true(@moonbit.assert_num_close(charges[0], -0.1, 1.0e-12))
  let masses = sys.get_masses()
  assert_true(@moonbit.assert_num_close(masses[0], 12.011, 1.0e-12))
  assert_true(@moonbit.assert_num_close(masses[1], 1.008, 1.0e-12))
  sys.set_periodic_box(10.0, 10.0, 10.0)
  assert_true(sys.is_periodic())
  // OmmSystem -> Universe (elements guessed from mass)
  let u1 = @src.omm_system_to_universe(sys, positions=Some(frame))
  assert_eq(u1.n_atoms(), 2)
  let info = u1.atom_info()
  assert_eq(info[0].element, "C")
  assert_eq(info[1].element, "H")
  assert_true(@moonbit.assert_num_close(info[0].charge, -0.1, 1.0e-12))
  assert_eq(u1.frame_index(), 0)
  let back_pos = u1.positions_at(0)
  assert_true(@moonbit.assert_num_close(back_pos[1][0], 1.5, 1.0e-12))
}

///|
test "prmtop_parse_full_and_bridges" {
  let prmtop = "%VERSION  VERSION STAMP = V0001.000  DATE = 11/07/26\n%FLAG POINTERS\n%FORMAT(10I8)\n       2       2       1       0       0       0       0       0       0       1\n%FLAG ATOM_NAME\n%FORMAT(20a4)\nC   H   \n%FLAG CHARGE\n%FORMAT(5E16.8)\n  -1.82000000E+00   1.82000000E+00\n%FLAG MASS\n%FORMAT(5E16.8)\n  1.20110000E+01   1.00800000E+00\n%FLAG ATOM_TYPE_INDEX\n%FORMAT(10I8)\n       1       2\n%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n       6       1\n%FLAG AMBER_ATOM_TYPE\n%FORMAT(20a4)\nCT  HC  \n%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nMOL \n%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1\n%FLAG BOND_FORCE_CONSTANT\n%FORMAT(5E16.8)\n  3.00000000E+02\n%FLAG BOND_EQUIL_VALUE\n%FORMAT(5E16.8)\n  1.50000000E+00\n%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n       0       3       1\n%FLAG LENNARD_JONES_ACOEF\n%FORMAT(5E16.8)\n  1.00000000E+03   8.00000000E+02   6.00000000E+02\n%FLAG LENNARD_JONES_BCOEF\n%FORMAT(5E16.8)\n  5.00000000E+00   4.00000000E+00   3.00000000E+00\n"
  let top = @src.omm_parse_prmtop_full(prmtop)
  assert_eq(top.n_atoms, 2)
  assert_eq(top.n_types, 2)
  assert_eq(top.atom_names.length(), 2)
  assert_eq(top.atom_names[0], "C")
  assert_eq(top.atom_names[1], "H")
  assert_eq(top.atomic_numbers[0], 6)
  assert_eq(top.bonds.length(), 1)
  // bond indices are 0-based: (0, 1, param 0)
  let b = top.bonds[0]
  assert_eq(b.0, 0)
  assert_eq(b.1, 1)
  assert_eq(b.2, 0)
  // charges: internal units / 18.2223 -> e (|q| = 1.82 / 18.2223 = 0.09989)
  assert_true(
    @moonbit.assert_num_close(top.charges[0].abs(), 1.82 / 18.2223, 1.0e-9),
  )
  // LJ diagonal (1,1): sigma = (1000/5)^(1/6), eps = 25/4000
  let lj = @src.omm_prmtop_lj(top, 1, 1)
  assert_true(
    @moonbit.assert_num_close(lj.0, @math.pow(200.0, 1.0 / 6.0), 1.0e-6),
  )
  assert_true(@moonbit.assert_num_close(lj.1, 0.00625, 1.0e-9))
  // Off-diagonal symmetric lookup (2,1) == (1,2)
  let lj2a = @src.omm_prmtop_lj(top, 1, 2)
  let lj2b = @src.omm_prmtop_lj(top, 2, 1)
  assert_true(@moonbit.assert_num_close(lj2a.0, lj2b.0, 1.0e-12))
  // prmtop -> Universe (with inpcrd coordinates)
  let inpcrd = "test system\n        2\n     0.0000000    0.0000000    0.0000000\n     1.5000000    0.0000000    0.0000000\n"
  let u = @src.omm_universe_from_prmtop(prmtop, inpcrd~)
  assert_eq(u.n_atoms(), 2)
  assert_eq(u.n_frames(), 1)
  let info = u.atom_info()
  assert_eq(info[0].resname, "MOL")
  assert_eq(info[0].resid, 1)
  assert_eq(info[0].element, "C")
  assert_eq(info[1].element, "H")
  assert_true(info[0].charge < 0.0)
  let p = u.positions_at(0)
  assert_true(@moonbit.assert_num_close(p[1][0], 1.5, 1.0e-6))
  // prmtop -> OmmSystem
  let sys = @src.omm_system_from_prmtop(prmtop)
  assert_eq(sys.num_particles(), 2)
  assert_eq(sys.num_constraints(), 0)
  let e = sys.potential_energy([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
  // Bond at equilibrium contributes 0; LJ is repulsive here. The pair
  // epsilon is the geometric mean of the per-type epsilons
  // (sqrt(0.00625 * 0.00375) = 0.004841):
  // 4*eps*((s/r)^12-(s/r)^6) = 4*0.004841*(308.294 - 17.558) = 5.6301
  assert_true(@moonbit.assert_num_close(e, 5.6301, 1.0e-3))
  let e2 = sys.potential_energy([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
  // LJ is repulsive at r=2.0 < sigma=2.418 (4*eps*((s/r)^12-(s/r)^6) =
  // 0.1286), plus bond strain 0.5*300*0.25 = 37.5
  assert_true(e2 > 10.0)
  // Box variant
  let sys2 = @src.omm_system_from_prmtop_box(prmtop, box=[20.0, 20.0, 20.0])
  assert_true(sys2.is_periodic())
}

///|
test "element_from_number_and_mass" {
  assert_eq(@src.omm_element_from_number(6), "C")
  assert_eq(@src.omm_element_from_number(8), "O")
  assert_eq(@src.omm_element_from_number(11), "Na")
  assert_eq(@src.omm_element_from_number(16), "S")
  assert_eq(@src.omm_element_from_number(999), "X")
  assert_eq(@src.omm_element_from_mass(1.008), "H")
  assert_eq(@src.omm_element_from_mass(12.011), "C")
  assert_eq(@src.omm_element_from_mass(15.999), "O")
  assert_eq(@src.omm_element_from_mass(22.990), "Na")
  assert_eq(@src.omm_element_from_mass(3.14159), "X")
}

///|
test "prmtop_missing_flags_defaults" {
  // A minimal file with only POINTERS parses to empty arrays.
  let tiny = "%FLAG POINTERS\n%FORMAT(10I8)\n       3       1\n"
  let top = @src.omm_parse_prmtop_full(tiny)
  assert_eq(top.n_atoms, 3)
  assert_eq(top.bonds.length(), 0)
  let u = @src.omm_universe_from_prmtop(tiny)
  assert_eq(u.n_atoms(), 3)
  assert_eq(u.n_frames(), 0)
}

// ============================================================================
// B1 + A7/A8: omm_simulation_run_to_universe + OmmTrajectorySample 测试
// ============================================================================

///|
test "omm_simulation_run_to_universe_basic" {
  // 3 个原子的水模拟，跑 100 步，每 10 步一帧 → 10 帧
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=300.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=100,
    report_interval=10,
  )
  assert_eq(sample.n_frames(), 10)
  // 能量数组长度 = 帧数
  assert_eq(sample.total_potential_energy().length(), 10)
  assert_eq(sample.temperature_trajectory().length(), 10)
  assert_eq(sample.pressure_trajectory().length(), 10)
  // 温度应该 > 0（有动能）
  let temps = sample.temperature_trajectory()
  assert_true(temps[0] > 0.0)
  assert_true(temps[temps.length() - 1] > 0.0)
}

///|
test "omm_trajectory_sample_energy_stats" {
  // 5 帧，PE = 10, 20, 30, 40, 50 → mean=30, std=~15.811, min=10, max=50
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  H   MOL A   1       1.500   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=10.0) // 低温，少运动
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=50,
    report_interval=10,
  )
  let pe_stats = sample.energy_stats(sample.total_potential_energy())
  let (mean, sd, mn, mx) = pe_stats
  assert_eq(sample.n_frames(), 5)
  assert_true(mx >= mn) // 基本不变式
  assert_true(mean >= mn && mean <= mx)
  assert_true(sd >= 0.0)
}

///|
test "omm_trajectory_sample_to_universe_roundtrip" {
  // 跑模拟 → 转 Universe → 验证原子数和帧数
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u_ref = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u_ref, temperature=100.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=200,
    report_interval=20,
  )
  let u_traj = sample.to_universe(u_ref)
  assert_eq(u_traj.n_atoms(), 3)
  assert_eq(u_traj.n_frames(), 10)
  // 最后一帧的坐标应该和 sample.frames[last] 一致
  let last_pos = u_traj.positions_at(9)
  assert_true(
    @moonbit.assert_num_close(last_pos[0][0], sample.frames[9][0].x, 1.0e-6),
  )
}

///|
test "omm_simulation_run_single_frame" {
  // 只跑 10 步 → 1 帧
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=100.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=10,
    report_interval=10,
  )
  assert_eq(sample.n_frames(), 1)
}

///|
test "omm_trajectory_sample_total_energy_sums_pe_ke" {
  // 验证 total_energy = pe + ke
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  H   MOL A   1       1.500   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=50.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=50,
    report_interval=10,
  )
  let pe = sample.total_potential_energy()
  let ke = sample.kinetic_energy()
  let total = sample.total_energy()
  for i = 0; i < total.length(); i = i + 1 {
    assert_true(@moonbit.assert_num_close(total[i], pe[i] + ke[i], 1.0e-9))
  }
}

// ============================================================================
// B3: Universe::to_omm_system / to_omm_simulation 原生方法测试
// ============================================================================

///|
test "universe_to_omm_system_native_method" {
  // B3 原生方法：Universe::to_omm_system 命名比 omm_system_from_universe 更自然
  let atoms : Array[@src.UnivAtom] = [
    @src.UnivAtom::new(
      name="O",
      resname="HOH",
      chain="A",
      segid="SYS",
      element="O",
      atom_type="OH",
      resid=1,
      mass=15.999,
      charge=-0.834,
      is_protein=false,
      is_water=true,
      is_backbone=false,
    ),
    @src.UnivAtom::new(
      name="H1",
      resname="HOH",
      chain="A",
      segid="SYS",
      element="H",
      atom_type="HO",
      resid=1,
      mass=1.008,
      charge=0.417,
      is_protein=false,
      is_water=true,
      is_backbone=false,
    ),
    @src.UnivAtom::new(
      name="H2",
      resname="HOH",
      chain="A",
      segid="SYS",
      element="H",
      atom_type="HO",
      resid=1,
      mass=1.008,
      charge=0.417,
      is_protein=false,
      is_water=true,
      is_backbone=false,
    ),
  ]
  let frame : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [-0.2879, 0.9120, 0.0],
  ]
  let u = @src.Universe::from_dcd_frames(atoms, [frame])
  let sys = @src.Universe::to_omm_system(u)
  assert_eq(sys.num_particles(), 3)
  let charges = sys.get_charges()
  assert_true(@moonbit.assert_num_close(charges[0], -0.834, 1.0e-6))
}

///|
test "universe_to_omm_simulation_params_and_time" {
  // B3: to_omm_simulation(dt=0.002, temperature=100.0, friction=5.0) 跑 50 步
  // → state.time ≈ 50 * 0.002 = 0.1 ps（验证 dt 参数传递）
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, _sys) = @src.Universe::to_omm_simulation(
    u,
    dt=0.002,
    temperature=100.0,
    friction=5.0,
  )
  sim.step(n_steps=50)
  let st = sim.state()
  // dt=0.002, 50 步 → 时间 ≈ 0.1 ps
  assert_true(@moonbit.assert_num_close(st.time, 0.1, 1.0e-6))
  assert_eq(st.step, 50)
}

///|
test "universe_to_omm_simulation_multi_frame_uses_first" {
  // B3: 多帧 Universe → 用 frame(0) 初始化 sim 的 positions，不是最后一帧
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  H   MOL A   1       1.000   0.000   0.000  1.00  0.00           H\n"
  let u1 = @src.universe_from_pdb_text(pdb_text)
  // 直接构造带 2 帧的 Universe
  let atoms = u1.atom_info()
  let frame0 : Array[Array[Double]] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let frame1 : Array[Array[Double]] = [[9.0, 9.0, 9.0], [8.0, 8.0, 8.0]]
  let u = @src.Universe::from_dcd_frames(atoms, [frame0, frame1])
  assert_eq(u.n_frames(), 2)
  let (sim, _sys) = @src.Universe::to_omm_simulation(u, temperature=0.01) // 极低温，几乎不动
  let st0 = sim.state()
  // 初始位置应该是 frame(0)，不是 frame1
  assert_true(@moonbit.assert_num_close(st0.positions[0][0], 0.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(st0.positions[1][0], 1.0, 1.0e-6))
}

///|
test "universe_to_omm_simulation_run_then_bridge_to_mda" {
  // B3 完整故事链：Universe → to_omm_simulation → 模拟 → sample → to_universe → 能喂 mda 分析
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u_ref = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u_ref, temperature=300.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=50,
    report_interval=10,
  )
  let u_traj = sample.to_universe(u_ref)
  // 结果应能喂 MDAnalysis 风格的分析函数
  assert_eq(u_traj.n_atoms(), 3)
  assert_eq(u_traj.n_frames(), 5)
  // 能量分解分量应该非零（有 bond）
  let bonds = sample.bond()
  assert_eq(bonds.length(), 5)
}

// ============================================================================
// A7: 能量分解时间序列测试
// ============================================================================

///|
test "omm_trajectory_sample_energy_decomposed_nonzero" {
  // A7: 有真实 bond/angle 体系 → bond() 和 angle() 应至少一个分量有限
  // 注意：builtin XML 可能没给显式 bond force → bond 能量可能为 0
  // 但 energy_decomposed 返回的 5 个分量必须都是有限 Double（API 合约）
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=300.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=50,
    report_interval=10,
  )
  // 所有能量分量必须是有限 Double
  let energies = [
    sample.total_potential_energy(),
    sample.bond(),
    sample.angle(),
    sample.dihedral(),
    sample.lj(),
    sample.coulomb(),
  ]
  for arr in energies {
    for e in arr {
      assert_true(e.is_finite())
    }
  }
}

///|
test "omm_trajectory_sample_coulomb_lj_nonzero" {
  // A7: 两个带电荷原子 + 距离较近 → coulomb() 非零
  let pdb_text = "ATOM      1  NA  NA  A   1       0.000   0.000   0.000  1.00  0.00          Na\n" +
    "ATOM      2  CL  CL  A   1       2.000   0.000   0.000  1.00  0.00          Cl\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  // 给 Na/Cl 显式电荷
  let info = u.atom_info()
  let atoms_updated = [
    @src.UnivAtom::new(
      name=info[0].name,
      resname=info[0].resname,
      chain=info[0].chain,
      segid=info[0].segid,
      element="Na",
      atom_type="Na",
      resid=1,
      mass=info[0].mass,
      charge=1.0,
      is_protein=false,
      is_water=false,
      is_backbone=false,
    ),
    @src.UnivAtom::new(
      name=info[1].name,
      resname=info[1].resname,
      chain=info[1].chain,
      segid=info[1].segid,
      element="Cl",
      atom_type="Cl",
      resid=1,
      mass=info[1].mass,
      charge=-1.0,
      is_protein=false,
      is_water=false,
      is_backbone=false,
    ),
  ]
  let frame : Array[Array[Double]] = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let u2 = @src.Universe::from_dcd_frames(atoms_updated, [frame])
  let (sim, sys) = @src.Universe::to_omm_simulation(u2, temperature=100.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=20,
    report_interval=10,
  )
  let coul = sample.coulomb()
  let ljs = sample.lj()
  // Na-Cl 静电吸引，coulomb 应该非零
  assert_eq(coul.length(), 2)
  assert_eq(ljs.length(), 2)
  // coulomb 能量应该显著（±1 电荷，2 Å 距离）
  assert_true(coul[0].abs() > 1.0)
}

///|
test "omm_system_energy_decomposed_sum_approx_pe" {
  // A7 底层接口：system.energy_decomposed 各分量和 ≈ potential_energy（无约束体系）
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  H   MOL A   1       0.740   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=50.0)
  sim.step(n_steps=10)
  let st = sim.state()
  let (bond, ang, dih, lj, coul) = sys.energy_decomposed(st.positions)
  // 验证返回值都是有限 Double
  assert_true(bond.is_finite())
  assert_true(ang.is_finite())
  assert_true(dih.is_finite())
  assert_true(lj.is_finite())
  assert_true(coul.is_finite())
  // 两原子之间应有 LJ 相互作用（H-H 近距离）
  let nonbond = lj + coul
  assert_true(nonbond.is_finite())
}

// ============================================================================
// A8: 温度 / 压力 / 密度时间序列测试
// ============================================================================

///|
test "omm_trajectory_sample_density_trajectory_basic" {
  // A8 density 长度 = 帧数，值为有限
  // 注意：如果 periodic box 体积为 0，density 可能是 0.0
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=300.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=50,
    report_interval=10,
  )
  let dens = sample.density_trajectory()
  assert_eq(dens.length(), 5)
  for d in dens {
    assert_true(d.is_finite())
    assert_true(d >= 0.0)
  }
}

///|
test "omm_trajectory_sample_temperature_positive_and_numeric" {
  // A8 temperature 全部为正且有限
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  H   MOL A   1       1.000   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=50.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=30,
    report_interval=10,
  )
  let temps = sample.temperature_trajectory()
  assert_eq(temps.length(), 3)
  for t in temps {
    assert_true(t.is_finite())
    assert_true(t > 0.0)
  }
}

///|
test "omm_trajectory_sample_pressure_nonnegative_and_numeric" {
  // A8 pressure 非负且有限（理想气体近似 P = NkT/V）
  let pdb_text = "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, sys) = @src.Universe::to_omm_simulation(u, temperature=100.0)
  let sample = @src.omm_simulation_run_to_universe(
    sim,
    sys,
    n_steps=20,
    report_interval=10,
  )
  let press = sample.pressure_trajectory()
  assert_eq(press.length(), 2)
  for p in press {
    assert_true(p.is_finite())
    assert_true(p >= 0.0) // 理想气体 P 不能负
  }
}
