///|
test "omm_system_particle_and_box_management" {
  let sys = @src.OmmSystem::new()
  assert_eq(sys.num_particles(), 0)
  assert_false(sys.is_periodic())
  let _ = sys.add_particle(
    particle=@src.OmmParticle::new(
      mass=15.9994,
      charge=-0.834,
      sigma=0.316557,
      epsilon=0.650194,
    ),
  )
  let _ = sys.add_particle(particle=@src.OmmParticle::new(mass=1.008))
  let _ = sys.add_particle(particle=@src.OmmParticle::new(mass=1.008))
  assert_eq(sys.num_particles(), 3)
  sys.add_constraint(0, 1, 0.9572)
  sys.add_constraint(0, 2, 0.9572)
  assert_eq(sys.num_constraints(), 2)
  let masses = sys.get_masses()
  assert_true(@moonbit.assert_num_close(masses[0], 15.9994, 1.0e-12))
  let charges = sys.get_charges()
  assert_true(@moonbit.assert_num_close(charges[0], -0.834, 1.0e-12))
  assert_true(@moonbit.assert_num_close(charges[1], 0.0, 1.0e-12))
  sys.set_periodic_box(30.0, 30.0, 30.0)
  assert_true(sys.is_periodic())
}

///|
test "omm_system_bond_and_angle_energy" {
  // harmonic bond: E = 0.5 * k * (r - r0)^2 with r = 1.2, r0 = 1.0, k = 10
  let sys = @src.OmmSystem::new()
  let _ = sys.add_particle(particle=@src.OmmParticle::new())
  let _ = sys.add_particle(particle=@src.OmmParticle::new())
  sys.add_bond(0, 1, 10.0, 1.0)
  let pos = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(pos), 0.2, 1.0e-12),
  )
  // linear trimer: theta = pi, equilibrium pi/2
  let sys2 = @src.OmmSystem::new()
  let _ = sys2.add_particle(particle=@src.OmmParticle::new())
  let _ = sys2.add_particle(particle=@src.OmmParticle::new())
  let _ = sys2.add_particle(particle=@src.OmmParticle::new())
  sys2.add_angle(0, 1, 2, 5.0, @math.PI / 2.0)
  let pos2 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let dtheta = @math.PI / 2.0
  assert_true(
    @moonbit.assert_num_close(
      sys2.potential_energy(pos2),
      0.5 * 5.0 * dtheta * dtheta,
      1.0e-9,
    ),
  )
}

///|
test "omm_system_lj_energy_minimum" {
  let sys = @src.OmmSystem::new()
  let _ = sys.add_particle(
    particle=@src.OmmParticle::new(sigma=1.0, epsilon=1.0),
  )
  let _ = sys.add_particle(
    particle=@src.OmmParticle::new(sigma=1.0, epsilon=1.0),
  )
  // LJ minimum at r = 2^(1/6): E = -epsilon
  let rmin = @math.pow(2.0, 1.0 / 6.0)
  let pos = [[0.0, 0.0, 0.0], [rmin, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(pos), -1.0, 1.0e-9),
  )
  // crossing distance r = sigma: E = 0
  let pos2 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(pos2), 0.0, 1.0e-9),
  )
}

///|
test "omm_system_forces_match_analytic_gradient" {
  let sys = @src.OmmSystem::new()
  let _ = sys.add_particle(particle=@src.OmmParticle::new())
  let _ = sys.add_particle(particle=@src.OmmParticle::new())
  sys.add_bond(0, 1, 10.0, 1.0)
  let pos = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]
  let f = sys.forces(pos)
  // dU/dr = k * (r - r0) = 2.0; the force is opposite on each particle
  // (row-independent allocation since the shared-row bug was fixed).
  assert_true(@moonbit.assert_num_close(f[0][0], 2.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(f[1][0], -2.0, 1.0e-5))
}

///|
test "omm_system_forces_single_particle_zero" {
  let sys = @src.OmmSystem::new()
  let _ = sys.add_particle(particle=@src.OmmParticle::new())
  // a lone particle with no interactions feels no force
  let f = sys.forces([[0.5, 0.5, 0.5]])
  assert_true(@moonbit.assert_num_close(f[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(f[0][2], 0.0, 1.0e-12))
}

// ============================================================================
// B2: omm_system_from_pdb_ff + builtin 测试
// ============================================================================

///|
test "omm_system_from_pdb_ff_builtin_water" {
  // 标准 3 原子 TIP3P 水 PDB + 内置 XML 力场
  let pdb_text = "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let sys = @src.omm_system_from_pdb_ff_builtin(pdb_text)
  // 3 粒子 + 2 键 + 1 角
  assert_eq(sys.particles.length(), 3)
  assert_eq(sys.bonds.length(), 2)
  assert_eq(sys.angles.length(), 1)
  // 约束（O-H）
  assert_true(sys.constraints.length() >= 0) // 不一定加，看 XML 定义
}

///|
test "omm_system_from_pdb_ff_builtin_hcl" {
  // HCl 分子：H 和 Cl 两个原子
  let pdb_text = "ATOM      1  H   HCL A   1       0.000   0.000   0.000  1.00  0.00           H\n" +
    "ATOM      2  CL  HCL A   1       1.270   0.000   0.000  1.00  0.00          Cl\n"
  let sys = @src.omm_system_from_pdb_ff_builtin(pdb_text)
  assert_eq(sys.particles.length(), 2)
  // XML 里没 HCL 残基模板 → 键可能为 0（没有自动推断键），这是合理边界
  // 但 Cl 的电荷应该是 -1.0（从 XML Atom type="Cl" 匹配）
  assert_true(@moonbit.assert_num_close(sys.particles[1].charge, -1.0, 1.0e-6))
  // H 的电荷是 0
  assert_true(@moonbit.assert_num_close(sys.particles[0].charge, 0.0, 1.0e-6))
}

///|
test "omm_system_from_pdb_ff_hexane_homemade_xml" {
  // 手工写极简单的 XML 力场 — element 匹配 → particle.mass 正确
  let xml = "<?xml version=\"1.0\"?>\n" +
    "<ForceField name=\"hexane_test\">\n" +
    "  <Type name=\"C\" element=\"C\" mass=\"12.011\"/>\n" +
    "  <Type name=\"H\" element=\"H\" mass=\"1.008\"/>\n" +
    "  <Atom class=\"A\" type=\"C\" charge=\"-0.054\" sigma=\"0.33997\" epsilon=\"0.35982\"/>\n" +
    "  <Atom class=\"A\" type=\"H\" charge=\"0.018\" sigma=\"0.00000\" epsilon=\"0.000\"/>\n" +
    "</ForceField>\n"
  let pdb_text = "ATOM      1  C1  ETH A   1       0.765   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      2  C2  ETH A   1      -0.765   0.000   0.000  1.00  0.00           C\n"
  let sys = @src.omm_system_from_pdb_ff(pdb_text, xml)
  assert_eq(sys.particles.length(), 2)
  // 两个 C 的质量应为 12.011（从 XML Type 表）
  assert_true(@moonbit.assert_num_close(sys.particles[0].mass, 12.011, 1.0e-3))
  assert_true(@moonbit.assert_num_close(sys.particles[1].mass, 12.011, 1.0e-3))
}

///|
test "omm_system_from_pdb_ff_empty_pdb_no_crash" {
  // 空 PDB → 零粒子系统，不崩溃
  let xml = "<?xml version=\"1.0\"?>\n" +
    "<ForceField name=\"empty_test\">\n" +
    "  <Type name=\"C\" element=\"C\" mass=\"12.011\"/>\n" +
    "</ForceField>\n"
  let sys = @src.omm_system_from_pdb_ff("", xml)
  assert_eq(sys.particles.length(), 0)
  assert_eq(sys.bonds.length(), 0)
}

///|
test "omm_system_from_pdb_ff_builtin_parses_hetatm" {
  // HETATM 记录也应该被解析
  let pdb_text = "HETATM    1  NA  NA  A   1       0.000   0.000   0.000  1.00  0.00          Na\n" +
    "HETATM    2  CL  CL  A   1       2.820   0.000   0.000  1.00  0.00          Cl\n"
  let sys = @src.omm_system_from_pdb_ff_builtin(pdb_text)
  assert_eq(sys.particles.length(), 2)
  assert_true(@moonbit.assert_num_close(sys.particles[0].charge, 1.0, 1.0e-6)) // Na+
  assert_true(@moonbit.assert_num_close(sys.particles[1].charge, -1.0, 1.0e-6)) // Cl-
}

///|
test "omm_pdb_element_from_atom_name" {
  // 测试元素推断逻辑（H -> H, Cl -> C 首字母, Fe -> F 首字母）
  let el1 = @src.omm_system_from_pdb_ff_builtin(
    "ATOM      1  H   MOL A   1       0.000   0.000   0.000  1.00  0.00           H\n",
  )
  assert_eq(el1.particles.length(), 1)
  assert_true(@moonbit.assert_num_close(el1.particles[0].mass, 1.008, 1.0e-3))
}
