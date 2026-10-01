test "openmm langevin step reduces energy for stretched bond" {
  @src.md_seed_rng(42L)
  let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let forces = @src.md_bond_forces(positions, bonds)
  let (pos_new, _v) = @src.md_langevin_step(
    positions, velocities, forces, masses, 0.001, 1.0, 10.0,
  )
  let dx = pos_new[1][0] - pos_new[0][0]
  assert_true(dx < 1.5)
}

test "openmd langevin finish updates velocities" {
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let forces = [[50.0, 0.0, 0.0], [-50.0, 0.0, 0.0]]
  let v_new = @src.md_langevin_finish(velocities, forces, masses, 0.01)
  assert_true(v_new[0][0] > 0.0)
  assert_true(v_new[1][0] < 0.0)
}

test "openmm brownian step moves toward equilibrium" {
  let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let forces = [[100.0, 0.0, 0.0], [-100.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let pos_new = @src.md_brownian_step(positions, forces, masses, 0.001, 300.0, 10.0)
  assert_true(pos_new[0][0] > 0.0)
  assert_true(pos_new[1][0] < 2.0)
}

test "openmm brownian simulation returns temperatures" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let temps = @src.md_simulate_brownian(positions, masses, bonds, 0.001, 300.0, 10.0, 10)
  assert_eq(temps.length(), 10)
  for i = 0; i < 10; i = i + 1 {
    assert_eq(temps[i], 300.0)
  }
}

test "openmm rbtorsion energy zero at default" {
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let torsions = [(0, 1, 2, 3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]
  let e = @src.md_rbtorsion_energy(positions, torsions)
  assert_true((e - 0.0).abs() < 1.0e-10)
}

test "openmm rbtorsion energy constant offset" {
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let torsions = [(0, 1, 2, 3, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0)]
  let e = @src.md_rbtorsion_energy(positions, torsions)
  assert_true((e - 5.0).abs() < 1.0e-10)
}

test "openmm rbtorsion forces sum to zero" {
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let torsions = [(0, 1, 2, 3, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)]
  let forces = @src.md_rbtorsion_forces(positions, torsions)
  let mut fx = 0.0
  let mut fy = 0.0
  let mut fz = 0.0
  for i = 0; i < 4; i = i + 1 {
    fx = fx + forces[i][0]
    fy = fy + forces[i][1]
    fz = fz + forces[i][2]
  }
  assert_true(fx.abs() < 1.0e-3)
  assert_true(fy.abs() < 1.0e-3)
  assert_true(fz.abs() < 1.0e-3)
}

test "openmm cmap energy bilinear interpolation" {
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [2.0, 1.0, 1.0],
  ]
  let grid = [0.0, 1.0, 2.0, 3.0]
  let cmaps = [(0, 1, 2, 3, 4, grid, 2)]
  let e = @src.md_cmap_energy(positions, cmaps)
  assert_true(e.abs() < 1.0e6)
}

test "openmm cmap energy zero grid" {
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [2.0, 1.0, 1.0],
  ]
  let grid = [0.0, 0.0, 0.0, 0.0]
  let cmaps = [(0, 1, 2, 3, 4, grid, 2)]
  let e = @src.md_cmap_energy(positions, cmaps)
  assert_true((e - 0.0).abs() < 1.0e-10)
}

test "openmm andersen thermostat returns correct size" {
  let velocities = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
  ]
  let masses = [1.0, 1.0, 1.0]
  let v_new = @src.md_andersen_thermostat(velocities, masses, 300.0, 1.0, 0.01)
  assert_eq(v_new.length(), 3)
  for i = 0; i < 3; i = i + 1 {
    assert_eq(v_new[i].length(), 3)
  }
}

test "openmm andersen simulation returns temperatures" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let temps = @src.md_simulate_andersen(
    positions, velocities, masses, bonds, 0.001, 300.0, 1.0, 5, 6,
  )
  assert_eq(temps.length(), 5)
}

test "openmm remove com motion zero total momentum" {
  let velocities = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
  ]
  let masses = [1.0, 1.0, 1.0]
  let v_new = @src.md_remove_com_motion(velocities, masses)
  let mut px = 0.0
  let mut py = 0.0
  let mut pz = 0.0
  for i = 0; i < 3; i = i + 1 {
    px = px + masses[i] * v_new[i][0]
    py = py + masses[i] * v_new[i][1]
    pz = pz + masses[i] * v_new[i][2]
  }
  assert_true(px.abs() < 1.0e-10)
  assert_true(py.abs() < 1.0e-10)
  assert_true(pz.abs() < 1.0e-10)
}

test "openmm remove com motion preserves relative velocity" {
  let velocities = [
    [1.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
  ]
  let masses = [1.0, 1.0]
  let v_new = @src.md_remove_com_motion(velocities, masses)
  let dv_before = velocities[1][0] - velocities[0][0]
  let dv_after = v_new[1][0] - v_new[0][0]
  assert_true((dv_before - dv_after).abs() < 1.0e-10)
}

test "openmm remove com angular removes rotation" {
  let positions = [
    [1.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
  ]
  let velocities = [
    [0.0, 1.0, 0.0],
    [0.0, -1.0, 0.0],
    [-1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
  ]
  let masses = [1.0, 1.0, 1.0, 1.0]
  let v_new = @src.md_remove_com_angular(positions, velocities, masses)
  for i = 0; i < 4; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      assert_true(v_new[i][d].abs() < 1.0e-6)
    }
  }
}

test "openmm pbc minimum image orthorhombic" {
  let dr = [6.0, 0.0, 0.0]
  let box = [10.0, 10.0, 10.0]
  let result = @src.md_pbc_minimum_image(dr, box)
  assert_true((result[0] - (-4.0)).abs() < 1.0e-10)
  assert_true(result[1].abs() < 1.0e-10)
  assert_true(result[2].abs() < 1.0e-10)
}

test "openmm pbc minimum image no wrap needed" {
  let dr = [3.0, -2.0, 1.0]
  let box = [10.0, 10.0, 10.0]
  let result = @src.md_pbc_minimum_image(dr, box)
  assert_true((result[0] - 3.0).abs() < 1.0e-10)
  assert_true((result[1] - (-2.0)).abs() < 1.0e-10)
  assert_true((result[2] - 1.0).abs() < 1.0e-10)
}

test "openmm pbc wrap position into unit cell" {
  let pos = [12.0, -3.0, 5.0]
  let box = [10.0, 10.0, 10.0]
  let result = @src.md_pbc_wrap_position(pos, box)
  assert_true((result[0] - 2.0).abs() < 1.0e-10)
  assert_true((result[1] - 7.0).abs() < 1.0e-10)
  assert_true((result[2] - 5.0).abs() < 1.0e-10)
}

test "openmm pbc wrap all positions" {
  let positions = [[12.0, 0.0, 0.0], [-3.0, 5.0, 15.0]]
  let box = [10.0, 10.0, 10.0]
  let result = @src.md_pbc_wrap_all(positions, box)
  assert_true((result[0][0] - 2.0).abs() < 1.0e-10)
  assert_true((result[1][0] - 7.0).abs() < 1.0e-10)
  assert_true((result[1][2] - 5.0).abs() < 1.0e-10)
}

test "openmm pbc distance" {
  let pos1 = [0.0, 0.0, 0.0]
  let pos2 = [9.0, 0.0, 0.0]
  let box = [10.0, 10.0, 10.0]
  let dist = @src.md_pbc_distance(pos1, pos2, box)
  assert_true((dist - 1.0).abs() < 1.0e-10)
}

test "openmm pbc no box no wrapping" {
  let dr = [100.0, 200.0, 300.0]
  let box = [0.0, 0.0, 0.0]
  let result = @src.md_pbc_minimum_image(dr, box)
  assert_true((result[0] - 100.0).abs() < 1.0e-10)
  assert_true((result[1] - 200.0).abs() < 1.0e-10)
  assert_true((result[2] - 300.0).abs() < 1.0e-10)
}

test "openmm gbsa obc energy two atoms" {
  let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
  let charges = [1.0, -1.0]
  let radii = [1.4, 1.4]
  let screen = [0.8, 0.8]
  let (e_gb, e_sa) = @src.md_gbsa_obc_energy(
    positions, charges, radii, screen, 78.5, 0.005,
  )
  // GB energy should be negative (attractive charges in dielectric)
  assert_true(e_gb < 0.0)
  // SA energy should be positive (surface area penalty)
  assert_true(e_sa >= 0.0)
}

test "openmm gbsa obc energy zero charges" {
  let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
  let charges = [0.0, 0.0]
  let radii = [1.4, 1.4]
  let screen = [0.8, 0.8]
  let (e_gb, _e_sa) = @src.md_gbsa_obc_energy(
    positions, charges, radii, screen, 78.5, 0.005,
  )
  // With zero charges, GB energy should be zero
  assert_true(e_gb.abs() < 1.0e-10)
}

test "openmm mc barostat step returns valid result" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let box = [10.0, 10.0, 10.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let energy_fn = fn(p : Array[Array[Double]]) -> Double {
    @src.md_bond_energy(p, bonds)
  }
  let (pos_new, box_new, _accepted) = @src.md_mc_barostat_step(
    positions, box, energy_fn, 1.0, 300.0, 2,
  )
  assert_eq(pos_new.length(), 2)
  assert_eq(box_new.length(), 3)
  // Box dimensions should be positive
  assert_true(box_new[0] > 0.0)
  assert_true(box_new[1] > 0.0)
  assert_true(box_new[2] > 0.0)
}

test "openmm seed rng produces reproducible sequence" {
  @src.md_seed_rng(42L)
  let v1 = @src.md_rand_double()
  @src.md_seed_rng(42L)
  let v2 = @src.md_rand_double()
  assert_true((v1 - v2).abs() < 1.0e-15)
}

test "openmm rand double in unit interval" {
  @src.md_seed_rng(123L)
  for _i = 0; _i < 100; _i = _i + 1 {
    let v = @src.md_rand_double()
    assert_true(v >= 0.0)
    assert_true(v < 1.0)
  }
}
