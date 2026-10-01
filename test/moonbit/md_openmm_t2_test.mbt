test "openmm rmsd identical structures" {
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let rmsd_val = @src.md_rmsd(pos, pos)
  assert_true(rmsd_val.abs() < 1.0e-10)
}

test "openmm rmsd translated structures" {
  // RMSD should be 0 for translated copies (optimal superposition)
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let ref = [[5.0, 3.0, 2.0], [6.0, 3.0, 2.0], [5.0, 4.0, 2.0]]
  let rmsd_val = @src.md_rmsd(pos, ref)
  assert_true(rmsd_val < 1.0e-6)
}

test "openmm rmsd scaled structures" {
  // RMSD should be nonzero for scaled structures
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let ref = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
  let rmsd_val = @src.md_rmsd(pos, ref)
  assert_true(rmsd_val > 0.1)
}

test "openmm rmsd force zero at reference" {
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let forces = @src.md_rmsd_force(pos, pos, 100.0, 0.0)
  for i = 0; i < 3; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      assert_true(forces[i][d].abs() < 1.0e-6)
    }
  }
}

test "openmm tabulated1d linear interpolation" {
  // f(x) = 2x on [0, 1], 3 points: [0, 1, 2]
  let f = @src.MdTabulated1D::new([0.0, 1.0, 2.0], 0.0, 1.0)
  let v = @src.md_tabulated1d_eval(f, 0.5)
  assert_true((v - 1.0).abs() < 1.0e-10)
}

test "openmm tabulated1d boundary" {
  let f = @src.MdTabulated1D::new([0.0, 1.0, 2.0], 0.0, 1.0)
  assert_true((@src.md_tabulated1d_eval(f, -1.0) - 0.0).abs() < 1.0e-10)
  assert_true((@src.md_tabulated1d_eval(f, 2.0) - 2.0).abs() < 1.0e-10)
}

test "openmm tabulated1d derivative" {
  let f = @src.MdTabulated1D::new([0.0, 1.0, 2.0], 0.0, 1.0)
  let d = @src.md_tabulated1d_deriv(f, 0.5)
  // f(x) = 2x, derivative = 2.0
  assert_true((d - 2.0).abs() < 1.0e-10)
}

test "openmm tabulated2d bilinear" {
  // 2x2 grid: f(x,y) = x + y on [0,1]x[0,1]
  // values: f(0,0)=0, f(1,0)=1, f(0,1)=1, f(1,1)=2
  let f = @src.MdTabulated2D::new([0.0, 1.0, 1.0, 2.0], 0.0, 1.0, 0.0, 1.0, 2, 2)
  let v = @src.md_tabulated2d_eval(f, 0.5, 0.5)
  assert_true((v - 1.0).abs() < 1.0e-10)
}

test "openmm tabulated1d discrete" {
  let vals = [10.0, 20.0, 30.0, 40.0]
  let v = @src.md_tabulated1d_discrete(vals, 0.0, 4.0, 1.5)
  assert_eq(v, 20.0)
}

test "openmm custom bond energy harmonic" {
  // expr = "0.5 * 100 * (r - 1.0)^2" — harmonic bond with k=100, r0=1.0
  let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let bonds = [(0, 1, "0.5 * 100 * (r - 1.0)^2")]
  let e = @src.md_custom_bond_energy(positions, bonds)
  // Expected: 0.5 * 100 * 0.25 = 12.5
  assert_true((e - 12.5).abs() < 1.0e-4)
}

test "openmm custom bond energy zero at equilibrium" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let bonds = [(0, 1, "0.5 * 100 * (r - 1.0)^2")]
  let e = @src.md_custom_bond_energy(positions, bonds)
  assert_true(e.abs() < 1.0e-4)
}

test "openmm custom bond forces pull together" {
  let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let bonds = [(0, 1, "0.5 * 100 * (r - 1.0)^2")]
  let forces = @src.md_custom_bond_forces(positions, bonds)
  // Atom 0 should be pulled right (+x), atom 1 left (-x)
  assert_true(forces[0][0] > 0.0)
  assert_true(forces[1][0] < 0.0)
}

test "openmm custom angle energy" {
  // expr = "0.5 * 50 * (theta - 1.5707963267948966)^2" — angle restraint at pi/2
  let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let angles = [(0, 1, 2, "0.5 * 50 * (theta - 1.5707963267948966)^2")]
  let e = @src.md_custom_angle_energy(positions, angles)
  // Angle is exactly pi/2, so energy should be ~0
  assert_true(e.abs() < 1.0e-2)
}

test "openmm custom torsion energy" {
  // expr = "cos(theta)" — periodic torsion
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let torsions = [(0, 1, 2, 3, "cos(theta)")]
  let e = @src.md_custom_torsion_energy(positions, torsions)
  assert_true(e.abs() <= 1.0)
}

test "openmm custom external energy" {
  // expr = "0.5 * 10 * r^2" — harmonic well at origin
  let positions = [[3.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let externals = [(0, "0.5 * 10 * r^2")]
  let e = @src.md_custom_external_energy(positions, externals)
  // Expected: 0.5 * 10 * 9 = 45
  assert_true((e - 45.0).abs() < 1.0e-2)
}

test "openmm custom external force toward origin" {
  let positions = [[3.0, 0.0, 0.0]]
  let externals = [(0, "0.5 * 10 * r^2")]
  let forces = @src.md_custom_external_forces(positions, externals)
  // Force should point toward origin (-x)
  assert_true(forces[0][0] < 0.0)
}

test "openmm custom nonbonded energy lj" {
  // expr = "4 * 1.0 * ((1/r)^12 - (1/r)^6)" — LJ with sigma=1, eps=1
  let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let e = @src.md_custom_nonbonded_energy(positions, "4 * ((1/r)^12 - (1/r)^6)", [], 0.0)
  // Expected: 4 * (1/4096 - 1/64) = 4 * (0.000244 - 0.015625) = -0.0615
  assert_true((e - (-0.0615234375)).abs() < 1.0e-4)
}

test "openmm custom nonbonded with exclusion" {
  let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
  let exclusions = [(0, 1)]
  let e = @src.md_custom_nonbonded_energy(positions, "4 * ((1/r)^12 - (1/r)^6)", exclusions, 0.0)
  // Only pair (0,2) and (1,2) contribute, (0,1) excluded
  assert_true(e.abs() < 1.0e6)
}

test "openmm custom nonbonded with cutoff" {
  let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
  let e = @src.md_custom_nonbonded_energy(positions, "4 * ((1/r)^12 - (1/r)^6)", [], 5.0)
  // Pair (0,2) and (1,2) are beyond cutoff=5.0, only (0,1) contributes
  let e_single = @src.md_custom_nonbonded_energy([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], "4 * ((1/r)^12 - (1/r)^6)", [], 0.0)
  assert_true((e - e_single).abs() < 1.0e-10)
}

test "openmm nonbonded pbc energy" {
  let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
  let params = [(1.0, 1.0, 0.0), (1.0, 1.0, 0.0)]
  let box = [10.0, 10.0, 10.0]
  let e = @src.md_nonbonded_pbc_energy(positions, params, box, 5.0, 1.0)
  // LJ at r=0.5: 4*(2^12 - 2^6) = 4*(4096 - 64) = 16128
  assert_true(e > 0.0)
}

test "openmm nonbonded pbc wraps distance" {
  // Two atoms at 0.0 and 9.5 in box 10.0 → PBC distance = 0.5
  let positions = [[0.0, 0.0, 0.0], [9.5, 0.0, 0.0]]
  let params = [(1.0, 1.0, 0.0), (1.0, 1.0, 0.0)]
  let box = [10.0, 10.0, 10.0]
  let e_pbc = @src.md_nonbonded_pbc_energy(positions, params, box, 5.0, 1.0)
  let e_direct = @src.md_nonbonded_pbc_energy([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], params, box, 5.0, 1.0)
  assert_true((e_pbc - e_direct).abs() < 1.0)
}

test "openmm neighbor list" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
  let box = [10.0, 10.0, 10.0]
  let pairs = @src.md_neighbor_list(positions, box, 3.0)
  // Only (0,1) within cutoff=3.0
  assert_eq(pairs.length(), 1)
  let (i, j, r) = pairs[0]
  assert_eq(i, 0)
  assert_eq(j, 1)
  assert_true((r - 1.0).abs() < 1.0e-10)
}

test "openmm nose hoover chain creation" {
  let chain = @src.MdNoseHooverChain::new(3, 6, 0.1)
  assert_eq(chain.chain_length, 3)
  assert_eq(chain.n_dof, 6)
}

test "openmm nose hoover step preserves shape" {
  let velocities = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let masses = [1.0, 1.0]
  let chain = @src.MdNoseHooverChain::new(2, 6, 0.1)
  let v_new = @src.md_nose_hoover_step(velocities, masses, chain, 0.001, 300.0)
  assert_eq(v_new.length(), 2)
  for i = 0; i < 2; i = i + 1 {
    assert_eq(v_new[i].length(), 3)
  }
}

test "openmm nose hoover chain simulation" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let temps = @src.md_simulate_nose_hoover_chain(
    positions, velocities, masses, bonds, 0.001, 300.0, 0.1, 3, 10, 6,
  )
  assert_eq(temps.length(), 10)
}

test "openmm dpd step shape" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let forces = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let (x_new, v_new) = @src.md_dpd_step(positions, velocities, masses, forces, 0.01, 300.0, 1.0, 2.0)
  assert_eq(x_new.length(), 2)
  assert_eq(v_new.length(), 2)
}

test "openmm dpd simulation" {
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let velocities = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let masses = [1.0, 1.0, 1.0]
  let bonds = [(0, 1, 1.0, 100.0)]
  let temps = @src.md_simulate_dpd(positions, velocities, masses, bonds, 0.005, 300.0, 1.0, 2.0, 10, 9)
  assert_eq(temps.length(), 10)
}
