///|
test "omm virial tensor single particle" {
  let coords = [(1.0, 2.0, 3.0)]
  let forces = [(4.0, 5.0, 6.0)]
  let v = @src.omm_virial_tensor(coords, forces)
  assert_eq(v.length(), 3)
  assert_true(@moonbit.assert_num_close(v[0][0], 4.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[0][1], 5.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[0][2], 6.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[1][0], 8.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[1][1], 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[1][2], 12.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[2][0], 12.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[2][1], 15.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[2][2], 18.0, 1.0e-12))
}

///|
test "omm virial tensor sums over particles" {
  let coords = [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
  let forces = [(1.0, 2.0, 0.0), (3.0, 4.0, 0.0)]
  let v = @src.omm_virial_tensor(coords, forces)
  // V_xx = 1*1 + 2*3 = 7, V_xy = 1*2 + 2*4 = 10
  assert_true(@moonbit.assert_num_close(v[0][0], 7.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[0][1], 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[2][2], 0.0, 1.0e-12))
}

///|
test "omm pressure tensor kinetic plus virial over volume" {
  // one particle m=2, v=(1,0,0), no forces, volume=4
  // P_xx = m*vx^2/V = 2/4 = 0.5, off-diagonals zero
  let masses = [2.0]
  let velocities = [(1.0, 0.0, 0.0)]
  let coords = [(0.0, 0.0, 0.0)]
  let forces = [(0.0, 0.0, 0.0)]
  let p = @src.omm_pressure_tensor(masses, velocities, coords, forces, 4.0)
  // implementation uses inv_vol = 1/(V + 1e-10) so allow small bias
  assert_true(@moonbit.assert_num_close(p[0][0], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(p[1][1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(p[0][1], 0.0, 1.0e-9))
  let scalar = @src.omm_scalar_pressure(p)
  assert_true(@moonbit.assert_num_close(scalar, 0.5 / 3.0, 1.0e-9))
}

///|
test "omm scalar pressure averages diagonal" {
  let p = [[3.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 9.0]]
  assert_true(
    @moonbit.assert_num_close(@src.omm_scalar_pressure(p), 6.0, 1.0e-12),
  )
}

///|
test "omm moment of inertia single particle on x axis" {
  // m=1 at (1,0,0): Ixx = 0, Iyy = Izz = 1, off-diagonal zero
  let masses = [1.0]
  let coords = [(1.0, 0.0, 0.0)]
  let moi = @src.omm_moment_of_inertia(masses, coords)
  assert_true(@moonbit.assert_num_close(moi[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[1][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[2][2], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[0][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[0][2], 0.0, 1.0e-12))
}

///|
test "omm moment of inertia diatomic along x" {
  // m=1 at (+-1,0,0): Ixx = 0, Iyy = Izz = 2
  let masses = [1.0, 1.0]
  let coords = [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)]
  let moi = @src.omm_moment_of_inertia(masses, coords)
  assert_true(@moonbit.assert_num_close(moi[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[1][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moi[2][2], 2.0, 1.0e-12))
}

///|
test "omm principal moments return diagonal elements" {
  let masses = [1.0]
  let coords = [(1.0, 0.0, 0.0)]
  let (i1, i2, i3) = @src.omm_principal_moments(masses, coords)
  assert_true(@moonbit.assert_num_close(i1, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(i2, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(i3, 1.0, 1.0e-12))
}

///|
test "omm decompose energy sums contributions" {
  let d = @src.omm_decompose_energy(1.0, 2.0, 3.0, 4.0, restraint_energy=5.0)
  assert_true(@moonbit.assert_num_close(d.total, 15.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.bond_energy, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.angle_energy, 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.torsion_energy, 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.nonbonded_energy, 4.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.restraint_energy, 5.0, 1.0e-12))
}

///|
test "omm decompose energy default restraint is zero" {
  let d = @src.omm_decompose_energy(1.0, 1.0, 1.0, 1.0)
  assert_true(@moonbit.assert_num_close(d.total, 4.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d.restraint_energy, 0.0, 1.0e-12))
}

///|
test "omm per atom energy coulomb pair symmetric" {
  // two atoms q=+1/-1 at r=2, no LJ -> per-atom contribution is half the pair energy
  let charges = [1.0, -1.0]
  let coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
  let epsilon = [0.0, 0.0]
  let sigma = [1.0, 1.0]
  let energies = @src.omm_per_atom_energy(charges, coords, epsilon, sigma)
  let expected = 0.5 * 138.935485 * 1.0 * -1.0 / 2.0
  assert_true(@moonbit.assert_num_close(energies[0], expected, 1.0e-9))
  assert_true(@moonbit.assert_num_close(energies[1], expected, 1.0e-9))
}

///|
test "omm per atom energy lj minimum value" {
  // r = 2^(1/6) sigma, epsilon=1 -> pair LJ energy = -1, per atom = -0.5
  let charges = [0.0, 0.0]
  let r = 1.122462048309373
  let coords = [(0.0, 0.0, 0.0), (r, 0.0, 0.0)]
  let epsilon = [1.0, 1.0]
  let sigma = [1.0, 1.0]
  let energies = @src.omm_per_atom_energy(charges, coords, epsilon, sigma)
  assert_true(@moonbit.assert_num_close(energies[0], -0.5, 1.0e-6))
  assert_true(@moonbit.assert_num_close(energies[1], -0.5, 1.0e-6))
}

///|
test "omm per atom energy self contribution zero" {
  let charges = [1.0, 1.0]
  let coords = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  let epsilon = [0.0, 0.0]
  let sigma = [1.0, 1.0]
  // coincident atoms are skipped (r ~ 0 guard)
  let energies = @src.omm_per_atom_energy(charges, coords, epsilon, sigma)
  assert_true(@moonbit.assert_num_close(energies[0], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(energies[1], 0.0, 1.0e-9))
}
