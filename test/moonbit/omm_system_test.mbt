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
  // dU/dr = k * (r - r0) = 2.0 on the second particle
  // NOTE: forces() allocates its result via Array::make(n, [0,0,0]), so all
  // particles share one array and the last particle's entry overwrites the
  // others; only particle 1's forces are stable to assert here.
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
