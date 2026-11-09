///|
test "simulation_minimize_lowers_energy" {
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=12.0)))
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  sys.add_bond(0, 1, 100.0, 1.0)
  let sim = @src.OmmSimulation::new(
    system=sys,
    dt=0.001,
    positions=Some([[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]]),
  )
  let e0 = sys.potential_energy([[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]])
  let dropped = sim.minimize(max_iter=300, tolerance=0.001)
  // Initial strain: 0.5 * 100 * 0.6^2 = 18 kcal/mol.
  assert_true(dropped > 10.0)
  let p = sim.state().positions
  let r = @src.md_dist(p[0], p[1])
  assert_true(@moonbit.assert_num_close(r, 1.0, 0.05))
  let e1 = sim.state().potential_energy
  assert_true(e1 < e0)
  // Velocities are reset by minimization.
  assert_eq(sim.state().kinetic_energy, 0.0)
}

///|
test "simulation_nve_small_drift" {
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=12.0)))
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  sys.add_bond(0, 1, 100.0, 1.0)
  let sim = @src.OmmSimulation::new(
    system=sys,
    dt=0.001,
    temperature=300.0,
    friction=0.0,
    positions=Some([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]),
  )
  sim.initialize_velocities(temperature=300.0, seed=7L)
  let e0 = sim.state().total_energy()
  sim.step(n_steps=100)
  let e1 = sim.state().total_energy()
  let scale = if e0.abs() > 1.0e-8 { e0.abs() } else { 1.0 }
  assert_true((e1 - e0).abs() / scale < 0.1)
  assert_eq(sim.state().step, 100)
  assert_true(@moonbit.assert_num_close(sim.state().time, 0.1, 1.0e-9))
}

///|
test "simulation_equilibrate_and_production" {
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=12.0)))
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  sys.add_bond(0, 1, 100.0, 1.0)
  let sim = @src.OmmSimulation::new(
    system=sys,
    dt=0.001,
    positions=Some([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]),
  )
  // Equilibration without reporting returns no rows but advances state.
  let eq = sim.equilibrate(n_steps=5, report_interval=0)
  assert_eq(eq.length(), 0)
  assert_eq(sim.state().step, 5)
  // Production reports every 2 steps -> steps 2, 4, 6.
  let prod = sim.production(n_steps=6, report_interval=2)
  assert_eq(prod.length(), 3)
  assert_eq(prod[0].step, 7)
  assert_eq(prod[2].step, 11)
  // The state snapshot is a deep copy: mutating it leaves the sim intact.
  let snap = sim.state()
  snap.positions[0][0] = 999.0
  assert_true(sim.state().positions[0][0] < 100.0)
}

///|
test "checkpoint_roundtrip" {
  let pos : Array[Array[Double]] = [[0.1, 0.2, 0.3], [1.1, 1.2, 1.3]]
  let vel : Array[Array[Double]] = [[0.01, 0.0, -0.02], [0.05, 0.1, 0.0]]
  let st = @src.OmmState::new(
    positions=pos,
    velocities=vel,
    box=[9.0, 9.0, 9.0],
    time=1.5,
    step=150,
    potential_energy=-3.5,
    kinetic_energy=4.25,
  )
  let text = @src.omm_checkpoint_save(st)
  let back = @src.omm_checkpoint_load(text)
  match back {
    Ok(s2) => {
      assert_eq(s2.n_atoms(), 2)
      assert_eq(s2.step, 150)
      assert_true(@moonbit.assert_num_close(s2.time, 1.5, 1.0e-9))
      assert_true(@moonbit.assert_num_close(s2.potential_energy, -3.5, 1.0e-9))
      assert_true(@moonbit.assert_num_close(s2.kinetic_energy, 4.25, 1.0e-9))
      for i = 0; i < 2; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          assert_true(
            @moonbit.assert_num_close(s2.positions[i][d], pos[i][d], 1.0e-9),
          )
          assert_true(
            @moonbit.assert_num_close(s2.velocities[i][d], vel[i][d], 1.0e-9),
          )
        }
      }
      assert_true(@moonbit.assert_num_close(s2.total_energy(), 0.75, 1.0e-9))
    }
    Err(msg) => fail("checkpoint roundtrip failed: " + msg)
  }
}

///|
test "checkpoint_load_rejects_garbage" {
  match @src.omm_checkpoint_load("not a checkpoint") {
    Ok(_) => fail("garbage should not parse")
    Err(_) => ()
  }
  // Header present but no NATOMS.
  match @src.omm_checkpoint_load("OMMCHECKPOINT 1\nBOX 1 1 1\n") {
    Ok(_) => fail("missing NATOMS should fail")
    Err(_) => ()
  }
  // Declared atom count mismatch.
  let bad = "OMMCHECKPOINT 1\nNATOMS 2\nBOX 0 0 0\nTIME 0\nSTEP 0\nPE 0\nKE 0\n0 0 0 0 0 0\n"
  match @src.omm_checkpoint_load(bad) {
    Ok(_) => fail("count mismatch should fail")
    Err(_) => ()
  }
}

///|
test "omm_state_kinetic_and_temperature" {
  // Single particle, mass 1, v = (1,0,0) -> KE = 0.5.
  let st = @src.OmmState::new(positions=[[0.0, 0.0, 0.0]], velocities=[
    [1.0, 0.0, 0.0],
  ])
  let ke = @src.omm_kinetic_energy(st.velocities, [1.0])
  assert_true(@moonbit.assert_num_close(ke, 0.5, 1.0e-12))
  // T = 2 KE / (N_dof kB); 3 DOF -> 2*0.5/(3*0.001987204) ~ 167.7 K
  let t = @src.omm_temperature_from_energy(ke, 3)
  assert_true(@moonbit.assert_num_close(t, 167.740, 0.01))
}
