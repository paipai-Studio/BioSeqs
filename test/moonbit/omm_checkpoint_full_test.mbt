// ============================================================================
// O1: CheckpointReporter（OmmFullCheckpointReporter + save/load/restore）单元测试
// ============================================================================

// --- OmmFullCheckpointReporter: 基本存/取 ---

test "omm_full_checkpoint_reporter_new_and_interval" {
  let rep = @src.OmmFullCheckpointReporter::new(interval=1000)
  assert_true(rep.should_report(0))
  assert_false(rep.should_report(500))
  assert_true(rep.should_report(1000))
  assert_true(rep.should_report(2000))
  assert_false(rep.should_report(2500))
  match rep.last() {
    Some(_) => fail("expected no checkpoint before reporting")
    None => ()
  }
}

test "omm_full_checkpoint_reporter_collects_single" {
  let rep = @src.OmmFullCheckpointReporter::new(interval=10)
  let pos = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let vel = [[0.1, 0.0, 0.0], [0.0, -0.1, 0.0]]
  let box = [50.0, 50.0, 50.0]
  rep.report(100, 0.5, box, pos, vel)
  match rep.last() {
    Some(cp) => {
      assert_eq(cp.step, 100)
      assert_true(@moonbit.assert_num_close(cp.time, 0.5, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.box[0], 50.0, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.positions[0][0], 1.0, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.positions[1][2], 6.0, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.velocities[0][0], 0.1, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.velocities[1][1], -0.1, 1.0e-10))
    }
    None => fail("expected a checkpoint")
  }
}

test "omm_full_checkpoint_reporter_stores_copy" {
  let rep = @src.OmmFullCheckpointReporter::new(interval=10)
  let pos = [[1.0, 0.0, 0.0]]
  let vel = [[0.5, 0.0, 0.0]]
  rep.report(1, 0.001, [10.0, 10.0, 10.0], pos, vel)
  pos[0][0] = 42.0
  vel[0][0] = -99.0
  match rep.last() {
    Some(cp) => {
      assert_true(@moonbit.assert_num_close(cp.positions[0][0], 1.0, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.velocities[0][0], 0.5, 1.0e-10))
    }
    None => fail("expected stored checkpoint copy")
  }
}

test "omm_full_checkpoint_reporter_multiple_frames" {
  let rep = @src.OmmFullCheckpointReporter::new(interval=10)
  rep.report(0, 0.0, [10.0, 10.0, 10.0], [[1.0, 0.0, 0.0]], [[0.1, 0.0, 0.0]])
  rep.report(10, 0.01, [10.0, 10.0, 10.0], [[2.0, 0.0, 0.0]], [[0.2, 0.0, 0.0]])
  rep.report(20, 0.02, [10.0, 10.0, 10.0], [[3.0, 0.0, 0.0]], [[0.3, 0.0, 0.0]])
  match rep.last() {
    Some(cp) => {
      assert_eq(cp.step, 20)
      assert_true(@moonbit.assert_num_close(cp.time, 0.02, 1.0e-10))
      assert_true(@moonbit.assert_num_close(cp.positions[0][0], 3.0, 1.0e-10))
    }
    None => fail("expected final checkpoint from multiple frames")
  }
}

// --- serialize_last / omm_checkpoint_save ---

test "omm_checkpoint_save_format" {
  let rep = @src.OmmFullCheckpointReporter::new(interval=1)
  rep.report(100, 0.5, [50.0, 50.0, 50.0],
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [[0.1, 0.0, 0.0], [0.0, -0.1, 0.0]])
  let txt = rep.serialize_last()
  assert_true(txt.length() > 0)
  assert_true(txt.has_prefix("OMMCHECKPOINT"))
  assert_true(txt.has_prefix("OMMCHECKPOINT 1\nNATOMS") || txt.has_prefix("OMMCHECKPOINT 1"))
  // checkpoint 文本包含 BOX/STEP 关键字
  let has_box = txt.contains("BOX")
  let has_step = txt.contains("STEP")
  assert_true(has_box)
  assert_true(has_step)
}

// --- omm_checkpoint_load / omm_checkpoint_restore roundtrip ---

test "omm_checkpoint_load_roundtrip" {
  let pdb_text =
    "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n" +
    "ATOM      2  H1  HOH A   1       0.9572  0.000   0.000  1.00  0.00           H\n" +
    "ATOM      3  H2  HOH A   1      -0.2879  0.9120  0.000  1.00  0.00           H\n"
  let u = @src.universe_from_pdb_text(pdb_text)
  let (sim, _sys) = @src.Universe::to_omm_simulation(u, temperature=300.0)
  let reporter = @src.OmmFullCheckpointReporter::new(interval=100)
  let mut step = 0
  while step < 500 {
    sim.step(n_steps=100)
    let st = sim.state()
    reporter.report(st.step, st.time, st.box, st.positions, st.velocities)
    step = step + 100
  }
  let last_ckpt = reporter.serialize_last()
  assert_true(last_ckpt.has_prefix("OMMCHECKPOINT"))
  let rep2 = @src.omm_checkpoint_restore(last_ckpt)
  match rep2.last() {
    Some(cp) => {
      assert_true(cp.time > 0.0)
      assert_eq(cp.positions.length(), 3)
      assert_true(@moonbit.assert_num_close(cp.box[0], sim.state().box[0], 1.0e-6))
    }
    None => fail("roundtrip should have one checkpoint")
  }
}

test "omm_checkpoint_load_garbage_input" {
  let rep = @src.omm_checkpoint_restore("not a checkpoint at all\n")
  match rep.last() {
    Some(_) => fail("garbage should produce no checkpoint")
    None => ()
  }
}

// --- omm_simulation_from_checkpoint ---

test "omm_simulation_from_checkpoint_returns_sim" {
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  let ckpt =
    "OMMCHECKPOINT 1\n" +
    "NATOMS 2\n" +
    "BOX 50.0 50.0 50.0\n" +
    "TIME 1.5\n" +
    "STEP 1500\n" +
    "PE -42.0\n" +
    "KE 1.23\n" +
    "1.0 2.0 3.0 0.1 0.2 0.3\n" +
    "4.0 5.0 6.0 0.4 0.5 0.6\n"
  let sim = @src.omm_simulation_from_checkpoint(sys, ckpt)
  let st = sim.state()
  assert_eq(st.positions.length(), 2)
  assert_true(@moonbit.assert_num_close(st.positions[0][0], 1.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(st.positions[1][2], 6.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(st.velocities[0][1], 0.2, 1.0e-10))
  assert_true(@moonbit.assert_num_close(st.velocities[1][0], 0.4, 1.0e-10))
}

test "omm_simulation_from_checkpoint_bogus_returns_fallback" {
  let sys = @src.OmmSystem::new()
  let sim = @src.omm_simulation_from_checkpoint(sys, "garbage\n")
  assert_eq(sim.state().positions.length(), 0)
}
