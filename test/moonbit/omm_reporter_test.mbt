///|
test "omm_state_reporter_collects_data" {
  let rep = @src.OmmStateReporter::new(interval=10)
  assert_true(rep.should_report(0))
  assert_false(rep.should_report(5))
  assert_true(rep.should_report(10))
  rep.report(0, 0.0, -100.0, 50.0, 300.0, volume=1.0)
  rep.report(10, 0.01, -80.0, 45.0, 295.0)
  assert_eq(rep.n_reports(), 2)
  let d = rep.data()
  assert_eq(d[0].step, 0)
  assert_true(@moonbit.assert_num_close(d[0].time, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[0].potential_energy, -100.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[0].kinetic_energy, 50.0, 1.0e-12))
  // total energy is derived as PE + KE
  assert_true(@moonbit.assert_num_close(d[0].total_energy, -50.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[0].temperature, 300.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[0].volume, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[1].total_energy, -35.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(d[1].temperature, 295.0, 1.0e-12))
}

///|
test "omm_format_state_table_layout" {
  let rep = @src.OmmStateReporter::new(interval=1)
  rep.report(0, 0.0, -100.0, 50.0, 300.0)
  rep.report(1, 0.001, -90.0, 40.0, 310.0)
  let table = @src.omm_format_state_table(rep.data())
  let lines = table.split("\n").to_array()
  assert_eq(
    lines[0],
    "Step    Time(ps)    PE(kJ/mol)    KE(kJ/mol)    Total(kJ/mol)    Temp(K)",
  )
  // one line per data record, each starting with its step number
  assert_true(lines[1].has_prefix("0    "))
  assert_true(lines[2].has_prefix("1    "))
}

///|
test "omm_dcd_reporter_frames" {
  let rep = @src.OmmDCDReporter::new(interval=5)
  assert_true(rep.should_report(0))
  assert_false(rep.should_report(3))
  assert_true(rep.should_report(5))
  let p1 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
  rep.report(p1)
  p1[0][0] = 99.0
  rep.report(p1)
  assert_eq(rep.n_frames(), 2)
  // frames are defensive copies
  assert_true(@moonbit.assert_num_close(rep.frames()[0][0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rep.frames()[1][0][0], 99.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rep.frames()[0][1][1], 2.0, 1.0e-12))
}

///|
test "omm_checkpoint_reporter_last_checkpoint" {
  let rep = @src.OmmCheckpointReporter::new(interval=100)
  assert_true(rep.should_report(0))
  assert_false(rep.should_report(50))
  match rep.last_checkpoint() {
    Some(_) => fail("expected no checkpoint before reporting")
    None => ()
  }
  let pos = [[1.0, 1.0, 1.0]]
  let vel = [[0.5, 0.0, 0.0]]
  rep.report(pos, vel)
  pos[0][0] = 42.0
  match rep.last_checkpoint() {
    Some((p, v)) => {
      // checkpoint stores a copy, not a reference
      assert_true(@moonbit.assert_num_close(p[0][0], 1.0, 1.0e-12))
      assert_true(@moonbit.assert_num_close(p[0][2], 1.0, 1.0e-12))
      assert_true(@moonbit.assert_num_close(v[0][0], 0.5, 1.0e-12))
    }
    None => fail("expected a stored checkpoint")
  }
}
