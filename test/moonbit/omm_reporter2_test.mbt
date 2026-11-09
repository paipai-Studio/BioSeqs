///|
test "omm_pdb_frame_fixed_columns" {
  let coords = [(1.0, 0.0, 0.5)]
  let text = @src.omm_pdb_frame(coords, ["O"], ["HOH"], [1], 0)
  let expected = "MODEL     0\n" +
    "ATOM  " +
    "    1" +
    " " +
    "   O" +
    " " +
    "HOH" +
    " A" +
    "   1" +
    "    " +
    "       1.000" +
    "       0.000" +
    "       0.500" +
    "  1.00  0.00\n" +
    "ENDMDL\n"
  assert_eq(text, expected)
}

///|
test "omm_pdb_frame_defaults_for_missing_labels" {
  let coords = [(0.0, 0.0, 0.0)]
  let text = @src.omm_pdb_frame(coords, [], [], [], 3)
  // missing atom name / residue name fall back to "X", residue id to 1
  let expected = "MODEL     3\n" +
    "ATOM  " +
    "    1" +
    " " +
    "   X" +
    " " +
    "  X" +
    " A" +
    "   1" +
    "    " +
    "       0.000" +
    "       0.000" +
    "       0.000" +
    "  1.00  0.00\n" +
    "ENDMDL\n"
  assert_eq(text, expected)
}

///|
test "omm_pdb_reporter_interval_and_frames" {
  let rep = @src.OmmPdbReporter::new(interval=2)
  assert_true(rep.should_report())
  let coords = [(0.0, 0.0, 0.0)]
  rep.report(coords, ["O"], ["HOH"], [1])
  assert_false(rep.should_report())
  rep.report(coords, ["O"], ["HOH"], [1])
  assert_true(rep.should_report())
  // frames stored with their internal step counter
  assert_eq(rep.frames.length(), 2)
  assert_true(rep.frames[0].has_prefix("MODEL     0\n"))
  assert_true(rep.frames[1].has_prefix("MODEL     1\n"))
  assert_eq(rep.step_count, 2)
}

///|
test "omm_xtc_reporter_frames" {
  let rep = @src.OmmXtcReporter::new(interval=2)
  assert_true(rep.should_report())
  let c0 = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  rep.report(c0, box=(10.0, 10.0, 10.0))
  assert_eq(rep.n_frames(), 1)
  assert_false(rep.should_report())
  let c1 = [(2.0, 0.0, 0.0)]
  rep.report(c1)
  assert_eq(rep.n_frames(), 2)
  assert_true(rep.should_report())
  match rep.get_frame(0) {
    Some((step, coords)) => {
      assert_eq(step, 0)
      assert_eq(coords.length(), 2)
      assert_true(@moonbit.assert_num_close(coords[1].0, 1.0, 1.0e-12))
    }
    None => fail("expected frame 0")
  }
  match rep.get_frame(1) {
    Some((step, _)) => assert_eq(step, 1)
    None => fail("expected frame 1")
  }
  match rep.get_frame(2) {
    Some(_) => fail("frame 2 should not exist")
    None => ()
  }
  match rep.get_frame(-1) {
    Some(_) => fail("negative index must be rejected")
    None => ()
  }
}
