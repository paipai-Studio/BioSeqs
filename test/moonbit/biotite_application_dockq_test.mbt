///|
/// Tests for biotite_application_dockq.mbt.
test "biotite_dockq_score" {
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_dockq_score(0.0, 0.0, 1.0),
      1.0,
      1.0e-12,
    ),
  )
  let mid = @src.biotite_dockq_score(8.5, 1.5, 0.0)
  assert_true(@moonbit.assert_num_close(mid, 1.0 / 3.0, 1.0e-12))
}

///|
test "biotite_dockq_class" {
  assert_eq(@src.biotite_dockq_class(0.9), "High")
  assert_eq(@src.biotite_dockq_class(0.5), "Medium")
  assert_eq(@src.biotite_dockq_class(0.3), "Acceptable")
  assert_eq(@src.biotite_dockq_class(0.1), "Incorrect")
}

///|
test "biotite_dockq_command_and_parse" {
  assert_eq(
    @src.biotite_dockq_command("model.pdb", "native.pdb"),
    "DockQ model.pdb native.pdb",
  )
  match @src.biotite_dockq_parse_output("Total DockQ 0.742\n") {
    Some(v) => assert_true(@moonbit.assert_num_close(v, 0.742, 1.0e-12))
    None => fail("expected DockQ value")
  }
  assert_true(@src.biotite_dockq_parse_output("no score\n") is None)
}
