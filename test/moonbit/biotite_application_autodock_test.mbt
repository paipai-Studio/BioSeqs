///|
/// Tests for biotite_application_autodock.mbt.
test "biotite_vina_config" {
  let config = @src.biotite_vina_config(
    "receptor.pdbqt", "ligand.pdbqt", 1.0, 2.0, 3.0, 20, 20, 20, "out.pdbqt",
  )
  assert_true(config.contains("receptor = receptor.pdbqt"))
  assert_true(config.contains("ligand = ligand.pdbqt"))
  assert_true(config.contains("center_x = 1"))
  assert_true(config.contains("size_x = 20"))
  assert_true(config.contains("out = out.pdbqt"))
}

///|
test "biotite_vina_command" {
  assert_eq(@src.biotite_vina_command("vina.conf"), "vina --config vina.conf")
  assert_eq(
    @src.biotite_vina_command("vina.conf", cpu=4),
    "vina --config vina.conf --cpu 4",
  )
}

///|
test "biotite_vina_parse_output" {
  let text = "mode | affinity | dist from best mode\n-----+------------+----------+----\n   1       -7.5      0.000      0.000\n   2       -7.0      1.200      2.100\n"
  let poses = @src.biotite_vina_parse_output(text)
  assert_eq(poses.length(), 2)
  assert_eq(poses[0].0, 1)
  assert_true(@moonbit.assert_num_close(poses[0].1, -7.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(poses[1].1, -7.0, 1.0e-12))
}

///|
test "biotite_vina_best_pose" {
  let poses : Array[(Int32, Double, Double, Double)] = [
    (1, -7.5, 0.0, 0.0),
    (2, -8.0, 1.0, 2.0),
  ]
  match @src.biotite_vina_best_pose(poses) {
    Some(best) => assert_eq(best.0, 2)
    None => fail("expected best pose")
  }
  assert_true(@src.biotite_vina_best_pose([]) is None)
}
