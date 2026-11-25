///|
/// Tests for biotite_structure_traj.mbt.
fn bt_traj_pdb(x1 : Double, x2 : Double) -> String {
  let a1 = "ATOM      1  CA  ALA A   1      " +
    x1.to_string() +
    "   0.000   0.000  1.00 20.00           C"
  let a2 = "ATOM      2  CB  ALA A   1      " +
    x2.to_string() +
    "   0.000   0.000  1.00 20.00           C"
  "MODEL        1\n" +
  a1 +
  "\n" +
  a2 +
  "\nENDMDL\nMODEL        2\n" +
  a1 +
  "\n" +
  a2 +
  "\nENDMDL\nEND\n"
}

///|
test "biotite_trajectory_from_pdb" {
  let traj = @src.biotite_trajectory_from_pdb(bt_traj_pdb(1.0, 2.0))
  assert_eq(@src.biotite_trajectory_n_models(traj), 2)
  assert_eq(@src.biotite_trajectory_get_model(traj, 0).len(), 2)
  assert_eq(@src.biotite_trajectory_get_model(traj, 9).len(), 0)
}

///|
test "biotite_trajectory_concatenate" {
  let traj = @src.biotite_trajectory_from_pdb(bt_traj_pdb(1.0, 2.0))
  let combined = @src.biotite_trajectory_concatenate([traj, traj])
  assert_eq(@src.biotite_trajectory_n_models(combined), 4)
}

///|
test "biotite_trajectory_to_pdb" {
  let traj = @src.biotite_trajectory_from_pdb(bt_traj_pdb(1.0, 2.0))
  let text = @src.biotite_trajectory_to_pdb(traj)
  assert_true(text.has_prefix("MODEL"))
  assert_true(text.contains("ENDMDL"))
}

///|
test "biotite_trajectory_single_model" {
  let pdb = "ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C\nEND\n"
  let traj = @src.biotite_trajectory_from_pdb(pdb)
  assert_eq(@src.biotite_trajectory_n_models(traj), 1)
  assert_eq(@src.biotite_trajectory_get_model(traj, 0).len(), 1)
}
