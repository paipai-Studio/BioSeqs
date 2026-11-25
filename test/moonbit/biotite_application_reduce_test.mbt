///|
/// Tests for biotite_application_reduce.mbt.
test "biotite_reduce_command" {
  let command = @src.biotite_reduce_command("in.pdb", "out.pdb")
  assert_true(command.contains("reduce"))
  assert_true(command.contains("-build"))
  assert_true(command.contains("-pdb in.pdb"))
  assert_true(command.contains("> out.pdb"))
}

///|
test "biotite_reduce_trim_command" {
  let command = @src.biotite_reduce_command("in.pdb", "out.pdb", build=false)
  assert_true(command.contains("-trim"))
  assert_true(!command.contains("-build"))
}

///|
test "biotite_reduce_hydrogen_name" {
  assert_eq(@src.biotite_reduce_hydrogen_name("N", 0), "HN")
  assert_eq(@src.biotite_reduce_hydrogen_name("N", 1), "H1N")
  assert_eq(@src.biotite_reduce_hydrogen_name("CA", 2), "H2CA")
}

///|
test "biotite_reduce_polar_only_command" {
  let command = @src.biotite_reduce_polar_only_command("in.pdb", "out.pdb")
  assert_true(command.contains("-noflip"))
  assert_true(command.contains("-build"))
}
