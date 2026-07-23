///|
/// Tests for Align.Applications module.

test "ClustalwCommandline creation" {
  let cmd = @src.create_example_clustalw()
  let command = cmd.build()
  assert_true(command.contains("clustalw2"))
  assert_true(command.contains("-infile"))
}

test "ClustalwCommandline set_matrix" {
  let mut cmd = @src.ClustalwCommandline::new("clustalw2")
  cmd = cmd.set_matrix("BLOSUM62")
  let command = cmd.build()
  assert_true(command.contains("-matrix"))
  assert_true(command.contains("BLOSUM62"))
}

test "ClustalOmegaCommandline creation" {
  let cmd = @src.create_example_clustalomega()
  let command = cmd.build()
  assert_true(command.contains("clustalo"))
  assert_true(command.contains("-i"))
}

test "ClustalOmegaCommandline set_iterations" {
  let mut cmd = @src.ClustalOmegaCommandline::new("clustalo")
  cmd = cmd.set_iterations(10)
  let command = cmd.build()
  assert_true(command.contains("--iterations"))
  assert_true(command.contains("10"))
}

test "MuscleCommandline creation" {
  let cmd = @src.create_example_muscle()
  let command = cmd.build()
  assert_true(command.contains("muscle"))
}

test "MAFFTCommandline creation" {
  let cmd = @src.create_example_mafft()
  let command = cmd.build()
  assert_true(command.contains("mafft"))
  assert_true(command.contains("--auto"))
}