///|
/// Tests for Application module.

test "AbstractCommandline creation" {
  let cmd = @src.create_example_commandline()
  assert_eq(cmd.executable, "blastn")
}

test "AbstractCommandline add_arg" {
  let mut cmd = @src.AbstractCommandline::new("blastn")
  cmd = cmd.add_arg("-query", Some("query.fasta"))
  assert_eq(cmd.arguments.length(), 2)
}

test "AbstractCommandline add_arg no value" {
  let mut cmd = @src.AbstractCommandline::new("blastn")
  cmd = cmd.add_arg("-help", None)
  assert_eq(cmd.arguments.length(), 1)
}

test "AbstractCommandline build_command" {
  let mut cmd = @src.AbstractCommandline::new("blastn")
  cmd = cmd.add_arg("-query", Some("query.fasta"))
  cmd = cmd.add_arg("-db", Some("nr"))
  let command = cmd.build_command()
  assert_true(command.contains("blastn"))
  assert_true(command.contains("-query"))
  assert_true(command.contains("-db"))
}

test "AbstractCommandline set_stdout" {
  let mut cmd = @src.AbstractCommandline::new("blastn")
  cmd = cmd.set_stdout("output.txt")
  assert_eq(cmd.stdout, "output.txt")
}