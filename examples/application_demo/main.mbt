///|
/// Application Demo - Demonstrates command-line tool wrappers.

fn main {
  println("=== Application Demo ===")
  
  let mut cmd = @bio.AbstractCommandline::new("blastn")
  cmd = cmd.add_arg("-query", Some("query.fasta"))
  cmd = cmd.add_arg("-db", Some("nr"))
  cmd = cmd.add_arg("-out", Some("results.txt"))
  cmd = cmd.add_arg("-evalue", Some("1e-5"))
  
  println("Command: \{cmd.build_command()}")
  println("Executable: \{cmd.executable}")
  println("Arguments count: \{cmd.arguments.length()}")
}