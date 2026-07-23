///|
/// Bio.Application - Wrapper for command-line tools.

///|
pub struct AbstractCommandline {
  executable : String
  arguments : Array[String]
  stdin : String
  stdout : String
  stderr : String
  env : Map[String, String]
}

///|
pub fn AbstractCommandline::new(executable : String) -> AbstractCommandline {
  AbstractCommandline::{
    executable,
    arguments: Array::new(),
    stdin: "",
    stdout: "",
    stderr: "",
    env: Map([], capacity=0)
  }
}

///|
pub fn AbstractCommandline::add_arg(self : AbstractCommandline, arg : String, value : String?) -> AbstractCommandline {
  let args = self.arguments.copy()
  args.push(arg)
  if value is Some(_) {
    args.push(value.unwrap())
  }
  AbstractCommandline::{
    executable: self.executable,
    arguments: args,
    stdin: self.stdin,
    stdout: self.stdout,
    stderr: self.stderr,
    env: self.env
  }
}

///|
pub fn AbstractCommandline::set_parameter(self : AbstractCommandline, param : String, value : String) -> AbstractCommandline {
  let args = self.arguments.copy()
  args.push(param)
  args.push(value)
  AbstractCommandline::{
    executable: self.executable,
    arguments: args,
    stdin: self.stdin,
    stdout: self.stdout,
    stderr: self.stderr,
    env: self.env
  }
}

///|
pub fn AbstractCommandline::set_stdout(self : AbstractCommandline, path : String) -> AbstractCommandline {
  AbstractCommandline::{
    executable: self.executable,
    arguments: self.arguments,
    stdin: self.stdin,
    stdout: path,
    stderr: self.stderr,
    env: self.env
  }
}

///|
pub fn AbstractCommandline::set_stderr(self : AbstractCommandline, path : String) -> AbstractCommandline {
  AbstractCommandline::{
    executable: self.executable,
    arguments: self.arguments,
    stdin: self.stdin,
    stdout: self.stdout,
    stderr: path,
    env: self.env
  }
}

///|
pub fn AbstractCommandline::set_stdin(self : AbstractCommandline, path : String) -> AbstractCommandline {
  AbstractCommandline::{
    executable: self.executable,
    arguments: self.arguments,
    stdin: path,
    stdout: self.stdout,
    stderr: self.stderr,
    env: self.env
  }
}

///|
pub fn AbstractCommandline::build_command(self : AbstractCommandline) -> String {
  let mut cmd = self.executable
  let mut i = 0
  while i < self.arguments.length() {
    cmd = cmd + " " + self.arguments[i]
    i = i + 1
  }
  cmd
}

///|
pub fn AbstractCommandline::run(self : AbstractCommandline) -> Int {
  0
}

///|
pub struct CommandlineError {
  message : String
  exit_code : Int
  command : String
}

///|
pub fn CommandlineError::new(message : String, exit_code : Int, command : String) -> CommandlineError {
  CommandlineError::{ message, exit_code, command }
}

///|
pub struct BioApplication {
  commandline : AbstractCommandline
}

///|
pub fn BioApplication::new(executable : String) -> BioApplication {
  BioApplication::{ commandline: AbstractCommandline::new(executable) }
}

///|
pub fn BioApplication::run(self : BioApplication) -> (String, String) {
  ("stdout output", "stderr output")
}

///|
pub fn create_example_commandline() -> AbstractCommandline {
  let mut cmd = AbstractCommandline::new("blastn")
  cmd = cmd.add_arg("-query", Some("query.fasta"))
  cmd = cmd.add_arg("-db", Some("nr"))
  cmd = cmd.add_arg("-out", Some("results.txt"))
  cmd = cmd.add_arg("-evalue", Some("1e-5"))
  cmd
}