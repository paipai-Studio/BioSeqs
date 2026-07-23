///|
/// Bio.Align.Applications - Wrappers for multiple sequence alignment tools.

///|
pub struct ClustalwCommandline {
  commandline : AbstractCommandline
}

///|
pub fn ClustalwCommandline::new(executable : String) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: AbstractCommandline::new(executable) }
}

///|
pub fn ClustalwCommandline::set_input(self : ClustalwCommandline, infile : String) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-infile", Some(infile)) }
}

///|
pub fn ClustalwCommandline::set_output(self : ClustalwCommandline, outfile : String) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-outfile", Some(outfile)) }
}

///|
pub fn ClustalwCommandline::set_output_format(self : ClustalwCommandline, format : String) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-output", Some(format)) }
}

///|
pub fn ClustalwCommandline::set_matrix(self : ClustalwCommandline, matrix : String) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-matrix", Some(matrix)) }
}

///|
pub fn ClustalwCommandline::set_gap_open(self : ClustalwCommandline, penalty : Double) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-gapopen", Some(penalty.to_string())) }
}

///|
pub fn ClustalwCommandline::set_gap_extend(self : ClustalwCommandline, penalty : Double) -> ClustalwCommandline {
  ClustalwCommandline::{ commandline: self.commandline.add_arg("-gapext", Some(penalty.to_string())) }
}

///|
pub fn ClustalwCommandline::build(self : ClustalwCommandline) -> String {
  self.commandline.build_command()
}

///|
pub struct ClustalOmegaCommandline {
  commandline : AbstractCommandline
}

///|
pub fn ClustalOmegaCommandline::new(executable : String) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: AbstractCommandline::new(executable) }
}

///|
pub fn ClustalOmegaCommandline::set_input(self : ClustalOmegaCommandline, infile : String) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: self.commandline.add_arg("-i", Some(infile)) }
}

///|
pub fn ClustalOmegaCommandline::set_output(self : ClustalOmegaCommandline, outfile : String) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: self.commandline.add_arg("-o", Some(outfile)) }
}

///|
pub fn ClustalOmegaCommandline::set_output_format(self : ClustalOmegaCommandline, format : String) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: self.commandline.add_arg("--outfmt", Some(format)) }
}

///|
pub fn ClustalOmegaCommandline::set_iterations(self : ClustalOmegaCommandline, n : Int) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: self.commandline.add_arg("--iterations", Some(n.to_string())) }
}

///|
pub fn ClustalOmegaCommandline::set_full_matrix(self : ClustalOmegaCommandline) -> ClustalOmegaCommandline {
  ClustalOmegaCommandline::{ commandline: self.commandline.add_arg("--full", None) }
}

///|
pub fn ClustalOmegaCommandline::build(self : ClustalOmegaCommandline) -> String {
  self.commandline.build_command()
}

///|
pub struct MuscleCommandline {
  commandline : AbstractCommandline
}

///|
pub fn MuscleCommandline::new(executable : String) -> MuscleCommandline {
  MuscleCommandline::{ commandline: AbstractCommandline::new(executable) }
}

///|
pub fn MuscleCommandline::set_input(self : MuscleCommandline, infile : String) -> MuscleCommandline {
  MuscleCommandline::{ commandline: self.commandline.add_arg("-in", Some(infile)) }
}

///|
pub fn MuscleCommandline::set_output(self : MuscleCommandline, outfile : String) -> MuscleCommandline {
  MuscleCommandline::{ commandline: self.commandline.add_arg("-out", Some(outfile)) }
}

///|
pub fn MuscleCommandline::set_diags(self : MuscleCommandline) -> MuscleCommandline {
  MuscleCommandline::{ commandline: self.commandline.add_arg("-diags", None) }
}

///|
pub fn MuscleCommandline::set_max_iterations(self : MuscleCommandline, n : Int) -> MuscleCommandline {
  MuscleCommandline::{ commandline: self.commandline.add_arg("-maxiters", Some(n.to_string())) }
}

///|
pub fn MuscleCommandline::build(self : MuscleCommandline) -> String {
  self.commandline.build_command()
}

///|
pub struct MAFFTCommandline {
  commandline : AbstractCommandline
}

///|
pub fn MAFFTCommandline::new(executable : String) -> MAFFTCommandline {
  MAFFTCommandline::{ commandline: AbstractCommandline::new(executable) }
}

///|
pub fn MAFFTCommandline::set_input(self : MAFFTCommandline, infile : String) -> MAFFTCommandline {
  MAFFTCommandline::{ commandline: self.commandline.add_arg("--input", Some(infile)) }
}

///|
pub fn MAFFTCommandline::set_output(self : MAFFTCommandline, outfile : String) -> MAFFTCommandline {
  MAFFTCommandline::{ commandline: self.commandline.add_arg("--output", Some(outfile)) }
}

///|
pub fn MAFFTCommandline::set_auto(self : MAFFTCommandline) -> MAFFTCommandline {
  MAFFTCommandline::{ commandline: self.commandline.add_arg("--auto", None) }
}

///|
pub fn MAFFTCommandline::set_threads(self : MAFFTCommandline, n : Int) -> MAFFTCommandline {
  MAFFTCommandline::{ commandline: self.commandline.add_arg("--thread", Some(n.to_string())) }
}

///|
pub fn MAFFTCommandline::build(self : MAFFTCommandline) -> String {
  self.commandline.build_command()
}

///|
pub fn create_example_clustalw() -> ClustalwCommandline {
  let cmd = ClustalwCommandline::new("clustalw2")
  let cmd = cmd.set_input("sequences.fasta")
  let cmd = cmd.set_output("aligned.fasta")
  let cmd = cmd.set_output_format("FASTA")
  let cmd = cmd.set_matrix("BLOSUM62")
  cmd
}

///|
pub fn create_example_clustalomega() -> ClustalOmegaCommandline {
  let cmd = ClustalOmegaCommandline::new("clustalo")
  let cmd = cmd.set_input("sequences.fasta")
  let cmd = cmd.set_output("aligned.fasta")
  let cmd = cmd.set_output_format("fasta")
  let cmd = cmd.set_iterations(10)
  cmd
}

///|
pub fn create_example_muscle() -> MuscleCommandline {
  let cmd = MuscleCommandline::new("muscle")
  let cmd = cmd.set_input("sequences.fasta")
  let cmd = cmd.set_output("aligned.fasta")
  let cmd = cmd.set_diags()
  let cmd = cmd.set_max_iterations(16)
  cmd
}

///|
pub fn create_example_mafft() -> MAFFTCommandline {
  let cmd = MAFFTCommandline::new("mafft")
  let cmd = cmd.set_input("sequences.fasta")
  let cmd = cmd.set_output("aligned.fasta")
  let cmd = cmd.set_auto()
  let cmd = cmd.set_threads(4)
  cmd
}