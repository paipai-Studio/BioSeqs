///|
/// biotite_application_reduce.mbt — Reduce (hydrogen placement) command
/// construction and naming helpers.

///|
/// Build a Reduce command line.
pub fn biotite_reduce_command(
  input_file : String,
  output_file : String,
  build~ : Bool = true,
) -> String {
  let flags = if build { "-build" } else { "-trim" }
  "reduce -quiet " + flags + " -pdb " + input_file + " > " + output_file
}

///|
/// Conventional hydrogen atom name for the n-th hydrogen on an atom.
pub fn biotite_reduce_hydrogen_name(
  atom_name : String,
  index : Int32,
) -> String {
  if index <= 0 {
    "H" + atom_name
  } else {
    "H" + index.to_string() + atom_name
  }
}

///|
/// Command to rebuild only polar hydrogens (common for docking).
pub fn biotite_reduce_polar_only_command(
  input_file : String,
  output_file : String,
) -> String {
  "reduce -quiet -build -noflip -pdb " + input_file + " > " + output_file
}
