///|
/// mmCIF format writer example.
///
/// Demonstrates writing a protein Structure to mmCIF format,
/// displaying the output, and showing structure statistics.
///
/// Run: moon run examples/mmcifio_demo/main.mbt
fn main {
  println("=== Bio.PDB.mmcifio Demo ===")
  println("Writes protein structures to mmCIF format.")
  println("")

  // === Create a sample structure ===
  let structure = @bio.mmcif_sample_structure()
  println("1. Sample structure created")
  println("   ID: \{structure.get_id()}")
  println("   Models: \{structure.get_num_models()}")
  println("   Chains: \{structure.get_num_chains()}")
  println("   Residues: \{structure.get_num_residues()}")
  println("   Atoms: \{structure.get_num_atoms()}")
  println("")

  // === Write header only ===
  println("2. mmCIF header (entry + entity categories)")
  println("-----")
  let header = @bio.write_mmcif_header(structure)
  println(header)
  println("-----")
  println("")

  // === Write atom_site category only ===
  println("3. mmCIF _atom_site loop (first 6 lines shown)")
  println("-----")
  let atom_site = @bio.write_mmcif_atom_site(structure)
  let lines = atom_site.split("\n").to_array()
  let limit = if lines.length() < 6 { lines.length() } else { 6 }
  for i = 0; i < limit; i = i + 1 {
    println(lines[i].to_owned())
  }
  println("... (\{lines.length()} lines total in atom_site)")
  println("-----")
  println("")

  // === Full mmCIF output ===
  println("4. Full mmCIF output")
  println("-----")
  let output = @bio.write_mmcif(structure)
  let out_lines = output.split("\n").to_array()
  for i = 0; i < out_lines.length(); i = i + 1 {
    println(out_lines[i].to_owned())
  }
  println("-----")
  println("")

  // === Statistics ===
  println("5. Output statistics")
  println("   Total output length: \{output.length()} characters")
  println("   Total lines: \{out_lines.length()}")
  let mut atom_lines = 0
  for line in out_lines {
    if line.contains("ATOM") {
      atom_lines = atom_lines + 1
    }
  }
  println("   ATOM data rows: \{atom_lines}")
  println("")

  // === Round-trip verification ===
  println("6. Round-trip: parse written output back")
  let block = @bio.parse_mmcif(output)
  println("   Parsed data block name: \{block.name}")
  let atoms = block.to_atoms()
  println("   Parsed back \{atoms.length()} atoms")
  if atoms.length() > 0 {
    let first = atoms[0]
    let coord = first.get_coord()
    println("   First atom: \{first.name} of \{first.resname} at (\{coord.x}, \{coord.y}, \{coord.z})")
  }
  println("")

  // === Single-atom formatting demo ===
  println("7. Single atom line formatting")
  let atom = @bio.Atom::new(
    name="CA",
    coord=@bio.Vector3::new(12.345, -6.789, 42.0),
    resname="GLY",
    chainid='B',
    resseq=7,
    occupancy=0.5,
    bfactor=35.75,
    element="C",
  )
  let line = @bio.mmcif_format_atom(atom, 1)
  println("   \{line}")
  println("")

  // === Value escaping demo ===
  println("8. mmCIF value escaping")
  let vals = ["ALA", "", "hello world", "it's a test"]
  for v in vals {
    let escaped = @bio.mmcif_escape_value(v)
    println("   '\{v}' -> \{escaped}")
  }
  println("")

  println("=== Done ===")
}
