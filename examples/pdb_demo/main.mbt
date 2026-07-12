///|
/// PDB structure parsing example.
///
/// Demonstrates parsing PDB files, accessing atoms/residues/chains,
/// and computing structural properties.
///
/// Run: moon run examples/pdb_demo/main.mbt
fn min_int(a : Int, b : Int) -> Int {
  if a < b {
    a
  } else {
    b
  }
}

///|
fn main {
  // Minimal PDB content (2 residues, 8 atoms)
  let pdb_content = "ATOM      1  N   ALA A   1      11.104  6.134  10.298  1.00 20.00           N\n" +
    "ATOM      2  CA  ALA A   1      11.639  6.168  10.426  1.00 20.00           C\n" +
    "ATOM      3  C   ALA A   1      12.221  6.201  10.542  1.00 20.00           C\n" +
    "ATOM      4  O   ALA A   1      12.803  6.234  10.658  1.00 20.00           O\n" +
    "ATOM      5  N   GLY A   2      13.385  6.267  10.774  1.00 20.00           N\n" +
    "ATOM      6  CA  GLY A   2      13.968  6.300  10.891  1.00 20.00           C\n" +
    "ATOM      7  C   GLY A   2      14.549  6.333  11.007  1.00 20.00           C\n" +
    "ATOM      8  O   GLY A   2      15.131  6.366  11.123  1.00 20.00           O\n" +
    "TER       9      ALA A   1\n" +
    "END\n"

  // === Parse PDB ===
  let structure = @bio.parse_pdb(pdb_content)
  println("=== PDB Structure ===")
  println("  ID: \{structure.get_id()}")
  println("  Models: \{structure.get_num_models()}")
  println("  Chains: \{structure.get_num_chains()}")
  println("  Residues: \{structure.get_num_residues()}")
  println("  Atoms: \{structure.get_num_atoms()}")

  // === Access chains ===
  let chains = structure.get_chains()
  println("\n=== Chains ===")
  for i = 0; i < chains.length(); i = i + 1 {
    let c = chains[i]
    println(
      "  Chain \{c.get_id()}: \{c.get_num_residues()} residues, \{c.get_atoms().length()} atoms",
    )
  }

  // === Access residues ===
  let residues = structure.get_residues()
  println("\n=== Residues ===")
  for i = 0; i < residues.length(); i = i + 1 {
    let r = residues[i]
    println(
      "  [\{i}] \{r.resname} \{r.chainid}\{r.resseq}: \{r.get_num_atoms()} atoms",
    )
    let is_het = r.is_het()
    println("       HETATM: \{is_het}")

    // Get alpha carbon
    let ca = r.get_ca()
    match ca {
      Some(atom) => {
        let coord = atom.get_coord()
        println("       CA: (\{coord.x}, \{coord.y}, \{coord.z})")
      }
      None => println("       CA: not found")
    }

    // Center of mass
    let com = r.center_of_mass()
    match com {
      Some(v) => println("       COM: (\{v.x}, \{v.y}, \{v.z})")
      None => println("       COM: N/A")
    }
  }

  // === Atom distance ===
  let atoms = structure.get_atoms()
  if atoms.length() >= 2 {
    let d = atoms[0].distance(atoms[1])
    println("\n=== Atom Distance ===")
    println("  \{atoms[0].name} -- \{atoms[1].name}: \{d}")
  }

  // === Write PDB ===
  let pdb_out = @bio.write_pdb(structure)
  println("\n=== Written PDB (first 5 lines) ===")
  let lines = pdb_out.split("\n").to_array()
  let limit = min_int(lines.length(), 5)
  for i = 0; i < limit; i = i + 1 {
    println("  \{lines[i].to_owned()}")
  }
}
