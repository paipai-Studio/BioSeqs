///|
/// PDB Dice Demo - PDB structure extraction and selection tools.
///
/// Run: moon run examples/pdb_dice_demo/main.mbt
fn main {
  println("=== PDB Dice Demo ===\n")

  // Create a simple PDB structure for demonstration
  let pdb_content =
    "ATOM      1  N   ALA A   1      11.104  6.134  10.298  1.00 20.00           N\n" +
    "ATOM      2  CA  ALA A   1      11.639  6.168  10.426  1.00 20.00           C\n" +
    "ATOM      3  C   ALA A   1      12.221  6.201  10.542  1.00 20.00           C\n" +
    "ATOM      4  O   ALA A   1      12.803  6.234  10.658  1.00 20.00           O\n" +
    "ATOM      5  CB  ALA A   1      12.000  7.000  11.000  1.00 20.00           C\n" +
    "ATOM      6  N   GLY A   2      13.385  6.267  10.774  1.00 25.00           N\n" +
    "ATOM      7  CA  GLY A   2      13.968  6.300  10.891  1.00 25.00           C\n" +
    "ATOM      8  C   GLY A   2      14.549  6.333  11.007  1.00 25.00           C\n" +
    "ATOM      9  O   GLY A   2      15.131  6.366  11.123  1.00 25.00           O\n" +
    "ATOM     10  N   VAL A   3      15.500  6.400  11.200  1.00 30.00           N\n" +
    "ATOM     11  CA  VAL A   3      16.000  6.450  11.300  1.00 30.00           C\n" +
    "ATOM     12  C   VAL A   3      16.500  6.500  11.400  1.00 30.00           C\n" +
    "ATOM     13  O   VAL A   3      17.000  6.550  11.500  1.00 30.00           O\n" +
    "ATOM     14  CB  VAL A   3      16.200  7.200  11.800  1.00 30.00           C\n" +
    "ATOM     15  CG1 VAL A   3      15.800  8.000  11.200  1.00 30.00           C\n" +
    "ATOM     16  CG2 VAL A   3      17.000  7.800  12.500  1.00 30.00           C\n" +
    "ATOM     17  N   LEU A   4      17.500  6.600  11.600  1.00 35.00           N\n" +
    "ATOM     18  CA  LEU A   4      18.000  6.650  11.700  1.00 35.00           C\n" +
    "ATOM     19  C   LEU A   4      18.500  6.700  11.800  1.00 35.00           C\n" +
    "ATOM     20  O   LEU A   4      19.000  6.750  11.900  1.00 35.00           O\n" +
    "TER      21      LEU A   4\n" +
    "END\n"

  let structure = @bio.parse_pdb(pdb_content)
  println("1. Original structure:")
  println("   Models: \{structure.get_num_models()}")
  println("   Chains: \{structure.get_num_chains()}")
  println("   Residues: \{structure.get_num_residues()}")
  println("   Atoms: \{structure.get_num_atoms()}")

  // 2. Extract chain A
  println("\n2. Extract chain A:")
  let chain_a = @bio.extract_chain(structure, 'A')
  println("   Chains extracted: \{chain_a.get_num_chains()}")
  println("   Residues: \{chain_a.get_num_residues()}")
  println("   Atoms: \{chain_a.get_num_atoms()}")

  // 3. Extract residue range
  println("\n3. Extract residues 2-3 from chain A:")
  let res_range = @bio.extract_residue_range(structure, 'A', 2, 3)
  println("   Residues in range: \{res_range.get_num_residues()}")
  println("   Atoms: \{res_range.get_num_atoms()}")

  // 4. Extract specific residues
  println("\n4. Extract specific residues (1, 4) from chain A:")
  let rid1 = @bio.ResidueId::new('A', 1)
  let rid2 = @bio.ResidueId::new('A', 4)
  let residue_ids = [rid1, rid2]
  let specific_res = @bio.extract_residues(structure, residue_ids)
  println("   Residues extracted: \{specific_res.get_num_residues()}")
  println("   Atoms: \{specific_res.get_num_atoms()}")

  // 5. Extract CA atoms
  println("\n5. Extract CA atoms only:")
  let ca_structure = @bio.extract_ca_atoms(structure)
  println("   Structure with CA atoms:")
  println("   Residues: \{ca_structure.get_num_residues()}")
  println("   Atoms: \{ca_structure.get_num_atoms()}")

  // 6. Filter by B-factor
  println("\n6. Filter atoms with B-factor between 25 and 40:")
  let filtered = @bio.filter_by_bfactor(structure, 25.0, 40.0)
  println("   Atoms after filtering: \{filtered.get_num_atoms()}")

  // 7. Select atoms near a point
  println("\n7. Select atoms near point (15.0, 6.5, 11.3) within 2.0Å:")
  let point = @bio.Vector3::new(15.0, 6.5, 11.3)
  let near_atoms = @bio.select_atoms_near_point(structure, point, 2.0)
  println("   Atoms near point: \{near_atoms.length()}")
  if near_atoms.length() > 0 {
    println("   First atom: \{near_atoms[0].name} (element: \{near_atoms[0].element})")
  }

  // 8. Structure statistics
  println("\n8. Structure statistics:")
  let stats = @bio.structure_stats(structure)
  println("   Structure ID: \{stats.structure_id}")
  println("   Atoms: \{stats.num_atoms}")
  println("   Residues: \{stats.num_residues}")
  println("   Chains: \{stats.num_chains}")
  println("   Models: \{stats.num_models}")

  // 9. Chain sequence
  println("\n9. Chain A sequence:")
  let seq = @bio.chain_sequence(structure, 'A')
  println("   Sequence: \{seq}")

  println("\n=== Demo Complete ===")
}
