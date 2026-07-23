///|
/// Bio.PDB.PSEA Demo - Protein secondary structure prediction from CA coordinates.
/// Inspired by Biopython's Bio.PDB.PSEA module.

fn main {
  println("=== Bio.PDB.PSEA Demo ===")

  // 1. Helix example
  println("")
  println("1. Predict secondary structure for an alpha-helix-like chain:")
  let helix_atoms = @bio.psea_create_example_helix()
  println("   CA atoms: " + helix_atoms.length().to_string())
  let helix_result = @bio.psea_run(helix_atoms)
  println("   SS string: " + helix_result.psea_ss_string())
  println("   " + helix_result.psea_summary())

  // 2. Strand example
  println("")
  println("2. Predict secondary structure for a beta-strand-like chain:")
  let strand_atoms = @bio.psea_create_example_strand()
  println("   CA atoms: " + strand_atoms.length().to_string())
  let strand_result = @bio.psea_run(strand_atoms)
  println("   SS string: " + strand_result.psea_ss_string())
  println("   " + strand_result.psea_summary())

  // 3. Mixed example
  println("")
  println("3. Predict secondary structure for a mixed chain:")
  let mixed_atoms = @bio.psea_create_example_mixed()
  println("   CA atoms: " + mixed_atoms.length().to_string())
  let mixed_result = @bio.psea_run(mixed_atoms)
  println("   SS string: " + mixed_result.psea_ss_string())
  println("   " + mixed_result.psea_summary())

  // 4. Geometry analysis
  println("")
  println("4. Geometric analysis of helix chain:")
  let geom = @bio.psea_compute_geometry(helix_atoms)
  println("   CA-CA distances: " + geom.ca_distances.length().to_string())
  println("   Bond angles: " + geom.bond_angles.length().to_string())
  println("   Dihedral angles: " + geom.dihedral_angles.length().to_string())
  if geom.dihedral_angles.length() > 0 {
    println("   First dihedral angle: " + geom.dihedral_angles[0].to_string() + " degrees")
  }

  // 5. Per-residue details
  println("")
  println("5. Per-residue details (first 6 residues of helix chain):")
  let mut i = 0
  while i < 6 && i < helix_result.psea_n_residues() {
    let res_num = helix_result.psea_get_res_num(i)
    let res_name = helix_result.psea_get_res_name(i)
    let ss = helix_result.psea_get_ss(i)
    let conf = helix_result.psea_get_confidence(i)
    println("   Residue " + res_num.to_string() + " (" + res_name + "): SS=" + ss + ", confidence=" + conf.to_string())
    i = i + 1
  }

  // 6. Three-to-eight state conversion
  println("")
  println("6. Convert 3-state to 8-state (DSSP-like) assignment:")
  let ss3 = ["C", "C", "H", "H", "H", "H", "H", "C", "E", "E", "E", "C"]
  let mut ss3_str = ""
  let mut idx = 0
  while idx < ss3.length() {
    ss3_str = ss3_str + ss3[idx]
    idx = idx + 1
  }
  let ss8 = @bio.psea_three_to_eight_state(ss3)
  let mut ss8_str = ""
  let mut j = 0
  while j < ss8.length() {
    ss8_str = ss8_str + ss8[j]
    j = j + 1
  }
  println("   3-state: " + ss3_str)
  println("   8-state: " + ss8_str)

  println("")
  println("=== Bio.PDB.PSEA Demo Completed ===")
}
