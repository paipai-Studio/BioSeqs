///|
fn main {
  println("=== Bio.PDB.DSSP Demo ===")
  println("Inspired by Biopython Bio.PDB.DSSP module")
  println("")

  println("1. Predict secondary structure")
  let sequence = "MKLILVLLVSLSL"
  let ss = @bio.predict_secondary_structure(sequence)
  println("  Sequence: \{sequence}")
  println("  Secondary structure: \{ss}")
  println("")

  println("2. Calculate solvent accessibility")
  let structure = "HHHHHHHHHHHH"
  let acc = @bio.calculate_solvent_accessibility(sequence, structure)
  println("  Structure: \{structure}")
  println("  Accessibility: \{acc.show()}")
  println("")

  println("3. Create example DSSP data")
  let data = @bio.create_example_dssp_data()
  println("  Created \{data.records.length()} DSSP records")
  let mut i = 0
  while i < data.records.length() {
    let r = data.records[i]
    println("    Chain: \{r.chain_id}, Residue: \{r.residue_number}")
    println("      AA: \{r.amino_acid}, SS: \{r.secondary_structure}")
    println("      Solvent accessibility: \{r.solvent_accessibility}")
    println("      Phi: \{r.phi}, Psi: \{r.psi}")
    i = i + 1
  }
  println("")

  println("4. Calculate Ramachandran plot points")
  let plot = @bio.calculate_ramachandran_plot(data.records)
  println("  Generated \{plot.length()} Ramachandran points")
  println("")

  println("5. Analyze structure composition")
  let comp = @bio.analyze_structure_composition(data.records)
  println("  Helix (H): \{comp.get("H").unwrap_or(0)}")
  println("  Sheet (E): \{comp.get("E").unwrap_or(0)}")
  println("  Coil (-): \{comp.get("-").unwrap_or(0)}")
}