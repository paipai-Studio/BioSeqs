///|
fn main {
  println("=== Bio.PDB.Polypeptide Demo ===")
  println("Inspired by Biopython Bio.PDB.Polypeptide module")
  println("")

  println("1. Calculate amino acid composition")
  let sequence = "MKLILVLLVSLSL"
  let comp = @bio.calculate_composition(sequence)
  println("  Sequence length: \{comp.total_residues}")
  println("  Alanine count: \{comp.amino_acids.get("A").unwrap_or(0)}")
  println("  Leucine count: \{comp.amino_acids.get("L").unwrap_or(0)}")
  println("")

  println("2. Calculate percent composition")
  let percent = @bio.calculate_percent_composition(comp)
  println("  Alanine: \{percent.get("A").unwrap_or(0.0)}%")
  println("  Leucine: \{percent.get("L").unwrap_or(0.0)}%")
  println("")

  println("3. Calculate hydrophobicity profile (Kyte-Doolittle)")
  let profile = @bio.calculate_hydrophobicity(sequence, "kyte-doolittle")
  println("  Generated \{profile.values.length()} hydrophobicity values")
  println("")

  println("4. Calculate hydrophobicity with sliding window")
  let window_profile = @bio.calculate_hydrophobicity_window(sequence, 5, "kyte-doolittle")
  println("  Window size: 5")
  println("  Generated \{window_profile.values.length()} windowed values")
  println("")

  println("5. Predict transmembrane regions")
  let tm_sequence = "MVLILVLLVSLSLILVLLV"
  let regions = @bio.predict_transmembrane_regions(tm_sequence, 7, 2.0)
  println("  Sequence: \{tm_sequence}")
  println("  Found \{regions.length()} transmembrane regions")
  let mut i = 0
  while i < regions.length() {
    let r = regions[i]
    println("    Region: \{r.start}-\{r.end}, Score: \{r.score}")
    i = i + 1
  }
  println("")

  println("6. Calculate isoelectric point")
  let pi = @bio.calculate_isoelectric_point(sequence)
  println("  Isoelectric point (pI): \{pi}")
  println("")

  println("7. Calculate aromaticity")
  let aromatic_seq = "FFYYWW"
  let aromaticity = @bio.calculate_aromaticity(aromatic_seq)
  println("  Sequence: \{aromatic_seq}")
  println("  Aromaticity: \{aromaticity}")
  println("")

  println("8. Create example polypeptide data")
  let (seq, comp_data, prof_data) = @bio.create_example_polypeptide_data()
  println("  Sequence length: \{seq.length()}")
  println("  Total residues: \{comp_data.total_residues}")
  println("  Hydrophobicity values: \{prof_data.values.length()}")
}