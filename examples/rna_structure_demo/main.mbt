///|
/// RNA Secondary Structure Prediction Demo
/// Demonstrates Bio.SeqUtils - Nussinov algorithm functionality
///
/// Run: moon run examples/rna_structure_demo/

fn main {
  println("=== RNA Secondary Structure Prediction Demo ===")

  // Create predictor with default parameters
  let predictor = @bio.rna_predictor_new()
  println("Predictor created with default parameters")
  println("  (loop_min_length=4, stacking_energy=-2.0, loop_penalty=0.5, mismatch_penalty=1.0)")

  // Test sequences
  let sequences = [
    "GCGCAAAAGCGC", // Simple hairpin
    "GGGAAACCCUUUGGG", // More complex
    "AUGCUAGCUAGCUAGC", // Random RNA
  ]

  for seq in sequences {
    println("")
    println("Sequence: \{seq}")
    println("Length: \{seq.length()} bases")

    // Check complementary
    println("Complementarity check:")
    let pairs_to_check : Array[(Char, Char)] = [('A', 'U'), ('G', 'C'), ('G', 'U'), ('A', 'C')]
    for pair in pairs_to_check {
      let (b1, b2) = pair
      let is_comp = @bio.is_complementary(b1.to_int().to_uint16(), b2.to_int().to_uint16())
      let comp_str = if is_comp { "true" } else { "false" }
      println("  \{b1.to_string()}-\{b2.to_string()}: \{comp_str}")
    }

    // Predict structure
    println("")
    println("Predicting secondary structure...")
    let result = @bio.rna_predict_secondary_structure(predictor, seq)

    // Display results
    let dot_bracket = @bio.structure_to_dot_bracket(result)
    println("Dot-bracket notation: \{dot_bracket}")
    println("MFE: \{@bio.get_mfe(result).to_string()} kcal/mol")

    // Display base pairs
    let base_pairs = @bio.get_base_pairs(result)
    println("Base pairs found: \{base_pairs.length().to_string()}")
    for bp in base_pairs {
      let (i, j) = bp
      println("  \{(i + 1).to_string()} - \{(j + 1).to_string()} (\{seq.unsafe_get(i).unsafe_to_char().to_string()}..\{seq.unsafe_get(j).unsafe_to_char().to_string()})")
    }

    // Visual representation
    println("  Sequence:  \{seq}")
    println("  Structure: \{dot_bracket}")
  }

  // Test custom parameters
  println("")
  println("=== Custom Predictor ===")
  let custom = @bio.rna_predictor_with_params(3, -3.0, 0.3, 0.5)
  let result2 = @bio.rna_predict_secondary_structure(custom, "GCGCAAAAGCGC")
  println("Custom predictor MFE: \{@bio.get_mfe(result2).to_string()} kcal/mol")
  println("Structure: \{@bio.structure_to_dot_bracket(result2)}")

  println("=== Demo complete ===")
}