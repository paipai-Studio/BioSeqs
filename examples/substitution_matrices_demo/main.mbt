///|
/// substitution_matrices demo - demonstrates the modern substitution matrix
/// infrastructure inspired by Biopython Bio.Align.substitution_matrices.
///
/// Run with: moon run examples/substitution_matrices_demo
fn main {
  println("=== Substitution Matrices Demo ===\n")

  // 1. Built-in matrices
  println("1. Built-in BLOSUM62 Matrix:")
  let blosum62 = @bio.subs_blosum62_matrix()
  println("   Name: " + blosum62.name())
  println("   Alphabet size: " + blosum62.n_letters().to_string())
  println("   A-A score: " + blosum62.get_score("A", "A").to_int().to_string())
  println("   R-R score: " + blosum62.get_score("R", "R").to_int().to_string())
  println("   A-R score: " + blosum62.get_score("A", "R").to_int().to_string())
  println("   W-W score: " + blosum62.get_score("W", "W").to_int().to_string())
  println("   C-C score: " + blosum62.get_score("C", "C").to_int().to_string())

  // 2. Matrix registry
  println("\n2. Matrix Registry:")
  @bio.subs_initialize_registry()
  let names = @bio.list_matrices()
  println("   Registered matrices: " + names.length().to_string())
  let mut i = 0
  while i < names.length() {
    println("   - " + names[i])
    i = i + 1
  }

  // Load a matrix by name
  let loaded = @bio.load_matrix("BLOSUM80")
  match loaded {
    Some(m) =>
      println(
        "   Loaded BLOSUM80: A-A = " +
        m.get_score("A", "A").to_int().to_string(),
      )
    None => println("   BLOSUM80 not found")
  }

  // 3. Nucleotide matrix
  println("\n3. Nucleotide Matrix NUC4.4:")
  let nuc = @bio.subs_nuc44_matrix()
  println("   Name: " + nuc.name())
  println("   A-A score: " + nuc.get_score("A", "A").to_int().to_string())
  println("   A-C score: " + nuc.get_score("A", "C").to_int().to_string())
  println("   G-T score: " + nuc.get_score("G", "T").to_int().to_string())

  // 4. Sub-matrix selection
  println("\n4. Sub-matrix Selection (A, C, D from BLOSUM62):")
  let sub = blosum62.select(["A", "C", "D"])
  println("   Selected alphabet size: " + sub.n_letters().to_string())
  println("   A-A: " + sub.get_score("A", "A").to_int().to_string())
  println("   C-C: " + sub.get_score("C", "C").to_int().to_string())
  println("   D-D: " + sub.get_score("D", "D").to_int().to_string())
  println("   A-C: " + sub.get_score("A", "C").to_int().to_string())

  // 5. Frequency matrix from alignment
  println("\n5. Frequency Matrix from Alignment:")
  let alignment = [
    "ACGTACGT",
    "ACGTACGT",
    "AGGTACGA",
    "ACGTTCGT",
  ]
  let alpha = @bio.subs_nucleotide_alphabet()
  let freq = @bio.calculate_frequency_matrix(alignment, alpha)
  println("   Alignment (4 sequences, 8 positions):")
  let mut si = 0
  while si < alignment.length() {
    println("     " + alignment[si])
    si = si + 1
  }
  println("   A-A frequency: " + freq.get(0, 0).to_string())
  println("   C-C frequency: " + freq.get(1, 1).to_string())
  println("   G-G frequency: " + freq.get(2, 2).to_string())
  println("   T-T frequency: " + freq.get(3, 3).to_string())
  println("   C-G frequency: " + freq.get(1, 2).to_string())

  // 6. Calculate substitution matrix from alignment
  println("\n6. Substitution Matrix from Alignment (log-odds):")
  let calc_matrix = @bio.calculate_substitution_matrix(alignment, alpha, scale=2.0)
  println("   A-A score: " + calc_matrix.get_score("A", "A").to_string())
  println("   C-C score: " + calc_matrix.get_score("C", "C").to_string())
  println("   G-G score: " + calc_matrix.get_score("G", "G").to_string())
  println("   T-T score: " + calc_matrix.get_score("T", "T").to_string())
  println("   A-C score: " + calc_matrix.get_score("A", "C").to_string())

  // 7. Shannon entropy
  println("\n7. Shannon Entropy of Frequency Matrix:")
  let entropy = @bio.subs_shannon_entropy(freq)
  println("   Entropy: " + entropy.to_string() + " bits")

  // 8. Relative entropy (KL divergence)
  println("\n8. Relative Entropy (KL Divergence):")
  let kl = @bio.subs_relative_entropy(freq)
  println("   KL divergence: " + kl.to_string() + " bits/substitution")

  // 9. Matrix correlation
  println("\n9. Matrix Correlation:")
  let m45 = @bio.subs_blosum45_matrix()
  let m62 = @bio.subs_blosum62_matrix()
  let corr = @bio.matrix_correlation(m45, m62)
  println(
    "   BLOSUM45 vs BLOSUM62 correlation: " + corr.to_string(),
  )

  // 10. NCBI matrix parsing
  println("\n10. NCBI Matrix File Parsing:")
  let ncbi_content =
    #|# Example NCBI matrix
    #|   A  C  G  T
    #|A  5 -4 -4 -4
    #|C -4  5 -4 -4
    #|G -4 -4  5 -4
    #|T -4 -4 -4  5
  let parsed = @bio.parse_ncbi_matrix(ncbi_content)
  match parsed {
    Some(m) => {
      println("   Parsed matrix successfully")
      println("   A-A: " + m.get_score("A", "A").to_int().to_string())
      println("   C-C: " + m.get_score("C", "C").to_int().to_string())
      println("   A-C: " + m.get_score("A", "C").to_int().to_string())
    }
    None => println("   Failed to parse matrix")
  }

  // 11. Matrix table display
  println("\n11. Matrix Table Display (BLOSUM62 subset):")
  let small_sub = blosum62.select(["A", "C", "D", "E", "F", "G", "H", "I"])
  println(small_sub.to_table_string())

  println("=== Demo Complete ===")
}
