///|
/// Tree Construction example: calculate distance matrices from DNA/protein
/// sequence alignments with Biopython-compatible models and build
/// phylogenetic trees.
///
/// Run: moon run examples/tree_construction_demo
fn main raise {
  println("=== Phylogenetic Tree Construction Demo ===")
  println("")

  // === Example 1: MoonBit extension models (JC69 / Kimura) ===
  println("--- 1. Jukes-Cantor / Kimura Corrected Distances (extension) ---")

  let jc_calc = @bio.DistanceCalculator::new(model="jukes_cantor")

  let seq1 = "ACGTACGTACGTACGTACGTACGTACGT"
  let seq2 = "ACGTACGTACGTACGTACGTACGTGCGA"
  let seq3 = "ACGTACGTACGTAAAGCTTGCATGCCTA"
  let seq4 = "TTTACGTACGTACGTACGTACGTACGTA"

  println("  Seq1: \{seq1}")
  println("  Seq2: \{seq2}")
  println("  JC69 distance (1 vs 2): \{fmt6(jc_calc.pairwise(seq1, seq2))}")
  println("  JC69 distance (1 vs 3): \{fmt6(jc_calc.pairwise(seq1, seq3))}")
  let kimura_calc = @bio.DistanceCalculator::new(model="kimura")
  println("  Kimura distance (1 vs 2): \{fmt6(kimura_calc.pairwise(seq1, seq2))}")
  println("  Kimura distance (1 vs 4): \{fmt6(kimura_calc.pairwise(seq1, seq4))}")
  println("")

  // === Example 2: identity model on the Biopython msa.phy alignment ===
  println("--- 2. Identity Model (Biopython TreeConstruction/msa.phy) ---")

  let msa = [
    ("Alpha", "AACGTGGCCACAT"),
    ("Beta", "AAGGTCGCCACAC"),
    ("Gamma", "CAGTTCGCCACAA"),
    ("Delta", "GAGATTTCCGCCT"),
    ("Epsilon", "GAGATCTCCGCCC"),
  ]
  let id_calc = @bio.DistanceCalculator::new(model="identity")
  print_distance_matrix(id_calc.get_distance(msa), ["Alpha", "Beta", "Gamma", "Delta", "Epsilon"])
  println("")

  // === Example 3: BLOSUM62 protein scoring-matrix model ===
  println("--- 3. BLOSUM62 Scoring-Matrix Model ---")

  let blosum_calc = @bio.DistanceCalculator::new(model="blosum62")
  print_distance_matrix(
    blosum_calc.get_distance(msa),
    ["Alpha", "Beta", "Gamma", "Delta", "Epsilon"],
  )
  println("")

  // === Example 4: DNA matrix models (blastn maps to the NUC.4.4 matrix) ===
  println("--- 4. DNA Matrix Models ---")

  let dna1 = "ACGTACGTAC"
  let dna2 = "ACGTACGAAC"
  let blastn_calc = @bio.DistanceCalculator::new(model="blastn")
  let hoxd70_calc = @bio.DistanceCalculator::new(model="hoxd70")
  println("  Seq1:      \{dna1}")
  println("  Seq2:      \{dna2}")
  println("  blastn distance:  \{fmt6(blastn_calc.pairwise(dna1, dna2))}")
  println("  hoxd70 distance:  \{fmt6(hoxd70_calc.pairwise(dna1, dna2))}")
  println("")

  // === Example 5: gaps, stop symbols and bad letters ===
  println("--- 5. Gap/Stop Handling and Error Reporting ---")

  // Columns with '-' or '*' are skipped by matrix models
  let gap_dist = blosum_calc.pairwise("AC-DE", "ACWDE")
  println("  BLOSUM62 distance 'AC-DE' vs 'ACWDE' (gap column skipped): \{fmt6(gap_dist)}")
  // The identity model skips nothing: gaps count as characters
  let id_gap_dist = id_calc.pairwise("A--A", "ACCA")
  println("  identity distance 'A--A' vs 'ACCA': \{fmt6(id_gap_dist)}")
  // Unknown models raise an error listing the available models
  try {
    ignore(@bio.DistanceCalculator::new(model="unknown_model"))
  } catch {
    SeqError(msg) => {
      let short =
        if msg.length() > 70 {
          msg.substring(start=0, end=70) + "..."
        } else {
          msg
        }
      println("  Unknown model raised: \{short}")
    }
  }
  println("")

  // === Example 6: available models ===
  println("--- 6. Available Distance Models ---")

  let models = @bio.distance_models()
  let line = StringBuilder()
  line.write_string("  (\{models.length()} models) ")
  let mut mi = 0
  while mi < models.length() {
    line.write_string(models[mi] + " ")
    mi = mi + 1
  }
  println(line.to_string())
  println("")

  // === Example 7: build trees from the BLOSUM62 distance matrix ===
  println("--- 7. Tree Construction from BLOSUM62 Distances ---")

  let matrix = blosum_calc.get_distance(msa)
  let names = ["Alpha", "Beta", "Gamma", "Delta", "Epsilon"]

  let upgma_cons = @bio.TreeConstructor::new(tree_method="upgma")
  let upgma_tree = upgma_cons.build_tree(matrix, names)
  println("  UPGMA root clades: \{upgma_tree.clades.length()}")

  let nj_cons = @bio.TreeConstructor::new(tree_method="nj")
  let nj_tree = nj_cons.build_tree(matrix, names)
  println("  NJ root clades: \{nj_tree.clades.length()}")

  let wpgma_cons = @bio.TreeConstructor::new(tree_method="wpgma")
  let wpgma_tree = wpgma_cons.build_tree(matrix, names)
  println("  WPGMA root clades: \{wpgma_tree.clades.length()}")
  println("")

  println("=== Demo Complete ===")
}

///|
/// Print a symmetric distance matrix with row/column labels.
fn print_distance_matrix(
  matrix : Array[Array[Double]],
  names : Array[String],
) -> Unit {
  let header = StringBuilder()
  header.write_string("        ")
  for name in names {
    header.write_string(name.pad_end(10, ' '))
  }
  println(header.to_string())
  for i in 0..<matrix.length() {
    let line = StringBuilder()
    line.write_string(names[i].pad_end(8, ' '))
    for j in 0..<(i + 1) {
      line.write_string(fmt6(matrix[i][j]).pad_end(10, ' '))
    }
    println(line.to_string())
  }
}

///|
/// Format a distance with six digits after the decimal point, matching the
/// Biopython distance-matrix presentation.
fn fmt6(v : Double) -> String {
  let neg = v < 0.0
  let av = if neg {
    -v
  } else {
    v
  }
  let scaled = (av * 1.0e6).round()
  let whole = (scaled / 1.0e6).to_int()
  let mut frac = (scaled - whole.to_double() * 1.0e6).to_int()
  let digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
  let buf = StringBuilder()
  if neg {
    buf.write_char('-')
  }
  buf.write_string(whole.to_string())
  buf.write_char('.')
  let mut factor = 100000
  while factor >= 1 {
    let d = frac / factor
    buf.write_string(digits[d])
    frac = frac - d * factor
    factor = factor / 10
  }
  buf.to_string()
}
