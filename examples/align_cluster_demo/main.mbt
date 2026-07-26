///|
/// Progressive Multiple Sequence Alignment Demo
/// Demonstrates Bio.Align.AlignClusterer functionality
///
/// Run: moon run examples/align_cluster_demo/

fn main {
  // Create aligner with default parameters
  let aligner = @bio.align_cluster_new()
  println("=== Progressive MSA Demo ===")
  println("Aligner created with default scoring parameters")
  println("  (match=2.0, mismatch=-1.0, open_gap=-2.0, extend_gap=-0.5)")

  // Prepare sequences (simulated homologous sequences)
  let sequences = [
    "ACGTACGTACGTACGT",
    "ACGTACGTTCGTACGT",
    "ACGTACGTACGTTCGT",
    "TCGTACGTACGTACGT",
  ]
  let ids = ["seq1", "seq2", "seq3", "seq4"]

  println("")
  println("Input sequences:")
  for i = 0; i < sequences.length(); i = i + 1 {
    println("  \{ids[i]}: \{sequences[i]}")
  }

  // Build distance matrix
  println("")
  println("Building distance matrix...")
  let dist_matrix = @bio.build_distance_matrix(aligner, sequences, ids)
  let matrix = dist_matrix.get_matrix()
  println("Distance matrix (\{matrix.length()}x\{matrix[0].length()}):")
  for i = 0; i < matrix.length(); i = i + 1 {
    let mut row = ""
    for j = 0; j < matrix[i].length(); j = j + 1 {
      row = row + matrix[i][j].to_string() + "  "
    }
    println("  \{row}")
  }

  // Build guide tree
  println("")
  println("Building UPGMA guide tree...")
  let tree = @bio.upgma_tree(dist_matrix)
  let nodes = tree.get_nodes()
  println("Guide tree has \{nodes.length()} nodes")
  println("Sequence count (expected leaves): \{sequences.length()}")
  println("Internal nodes: \{nodes.length() - sequences.length()}")

  // Perform progressive alignment
  println("")
  println("Performing progressive alignment...")
  let msa = @bio.progressive_alignment(aligner, sequences, ids)
  println("Alignment complete!")
  println("Number of sequences: \{msa.num_records()}")
  println("Alignment length: \{msa.get_alignment_length()}")

  // Display aligned sequences using get_aligned_sequences
  println("")
  println("Aligned sequences (via get_aligned_sequences):")
  let aligned = @bio.get_aligned_sequences(aligner, sequences)
  for i = 0; i < aligned.length(); i = i + 1 {
    println("  \{ids[i]}: \{aligned[i]}")
  }

  // Display MSA representation
  println("")
  println("MSA representation:")
  println(msa.repr())

  println("=== Demo complete ===")
}