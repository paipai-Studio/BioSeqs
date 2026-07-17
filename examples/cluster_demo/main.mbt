///|
/// Sequence clustering example inspired by Biopython Bio.Cluster.
///
/// Run: moon run examples/cluster_demo/main.mbt
fn main {
  println("=== Bio.Cluster Demo ===")

  let sequences = [
    "ACGTACGTACGT", "ACGTACATACGT", "ACATACGTACGT", "TTTTTTTTTTTT", "TTTTTTATTTTT",
    "AAAAAAAAAAAA",
  ]
  let names = ["Seq1", "Seq2", "Seq3", "Seq4", "Seq5", "Seq6"]

  println("\n1. Sequences to cluster:")
  for i = 0; i < sequences.length(); i = i + 1 {
    println("   \{names[i]}: \{sequences[i]}")
  }

  println("\n2. Distance matrix:")
  let distances = @bio.distance_matrix(sequences)
  for i = 0; i < distances.length(); i = i + 1 {
    let row = distances[i].map(fn(d) { d.to_string() }).join("\t")
    println("   \{names[i]}:\t\{row}")
  }

  println("\n3. Single linkage clustering (threshold=0.2):")
  let single_result = @bio.single_linkage_clustering(distances, 0.2)
  println("   Number of clusters: \{@bio.total_clusters(single_result)}")

  println("\n4. Complete linkage clustering (threshold=0.2):")
  let complete_result = @bio.complete_linkage_clustering(distances, 0.2)
  println("   Number of clusters: \{@bio.total_clusters(complete_result)}")

  println("\n5. Average linkage clustering (threshold=0.2):")
  let avg_result = @bio.average_linkage_clustering(distances, 0.2)
  println("   Number of clusters: \{@bio.total_clusters(avg_result)}")

  println("\n6. Convert to Newick format:")
  let newick = @bio.hclust_to_newick(single_result, names)
  println("   Newick: \{newick}")

  println("\n7. Silhouette score (manual labels):")
  let labels = [0, 0, 0, 1, 1, 2]
  let score = @bio.seq_silhouette_score(sequences, labels)
  println("   Score: \{score}")

  println("\n8. Single linkage with higher threshold (0.5):")
  let high_result = @bio.single_linkage_clustering(distances, 0.5)
  println("   Number of clusters: \{@bio.total_clusters(high_result)}")
}
