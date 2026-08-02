///|
/// clusterExperiment RSEC workflow demo.
///
/// Demonstrates consensus and sequential clustering on synthetic
/// high-dimensional biological data using the RSEC pipeline.

fn main() -> Unit {
  println("=== clusterExperiment RSEC Demo ===\n")

  // Create synthetic data: 3 well-separated clusters in 5D space
  let data : Array[Array[Double]] = Array::new()

  // Cluster 0: 10 samples centered at (1, 1, 1, 1, 1)
  for i in 0..<10 {
    let pt = Array::new()
    for j in 0..<5 {
      pt.push(1.0 + (i.to_double() - 5.0) * 0.05)
    }
    data.push(pt)
  }

  // Cluster 1: 10 samples centered at (5, 5, 5, 5, 5)
  for i in 0..<10 {
    let pt = Array::new()
    for j in 0..<5 {
      pt.push(5.0 + (i.to_double() - 5.0) * 0.05)
    }
    data.push(pt)
  }

  // Cluster 2: 8 samples centered at (9, 1, 9, 1, 9)
  for i in 0..<8 {
    let pt = Array::new()
    pt.push(9.0 + (i.to_double() - 4.0) * 0.05)
    pt.push(1.0 + (i.to_double() - 4.0) * 0.05)
    pt.push(9.0 + (i.to_double() - 4.0) * 0.05)
    pt.push(1.0 + (i.to_double() - 4.0) * 0.05)
    pt.push(9.0 + (i.to_double() - 4.0) * 0.05)
    data.push(pt)
  }

  println("Input: " + data.length().to_string() + " samples, 5 features")
  println("Expected: 3 clusters (10 + 10 + 8 samples)\n")

  // Get default RSEC parameters
  let params = @bio.ce_default_params()
  println("RSEC Parameters:")
  println("  k values: [2, 3, 4]")
  println("  reduced dims: [0, 5]")
  println("  subsamples: [1.0, 0.8]")
  println("  merge cutoff: " + params.merge_cutoff.to_string())
  println("  sequential alpha: " + params.sequential_alpha.to_string())
  println("")

  // Step 1: clusterMany - ensemble of clusterings
  let ensemble = @bio.ce_cluster_many(data, params)
  println("clusterMany: generated " + ensemble.length().to_string() + " clusterings")

  // Step 2: Co-clustering matrix
  let cocluster = @bio.ce_cocluster_matrix(ensemble, data.length())
  println("Co-clustering matrix: " + cocluster.length().to_string() + "x" + cocluster.length().to_string())

  // Print a small sample of the co-clustering matrix
  println("\nCo-clustering matrix (first 5x5 block):")
  for i in 0..<5 {
    let mut row = ""
    for j in 0..<5 {
      row = row + cocluster[i][j].to_string() + "  "
    }
    println("  " + row)
  }

  // Step 3: Consensus clustering
  let consensus_labels = @bio.ce_make_consensus(cocluster, 3)
  println("\nConsensus clustering (k=3):")
  let consensus_sizes = @bio.ce_cluster_sizes(consensus_labels)
  let mut ci = 0
  while ci < consensus_sizes.length() {
    let (id, size) = consensus_sizes[ci]
    println("  cluster " + id.to_string() + ": " + size.to_string() + " samples")
    ci = ci + 1
  }

  // Step 4: Dendrogram
  let dendro = @bio.ce_make_dendrogram(data, consensus_labels)
  println("\nDendrogram built: " + (if dendro is Some(_) { "yes" } else { "no" }))

  // Step 5: Merge clusters by significance
  let (merged_labels, merge_info) = @bio.ce_merge_clusters(data, consensus_labels, dendro, params.merge_cutoff)
  println("Merge clusters: " + merge_info.length().to_string() + " merge steps evaluated")
  let mut merged_count = 0
  for info in merge_info {
    if info.merged {
      merged_count = merged_count + 1
    }
  }
  println("  " + merged_count.to_string() + " merges applied")

  // Step 6: Sequential clustering
  let seq_labels = @bio.ce_sequential_cluster(data, merged_labels, params.sequential_alpha, params.sequential_beta)
  let n_seq = @bio.ce_count_clusters(seq_labels)
  println("\nSequential clustering: " + n_seq.to_string() + " final clusters")

  // Full RSEC pipeline
  println("\n--- Full RSEC Pipeline ---")
  let result = @bio.ce_rsec(data, params)
  println(@bio.ce_summary(result))

  // Show per-sample assignments
  println("Per-sample cluster assignments:")
  let mut i = 0
  while i < result.labels.length() {
    let label_str = result.labels[i].to_string()
    let mut line = "  sample " + i.to_string() + " -> cluster " + label_str
    if i < 10 {
      line = line + " (cluster 0 group)"
    } else if i < 20 {
      line = line + " (cluster 1 group)"
    } else {
      line = line + " (cluster 2 group)"
    }
    println(line)
    i = i + 1
  }

  println("\n=== Demo Complete ===")
}
