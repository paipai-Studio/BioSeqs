///|
/// Demo for rrvgo module: GO term semantic redundancy reduction.

fn main {
  println("=== rrvgo Demo: GO Term Redundancy Reduction ===")
  println("")

  // Example 1: Build a small GO graph using gosemsim demo graph
  println("--- Step 1: Build GO Graph ---")
  let graph = @src.gosemsim_example_graph()
  let ids = graph.term_ids()
  println("GO terms count: " + ids.length().to_string())
  println("")

  // Example 2: Compute semantic similarity matrix (Wang method)
  println("--- Step 2: Semantic Similarity Matrix ---")
  let sims = @src.rvgo_sim_matrix(graph, ids)
  println("Similarity matrix: " + sims.length().to_string() + " x " + sims[0].length().to_string())
  println("GO:0007165 vs GO:0023046 sim = " + sims[5][4].to_string())
  println("")

  // Example 3: Distance matrix
  println("--- Step 3: Distance Matrix ---")
  let dist = @src.rvgo_dist_matrix(sims)
  println("Dist range: [0.0, " + dist[5][4].to_string() + "]")
  println("")

  // Example 4: Hierarchical clustering
  println("--- Step 4: Hierarchical Clustering ---")
  let merges = @src.rvgo_hclust(dist)
  println("Merges: " + merges.length().to_string())
  println("")

  // Example 5: Cut into clusters
  println("--- Step 5: Cut Tree into Clusters ---")
  let clusters = @src.rvgo_cutree(merges, ids.length(), 0.7)
  println("Number of clusters: " + clusters.length().to_string())
  let cluster_sizes : Map[Int, Int] = Map([], capacity=clusters.length())
  for i in 0..<ids.length() {
    let c = clusters[i]
    let old = if cluster_sizes.contains(c) { cluster_sizes.get(c).unwrap_or(0) } else { 0 }
    cluster_sizes.set(c, old + 1)
  }
  for key in cluster_sizes.keys() {
    println("  Cluster " + key.to_string() + ": " + cluster_sizes.get(key).unwrap_or(0).to_string() + " terms")
  }
  println("")

  // Example 6: Get representative terms
  println("--- Step 6: Get Representative Terms ---")
  let term_ids = ["GO:0008150", "GO:0050896", "GO:0009987", "GO:0006810",
                  "GO:0023046", "GO:0007165", "GO:0009058", "GO:0010467"]
  let scores_arr = [0.01, 0.05, 0.02, 0.1, 0.03, 0.001, 0.08, 0.04]
  let reduced = @src.rvgo_reduce(graph, term_ids, scores_arr, threshold=0.7)
  println("Representative terms: " + reduced.representatives.length().to_string())
  for rep in reduced.representatives {
    println("  " + rep)
  }
  println("")

  // Example 7: PCoA 2D visualization coordinates
  println("--- Step 7: PCoA 2D Coordinates ---")
  let coords = @src.rvgo_mds_2d(dist)
  println("Coordinates: " + coords.length().to_string() + " points")
  for i in 0..<coords.length() {
    let (x, y) = coords[i]
    println("  " + ids[i] + " -> (" + x.to_string() + ", " + y.to_string() + ")")
  }
  println("")

  println("=== rrvgo Demo Complete ===")
}
