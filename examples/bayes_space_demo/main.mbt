///|
/// BayesSpace module example - Spatial transcriptomics clustering.
///
/// Run: moon run examples/bayes_space_demo
fn main {
  println("=== Bioconductor BayesSpace Module Demo ===")
  println("Spatial Transcriptomics Clustering with Spatial Prior")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Spatial Data ===")
  let (expr, spots) = @src.bayes_space_sample_data()
  println("Spots: " + spots.length().to_string())
  println("Features per spot: " + expr[0].length().to_string())
  println("Grid: 4x4")
  println("")

  // ---- Section 2: Build neighbor list ----
  println("=== 2. Build Neighbor List (hex grid) ===")
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let mut i = 0
  while i < spots.length() {
    println("  " + spots[i].spot_id + " (row=" + spots[i].row.to_string() + ", col=" + spots[i].col.to_string() + ") -> " + neighbors[i].length().to_string() + " neighbors")
    i = i + 1
  }
  println("")

  // ---- Section 3: Run BayesSpace clustering ----
  println("=== 3. Run BayesSpace (q=2 clusters) ===")
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=50, seed=42)
  println("Iterations: " + result.n_iterations.to_string())
  println("Converged: " + result.converged.to_string())
  println("Log-likelihood: " + result.log_likelihood.to_string())
  println("Number of clusters: " + result.q.to_string())
  println("")

  // ---- Section 4: Cluster assignments ----
  println("=== 4. Cluster Assignments ===")
  let mut si = 0
  while si < result.spot_ids.length() {
    println("  " + result.spot_ids[si] + " -> cluster " + result.clusters[si].to_string())
    si = si + 1
  }
  println("")

  // ---- Section 5: Cluster counts ----
  println("=== 5. Cluster Sizes ===")
  let counts = @src.bayes_space_cluster_counts(result)
  let mut ci = 0
  while ci < counts.length() {
    println("  Cluster " + ci.to_string() + ": " + counts[ci].to_string() + " spots")
    ci = ci + 1
  }
  println("")

  // ---- Section 6: Spots in each cluster ----
  println("=== 6. Spots in Cluster 0 ===")
  let c0 = @src.bayes_space_spots_in_cluster(result, 0)
  let mut c0i = 0
  while c0i < c0.length() {
    println("  " + c0[c0i])
    c0i = c0i + 1
  }
  println("")

  println("=== 7. Spots in Cluster 1 ===")
  let c1 = @src.bayes_space_spots_in_cluster(result, 1)
  let mut c1i = 0
  while c1i < c1.length() {
    println("  " + c1[c1i])
    c1i = c1i + 1
  }
  println("")

  // ---- Section 8: Render cluster map ----
  println("=== 8. Cluster Map ===")
  let map = @src.bayes_space_render_clusters(result, spots)
  println(map)
  println("")

  // ---- Section 9: Get cluster by spot ID ----
  println("=== 9. Get Cluster by Spot ID ===")
  let cl = @src.bayes_space_get_cluster(result, "spot_0_0")
  println("  spot_0_0 -> cluster " + cl.to_string())
  let cl2 = @src.bayes_space_get_cluster(result, "spot_3_3")
  println("  spot_3_3 -> cluster " + cl2.to_string())
  println("")

  // ---- Section 10: Square grid neighbors ----
  println("=== 10. Square Grid Neighbors ===")
  let sq_neighbors = @src.bayes_space_square_neighbors(spots)
  let mut sqi = 0
  while sqi < 4 {
    println("  " + spots[sqi].spot_id + " -> " + sq_neighbors[sqi].length().to_string() + " neighbors")
    sqi = sqi + 1
  }
  println("")

  println("=== Done ===")
}
