///|
/// BiocNeighbors example: Nearest neighbor search for single-cell analysis.
/// Demonstrates KMKNN, Annoy, and BruteForce methods with various distance metrics.
///
/// Run: moon run examples/bioc_neighbors_demo/main.mbt
fn main {
  println("=== BiocNeighbors Demo ===")
  println("")

  // === Create sample data ===
  println("--- 1. Create Sample Data (2D points) ---")

  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 2.0],
    [2.0, 3.0],
    [3.0, 3.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let n_points = data.length()
  let n_dim = 2

  println("  Created \{n_points} points in \{n_dim}D space")
  for i = 0; i < n_points; i = i + 1 {
    println("    Point \{i}: (\{data[i][0]}, \{data[i][1]})")
  }
  println("")

  // === Brute Force KNN ===
  println("--- 2. Brute Force KNN ---")
  let query : Array[Array[Double]] = [[0.0, 0.0]]
  let result1 = @bio.knn_brute_force(data, n_points, n_dim, query, 1, 3, "euclidean")
  println("  Query: (0.0, 0.0)")
  println("  k=3, euclidean distance")
  println("  Found neighbors:")
  for i = 0; i < result1.indices[0].length(); i = i + 1 {
    let idx = result1.indices[0][i]
    let dist = result1.distances[0][i]
    println("    Neighbor \{i}: index=\{idx}, point=(\{data[idx][0]}, \{data[idx][1]}), dist=\{dist}")
  }
  println("")

  // === Build KMKNN Index ===
  println("--- 3. KMKNN (K-Means based KNN) ---")
  let kmknn_opts = @bio.knn_options_full(@bio.kmknn_method(), 3, "euclidean", 50, -1)
  let kmknn_index = @bio.build_kmknn_index(data, n_points, n_dim, kmknn_opts)
  println("  Index built: \{kmknn_index.centroids.length()} centroids")
  for c = 0; c < kmknn_index.centroids.length(); c = c + 1 {
    println("    Centroid \{c}: (\{kmknn_index.centroids[c][0]}, \{kmknn_index.centroids[c][1]})")
  }

  let query2 : Array[Array[Double]] = [[0.0, 0.0], [10.0, 10.0]]
  let kmknn_result = @bio.knn_kmknn(kmknn_index, query2, 2, kmknn_opts)
  println("  Query 1: (0.0, 0.0)")
  for i = 0; i < kmknn_result.indices[0].length(); i = i + 1 {
    let idx = kmknn_result.indices[0][i]
    let dist = kmknn_result.distances[0][i]
    println("    Neighbor \{i}: index=\{idx}, dist=\{dist}")
  }
  println("  Query 2: (10.0, 10.0)")
  for i = 0; i < kmknn_result.indices[1].length(); i = i + 1 {
    let idx = kmknn_result.indices[1][i]
    let dist = kmknn_result.distances[1][i]
    println("    Neighbor \{i}: index=\{idx}, dist=\{dist}")
  }
  println("")

  // === Build Annoy Index ===
  println("--- 4. Annoy (Approximate Nearest Neighbors) ---")
  let annoy_opts = @bio.knn_options_full(@bio.annoy_method(), 3, "euclidean", 5, -1)
  let annoy_index = @bio.build_annoy_index(data, n_points, n_dim, annoy_opts)
  println("  Index built: \{annoy_index.tree_nodes.length()} tree nodes (5 trees)")

  let annoy_result = @bio.knn_annoy(annoy_index, query2, 2, annoy_opts)
  println("  Query 1: (0.0, 0.0)")
  for i = 0; i < annoy_result.indices[0].length(); i = i + 1 {
    let idx = annoy_result.indices[0][i]
    let dist = annoy_result.distances[0][i]
    println("    Neighbor \{i}: index=\{idx}, dist=\{dist}")
  }
  println("  Query 2: (10.0, 10.0)")
  for i = 0; i < annoy_result.indices[1].length(); i = i + 1 {
    let idx = annoy_result.indices[1][i]
    let dist = annoy_result.distances[1][i]
    println("    Neighbor \{i}: index=\{idx}, dist=\{dist}")
  }
  println("")

  // === Distance Metrics ===
  println("--- 5. Distance Metrics ---")
  let p1 = [1.0, 2.0, 3.0]
  let p2 = [4.0, 6.0, 3.0]
  println("  Point 1: (1.0, 2.0, 3.0)")
  println("  Point 2: (4.0, 6.0, 3.0)")
  println("    Euclidean: \{@bio.euclidean_distance(p1, p2)}")
  println("    Manhattan: \{@bio.manhattan_distance(p1, p2)}")
  println("    Cosine:    \{@bio.cosine_distance(p1, p2)}")
  println("")

  // === Run KNN with Index ===
  println("--- 6. Unified run_knn Interface ---")
  let brute_opts = @bio.knn_options_full(@bio.brute_force_method(), 2, "euclidean", 50, -1)
  let brute_index = @bio.build_knn_index(data, n_points, n_dim, brute_opts)
  let query3 : Array[Array[Double]] = [[0.5, 0.5]]
  let final_result = @bio.run_knn(brute_index, query3, 1, brute_opts)
  println("  Query: (0.5, 0.5)")
  println("  k=2, euclidean distance")
  for i = 0; i < final_result.indices[0].length(); i = i + 1 {
    let idx = final_result.indices[0][i]
    let dist = final_result.distances[0][i]
    println("    Neighbor \{i}: index=\{idx}, point=(\{data[idx][0]}, \{data[idx][1]}), dist=\{dist}")
  }
  println("")

  // === Edge Cases ===
  println("--- 7. Edge Cases ---")

  // Empty data
  let empty_data : Array[Array[Double]] = []
  let empty_result = @bio.knn_brute_force(empty_data, 0, 2, query, 1, 2, "euclidean")
  println("  Empty data: n_query=\{empty_result.n_query}, k=\{empty_result.k}")

  // k > n_points
  let big_k_result = @bio.knn_brute_force(data, n_points, n_dim, query, 1, 100, "euclidean")
  println("  k > n_points: k=\{big_k_result.k}, indices length=\{big_k_result.indices[0].length()}")

  // Single point
  let single_data = [[5.0, 5.0]]
  let single_result = @bio.knn_brute_force(single_data, 1, 2, query, 1, 1, "euclidean")
  println("  Single point: k=\{single_result.k}, indices=\{single_result.indices[0][0]}")

  println("")
  println("=== Demo Complete ===")
}