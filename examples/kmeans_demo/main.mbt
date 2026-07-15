///|
/// K-means Clustering Demo - Biological Data Analysis
///
/// This demo shows how to use K-means clustering for analyzing
/// gene expression data and other biological data.
fn main {
  println("=== K-means Clustering Demo ===")
  println("")

  println("1. Distance calculations:")
  let a = [1.0, 2.0]
  let b = [4.0, 6.0]
  println("   Point A: \{a}")
  println("   Point B: \{b}")
  println("   Euclidean distance: \{@src.euclidean_distance(a, b)}")
  println("   Manhattan distance: \{@src.manhattan_distance(a, b)}")
  println("")

  println("2. Simple clustering example:")
  let data = [
    [1.0, 2.0], [2.0, 1.0], [2.0, 3.0],
    [8.0, 7.0], [9.0, 8.0], [7.0, 9.0],
    [4.0, 5.0], [5.0, 4.0], [5.0, 6.0]
  ]
  let kmeans = @src.KMeans::new(3, 100, 0.000001).fit(data)
  println("   Data points: \{data.length()}")
  println("   Labels: \{kmeans.labels}")
  println("   Centroids: \{kmeans.centroids}")
  println("")

  println("3. Gene expression clustering:")
  let gene_expr = [
    [1.2, 3.4, 2.1], [1.5, 3.6, 2.3], [1.3, 3.2, 2.0],
    [8.5, 7.2, 6.8], [8.8, 7.5, 7.1], [8.3, 7.0, 6.5],
    [4.2, 5.1, 4.8], [4.5, 5.3, 5.0], [4.3, 5.0, 4.7]
  ]
  let gene_kmeans = @src.cluster_gene_expression(gene_expr, 3)
  println("   Genes: \{gene_expr.length()}")
  println("   Clusters: \{gene_kmeans.labels}")
  println("")

  println("4. Prediction on new data:")
  let new_points = [[1.4, 2.5], [8.6, 7.3], [4.4, 5.2]]
  let predictions = gene_kmeans.predict(new_points)
  println("   New points: \{new_points}")
  println("   Predictions: \{predictions}")
  println("")

  println("5. Single point prediction:")
  let single_pred = gene_kmeans.predict_single([5.0, 5.5, 5.2])
  println("   Point: [5.0, 5.5, 5.2]")
  println("   Cluster: \{single_pred}")
  println("")

  println("6. Inertia (sum of squared distances):")
  let inertia = gene_kmeans.inertia(gene_expr)
  println("   Inertia: \{inertia}")
  println("")

  println("7. Silhouette score:")
  let simple_data = [
    [1.0, 2.0], [2.0, 1.0], [2.0, 3.0],
    [8.0, 7.0], [9.0, 8.0], [7.0, 9.0]
  ]
  let simple_labels = [0, 0, 0, 1, 1, 1]
  let score = @src.silhouette_score(simple_data, simple_labels)
  println("   Score: \{score}")
  println("")

  println("8. Find optimal k:")
  let opt_data = [
    [1.0, 2.0], [2.0, 1.0], [2.0, 3.0],
    [8.0, 7.0], [9.0, 8.0], [7.0, 9.0],
    [4.0, 5.0], [5.0, 4.0], [5.0, 6.0]
  ]
  let optimal_k = @src.find_optimal_k(opt_data, 2, 5)
  println("   Optimal k: \{optimal_k}")
  println("")

  println("=== Demo Complete ===")
}
