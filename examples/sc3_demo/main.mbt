///|
/// SC3 demo.
///
/// Demonstrates single-cell consensus clustering.

pub fn main() {
  println("=== SC3 Demo ===")
  println("")
  
  let data = [
    [100.0, 200.0, 30.0, 40.0, 50.0],
    [110.0, 210.0, 35.0, 45.0, 55.0],
    [95.0, 195.0, 28.0, 38.0, 48.0],
    [10.0, 20.0, 300.0, 400.0, 500.0],
    [12.0, 22.0, 310.0, 410.0, 510.0],
    [8.0, 18.0, 290.0, 390.0, 490.0],
    [500.0, 510.0, 520.0, 10.0, 20.0],
    [510.0, 520.0, 530.0, 12.0, 22.0],
  ]
  
  println("Input data: 8 cells, 5 genes")
  println("")
  
  let processed = bio_sc3_preprocess(data)
  println("Preprocessing complete")
  println("")
  
  let pcs = bio_sc3_pca(data, 3)
  println("PCA complete (3 components)")
  println("  Eigenvalues: {:?}", pcs.eigenvalues)
  println("  Variance explained: {:?}", pcs.var_explained)
  println("")
  
  let result = bio_sc3_cluster(data, 3)
  
  println("SC3 clustering result (k={}):", result.k)
  println("  Optimal k: {}", result.optimal_k)
  println("  Cluster labels: {:?}", result.cluster_labels)
  println("")
  
  let mut avg_silhouette = 0.0
  let mut i = 0
  while i < result.silhouette_scores.length() {
    avg_silhouette = avg_silhouette + result.silhouette_scores[i]
    i = i + 1
  }
  avg_silhouette = avg_silhouette / result.silhouette_scores.length().to_double()
  println("  Average silhouette score: {:.4f}", avg_silhouette)
  println("")
  
  println("Consensus matrix:")
  let mut j = 0
  while j < result.consensus_matrix.length() {
    let mut row_str = ""
    let mut k = 0
    while k < result.consensus_matrix[j].length() {
      row_str = row_str.concat("{:.2f} ".format(result.consensus_matrix[j][k]))
      k = k + 1
    }
    println("  {}", row_str)
    j = j + 1
  }
  
  println("")
  println("=== Demo Complete ===")
}