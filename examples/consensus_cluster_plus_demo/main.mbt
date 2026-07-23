///|
/// ConsensusClusterPlus demo.
///
/// Demonstrates consensus clustering for cancer subtyping.

pub fn main() {
  println("=== ConsensusClusterPlus Demo ===")
  println("")
  
  let data = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [10.0, 11.0, 12.0, 13.0],
    [11.0, 12.0, 13.0, 14.0],
    [12.0, 13.0, 14.0, 15.0],
    [20.0, 21.0, 22.0, 23.0],
    [21.0, 22.0, 23.0, 24.0],
  ]
  
  println("Input data: 8 samples, 4 features")
  println("")
  
  let (best_k, scores) = bio_consensus_cluster_find_optimal_k(data, 2, 4)
  
  println("Optimal k search results:")
  let mut k = 2
  let mut idx = 0
  while idx < scores.length() {
    println("  k={}: consensus_score={:.4f}", k, scores[idx])
    k = k + 1
    idx = idx + 1
  }
  println("  Best k: {}", best_k)
  println("")
  
  let result = bio_consensus_cluster(data, best_k)
  
  println("Consensus clustering result (k={}):", best_k)
  println("  Consensus score: {:.4f}", result.consensus_score)
  println("  Cluster labels: {:?}", result.cluster_labels)
  println("  Stability scores: {:?}", result.stability_scores)
  println("")
  
  println("Consensus matrix:")
  let mut i = 0
  while i < result.consensus_matrix.length() {
    let mut row_str = ""
    let mut j = 0
    while j < result.consensus_matrix[i].length() {
      row_str = row_str.concat("{:.2f} ".format(result.consensus_matrix[i][j]))
      j = j + 1
    }
    println("  {}", row_str)
    i = i + 1
  }
  
  println("")
  println("=== Demo Complete ===")
}