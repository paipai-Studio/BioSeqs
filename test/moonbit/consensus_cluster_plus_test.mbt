///|
/// Tests for ConsensusClusterPlus module.

///|
test "consensus_cluster_plus_consensus_matrix" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [10.0, 11.0, 12.0],
    [11.0, 12.0, 13.0],
  ]
  
  let consensus = @src.ccp_calculate_consensus_matrix(data, 2, 10, 0.8)
  
  assert_true(consensus.length() == 4)
  assert_true(consensus[0].length() == 4)
}

///|
test "consensus_cluster_plus_kmeans" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
    [11.0, 12.0],
  ]
  
  let labels = @src.ccp_kmeans_cluster(data, 2)
  
  assert_true(labels.length() == 4)
}

///|
test "consensus_cluster_plus_consensus_score" {
  let consensus = [
    [1.0, 0.9, 0.8, 0.7],
    [0.9, 1.0, 0.8, 0.7],
    [0.8, 0.8, 1.0, 0.9],
    [0.7, 0.7, 0.9, 1.0],
  ]
  
  let score = @src.ccp_calculate_consensus_score(consensus)
  
  assert_true(score > 0.5)
}

///|
test "consensus_cluster_plus_find_optimal_k" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [10.0, 11.0, 12.0],
    [11.0, 12.0, 13.0],
    [20.0, 21.0, 22.0],
    [21.0, 22.0, 23.0],
  ]
  
  let (best_k, scores) = @src.ccp_find_optimal_k(data, 2, 4)
  
  assert_true(best_k >= 2 && best_k <= 4)
  assert_true(scores.length() == 3)
}

///|
test "consensus_cluster_plus_main" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [10.0, 11.0, 12.0],
    [11.0, 12.0, 13.0],
  ]
  
  let result = @src.bio_consensus_cluster(data, 2)
  
  assert_true(result.k == 2)
  assert_true(result.cluster_labels.length() == 4)
  assert_true(result.consensus_matrix.length() == 4)
}

///|
test "consensus_cluster_plus_bio_api" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
  ]
  
  let (best_k, scores) = @src.bio_consensus_cluster_find_optimal_k(data, 2, 3)
  let consensus = @src.ccp_calculate_consensus_matrix(data, 2, 10, 0.8)
  let score = @src.bio_consensus_cluster_consensus_score(consensus)
  
  assert_true(best_k >= 2)
  assert_true(scores.length() > 0)
  assert_true(score >= 0.0)
}