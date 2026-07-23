///|
/// Tests for SC3 module.

///|
test "sc3_preprocess" {
  let data = [
    [0.0, 1.0, 10.0],
    [0.0, 2.0, 20.0],
    [0.0, 3.0, 30.0],
  ]
  
  let processed = @src.sc3_preprocess(data)
  
  assert_true(processed.length() == 3)
  assert_true(processed[0].length() == 3)
}

///|
test "sc3_pca" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [10.0, 11.0, 12.0],
  ]
  
  let pcs = @src.sc3_pca(data, 2)
  
  assert_true(pcs.length() == 4)
  assert_true(pcs[0].length() == 2)
}

///|
test "sc3_kmeans" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
    [11.0, 12.0],
  ]
  
  let labels = @src.sc3_kmeans(data, 2)
  
  assert_true(labels.length() == 4)
}

///|
test "sc3_silhouette" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
    [11.0, 12.0],
  ]
  let labels = [0, 0, 1, 1]
  
  let scores = @src.sc3_calculate_silhouette(data, labels, 2)
  
  assert_true(scores.length() == 4)
}

///|
test "sc3_gap_statistics" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
    [11.0, 12.0],
  ]
  
  let gaps = @src.sc3_calculate_gap_statistics(data, 3)
  
  assert_true(gaps.length() == 3)
}

///|
test "sc3_consensus_cluster" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [10.0, 11.0, 12.0],
    [11.0, 12.0, 13.0],
  ]
  
  let result = @src.bio_sc3_cluster(data, 2)
  
  assert_true(result.k == 2)
  assert_true(result.cluster_labels.length() == 4)
  assert_true(result.consensus_matrix.length() == 4)
}

///|
test "sc3_bio_api" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [10.0, 11.0],
  ]
  
  let pcs = @src.bio_sc3_pca(data, 2)
  let labels = @src.sc3_kmeans(pcs, 2)
  let scores = @src.bio_sc3_silhouette(pcs, labels, 2)
  
  assert_true(pcs.length() == 3)
  assert_true(scores.length() == 3)
}