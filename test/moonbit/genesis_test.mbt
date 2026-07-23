///|
/// Tests for GENESIS module.

///|
test "genesis_estimate_kinship" {
  let genotypes = [
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [1.0, 2.0, 0.0],
    [2.0, 0.0, 1.0],
  ]
  
  let result = @src.bio_estimate_kinship(genotypes)
  
  assert_true(result.n_samples == 4)
  assert_true(result.kinship_matrix.length() == 4)
  assert_true(result.sample_ids.length() == 4)
}

///|
test "genesis_pca" {
  let genotypes = [
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [1.0, 2.0, 0.0],
    [2.0, 0.0, 1.0],
  ]
  
  let result = @src.bio_pca(genotypes, 2)
  
  assert_true(result.eigenvalues.length() == 2)
  assert_true(result.eigenvectors.length() > 0)
  assert_true(result.var_explained.length() == 2)
}

///|
test "genesis_genetic_distance" {
  let genotypes = [
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [1.0, 2.0, 0.0],
  ]
  
  let distance = @src.bio_genetic_distance(genotypes)
  
  assert_true(distance.distance_matrix.length() == 3)
  assert_true(distance.sample_ids.length() == 3)
}

///|
test "genesis_euclidean_distance" {
  let v1 = [0.0, 1.0, 2.0]
  let v2 = [0.0, 1.0, 2.0]
  
  let dist = @src.gs_euclidean_distance(v1, v2)
  
  assert_true(dist == 0.0)
}

///|
test "genesis_manhattan_distance" {
  let v1 = [0.0, 1.0, 2.0]
  let v2 = [0.0, 1.0, 2.0]
  
  let dist = @src.gs_manhattan_distance(v1, v2)
  
  assert_true(dist == 0.0)
}

///|
test "genesis_ibs_distance" {
  let v1 = [0.0, 1.0, 2.0]
  let v2 = [0.0, 1.0, 2.0]
  
  let dist = @src.gs_ibs_distance(v1, v2)
  
  assert_true(dist >= 0.0)
}
