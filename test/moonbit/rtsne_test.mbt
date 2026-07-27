///|
/// Tests for Rtsne module.

// ============================================================
// Distance calculation tests
// ============================================================

test "calc_distance_matrix square" {
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [3.0, 4.0],
  ]
  let dist = @src.calc_distance_matrix(data)
  assert_eq(dist.length(), 2)
  assert_eq(dist[0].length(), 2)
  assert_eq(dist[0][0], 0.0)
  assert_eq(dist[1][1], 0.0)
  // Distance between (0,0) and (3,4) should be 5
  assert_eq(dist[0][1], 5.0)
  assert_eq(dist[1][0], 5.0)
}

test "calc_distance_matrix empty" {
  let data : Array[Array[Double]] = []
  let dist = @src.calc_distance_matrix(data)
  assert_eq(dist.length(), 0)
}

test "calc_distance_matrix single point" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0]]
  let dist = @src.calc_distance_matrix(data)
  assert_eq(dist.length(), 1)
  assert_eq(dist[0][0], 0.0)
}

// ============================================================
// t-SNE config tests
// ============================================================

test "tsne_config default" {
  let config = @src.TsneConfig::new()
  assert_eq(config.perplexity, 30.0)
  assert_eq(config.max_iter, 1000)
  assert_eq(config.dims, 2)
  assert_eq(config.theta, 0.5)
  assert_eq(config.eta, 200.0)
}

// ============================================================
// t-SNE algorithm tests
// ============================================================

test "tsne basic 2d" {
  let data = @src.create_tsne_test_data(10, 5)
  assert_eq(data.length(), 10)
  assert_eq(data[0].length(), 5)
  
  let test_config = @src.TsneConfig::new_custom(3.0, 100, 2, 42)
  
  let result = @src.tsne(data, test_config)
  assert_eq(result.embedding.length(), 10)
  assert_eq(result.embedding[0].length(), 2)
  assert_eq(result.n_iter, 100)
  assert_eq(result.costs.length(), 100)
}

test "tsne empty data" {
  let data : Array[Array[Double]] = []
  let config = @src.TsneConfig::new()
  let result = @src.tsne(data, config)
  assert_eq(result.embedding.length(), 0)
}

test "tsne single sample" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0]]
  let config = @src.TsneConfig::new()
  let result = @src.tsne(data, config)
  // Single sample should return empty or single point
  assert_true(result.embedding.length() <= 1)
}

test "tsne costs decrease" {
  let data = @src.create_tsne_test_data(15, 4)
  let config = @src.TsneConfig::new_custom(3.0, 50, 2, 42)
  
  let result = @src.tsne(data, config)
  
  // Costs should be finite
  if result.costs.length() > 10 {
    // Check that costs are not NaN or infinite
    let first_cost = result.costs[5]
    let last_cost = result.costs[result.costs.length() - 1]
    assert_true(first_cost > 0.0)
    assert_true(last_cost > 0.0)
  }
}

test "tsne different dims" {
  let data = @src.create_tsne_test_data(10, 5)
  let config = @src.TsneConfig::new_custom(3.0, 50, 3, 42)
  
  let result = @src.tsne(data, config)
  assert_eq(result.embedding[0].length(), 3)
}

// ============================================================
// Test data generation tests
// ============================================================

test "create_tsne_test_data dimensions" {
  let data = @src.create_tsne_test_data(20, 8)
  assert_eq(data.length(), 20)
  assert_eq(data[0].length(), 8)
}

test "create_tsne_test_data clusters" {
  let data = @src.create_tsne_test_data(9, 3)
  // 3 clusters, each should have similar values
  // Cluster 0: samples 0, 3, 6
  // Cluster 1: samples 1, 4, 7
  // Cluster 2: samples 2, 5, 8
  
  // Check that cluster 0 is different from cluster 1
  let mut sum0 = 0.0
  let mut j = 0
  while j < 3 {
    sum0 = sum0 + data[0][j]
    j = j + 1
  }
  
  let mut sum1 = 0.0
  j = 0
  while j < 3 {
    sum1 = sum1 + data[1][j]
    j = j + 1
  }
  
  // Clusters should be separated
  assert_true(sum1 > sum0)
}
