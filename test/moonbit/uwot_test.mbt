///|
/// Tests for uwot (UMAP) module.

// ============================================================
// Distance calculation tests
// ============================================================

test "uwot_distance_matrix square" {
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [3.0, 4.0],
  ]
  let dist = @src.uwot_distance_matrix(data)
  assert_eq(dist.length(), 2)
  assert_eq(dist[0].length(), 2)
  assert_eq(dist[0][0], 0.0)
  assert_eq(dist[1][1], 0.0)
  // Distance between (0,0) and (3,4) should be 5
  assert_eq(dist[0][1], 5.0)
  assert_eq(dist[1][0], 5.0)
}

test "uwot_distance_matrix empty" {
  let data : Array[Array[Double]] = []
  let dist = @src.uwot_distance_matrix(data)
  assert_eq(dist.length(), 0)
}

// ============================================================
// UMAP config tests
// ============================================================

test "umap_config default" {
  let config = @src.UmapConfig::new()
  assert_eq(config.n_neighbors, 15)
  assert_eq(config.n_components, 2)
  assert_eq(config.n_epochs, 200)
  assert_eq(config.min_dist, 0.1)
  assert_eq(config.spread, 1.0)
  assert_eq(config.learning_rate, 1.0)
}

// ============================================================
// UMAP algorithm tests
// ============================================================

test "umap basic 2d" {
  let data = @src.create_umap_test_data(10, 5)
  assert_eq(data.length(), 10)
  assert_eq(data[0].length(), 5)
  
  let config = @src.UmapConfig::new()
  
  let result = @src.uwot_umap(data, config)
  assert_eq(result.embedding.length(), 10)
  assert_eq(result.embedding[0].length(), 2)
  assert_eq(result.n_epochs, 200)
}

test "umap empty data" {
  let data : Array[Array[Double]] = []
  let config = @src.UmapConfig::new()
  let result = @src.uwot_umap(data, config)
  assert_eq(result.embedding.length(), 0)
}

test "umap single sample" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0]]
  let config = @src.UmapConfig::new()
  let result = @src.uwot_umap(data, config)
  assert_true(result.embedding.length() <= 1)
}

test "umap different n_components" {
  let data = @src.create_umap_test_data(10, 5)
  let custom_config = @src.UmapConfig::new_custom(3, 3, 50, 42)
  
  let result = @src.uwot_umap(data, custom_config)
  assert_eq(result.embedding[0].length(), 3)
}

test "umap produces finite values" {
  let data = @src.create_umap_test_data(12, 4)
  let custom_config = @src.UmapConfig::new_custom(3, 2, 50, 42)
  
  let result = @src.uwot_umap(data, custom_config)
  
  // Check that all embedding values are finite
  let mut i = 0
  while i < result.embedding.length() {
    let mut j = 0
    while j < result.embedding[i].length() {
      let val = result.embedding[i][j]
      assert_true(val == val)  // NaN check
      assert_true(val < 1.0e300)  // infinity check
      assert_true(val > -1.0e300)
      j = j + 1
    }
    i = i + 1
  }
}

// ============================================================
// Test data generation tests
// ============================================================

test "create_umap_test_data dimensions" {
  let data = @src.create_umap_test_data(20, 8)
  assert_eq(data.length(), 20)
  assert_eq(data[0].length(), 8)
}

test "create_umap_test_data clusters" {
  let data = @src.create_umap_test_data(8, 3)
  // 4 clusters, each should have similar values
  // Cluster 0: samples 0, 4
  // Cluster 1: samples 1, 5
  // Cluster 2: samples 2, 6
  // Cluster 3: samples 3, 7
  
  // Check that cluster 0 is different from cluster 3
  let mut sum0 = 0.0
  let mut j = 0
  while j < 3 {
    sum0 = sum0 + data[0][j]
    j = j + 1
  }
  
  let mut sum3 = 0.0
  j = 0
  while j < 3 {
    sum3 = sum3 + data[3][j]
    j = j + 1
  }
  
  // Clusters should be separated
  assert_true(sum3 > sum0)
}
