///|
/// Tests for Bioconductor BayesSpace module - Spatial transcriptomics clustering.

test "spot_coord_create" {
  let s = @src.spot_coord("spot_1_2", 1, 2, 2.0, 1.0)
  assert_eq(s.spot_id, "spot_1_2")
  assert_eq(s.row, 1)
  assert_eq(s.col, 2)
  assert_true((s.x - 2.0).abs() < 1.0e-9)
  assert_true((s.y - 1.0).abs() < 1.0e-9)
}

test "bayes_space_hex_neighbors_basic" {
  // 3x3 grid of spots
  let spots : Array[@src.SpotCoord] = []
  let mut r = 0
  while r < 3 {
    let mut c = 0
    while c < 3 {
      spots.push(@src.spot_coord("s_" + r.to_string() + "_" + c.to_string(), r, c, c.to_double(), r.to_double()))
      c = c + 1
    }
    r = r + 1
  }
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  // Center spot (1, 1) should have several neighbors
  let center_idx = 4 // (1, 1) in 3x3 grid
  assert_true(neighbors[center_idx].length() > 0)
  // Corner spot (0, 0) should have fewer neighbors
  let corner_idx = 0
  assert_true(neighbors[corner_idx].length() < neighbors[center_idx].length())
}

test "bayes_space_square_neighbors_basic" {
  let spots : Array[@src.SpotCoord] = []
  let mut r = 0
  while r < 3 {
    let mut c = 0
    while c < 3 {
      spots.push(@src.spot_coord("s_" + r.to_string() + "_" + c.to_string(), r, c, c.to_double(), r.to_double()))
      c = c + 1
    }
    r = r + 1
  }
  let neighbors = @src.bayes_space_square_neighbors(spots)
  // Center spot (1, 1) should have 4 neighbors in square grid
  let center_idx = 4
  assert_eq(neighbors[center_idx].length(), 4)
}

test "bayes_space_sample_data_shape" {
  let (expr, spots) = @src.bayes_space_sample_data()
  assert_eq(expr.length(), 16) // 4x4 grid
  assert_eq(spots.length(), 16)
  assert_eq(expr[0].length(), 3) // 3 features per spot
}

test "bayes_space_run_basic" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=20, seed=42)
  assert_eq(result.q, 2)
  assert_eq(result.spot_ids.length(), 16)
  assert_eq(result.clusters.length(), 16)
  assert_eq(result.cluster_centers.length(), 2)
  // Each cluster center should have 3 features
  assert_eq(result.cluster_centers[0].length(), 3)
  assert_eq(result.cluster_centers[1].length(), 3)
  // Each spot's responsibilities sum to ~1
  let mut i = 0
  while i < 16 {
    let mut sum = 0.0
    let mut k = 0
    while k < 2 {
      sum = sum + result.responsibilities[i][k]
      k = k + 1
    }
    assert_true((sum - 1.0).abs() < 1.0e-6)
    i = i + 1
  }
}

test "bayes_space_run_convergence_or_max_iters" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=50, seed=42)
  // Should either converge or hit max iters
  assert_true(result.n_iterations > 0)
  assert_true(result.n_iterations <= 50)
}

test "bayes_space_two_clusters_separated" {
  // The sample data has top-left (low) and bottom-right (high) clusters.
  // After clustering with q=2, the corner spots should be in different clusters.
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=30, seed=42)
  // Find indices of corner spots
  let mut tl_idx = -1
  let mut br_idx = -1
  let mut i = 0
  while i < spots.length() {
    if spots[i].row == 0 && spots[i].col == 0 {
      tl_idx = i
    }
    if spots[i].row == 3 && spots[i].col == 3 {
      br_idx = i
    }
    i = i + 1
  }
  assert_true(tl_idx >= 0)
  assert_true(br_idx >= 0)
  // They should be in different clusters (top-left is low, bottom-right is high)
  assert_not_eq(result.clusters[tl_idx], result.clusters[br_idx])
}

test "bayes_space_get_cluster_by_id" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=20, seed=42)
  let c = @src.bayes_space_get_cluster(result, "spot_0_0")
  assert_true(c >= 0 && c < 2)
}

test "bayes_space_get_cluster_unknown_id" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=10, seed=42)
  let c = @src.bayes_space_get_cluster(result, "UNKNOWN_SPOT")
  assert_eq(c, -1)
}

test "bayes_space_spots_in_cluster" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=20, seed=42)
  let c0 = @src.bayes_space_spots_in_cluster(result, 0)
  let c1 = @src.bayes_space_spots_in_cluster(result, 1)
  // Both clusters should have at least one spot
  assert_true(c0.length() > 0)
  assert_true(c1.length() > 0)
  // Total should equal total spots
  assert_eq(c0.length() + c1.length(), 16)
}

test "bayes_space_cluster_counts" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=20, seed=42)
  let counts = @src.bayes_space_cluster_counts(result)
  assert_eq(counts.length(), 2)
  // Sum of counts should equal total spots
  assert_eq(counts[0] + counts[1], 16)
}

test "bayes_space_render_clusters" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=10, seed=42)
  let render = @src.bayes_space_render_clusters(result, spots)
  // Should contain 4 rows separated by newlines
  assert_true(render.length() > 0)
  let mut nl_count = 0
  let mut i = 0
  while i < render.length() {
    if render[i] == '\n' {
      nl_count = nl_count + 1
    }
    i = i + 1
  }
  assert_eq(nl_count, 4)
}

test "bayes_space_q_equals_3" {
  // Test with q=3 clusters
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 3, neighbors, max_iters=15, seed=42)
  assert_eq(result.q, 3)
  assert_eq(result.cluster_centers.length(), 3)
  // Each spot's cluster should be in [0, 3)
  let mut i = 0
  while i < result.clusters.length() {
    assert_true(result.clusters[i] >= 0 && result.clusters[i] < 3)
    i = i + 1
  }
}

test "bayes_space_empty_grid_neighbors" {
  // Empty spots list → empty neighbor list
  let neighbors = @src.bayes_space_hex_neighbors([])
  assert_eq(neighbors.length(), 0)
}

test "bayes_space_single_spot" {
  // Single spot → no neighbors
  let spots = [@src.spot_coord("only", 0, 0, 0.0, 0.0)]
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  assert_eq(neighbors.length(), 1)
  assert_eq(neighbors[0].length(), 0)
}

test "bayes_space_log_likelihood_finite" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let result = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=20, seed=42)
  // Log-likelihood should be a finite number (not NaN/Inf)
  assert_false(result.log_likelihood.is_nan())
  assert_false(result.log_likelihood.abs() > 1.0e300)
}

test "bayes_space_reproducible_with_same_seed" {
  let (expr, spots) = @src.bayes_space_sample_data()
  let neighbors = @src.bayes_space_hex_neighbors(spots)
  let r1 = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=15, seed=123)
  let r2 = @src.bayes_space_run(expr, spots, 2, neighbors, max_iters=15, seed=123)
  // Same seed → same cluster assignments
  let mut i = 0
  while i < r1.clusters.length() {
    assert_eq(r1.clusters[i], r2.clusters[i])
    i = i + 1
  }
}
