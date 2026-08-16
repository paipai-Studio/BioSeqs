// Tests for the Bioconductor Voyager-inspired spatial autocorrelation module.

///|
fn voy_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "Voyager value " +
      actual.to_string() +
      " differs from " +
      expected.to_string() +
      " (tolerance " +
      tolerance.to_string() +
      ")",
    )
  }
}

///|
fn voy_test_grid_coords() -> Array[Array[Double]] {
  // 3x3 grid of spots, coordinates [x, y].
  let coords : Array[Array[Double]] = []
  for row in 0..<3 {
    for col in 0..<3 {
      coords.push([col.to_double(), row.to_double()])
    }
  }
  coords
}

///|
fn voy_test_chain_weights() -> @src.VoyagerWeights {
  // 4 spots in a chain 0-1-2-3 with binary (style B) symmetric weights.
  @src.VoyagerWeights::new(
    [[1], [0, 2], [1, 3], [2]],
    [[1.0], [1.0, 1.0], [1.0, 1.0], [1.0]],
    style="B",
  ) catch {
    _ => abort("valid chain weights should build")
  }
}

///|
fn voy_test_gradient_values() -> Array[Double] {
  // Values increasing with the x coordinate on the 3x3 grid: strong positive
  // spatial autocorrelation.
  [
    0.0, 1.0, 2.0, // row 0
     0.0, 1.0, 2.0, // row 1
     0.0, 1.0, 2.0, // row 2
  ]
}

///|
fn voy_test_checkerboard_values() -> Array[Double] {
  // Alternating high/low: negative spatial autocorrelation on the grid.
  [0.0, 9.0, 0.0, 9.0, 0.0, 9.0, 0.0, 9.0, 0.0]
}

///|
fn voy_test_hotspot_values() -> Array[Double] {
  // 3x3 grid with a hotspot in the centre.
  [0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0]
}

// ===========================================================================
// Weights construction
// ===========================================================================

///|
test "Voyager weights kNN builds row-standardized neighbours" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  assert_eq(weights.n, 9)
  assert_eq(weights.style, "W")
  // Each spot on a 3x3 grid has 4 neighbours (corners: 2, edges: 3, centre: 4
  // are clamped to k=4 but corners/edges have fewer available neighbours).
  assert_true(weights.neighbors[4].length() == 4) // centre
  // Row-standardized: each row sums to 1.
  for index in 0..<weights.n {
    let mut row_sum = 0.0
    for value in weights.weights[index] {
      row_sum = row_sum + value
    }
    voy_test_close(row_sum, 1.0, 1.0e-12)
  }
}

///|
test "Voyager weights kNN rejects k below one" {
  let coords = voy_test_grid_coords()
  let failed = try {
    ignore(@src.voyager_weights_knn(coords, 0))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager weights kNN clamps k above n-1" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 100) catch {
    _ => abort("kNN weights with large k should clamp")
  }
  // With k clamped to n-1, every spot links to all other spots.
  assert_true(weights.neighbors[0].length() == 8)
}

///|
test "Voyager weights distance band links immediate grid neighbours" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_distance_band(coords, 1.0) catch {
    _ => abort("distance-band weights should build")
  }
  assert_eq(weights.n, 9)
  // Centre spot (index 4) has 4 immediate neighbours at distance 1.
  assert_eq(weights.neighbors[4].length(), 4)
  // Corner spot (index 0) has 2 immediate neighbours.
  assert_eq(weights.neighbors[0].length(), 2)
}

///|
test "Voyager weights inverse distance produces positive finite weights" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_inverse_distance(coords, 1.0) catch {
    _ => abort("inverse-distance weights should build")
  }
  assert_eq(weights.n, 9)
  for index in 0..<weights.n {
    for value in weights.weights[index] {
      assert_true(value > 0.0)
      assert_true(value <= 1.0e300)
    }
  }
}

///|
test "Voyager weights style C sums to one globally" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4, style="C") catch {
    _ => abort("style C weights should build")
  }
  let mut total = 0.0
  for index in 0..<weights.n {
    for value in weights.weights[index] {
      total = total + value
    }
  }
  voy_test_close(total, 1.0, 1.0e-12)
}

///|
test "Voyager weights style B uses binary ones" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4, style="B") catch {
    _ => abort("style B weights should build")
  }
  for index in 0..<weights.n {
    for value in weights.weights[index] {
      voy_test_close(value, 1.0, 0.0)
    }
  }
}

///|
test "Voyager weights S0 S1 S2 match hand-computed chain values" {
  let weights = voy_test_chain_weights()
  voy_test_close(@src.voyager_weights_s0(weights), 6.0, 1.0e-12)
  voy_test_close(@src.voyager_weights_s1(weights), 12.0, 1.0e-12)
  voy_test_close(@src.voyager_weights_s2(weights), 40.0, 1.0e-12)
}

///|
test "Voyager weights constructor rejects mismatched rows" {
  let failed = try {
    ignore(@src.VoyagerWeights::new([[1]], [[1.0], [1.0]], style="B"))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Global Moran's I
// ===========================================================================

///|
test "Voyager global Moran's I matches hand-computed chain example" {
  let weights = voy_test_chain_weights()
  let result = @src.voyager_global_morans_i([1.0, 1.0, 9.0, 9.0], weights) catch {
    _ => abort("global Moran's I should run")
  }
  voy_test_close(result.estimate, 1.0 / 3.0, 1.0e-10)
  voy_test_close(result.expectation, -1.0 / 3.0, 1.0e-12)
  voy_test_close(result.s0, 6.0, 1.0e-12)
  voy_test_close(result.s1, 12.0, 1.0e-12)
  voy_test_close(result.s2, 40.0, 1.0e-12)
  // b2 = 1, n = 4 -> variance = 88/216 - 1/9 = 0.296296...
  voy_test_close(result.variance, 88.0 / 216.0 - 1.0 / 9.0, 1.0e-10)
  assert_eq(result.n, 4)
}

///|
test "Voyager global Moran's I is positive for a spatial gradient" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_global_morans_i(voy_test_gradient_values(), weights) catch {
    _ => abort("global Moran's I should run")
  }
  assert_true(result.estimate > 0.0)
  assert_true(result.p_value < 0.05)
}

///|
test "Voyager global Moran's I is negative for a checkerboard" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_global_morans_i(
    voy_test_checkerboard_values(),
    weights,
  ) catch {
    _ => abort("global Moran's I should run")
  }
  assert_true(result.estimate < 0.0)
}

///|
test "Voyager global Moran's I rejects non-finite values" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let failed = try {
    ignore(@src.voyager_global_morans_i([1.0, 2.0, 1.0e400, 4.0], weights))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager global Moran's I rejects length mismatch" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let failed = try {
    ignore(@src.voyager_global_morans_i([1.0, 2.0], weights))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Global Geary's c
// ===========================================================================

///|
test "Voyager global Geary's c matches hand-computed chain example" {
  let weights = voy_test_chain_weights()
  let result = @src.voyager_global_gearys_c([1.0, 1.0, 9.0, 9.0], weights) catch {
    _ => abort("global Geary's c should run")
  }
  voy_test_close(result.estimate, 0.5, 1.0e-10)
  voy_test_close(result.expectation, 1.0, 1.0e-12)
  assert_eq(result.n, 4)
}

///|
test "Voyager global Geary's c is below one for a gradient" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_global_gearys_c(voy_test_gradient_values(), weights) catch {
    _ => abort("global Geary's c should run")
  }
  assert_true(result.estimate < 1.0)
}

///|
test "Voyager global Geary's c is above one for a checkerboard" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_global_gearys_c(
    voy_test_checkerboard_values(),
    weights,
  ) catch {
    _ => abort("global Geary's c should run")
  }
  assert_true(result.estimate > 1.0)
}

// ===========================================================================
// Local Moran's I (LISA)
// ===========================================================================

///|
test "Voyager local Moran's I returns one statistic per spot" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=99,
    seed=123,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  assert_eq(result.n, 9)
  assert_eq(result.local_i.length(), 9)
  assert_eq(result.quadrants.length(), 9)
  assert_eq(result.fdr.length(), 9)
  assert_eq(result.permutations, 99)
}

///|
test "Voyager local Moran's I permutation is deterministic for a fixed seed" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let first = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=50,
    seed=777,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  let second = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=50,
    seed=777,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  for index in 0..<first.n {
    voy_test_close(first.p_values[index], second.p_values[index], 0.0)
    voy_test_close(first.z_scores[index], second.z_scores[index], 0.0)
  }
}

///|
test "Voyager local Moran's I quadrant labels are valid" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=0,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  for label in result.quadrants {
    assert_true(
      label == "HH" ||
      label == "LL" ||
      label == "HL" ||
      label == "LH" ||
      label == "not significant",
    )
  }
}

///|
test "Voyager local Moran's I BH-FDR is monotonic non-decreasing in rank" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=99,
    seed=42,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  let indexed : Array[(Double, Int)] = []
  for index in 0..<result.n {
    indexed.push((result.p_values[index], index))
  }
  indexed.sort_by(fn(a : (Double, Int), b : (Double, Int)) -> Int {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  let mut previous = 0.0
  for entry in indexed {
    let adjusted = result.fdr[entry.1]
    assert_true(adjusted >= previous - 1.0e-12)
    previous = adjusted
  }
}

///|
test "Voyager local Moran's I rejects zero-variance input" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let failed = try {
    ignore(
      @src.voyager_local_morans_i(
        [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        weights,
      ),
    )
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Local Geary's c
// ===========================================================================

///|
test "Voyager local Geary's c returns per-spot statistics and classifications" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_gearys_c(
    voy_test_gradient_values(),
    weights,
    permutations=99,
    seed=9,
  ) catch {
    _ => abort("local Geary's c should run")
  }
  assert_eq(result.n, 9)
  assert_eq(result.local_c.length(), 9)
  assert_eq(result.classifications.length(), 9)
  for label in result.classifications {
    assert_true(
      label == "similar" || label == "dissimilar" || label == "not significant",
    )
  }
}

///|
test "Voyager local Geary's c is non-negative everywhere" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_gearys_c(voy_test_gradient_values(), weights) catch {
    _ => abort("local Geary's c should run")
  }
  for value in result.local_c {
    assert_true(value >= 0.0)
  }
}

// ===========================================================================
// Local Getis–Ord Gi / Gi*
// ===========================================================================

///|
test "Voyager local Getis-Ord Gi* flags the centre of a hotspot" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_getis_ord(
    voy_test_hotspot_values(),
    weights,
    star=true,
    fdr_threshold=0.2,
  ) catch {
    _ => abort("Getis–Ord should run")
  }
  assert_eq(result.n, 9)
  assert_eq(result.star, true)
  // Centre spot (index 4) is the hotspot.
  assert_true(result.z_scores[4] > 0.0)
  assert_true(result.classifications[4] == "hotspot")
}

///|
test "Voyager local Getis-Ord Gi excludes self from the weighted sum" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let gi = @src.voyager_local_getis_ord(
    voy_test_hotspot_values(),
    weights,
    star=false,
  ) catch {
    _ => abort("Getis–Ord should run")
  }
  let gi_star = @src.voyager_local_getis_ord(
    voy_test_hotspot_values(),
    weights,
    star=true,
  ) catch {
    _ => abort("Getis–Ord should run")
  }
  assert_eq(gi.star, false)
  // The two statistics differ for the hotspot centre because Gi* includes self.
  assert_true((gi.statistic[4] - gi_star.statistic[4]).abs() > 1.0e-9)
}

///|
test "Voyager local Getis-Ord rejects negative values" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let failed = try {
    ignore(
      @src.voyager_local_getis_ord(
        [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        weights,
      ),
    )
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager local Getis-Ord z-scores are finite" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_getis_ord(
    voy_test_gradient_values(),
    weights,
    star=true,
  ) catch {
    _ => abort("Getis–Ord should run")
  }
  for value in result.z_scores {
    assert_true(value == value)
    assert_true(value.abs() <= 1.0e300)
  }
}

// ===========================================================================
// Lee's L
// ===========================================================================

///|
test "Voyager Lee's L is positive for co-varying features" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_lees_l(
    voy_test_gradient_values(),
    voy_test_gradient_values(),
    weights,
  ) catch {
    _ => abort("Lee's L should run")
  }
  assert_eq(result.n, 9)
  // A feature with itself: Lee's L is positive and bounded by 1.
  assert_true(result.global_l > 0.0)
  assert_true(result.global_l <= 1.0 + 1.0e-9)
}

///|
test "Voyager Lee's L is negative for anti-correlated features" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let anti : Array[Double] = []
  for value in voy_test_gradient_values() {
    anti.push(10.0 - value)
  }
  let result = @src.voyager_lees_l(voy_test_gradient_values(), anti, weights) catch {
    _ => abort("Lee's L should run")
  }
  assert_true(result.global_l < 0.0)
}

///|
test "Voyager Lee's L local values length matches spot count" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_lees_l(
    voy_test_gradient_values(),
    voy_test_hotspot_values(),
    weights,
    permutations=49,
    seed=5,
  ) catch {
    _ => abort("Lee's L should run")
  }
  assert_eq(result.local_l.length(), 9)
  assert_eq(result.fdr.length(), 9)
}

///|
test "Voyager Lee's L rejects unequal feature lengths" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let failed = try {
    ignore(
      @src.voyager_lees_l([1.0, 2.0, 3.0], voy_test_gradient_values(), weights),
    )
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Multivariate local Geary
// ===========================================================================

///|
test "Voyager multivariate local Geary combines features into one statistic" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let matrix : Array[Array[Double]] = [
    voy_test_gradient_values(),
    voy_test_hotspot_values(),
  ]
  let result = @src.voyager_multivariate_local_geary(
    matrix,
    weights,
    permutations=49,
    seed=3,
  ) catch {
    _ => abort("multivariate local Geary should run")
  }
  assert_eq(result.n, 9)
  assert_eq(result.n_features, 2)
  assert_eq(result.local_statistic.length(), 9)
  assert_eq(result.classifications.length(), 9)
  assert_eq(result.feature_set, ["feature1", "feature2"])
}

///|
test "Voyager multivariate local Geary uses supplied feature names" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let matrix : Array[Array[Double]] = [
    voy_test_gradient_values(),
    voy_test_hotspot_values(),
  ]
  let result = @src.voyager_multivariate_local_geary(matrix, weights, feature_names=[
    "gradient", "hotspot",
  ]) catch {
    _ => abort("multivariate local Geary should run")
  }
  assert_eq(result.feature_set, ["gradient", "hotspot"])
}

///|
test "Voyager multivariate local Geary is non-negative" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let matrix : Array[Array[Double]] = [voy_test_gradient_values()]
  let result = @src.voyager_multivariate_local_geary(matrix, weights) catch {
    _ => abort("multivariate local Geary should run")
  }
  for value in result.local_statistic {
    assert_true(value >= 0.0)
  }
}

// ===========================================================================
// Empirical variogram and model fitting
// ===========================================================================

///|
test "Voyager empirical variogram bins pairs and computes semivariance" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=5,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  assert_true(empirical.length() >= 2)
  for point in empirical {
    assert_true(point.npairs > 0)
    assert_true(point.gamma >= 0.0)
    assert_true(point.lag > 0.0)
  }
}

///|
test "Voyager variogram fit recovers a positive range and bounded SSE" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=8,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  let model = @src.voyager_fit_variogram(empirical, model_type="spherical") catch {
    _ => abort("spherical fit should run")
  }
  assert_eq(model.model_type, "spherical")
  assert_true(model.range > 0.0)
  assert_true(model.sill >= 0.0)
  assert_true(model.nugget >= 0.0)
  assert_true(model.fitted_sse >= 0.0)
}

///|
test "Voyager variogram predict increases with distance up to the range" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=6,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  let model = @src.voyager_fit_variogram(empirical, model_type="exponential") catch {
    _ => abort("exponential fit should run")
  }
  let small = @src.voyager_variogram_predict(model, 0.01)
  let large = @src.voyager_variogram_predict(model, model.range * 3.0)
  assert_true(large >= small)
}

///|
test "Voyager variogram fit supports gaussian model" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=6,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  let model = @src.voyager_fit_variogram(empirical, model_type="gaussian") catch {
    _ => abort("gaussian fit should run")
  }
  assert_eq(model.model_type, "gaussian")
}

///|
test "Voyager variogram fit rejects an unknown model type" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=4,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  let failed = try {
    ignore(@src.voyager_fit_variogram(empirical, model_type="cubic"))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager variogram fit respects a fixed nugget" {
  let coords = voy_test_grid_coords()
  let empirical = @src.voyager_empirical_variogram(
    coords,
    voy_test_gradient_values(),
    n_lags=6,
  ) catch {
    _ => abort("empirical variogram should run")
  }
  let model = @src.voyager_fit_variogram(
    empirical,
    model_type="spherical",
    nugget=0.25,
    fix_nugget=true,
  ) catch {
    _ => abort("fixed-nugget fit should run")
  }
  voy_test_close(model.nugget, 0.25, 1.0e-12)
}

// ===========================================================================
// Correlogram
// ===========================================================================

///|
test "Voyager correlogram returns one point per non-empty lag" {
  let coords = voy_test_grid_coords()
  let points = @src.voyager_correlogram(
    coords,
    voy_test_gradient_values(),
    n_lags=4,
  ) catch {
    _ => abort("correlogram should run")
  }
  // Empty distance bins are skipped; every returned point has observed pairs.
  assert_true(points.length() >= 1)
  assert_true(points.length() <= 4)
  let mut previous_lag = -1.0
  for point in points {
    assert_true(point.npairs > 0)
    assert_true(point.lag > previous_lag)
    previous_lag = point.lag
    voy_test_close(point.expectation, -1.0 / 8.0, 1.0e-12)
  }
}

///|
test "Voyager correlogram gradient shows positive Moran at short lags" {
  let coords = voy_test_grid_coords()
  let points = @src.voyager_correlogram(
    coords,
    voy_test_gradient_values(),
    n_lags=3,
  ) catch {
    _ => abort("correlogram should run")
  }
  assert_true(points[0].morans_i > 0.0)
}

// ===========================================================================
// SpatialExperiment integration
// ===========================================================================

///|
test "Voyager example SpatialExperiment has the expected shape" {
  let se = @src.voyager_example_spatial_experiment()
  assert_eq(se.col_data.length(), 36)
  assert_eq(se.spatial_coords.length(), 36)
  let assay = se.assay["logcounts"]
  assert_eq(assay.length(), 4)
  assert_eq(assay[0].length(), 36)
}

///|
test "Voyager univariate SFE Moran writes back local results without mutating input" {
  let se = @src.voyager_example_spatial_experiment()
  let original_keys = se.col_data[0].keys().length()
  let output = @src.voyager_run_univariate_sfe(
    se,
    [0],
    assay_name="logcounts",
    stat_method="moran",
    permutations=19,
    seed=2024,
    output_prefix="voy",
  ) catch {
    _ => abort("Voyager SFE Moran integration should run")
  }
  // Input SpatialExperiment unchanged.
  assert_eq(se.col_data[0].keys().length(), original_keys)
  assert_true(!se.col_data[0].contains("voy.moran.local.gene_gradient"))
  // Output SpatialExperiment carries the write-back.
  assert_true(
    output.experiment.col_data[0].contains("voy.moran.local.gene_gradient"),
  )
  assert_true(
    output.experiment.row_data[0].contains("voy.moran.I.gene_gradient"),
  )
  assert_eq(output.results.length(), 1)
  assert_eq(output.results[0].stat_method, "moran")
  assert_eq(output.weights.n, 36)
}

///|
test "Voyager univariate SFE Geary writes classification column" {
  let se = @src.voyager_example_spatial_experiment()
  let output = @src.voyager_run_univariate_sfe(
    se,
    [0, 1],
    assay_name="logcounts",
    stat_method="geary",
    permutations=19,
    seed=7,
    output_prefix="vg",
  ) catch {
    _ => abort("Voyager SFE Geary integration should run")
  }
  assert_eq(output.results.length(), 2)
  assert_true(
    output.experiment.col_data[0].contains("vg.geary.class.gene_gradient"),
  )
  assert_true(output.experiment.row_data[1].contains("vg.geary.C.gene_hotspot"))
}

///|
test "Voyager univariate SFE Getis writes hotspot classification" {
  let se = @src.voyager_example_spatial_experiment()
  let output = @src.voyager_run_univariate_sfe(
    se,
    [1],
    assay_name="logcounts",
    stat_method="getis",
    output_prefix="g",
  ) catch {
    _ => abort("Voyager SFE Getis integration should run")
  }
  assert_eq(output.results.length(), 1)
  assert_eq(output.results[0].stat_method, "getis")
  assert_true(
    output.experiment.col_data[0].contains("g.getis.class.gene_hotspot"),
  )
  assert_true(output.experiment.col_data[0].contains("g.getis.z.gene_hotspot"))
}

///|
test "Voyager univariate SFE can use a distance-band bandwidth" {
  let se = @src.voyager_example_spatial_experiment()
  let output = @src.voyager_run_univariate_sfe(
    se,
    [0],
    assay_name="logcounts",
    stat_method="moran",
    bandwidth=1.0,
    output_prefix="db",
  ) catch {
    _ => abort("Voyager SFE distance-band integration should run")
  }
  assert_eq(output.weights.n, 36)
  // Distance band of 1.0 on a unit grid links only orthogonal immediate
  // neighbours (excluding diagonals at sqrt(2)), so centre spots have 4.
  let mut has_four = false
  for index in 0..<output.weights.n {
    if output.weights.neighbors[index].length() == 4 {
      has_four = true
    }
  }
  assert_true(has_four)
}

///|
test "Voyager univariate SFE rejects an out-of-range feature index" {
  let se = @src.voyager_example_spatial_experiment()
  let failed = try {
    ignore(@src.voyager_run_univariate_sfe(se, [99], assay_name="logcounts"))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager univariate SFE rejects an unknown method" {
  let se = @src.voyager_example_spatial_experiment()
  let failed = try {
    ignore(
      @src.voyager_run_univariate_sfe(
        se,
        [0],
        assay_name="logcounts",
        stat_method="unknown",
      ),
    )
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager univariate SFE rejects a missing assay" {
  let se = @src.voyager_example_spatial_experiment()
  let failed = try {
    ignore(@src.voyager_run_univariate_sfe(se, [0], assay_name="counts"))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager univariate SFE metadata records method and feature count" {
  let se = @src.voyager_example_spatial_experiment()
  let output = @src.voyager_run_univariate_sfe(
    se,
    [0, 1, 2],
    assay_name="logcounts",
    stat_method="moran",
    output_prefix="meta",
  ) catch {
    _ => abort("Voyager SFE integration should run")
  }
  assert_eq(output.experiment.metadata["meta.method"], "moran")
  assert_eq(output.experiment.metadata["meta.n_features"], "3")
  assert_eq(output.experiment.metadata["meta.weights_style"], "W")
}

// ===========================================================================
// Determinism and edge cases
// ===========================================================================

///|
test "Voyager PRNG produces reproducible shuffles for the same seed" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let a = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=30,
    seed=20240501,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  let b = @src.voyager_local_morans_i(
    voy_test_gradient_values(),
    weights,
    permutations=30,
    seed=20240501,
  ) catch {
    _ => abort("local Moran's I should run")
  }
  voy_test_close(a.p_values[0], b.p_values[0], 0.0)
  voy_test_close(a.z_scores[4], b.z_scores[4], 0.0)
}

///|
test "Voyager global Moran's I expectation equals negative one over n minus one" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_global_morans_i(voy_test_gradient_values(), weights) catch {
    _ => abort("global Moran's I should run")
  }
  voy_test_close(result.expectation, -1.0 / 8.0, 1.0e-12)
}

///|
test "Voyager weights distance band rejects non-positive bandwidth" {
  let coords = voy_test_grid_coords()
  let failed = try {
    ignore(@src.voyager_weights_distance_band(coords, 0.0))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager weights inverse distance rejects non-positive power" {
  let coords = voy_test_grid_coords()
  let failed = try {
    ignore(@src.voyager_weights_inverse_distance(coords, 0.0))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager empirical variogram rejects coincident coordinates" {
  let coords : Array[Array[Double]] = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
  let failed = try {
    ignore(@src.voyager_empirical_variogram(coords, [1.0, 2.0, 3.0]))
    false
  } catch {
    VoyagerError(_) => true
  }
  assert_true(failed)
}

///|
test "Voyager local Getis-Ord FDR is bounded in [0, 1]" {
  let coords = voy_test_grid_coords()
  let weights = @src.voyager_weights_knn(coords, 4) catch {
    _ => abort("kNN weights should build")
  }
  let result = @src.voyager_local_getis_ord(
    voy_test_hotspot_values(),
    weights,
    star=true,
  ) catch {
    _ => abort("Getis–Ord should run")
  }
  for value in result.fdr {
    assert_true(value >= 0.0 && value <= 1.0)
  }
}
