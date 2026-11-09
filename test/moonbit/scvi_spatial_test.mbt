///|
test "scvi_distribution_sampling_shapes" {
  let nb = @src.ScviNegativeBinomial::new([2.0, 3.0, 4.0], [1.0, 1.0, 1.0])
  let sample = nb.sample(11)
  assert_eq(sample.length(), 3)
  for i in 0..<sample.length() {
    assert_true(sample[i] >= 0.0)
  }
  let pois = @src.ScviPoisson::new([1.0, 2.0])
  assert_eq(pois.sample(3).length(), 2)
  let bern = @src.ScviBernoulli::new([0.5, 0.5])
  assert_eq(bern.sample(4).length(), 2)
  let zinb = @src.ScviZeroInflatedNegativeBinomial::new(
    [3.0, 3.0],
    [1.0, 1.0],
    [0.0, 0.0],
  )
  assert_eq(zinb.mean().length(), 2)
  assert_eq(zinb.sample(5).length(), 2)
}

///|
test "scvi_stereoscope_fit_and_deconvolve" {
  let reference = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
  let labels = [0, 1, 2]
  let config = @src.ScviStereoscopeConfig::default()
  let model = @src.scvi_stereoscope_fit(reference, labels, 3, config)
  assert_eq(model.n_cell_types(), 3)
  assert_eq(model.n_genes(), 3)
  let proportions = @src.scvi_stereoscope_deconvolve(
    model,
    [10.0, 1.0, 0.0],
    config,
  )
  assert_eq(proportions.length(), 3)
  let mut total = 0.0
  for i in 0..<proportions.length() {
    total += proportions[i]
    assert_true(proportions[i] >= 0.0)
  }
  assert_true((total - 1.0).abs() < 1.0e-6)
  assert_true(proportions[0] > proportions[1])
  assert_true(proportions[0] > proportions[2])
  let spatial = @src.scvi_stereoscope_spatial_deconvolve(
    model,
    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0]],
    config,
  )
  assert_eq(spatial.length(), 2)
  let dominant = @src.scvi_stereoscope_dominant_type(spatial)
  assert_eq(dominant[0], 0)
  assert_eq(dominant[1], 1)
}

///|
test "scvi_stereoscope_aliases" {
  let reference = [[5.0, 1.0], [1.0, 5.0]]
  let labels = [0, 1]
  let config = @src.ScviStereoscopeConfig::default()
  let model = @src.scvi_rnastereoscope_fit(reference, labels, 2, config)
  let w = @src.scvi_rnastereoscope_deconvolve(model, [5.0, 1.0], config)
  assert_eq(w.length(), 2)
  let spatial_model = @src.scvi_spatialstereoscope_fit(
    reference, labels, 2, config,
  )
  let spots = @src.scvi_spatialstereoscope_deconvolve(
    spatial_model,
    [[5.0, 1.0]],
    config,
  )
  assert_eq(spots.length(), 1)
}

///|
test "scvi_spatial_knn_weights" {
  let coords = [[0.0, 0.0], [0.0, 1.0], [0.0, 3.0], [9.0, 9.0]]
  let weights = @src.scvi_spatial_knn_weights(coords, 1)
  assert_eq(weights.length(), 4)
  for i in 0..<weights.length() {
    let mut row_sum = 0.0
    for j in 0..<weights[i].length() {
      row_sum += weights[i][j]
    }
    assert_true((row_sum - 1.0).abs() < 1.0e-9)
    assert_true(weights[i][i] == 0.0)
  }
}

///|
test "scvi_spatial_autocorrelation_detects_clusters" {
  let coords = [
    [0.0, 0.0],
    [0.0, 1.0],
    [0.0, 2.0],
    [0.0, 3.0],
    [10.0, 0.0],
    [10.0, 1.0],
    [10.0, 2.0],
    [10.0, 3.0],
  ]
  let values = [1.0, 2.0, 3.0, 4.0, 100.0, 200.0, 300.0, 400.0]
  let (moran, geary) = @src.scvi_spatial_autocorrelation(values, coords, 2)
  assert_true(moran > 0.0)
  assert_true(geary < 1.0)
  let weights = @src.scvi_spatial_knn_weights(coords, 2)
  assert_true((@src.scvi_morans_i(values, weights) - moran).abs() < 1.0e-12)
  assert_true((@src.scvi_gearys_c(values, weights) - geary).abs() < 1.0e-12)
}

///|
test "scvi_spatial_permutation_test_range" {
  let coords = [
    [0.0, 0.0],
    [0.0, 1.0],
    [0.0, 2.0],
    [10.0, 0.0],
    [10.0, 1.0],
    [10.0, 2.0],
  ]
  let values = [1.0, 2.0, 3.0, 50.0, 60.0, 70.0]
  let p = @src.scvi_spatial_permutation_test(values, coords, 1, 99, 42)
  assert_true(p > 0.0 && p <= 1.0)
}

///|
test "scvi_spatial_neighbor_enrichment" {
  let coords = [
    [0.0, 0.0],
    [0.0, 1.0],
    [0.0, 2.0],
    [10.0, 0.0],
    [10.0, 1.0],
    [10.0, 2.0],
  ]
  let labels = [0, 0, 0, 1, 1, 1]
  let same = @src.scvi_spatial_neighbor_enrichment(labels, coords, 1)
  assert_true(same > 0.5)
  assert_true(
    (@src.scvi_spatial_label_fraction(labels, 1) - 0.5).abs() < 1.0e-12,
  )
}
