// Black-box tests for the Bioconductor SpatialDecon-inspired workflow.

///|
fn sd_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "SpatialDecon value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn sd_test_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn sd_test_config(
  rescale_profile? : Bool = false,
  refit_outliers? : Bool = false,
  residual_threshold? : Double = 3.0,
) -> @src.SpatialDeconConfig {
  @src.SpatialDeconConfig::create(
    residual_threshold~,
    rescale_profile~,
    refit_outliers~,
    tolerance=1.0e-9,
  ) catch {
    _ => abort("SpatialDecon test configuration should be valid")
  }
}

///|
fn sd_test_exact_profile() -> @src.SpatialDeconProfile {
  @src.SpatialDeconProfile::create(["g1", "g2", "g3"], ["Type_A"], [
    [1.0],
    [2.0],
    [3.0],
  ]) catch {
    _ => abort("SpatialDecon exact profile should be valid")
  }
}

///|
fn sd_test_exact_data() -> @src.SpatialDeconData {
  @src.SpatialDeconData::create(
    ["g1", "g2", "g3"],
    ["spot_1", "spot_2"],
    [[7.0, 9.0], [9.0, 13.0], [11.0, 17.0]],
    background=[[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]],
  ) catch {
    _ => abort("SpatialDecon exact data should be valid")
  }
}

///|
fn sd_test_exact_result(
  nuclei_counts? : Array[Double] = [],
) -> @src.SpatialDeconResult {
  @src.spatial_decon(
    sd_test_exact_data(),
    sd_test_exact_profile(),
    nuclei_counts~,
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon exact fit should succeed")
  }
}

///|
fn sd_test_example_result() -> @src.SpatialDeconResult {
  let (data, profile, nuclei_counts) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon example data should be valid")
  }
  @src.spatial_decon(
    data,
    profile,
    nuclei_counts~,
    config=sd_test_config(refit_outliers=true),
  ) catch {
    _ => abort("SpatialDecon example fit should succeed")
  }
}

///|
fn sd_test_experiment() -> (
  @src.SpatialExperiment,
  @src.SpatialDeconProfile,
  Array[Double],
) {
  let (data, profile, nuclei_counts) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon example data should be valid")
  }
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "counts", data.values))
  ignore(@src.se_add_assay(experiment, "background", data.background))
  ignore(@src.se_add_assay(experiment, "weights", data.weights))
  for gene in data.gene_names {
    ignore(
      @src.se_add_row(
        experiment,
        Map([("gene_id", gene), ("feature_type", "target")]),
      ),
    )
  }
  for spot in 0..<data.spot_ids.length() {
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", data.spot_ids[spot]), ("sample", "sample_1")]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(spot.to_double(), 1.0),
      ),
    )
  }
  experiment.metadata["owner"] = "input"
  (experiment, profile, nuclei_counts)
}

///|
test "SpatialDecon: default configuration follows upstream controls" {
  let config = @src.SpatialDeconConfig::default()
  assert_eq(config.residual_threshold, 3.0)
  assert_eq(config.lower_threshold, 0.5)
  assert_eq(config.signal_floor, 1.0e-4)
  assert_eq(config.rescale_profile, true)
  assert_eq(config.profile_quantile, 0.99)
  assert_eq(config.profile_target, 2.0)
  assert_eq(config.refit_outliers, true)
  assert_eq(config.max_iterations, 1000)
  assert_eq(config.line_search_steps, 30)
  assert_eq(config.tolerance, 1.0e-8)
  assert_eq(config.ridge, 1.0e-8)
}

///|
test "SpatialDecon: custom configuration preserves optimization controls" {
  let config = @src.SpatialDeconConfig::create(
    residual_threshold=2.0,
    lower_threshold=0.25,
    signal_floor=1.0e-5,
    rescale_profile=false,
    profile_quantile=0.9,
    profile_target=3.0,
    refit_outliers=false,
    max_iterations=200,
    line_search_steps=12,
    tolerance=1.0e-6,
    ridge=1.0e-6,
  ) catch {
    _ => abort("SpatialDecon custom configuration should be valid")
  }
  assert_eq(config.residual_threshold, 2.0)
  assert_eq(config.lower_threshold, 0.25)
  assert_eq(config.rescale_profile, false)
  assert_eq(config.profile_target, 3.0)
  assert_eq(config.refit_outliers, false)
  assert_eq(config.max_iterations, 200)
}

///|
test "SpatialDecon: configuration rejects invalid thresholds" {
  let residual_failed = try {
    ignore(@src.SpatialDeconConfig::create(residual_threshold=0.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let lower_failed = try {
    ignore(@src.SpatialDeconConfig::create(lower_threshold=0.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let floor_failed = try {
    ignore(@src.SpatialDeconConfig::create(signal_floor=0.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(residual_failed && lower_failed && floor_failed)
}

///|
test "SpatialDecon: configuration rejects invalid profile scaling" {
  let quantile_failed = try {
    ignore(@src.SpatialDeconConfig::create(profile_quantile=1.1))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let target_failed = try {
    ignore(@src.SpatialDeconConfig::create(profile_target=0.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(quantile_failed && target_failed)
}

///|
test "SpatialDecon: configuration rejects invalid solver controls" {
  let iterations_failed = try {
    ignore(@src.SpatialDeconConfig::create(max_iterations=0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let search_failed = try {
    ignore(@src.SpatialDeconConfig::create(line_search_steps=0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let tolerance_failed = try {
    ignore(@src.SpatialDeconConfig::create(tolerance=1.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let ridge_failed = try {
    ignore(@src.SpatialDeconConfig::create(ridge=0.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(
    iterations_failed && search_failed && tolerance_failed && ridge_failed,
  )
}

///|
test "SpatialDecon: profile constructor exposes dimensions" {
  let profile = sd_test_exact_profile()
  assert_eq(profile.n_genes(), 3)
  assert_eq(profile.n_cell_types(), 1)
  assert_eq(profile.gene_names, ["g1", "g2", "g3"])
  assert_eq(profile.cell_types, ["Type_A"])
}

///|
test "SpatialDecon: profile constructor defensively copies inputs" {
  let genes = ["g1", "g2"]
  let cell_types = ["A"]
  let values = [[1.0], [2.0]]
  let profile = @src.SpatialDeconProfile::create(genes, cell_types, values) catch {
    _ => abort("SpatialDecon profile should be valid")
  }
  genes[0] = "changed"
  cell_types[0] = "changed"
  values[0][0] = 99.0
  assert_eq(profile.gene_names, ["g1", "g2"])
  assert_eq(profile.cell_types, ["A"])
  assert_eq(profile.values, [[1.0], [2.0]])
}

///|
test "SpatialDecon: profile rejects duplicate and empty names" {
  let duplicate_failed = try {
    ignore(
      @src.SpatialDeconProfile::create(["g1", "g1"], ["A"], [[1.0], [2.0]]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let empty_failed = try {
    ignore(@src.SpatialDeconProfile::create(["g1"], [""], [[1.0]]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(duplicate_failed && empty_failed)
}

///|
test "SpatialDecon: profile rejects malformed and negative matrices" {
  let shape_failed = try {
    ignore(
      @src.SpatialDeconProfile::create(["g1", "g2"], ["A", "B"], [
        [1.0, 2.0],
        [3.0],
      ]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let negative_failed = try {
    ignore(@src.SpatialDeconProfile::create(["g1"], ["A"], [[-1.0]]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(shape_failed && negative_failed)
}

///|
test "SpatialDecon: profile rejects non-finite and empty cell-type signals" {
  let finite_failed = try {
    ignore(
      @src.SpatialDeconProfile::create(["g1"], ["A"], [[@double.not_a_number]]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let signal_failed = try {
    ignore(
      @src.SpatialDeconProfile::create(["g1", "g2"], ["A", "B"], [
        [1.0, 0.0],
        [2.0, 0.0],
      ]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(finite_failed && signal_failed)
}

///|
test "SpatialDecon: data constructor supplies background and weights" {
  let data = @src.SpatialDeconData::create(["g1", "g2"], ["s1"], [[2.0], [3.0]]) catch {
    _ => abort("SpatialDecon data should be valid")
  }
  assert_eq(data.n_genes(), 2)
  assert_eq(data.n_spots(), 1)
  assert_eq(data.background, [[0.0], [0.0]])
  assert_eq(data.weights, [[1.0], [1.0]])
}

///|
test "SpatialDecon: data constructor defensively copies all matrices" {
  let genes = ["g1", "g2"]
  let spots = ["s1"]
  let values = [[2.0], [3.0]]
  let background = [[0.5], [0.6]]
  let weights = [[1.0], [2.0]]
  let data = @src.SpatialDeconData::create(
    genes,
    spots,
    values,
    background~,
    weights~,
  ) catch {
    _ => abort("SpatialDecon data should be valid")
  }
  genes[0] = "changed"
  spots[0] = "changed"
  values[0][0] = 99.0
  background[0][0] = 99.0
  weights[0][0] = 99.0
  assert_eq(data.gene_names, ["g1", "g2"])
  assert_eq(data.spot_ids, ["s1"])
  assert_eq(data.values[0][0], 2.0)
  assert_eq(data.background[0][0], 0.5)
  assert_eq(data.weights[0][0], 1.0)
}

///|
test "SpatialDecon: scalar background fills every observation" {
  let data = @src.SpatialDeconData::with_scalar_background(
    ["g1", "g2"],
    ["s1", "s2"],
    [[2.0, 3.0], [4.0, 5.0]],
    0.75,
  ) catch {
    _ => abort("SpatialDecon scalar background should be valid")
  }
  assert_eq(data.background, [[0.75, 0.75], [0.75, 0.75]])
}

///|
test "SpatialDecon: data rejects duplicate spots and invalid shapes" {
  let duplicate_failed = try {
    ignore(@src.SpatialDeconData::create(["g1"], ["s1", "s1"], [[1.0, 2.0]]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let shape_failed = try {
    ignore(
      @src.SpatialDeconData::create(["g1", "g2"], ["s1"], [[1.0], [2.0, 3.0]]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(duplicate_failed && shape_failed)
}

///|
test "SpatialDecon: data rejects invalid expression and background" {
  let expression_failed = try {
    ignore(@src.SpatialDeconData::create(["g1"], ["s1"], [[-1.0]]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let background_failed = try {
    ignore(
      @src.SpatialDeconData::create(["g1"], ["s1"], [[1.0]], background=[
        [@double.not_a_number],
      ]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let scalar_failed = try {
    ignore(
      @src.SpatialDeconData::with_scalar_background(
        ["g1"],
        ["s1"],
        [[1.0]],
        -1.0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(expression_failed && background_failed && scalar_failed)
}

///|
test "SpatialDecon: data requires strictly positive finite weights" {
  let zero_failed = try {
    ignore(
      @src.SpatialDeconData::create(["g1"], ["s1"], [[1.0]], weights=[[0.0]]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let finite_failed = try {
    ignore(
      @src.SpatialDeconData::create(["g1"], ["s1"], [[1.0]], weights=[
        [@double.not_a_number],
      ]),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(zero_failed && finite_failed)
}

///|
test "SpatialDecon: exact single-cell-type mixture is recovered" {
  let result = sd_test_exact_result()
  sd_test_close(result.beta[0][0], 2.0, 1.0e-6)
  sd_test_close(result.beta[0][1], 4.0, 1.0e-6)
  assert_eq(result.converged, [true, true])
  assert_true(result.objectives[0] < 1.0e-12)
  assert_true(result.objectives[1] < 1.0e-12)
}

///|
test "SpatialDecon: background-aware fit differs from uncorrected fit" {
  let corrected = sd_test_exact_result()
  let uncorrected_data = @src.SpatialDeconData::create(
    ["g1", "g2", "g3"],
    ["spot_1", "spot_2"],
    [[7.0, 9.0], [9.0, 13.0], [11.0, 17.0]],
  ) catch {
    _ => abort("SpatialDecon uncorrected data should be valid")
  }
  let uncorrected = @src.spatial_decon(
    uncorrected_data,
    sd_test_exact_profile(),
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon uncorrected fit should succeed")
  }
  sd_test_close(corrected.beta[0][0], 2.0, 1.0e-6)
  assert_true(uncorrected.beta[0][0] > 3.0)
}

///|
test "SpatialDecon: profile rescaling inversely rescales abundance" {
  let profile = @src.SpatialDeconProfile::create(["g1", "g2"], ["A"], [
    [2.0],
    [4.0],
  ]) catch {
    _ => abort("SpatialDecon scaling profile should be valid")
  }
  let data = @src.SpatialDeconData::create(["g1", "g2"], ["s1"], [[6.0], [12.0]]) catch {
    _ => abort("SpatialDecon scaling data should be valid")
  }
  let plain = @src.spatial_decon(data, profile, config=sd_test_config()) catch {
    _ => abort("SpatialDecon unscaled fit should succeed")
  }
  let scaled_config = @src.SpatialDeconConfig::create(
    rescale_profile=true,
    profile_quantile=1.0,
    profile_target=2.0,
    refit_outliers=false,
    tolerance=1.0e-9,
  ) catch {
    _ => abort("SpatialDecon scaling configuration should be valid")
  }
  let scaled = @src.spatial_decon(data, profile, config=scaled_config) catch {
    _ => abort("SpatialDecon scaled fit should succeed")
  }
  sd_test_close(plain.beta[0][0], 3.0, 1.0e-6)
  sd_test_close(scaled.beta[0][0], 6.0, 1.0e-6)
}

///|
test "SpatialDecon: zero profile quantile is rejected during rescaling" {
  let profile = @src.SpatialDeconProfile::create(["g1", "g2", "g3"], ["A"], [
    [0.0],
    [0.0],
    [1.0],
  ]) catch {
    _ => abort("SpatialDecon sparse profile should be valid")
  }
  let data = @src.SpatialDeconData::create(["g1", "g2", "g3"], ["s1"], [
    [1.0],
    [1.0],
    [1.0],
  ]) catch {
    _ => abort("SpatialDecon sparse data should be valid")
  }
  let config = @src.SpatialDeconConfig::create(profile_quantile=0.5) catch {
    _ => abort("SpatialDecon sparse configuration should be valid")
  }
  let failed = try {
    ignore(@src.spatial_decon(data, profile, config~))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(failed)
}

///|
test "SpatialDecon: precision weights favor high-confidence genes" {
  let profile = @src.SpatialDeconProfile::create(["g1", "g2"], ["A"], [
    [1.0],
    [1.0],
  ]) catch {
    _ => abort("SpatialDecon weighted profile should be valid")
  }
  let plain_data = @src.SpatialDeconData::create(["g1", "g2"], ["s1"], [
    [2.0],
    [8.0],
  ]) catch {
    _ => abort("SpatialDecon unweighted data should be valid")
  }
  let weighted_data = @src.SpatialDeconData::create(
    ["g1", "g2"],
    ["s1"],
    [[2.0], [8.0]],
    weights=[[100.0], [1.0]],
  ) catch {
    _ => abort("SpatialDecon weighted data should be valid")
  }
  let plain = @src.spatial_decon(plain_data, profile, config=sd_test_config()) catch {
    _ => abort("SpatialDecon unweighted fit should succeed")
  }
  let weighted = @src.spatial_decon(
    weighted_data,
    profile,
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon weighted fit should succeed")
  }
  sd_test_close(plain.beta[0][0], 4.0, 1.0e-5)
  assert_true(weighted.beta[0][0] < 2.1)
  assert_true(weighted.beta[0][0] > 2.0)
}

///|
test "SpatialDecon: two-stage fitting flags and removes a log outlier" {
  let profile = @src.SpatialDeconProfile::create(
    ["g1", "g2", "g3", "g4", "g5"],
    ["A"],
    [[1.0], [1.0], [1.0], [1.0], [1.0]],
  ) catch {
    _ => abort("SpatialDecon outlier profile should be valid")
  }
  let data = @src.SpatialDeconData::create(
    ["g1", "g2", "g3", "g4", "g5"],
    ["s1"],
    [[2.0], [2.0], [2.0], [2.0], [1024.0]],
  ) catch {
    _ => abort("SpatialDecon outlier data should be valid")
  }
  let result = @src.spatial_decon(
    data,
    profile,
    config=sd_test_config(refit_outliers=true),
  ) catch {
    _ => abort("SpatialDecon outlier fit should succeed")
  }
  sd_test_close(result.beta[0][0], 2.0, 1.0e-6)
  assert_eq(result.outliers[4][0], true)
  assert_eq(result.outliers[0][0], false)
  assert_true(result.residuals[4][0] != result.residuals[4][0])
}

///|
test "SpatialDecon: disabling refit preserves every observation" {
  let profile = @src.SpatialDeconProfile::create(
    ["g1", "g2", "g3", "g4", "g5"],
    ["A"],
    [[1.0], [1.0], [1.0], [1.0], [1.0]],
  ) catch {
    _ => abort("SpatialDecon outlier profile should be valid")
  }
  let data = @src.SpatialDeconData::create(
    ["g1", "g2", "g3", "g4", "g5"],
    ["s1"],
    [[2.0], [2.0], [2.0], [2.0], [1024.0]],
  ) catch {
    _ => abort("SpatialDecon outlier data should be valid")
  }
  let result = @src.spatial_decon(data, profile, config=sd_test_config()) catch {
    _ => abort("SpatialDecon non-refit fit should succeed")
  }
  for row in result.outliers {
    assert_eq(row[0], false)
  }
  assert_true(result.beta[0][0] > 2.0)
}

///|
test "SpatialDecon: synthetic multi-cell mixtures recover known abundances" {
  let result = sd_test_example_result()
  let expected = [
    [4.0, 1.0, 2.0, 0.5, 3.0],
    [1.0, 4.0, 2.0, 1.0, 0.5],
    [0.5, 1.0, 2.0, 4.0, 3.0],
  ]
  for cell_type in 0..<expected.length() {
    for spot in 0..<expected[cell_type].length() {
      sd_test_close(
        result.beta[cell_type][spot],
        expected[cell_type][spot],
        0.04,
      )
    }
  }
}

///|
test "SpatialDecon: fit diagnostics and Hessian uncertainty are finite" {
  let result = sd_test_example_result()
  for spot in 0..<result.n_spots() {
    assert_true(result.converged[spot])
    assert_true(result.objectives[spot] >= 0.0)
    assert_true(result.rmse[spot] < 0.03)
    assert_true(result.correlations[spot] > 0.99)
    for cell_type in 0..<result.n_cell_types() {
      assert_true(result.standard_errors[cell_type][spot] > 0.0)
      assert_true(sd_test_finite(result.t_statistics[cell_type][spot]))
      assert_true(result.p_values[cell_type][spot] >= 0.0)
      assert_true(result.p_values[cell_type][spot] <= 1.0)
      assert_true(result.covariances[spot][cell_type][cell_type] > 0.0)
    }
  }
}

///|
test "SpatialDecon: proportions sum to one per spot" {
  let result = sd_test_example_result()
  for spot in 0..<result.n_spots() {
    let mut total = 0.0
    for cell_type in 0..<result.n_cell_types() {
      total = total + result.proportions[cell_type][spot]
    }
    sd_test_close(total, 1.0, 1.0e-12)
  }
}

///|
test "SpatialDecon: cells-per-100 uses the maximum spot total" {
  let result = sd_test_example_result()
  let mut maximum = 0.0
  for spot in 0..<result.n_spots() {
    let mut total = 0.0
    for cell_type in 0..<result.n_cell_types() {
      total = total + result.cells_per_100[cell_type][spot]
    }
    if total > maximum {
      maximum = total
    }
  }
  sd_test_close(maximum, 100.0, 1.0e-10)
}

///|
test "SpatialDecon: nuclei counts scale cells-per-100 by spot" {
  let result = sd_test_exact_result(nuclei_counts=[100.0, 80.0])
  sd_test_close(result.cells_per_100[0][0], 50.0, 1.0e-5)
  sd_test_close(result.cells_per_100[0][1], 100.0, 1.0e-5)
  sd_test_close(result.cell_counts[0][0], 50.0, 1.0e-5)
  sd_test_close(result.cell_counts[0][1], 80.0, 1.0e-5)
}

///|
test "SpatialDecon: omitted nuclei counts leave cell counts empty" {
  let result = sd_test_exact_result()
  assert_eq(result.cell_counts.length(), 0)
}

///|
test "SpatialDecon: nuclei counts are strictly validated" {
  let length_failed = try {
    ignore(
      @src.spatial_decon(
        sd_test_exact_data(),
        sd_test_exact_profile(),
        nuclei_counts=[1.0],
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let value_failed = try {
    ignore(
      @src.spatial_decon(
        sd_test_exact_data(),
        sd_test_exact_profile(),
        nuclei_counts=[1.0, -1.0],
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(length_failed && value_failed)
}

///|
test "SpatialDecon: data and profile align by gene name and data order" {
  let profile = @src.SpatialDeconProfile::create(["extra", "g2", "g1"], ["A"], [
    [5.0],
    [2.0],
    [1.0],
  ]) catch {
    _ => abort("SpatialDecon reordered profile should be valid")
  }
  let data = @src.SpatialDeconData::create(["g1", "g2"], ["s1"], [[3.0], [6.0]]) catch {
    _ => abort("SpatialDecon reordered data should be valid")
  }
  let result = @src.spatial_decon(data, profile, config=sd_test_config()) catch {
    _ => abort("SpatialDecon reordered fit should succeed")
  }
  assert_eq(result.gene_names, ["g1", "g2"])
  sd_test_close(result.beta[0][0], 3.0, 1.0e-6)
}

///|
test "SpatialDecon: fitting rejects insufficient shared genes" {
  let no_shared_profile = @src.SpatialDeconProfile::create(["other"], ["A"], [
    [1.0],
  ]) catch {
    _ => abort("SpatialDecon disjoint profile should be valid")
  }
  let no_shared_failed = try {
    ignore(
      @src.spatial_decon(
        sd_test_exact_data(),
        no_shared_profile,
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let underdetermined_profile = @src.SpatialDeconProfile::create(
    ["g1", "other"],
    ["A", "B"],
    [[1.0, 1.0], [1.0, 2.0]],
  ) catch {
    _ => abort("SpatialDecon underdetermined profile should be valid")
  }
  let underdetermined_failed = try {
    ignore(
      @src.spatial_decon(
        sd_test_exact_data(),
        underdetermined_profile,
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(no_shared_failed && underdetermined_failed)
}

///|
test "SpatialDecon: fitting rejects all-zero expression" {
  let data = @src.SpatialDeconData::create(["g1", "g2", "g3"], ["s1"], [
    [0.0],
    [0.0],
    [0.0],
  ]) catch {
    _ => abort("SpatialDecon zero data constructor should succeed")
  }
  let failed = try {
    ignore(
      @src.spatial_decon(data, sd_test_exact_profile(), config=sd_test_config()),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(failed)
}

///|
test "SpatialDecon: abundance query and ranking expose fitted values" {
  let result = sd_test_example_result()
  match result.abundance_for("T_cell", "region_1") {
    Some(value) => sd_test_close(value, result.beta[0][0], 0.0)
    None => abort("SpatialDecon known abundance should exist")
  }
  assert_eq(result.abundance_for("missing", "region_1"), None)
  let ranked = result.top_cell_types("region_1", limit=2) catch {
    _ => abort("SpatialDecon ranking should succeed")
  }
  assert_eq(ranked.length(), 2)
  assert_eq(ranked[0].cell_type, "T_cell")
  assert_true(ranked[0].abundance >= ranked[1].abundance)
}

///|
test "SpatialDecon: ranking validates spot and limit" {
  let result = sd_test_example_result()
  let spot_failed = try {
    ignore(result.top_cell_types("missing"))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let limit_failed = try {
    ignore(result.top_cell_types("region_1", limit=0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(spot_failed && limit_failed)
}

///|
test "SpatialDecon: summary reports dimensions and convergence" {
  let summary = sd_test_example_result().summary()
  assert_true(summary.contains("genes=9"))
  assert_true(summary.contains("spots=5"))
  assert_true(summary.contains("cell_types=3"))
  assert_true(summary.contains("converged=5/5"))
}

///|
test "SpatialDecon: cell merge constructor defensively copies members" {
  let members = ["A", "B"]
  let merge = @src.SpatialDeconCellMerge::create("AB", members) catch {
    _ => abort("SpatialDecon merge should be valid")
  }
  members[0] = "changed"
  assert_eq(merge.name, "AB")
  assert_eq(merge.members, ["A", "B"])
}

///|
test "SpatialDecon: collapse sums abundance proportions and counts" {
  let result = sd_test_example_result()
  let lymphoid = @src.SpatialDeconCellMerge::create("Lymphoid", [
    "T_cell", "B_cell",
  ]) catch {
    _ => abort("SpatialDecon lymphoid merge should be valid")
  }
  let myeloid = @src.SpatialDeconCellMerge::create("Myeloid", ["Myeloid"]) catch {
    _ => abort("SpatialDecon myeloid merge should be valid")
  }
  let collapsed = @src.collapse_spatial_decon(result, [lymphoid, myeloid]) catch {
    _ => abort("SpatialDecon collapse should succeed")
  }
  assert_eq(collapsed.cell_types, ["Lymphoid", "Myeloid"])
  for spot in 0..<result.n_spots() {
    sd_test_close(
      collapsed.beta[0][spot],
      result.beta[0][spot] + result.beta[1][spot],
      1.0e-12,
    )
    sd_test_close(
      collapsed.cell_counts[0][spot],
      result.cell_counts[0][spot] + result.cell_counts[1][spot],
      1.0e-12,
    )
    sd_test_close(
      collapsed.proportions[0][spot],
      result.proportions[0][spot] + result.proportions[1][spot],
      1.0e-12,
    )
  }
}

///|
test "SpatialDecon: collapse applies A covariance A transpose" {
  let result = sd_test_example_result()
  let lymphoid = @src.SpatialDeconCellMerge::create("Lymphoid", [
    "T_cell", "B_cell",
  ]) catch {
    _ => abort("SpatialDecon lymphoid merge should be valid")
  }
  let collapsed = @src.collapse_spatial_decon(result, [lymphoid]) catch {
    _ => abort("SpatialDecon covariance collapse should succeed")
  }
  let expected = result.covariances[0][0][0] +
    result.covariances[0][0][1] +
    result.covariances[0][1][0] +
    result.covariances[0][1][1]
  sd_test_close(collapsed.covariances[0][0][0], expected, 1.0e-10)
  sd_test_close(collapsed.standard_errors[0][0], expected.sqrt(), 1.0e-10)
}

///|
test "SpatialDecon: collapse preserves gene-level fit diagnostics" {
  let result = sd_test_example_result()
  let all = @src.SpatialDeconCellMerge::create("All", [
    "T_cell", "B_cell", "Myeloid",
  ]) catch {
    _ => abort("SpatialDecon complete merge should be valid")
  }
  let collapsed = @src.collapse_spatial_decon(result, [all]) catch {
    _ => abort("SpatialDecon complete collapse should succeed")
  }
  assert_eq(collapsed.fitted, result.fitted)
  assert_eq(collapsed.residuals, result.residuals)
  assert_eq(collapsed.objectives, result.objectives)
  assert_eq(collapsed.converged, result.converged)
}

///|
test "SpatialDecon: collapse rejects empty unknown and repeated members" {
  let result = sd_test_example_result()
  let empty_failed = try {
    ignore(@src.collapse_spatial_decon(result, []))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let unknown = @src.SpatialDeconCellMerge::create("Unknown", ["missing"]) catch {
    _ => abort("SpatialDecon unknown merge constructor should succeed")
  }
  let unknown_failed = try {
    ignore(@src.collapse_spatial_decon(result, [unknown]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let first = @src.SpatialDeconCellMerge::create("First", ["T_cell"]) catch {
    _ => abort("SpatialDecon first merge should be valid")
  }
  let second = @src.SpatialDeconCellMerge::create("Second", ["T_cell"]) catch {
    _ => abort("SpatialDecon second merge should be valid")
  }
  let repeated_failed = try {
    ignore(@src.collapse_spatial_decon(result, [first, second]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(empty_failed && unknown_failed && repeated_failed)
}

///|
test "SpatialDecon: merge constructor rejects invalid names and members" {
  let name_failed = try {
    ignore(@src.SpatialDeconCellMerge::create("", ["A"]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let members_failed = try {
    ignore(@src.SpatialDeconCellMerge::create("AB", ["A", "A"]))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(name_failed && members_failed)
}

///|
test "SpatialDecon: probe-pool background averages negative probes" {
  let background = @src.spatial_decon_background(
    ["Neg_A1", "Neg_A2", "Marker_A", "Neg_B", "Marker_B"],
    ["s1", "s2"],
    [[1.0, 3.0], [3.0, 5.0], [10.0, 12.0], [2.0, 4.0], [20.0, 24.0]],
    ["A", "A", "A", "B", "B"],
    ["Neg_A1", "Neg_A2", "Neg_B"],
  ) catch {
    _ => abort("SpatialDecon background estimation should succeed")
  }
  assert_eq(background[0], [2.0, 4.0])
  assert_eq(background[2], [2.0, 4.0])
  assert_eq(background[3], [2.0, 4.0])
  assert_eq(background[4], [2.0, 4.0])
}

///|
test "SpatialDecon: probe-pool background validates pool coverage" {
  let length_failed = try {
    ignore(
      @src.spatial_decon_background(
        ["Neg", "Marker"],
        ["s1"],
        [[1.0], [2.0]],
        ["A"],
        ["Neg"],
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let coverage_failed = try {
    ignore(
      @src.spatial_decon_background(
        ["Neg", "Marker"],
        ["s1"],
        [[1.0], [2.0]],
        ["A", "B"],
        ["Neg"],
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(length_failed && coverage_failed)
}

///|
test "SpatialDecon: single-cell counts create mean cell-type profiles" {
  let profile = @src.create_spatial_decon_profile(
    ["g1", "g2"],
    ["a1", "a2", "b1", "b2"],
    ["A", "A", "B", "B"],
    [[10.0, 8.0, 1.0, 2.0], [1.0, 2.0, 9.0, 11.0]],
    scaling_factor=2.0,
    min_cells=1,
    min_genes=0,
  ) catch {
    _ => abort("SpatialDecon profile construction should succeed")
  }
  assert_eq(profile.gene_names, ["g1", "g2"])
  assert_eq(profile.cell_types, ["A", "B"])
  assert_eq(profile.values, [[18.0, 3.0], [3.0, 20.0]])
}

///|
test "SpatialDecon: single-cell profile supports normalization and gene filter" {
  let profile = @src.create_spatial_decon_profile(
    ["g1", "g2", "g3"],
    ["a1", "a2", "b1", "b2"],
    ["A", "A", "B", "B"],
    [[10.0, 8.0, 1.0, 2.0], [1.0, 2.0, 9.0, 11.0], [2.0, 2.0, 2.0, 2.0]],
    normalize=true,
    min_cells=1,
    min_genes=0,
    gene_filter=["g1", "g2"],
  ) catch {
    _ => abort("SpatialDecon normalized profile should succeed")
  }
  assert_eq(profile.gene_names, ["g1", "g2"])
  for row in profile.values {
    for value in row {
      assert_true(sd_test_finite(value) && value > 0.0)
    }
  }
}

///|
test "SpatialDecon: single-cell profile filters undersized cell types" {
  let failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1", "g2"],
        ["a1", "a2"],
        ["A", "A"],
        [[1.0, 2.0], [2.0, 1.0]],
        min_cells=2,
        min_genes=0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(failed)
}

///|
test "SpatialDecon: single-cell profile rejects invalid metadata and controls" {
  let metadata_failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1"],
        ["c1", "c2"],
        ["A"],
        [[1.0, 2.0]],
        min_cells=0,
        min_genes=0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let scaling_failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1"],
        ["c1"],
        ["A"],
        [[1.0]],
        scaling_factor=0.0,
        min_cells=0,
        min_genes=0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let threshold_failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1"],
        ["c1"],
        ["A"],
        [[1.0]],
        min_cells=-1,
        min_genes=0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(metadata_failed && scaling_failed && threshold_failed)
}

///|
test "SpatialDecon: single-cell profile rejects empty libraries and gene filters" {
  let library_failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1"],
        ["c1"],
        ["A"],
        [[0.0]],
        min_cells=0,
        min_genes=0,
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let gene_failed = try {
    ignore(
      @src.create_spatial_decon_profile(
        ["g1"],
        ["c1"],
        ["A"],
        [[1.0]],
        min_cells=0,
        min_genes=0,
        gene_filter=["missing"],
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(library_failed && gene_failed)
}

///|
test "SpatialDecon: reverse deconvolution fits genes from varying scores" {
  let (data, _, _) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon reverse data should be valid")
  }
  let result = sd_test_example_result()
  let reverse = @src.reverse_spatial_decon(data, result) catch {
    _ => abort("SpatialDecon reverse fit should succeed")
  }
  assert_eq(reverse.gene_names, data.gene_names)
  assert_eq(reverse.spot_ids, data.spot_ids)
  assert_eq(reverse.cell_types, result.cell_types)
  assert_eq(reverse.coefficients.length(), data.n_genes())
  assert_eq(reverse.coefficients[0].length(), result.n_cell_types() + 1)
  for row in reverse.coefficients {
    for value in row {
      assert_true(value >= 0.0 && sd_test_finite(value))
    }
  }
}

///|
test "SpatialDecon: reverse deconvolution returns residual diagnostics" {
  let (data, _, _) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon reverse data should be valid")
  }
  let reverse = @src.reverse_spatial_decon(data, sd_test_example_result()) catch {
    _ => abort("SpatialDecon reverse diagnostics should succeed")
  }
  assert_eq(reverse.fitted.length(), data.n_genes())
  assert_eq(reverse.fitted[0].length(), data.n_spots())
  for gene in 0..<data.n_genes() {
    assert_true(reverse.residual_sds[gene] >= 0.0)
    assert_true(sd_test_finite(reverse.correlations[gene]))
    assert_true(reverse.objectives[gene] >= 0.0)
  }
}

///|
test "SpatialDecon: reverse deconvolution validates spot identity and epsilon" {
  let result = sd_test_example_result()
  let (source, _, _) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon reverse source should be valid")
  }
  let mismatched = @src.SpatialDeconData::create(
    source.gene_names,
    ["x1", "x2", "x3", "x4", "x5"],
    source.values,
    background=source.background,
  ) catch {
    _ => abort("SpatialDecon mismatched reverse data should be valid")
  }
  let spots_failed = try {
    ignore(@src.reverse_spatial_decon(mismatched, result))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let epsilon_failed = try {
    ignore(@src.reverse_spatial_decon(source, result, epsilon=-1.0))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(spots_failed && epsilon_failed)
}

///|
test "SpatialDecon: reverse deconvolution requires varying cell scores" {
  let result = sd_test_exact_result()
  let constant_data = @src.SpatialDeconData::create(
    ["g1", "g2", "g3"],
    ["spot_1", "spot_2"],
    [[7.0, 7.0], [9.0, 9.0], [11.0, 11.0]],
    background=[[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]],
  ) catch {
    _ => abort("SpatialDecon constant data should be valid")
  }
  let constant_result = @src.spatial_decon(
    constant_data,
    sd_test_exact_profile(),
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon constant fit should succeed")
  }
  let failed = try {
    ignore(@src.reverse_spatial_decon(constant_data, constant_result))
    false
  } catch {
    SpatialDeconError(_) => true
  }
  ignore(result)
  assert_true(failed)
}

///|
test "SpatialDecon: SpatialExperiment wrapper writes scores and metadata" {
  let (experiment, profile, nuclei_counts) = sd_test_experiment()
  let output = @src.spatial_decon_spatial_experiment(
    experiment,
    profile,
    background_assay_name="background",
    weight_assay_name="weights",
    nuclei_counts~,
    config=sd_test_config(refit_outliers=true),
  ) catch {
    _ => abort("SpatialDecon SpatialExperiment wrapper should succeed")
  }
  assert_eq(output.experiment.metadata["spatialdecon_genes"], "9")
  assert_eq(output.experiment.metadata["spatialdecon_spots"], "5")
  assert_eq(output.experiment.metadata["spatialdecon_cell_types"], "3")
  assert_true(output.experiment.col_data[0].contains("SpatialDecon:T_cell"))
  assert_true(
    output.experiment.col_data[0].contains("SpatialDeconProp:Myeloid"),
  )
  assert_eq(output.result.cell_counts.length(), 3)
}

///|
test "SpatialDecon: SpatialExperiment wrapper does not mutate input" {
  let (experiment, profile, nuclei_counts) = sd_test_experiment()
  let original_value = experiment.assay["counts"][0][0]
  let output = @src.spatial_decon_spatial_experiment(
    experiment,
    profile,
    background_assay_name="background",
    nuclei_counts~,
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon immutable wrapper should succeed")
  }
  assert_true(!experiment.col_data[0].contains("SpatialDecon:T_cell"))
  assert_true(!experiment.metadata.contains("spatialdecon_genes"))
  output.experiment.assay["counts"][0][0] = 999.0
  output.experiment.col_data[0]["sample"] = "changed"
  output.experiment.metadata["owner"] = "changed"
  assert_eq(experiment.assay["counts"][0][0], original_value)
  assert_eq(experiment.col_data[0]["sample"], "sample_1")
  assert_eq(experiment.metadata["owner"], "input")
}

///|
test "SpatialDecon: SpatialExperiment wrapper generates missing identifiers" {
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "counts", [[2.0, 4.0], [4.0, 8.0]]))
  let profile = @src.SpatialDeconProfile::create(["gene_1", "gene_2"], ["A"], [
    [1.0],
    [2.0],
  ]) catch {
    _ => abort("SpatialDecon generated-name profile should be valid")
  }
  let output = @src.spatial_decon_spatial_experiment(
    experiment,
    profile,
    config=sd_test_config(),
  ) catch {
    _ => abort("SpatialDecon generated identifiers should succeed")
  }
  assert_eq(output.result.gene_names, ["gene_1", "gene_2"])
  assert_eq(output.result.spot_ids, ["spot_1", "spot_2"])
  assert_eq(output.experiment.col_data.length(), 2)
}

///|
test "SpatialDecon: SpatialExperiment wrapper validates assay names" {
  let (experiment, profile, _) = sd_test_experiment()
  let assay_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        experiment,
        profile,
        assay_name="missing",
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let background_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        experiment,
        profile,
        background_assay_name="missing",
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let weights_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        experiment,
        profile,
        weight_assay_name="missing",
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(assay_failed && background_failed && weights_failed)
}

///|
test "SpatialDecon: SpatialExperiment wrapper validates annotations" {
  let (experiment, profile, _) = sd_test_experiment()
  let broken_rows = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(broken_rows, "counts", experiment.assay["counts"]))
  ignore(@src.se_add_row(broken_rows, Map([("gene_id", "g1")])))
  let rows_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        broken_rows,
        profile,
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let missing_key = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(missing_key, "counts", experiment.assay["counts"]))
  for _ in 0..<9 {
    ignore(@src.se_add_row(missing_key, Map([("symbol", "x")])))
  }
  let key_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        missing_key,
        profile,
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(rows_failed && key_failed)
}

///|
test "SpatialDecon: SpatialExperiment wrapper validates prefixes and background" {
  let (experiment, profile, _) = sd_test_experiment()
  let prefix_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        experiment,
        profile,
        abundance_prefix="",
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  let background_failed = try {
    ignore(
      @src.spatial_decon_spatial_experiment(
        experiment,
        profile,
        scalar_background=-1.0,
        config=sd_test_config(),
      ),
    )
    false
  } catch {
    SpatialDeconError(_) => true
  }
  assert_true(prefix_failed && background_failed)
}

///|
test "SpatialDecon: example fixture contains deterministic spatial mixtures" {
  let (data, profile, nuclei_counts) = @src.spatial_decon_example_data() catch {
    _ => abort("SpatialDecon example fixture should succeed")
  }
  assert_eq(data.n_genes(), 9)
  assert_eq(data.n_spots(), 5)
  assert_eq(profile.n_cell_types(), 3)
  assert_eq(profile.cell_types, ["T_cell", "B_cell", "Myeloid"])
  assert_eq(nuclei_counts, [100.0, 120.0, 110.0, 90.0, 105.0])
  assert_true(data.background[0][4] > data.background[0][0])
}
