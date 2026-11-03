// Tests for the Bioconductor Banksy-inspired spatial clustering workflow.

///|
fn banksy_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "BANKSY value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn banksy_test_config(
  mode : @src.BanksySpatialMode,
) -> @src.BanksyComputeConfig {
  @src.BanksyComputeConfig::create(
    max_harmonic=1,
    k_geom=[4, 4],
    spatial_mode=mode,
    distance_power=2.0,
    sigma=1.5,
    alpha=0.05,
    k_spatial=4,
    seed=7,
  ) catch {
    _ => abort("BANKSY test configuration should be valid")
  }
}

///|
fn banksy_test_model() -> @src.BanksyModel {
  let (expression, coordinates, gene_names, spot_names) = @src.banksy_example_data()
  @src.banksy_compute_with_config(
    expression,
    coordinates,
    banksy_test_config(@src.banksy_knn_median()),
    gene_names~,
    spot_names~,
  ) catch {
    _ => abort("BANKSY example model should compute")
  }
}

///|
fn banksy_test_matrix() -> @src.BanksyMatrix {
  @src.banksy_get_matrix(
    banksy_test_model(),
    lambda=0.2,
    max_harmonic=1,
    scale=true,
  ) catch {
    _ => abort("BANKSY example matrix should compute")
  }
}

///|
fn banksy_test_pca() -> @src.BanksyPcaResult {
  @src.banksy_run_pca(banksy_test_matrix(), n_components=3) catch {
    _ => abort("BANKSY example PCA should compute")
  }
}

///|
fn banksy_test_cluster_config() -> @src.BanksyClusterConfig {
  @src.BanksyClusterConfig::create(
    n_clusters=2,
    n_starts=5,
    max_iterations=100,
    tolerance=1.0e-8,
    seed=3,
  ) catch {
    _ => abort("BANKSY cluster configuration should be valid")
  }
}

///|
fn banksy_test_experiment() -> @src.SpatialExperiment {
  let (expression, coordinates, gene_names, spot_names) = @src.banksy_example_data()
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", expression))
  for gene_name in gene_names {
    ignore(@src.se_add_row(experiment, Map([("gene_name", gene_name)])))
  }
  for spot in 0..<coordinates.length() {
    let section = if spot < coordinates.length() / 2 { "A" } else { "B" }
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", spot_names[spot]), ("section", section)]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(coordinates[spot][0], coordinates[spot][1]),
      ),
    )
  }
  experiment
}

///|
test "BANKSY: spatial mode helpers expose all kernels" {
  assert_true(
    @src.banksy_knn_median() is @src.BanksySpatialMode::BanksyKnnMedian,
  )
  assert_true(
    @src.banksy_knn_inverse_distance()
    is @src.BanksySpatialMode::BanksyKnnInverseDistance,
  )
  assert_true(
    @src.banksy_knn_inverse_power()
    is @src.BanksySpatialMode::BanksyKnnInversePower,
  )
  assert_true(@src.banksy_knn_rank() is @src.BanksySpatialMode::BanksyKnnRank)
  assert_true(
    @src.banksy_knn_uniform() is @src.BanksySpatialMode::BanksyKnnUniform,
  )
  assert_true(
    @src.banksy_radius_gaussian()
    is @src.BanksySpatialMode::BanksyRadiusGaussian,
  )
}

///|
test "BANKSY: default compute configuration matches upstream H0 behavior" {
  let config = @src.BanksyComputeConfig::default()
  assert_eq(config.max_harmonic, 0)
  assert_eq(config.k_geom, [15])
  assert_true(config.spatial_mode is @src.BanksySpatialMode::BanksyKnnMedian)
  assert_eq(config.center_harmonics, true)
}

///|
test "BANKSY: custom compute configuration preserves controls" {
  let config = @src.BanksyComputeConfig::create(
    max_harmonic=2,
    k_geom=[3, 4, 5],
    spatial_mode=@src.banksy_knn_inverse_power(),
    distance_power=3.0,
    sample_size=2,
    sample_renormalize=false,
    seed=9,
    center_harmonics=false,
  ) catch {
    _ => abort("custom BANKSY configuration should be valid")
  }
  assert_eq(config.max_harmonic, 2)
  assert_eq(config.k_geom, [3, 4, 5])
  assert_eq(config.distance_power, 3.0)
  assert_eq(config.sample_size, 2)
  assert_eq(config.sample_renormalize, false)
  assert_eq(config.seed, 9)
  assert_eq(config.center_harmonics, false)
}

///|
test "BANKSY: compute configuration copies neighbor counts" {
  let counts = [2, 3]
  let config = @src.BanksyComputeConfig::create(max_harmonic=1, k_geom=counts) catch {
    _ => abort("BANKSY configuration should copy k_geom")
  }
  counts[0] = 9
  assert_eq(config.k_geom, [2, 3])
}

///|
test "BANKSY: compute configuration rejects negative harmonic" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(max_harmonic=-1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects wrong k_geom length" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(max_harmonic=2, k_geom=[2, 3]))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects non-positive k_geom" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(k_geom=[0]))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects invalid distance power" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(distance_power=0.0))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects invalid Gaussian sigma" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(sigma=-1.0))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects radial alpha endpoints" {
  let zero = try {
    ignore(@src.BanksyComputeConfig::create(alpha=0.0))
    false
  } catch {
    BanksyError(_) => true
  }
  let one = try {
    ignore(@src.BanksyComputeConfig::create(alpha=1.0))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(zero)
  assert_true(one)
}

///|
test "BANKSY: compute configuration rejects radial neighbor count" {
  let failed = try {
    ignore(@src.BanksyComputeConfig::create(k_spatial=0))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute configuration rejects negative sampling controls" {
  let sample = try {
    ignore(@src.BanksyComputeConfig::create(sample_size=-1))
    false
  } catch {
    BanksyError(_) => true
  }
  let seed = try {
    ignore(@src.BanksyComputeConfig::create(seed=-1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(sample)
  assert_true(seed)
}

///|
test "BANKSY: example data has gene by spot orientation" {
  let (expression, coordinates, gene_names, spot_names) = @src.banksy_example_data()
  assert_eq(expression.length(), 5)
  assert_eq(expression[0].length(), 24)
  assert_eq(coordinates.length(), 24)
  assert_eq(coordinates[0].length(), 2)
  assert_eq(gene_names.length(), 5)
  assert_eq(spot_names.length(), 24)
}

///|
test "BANKSY: simple compute defaults to H0 and generated names" {
  let model = @src.banksy_compute(
    [[1.0, 2.0, 3.0, 4.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
    k_geom=2,
  ) catch {
    _ => abort("simple BANKSY compute should succeed")
  }
  assert_eq(model.n_genes(), 1)
  assert_eq(model.n_spots(), 4)
  assert_eq(model.max_harmonic(), 0)
  assert_eq(model.gene_names, ["gene_1"])
  assert_eq(model.spot_names, ["spot_1", "spot_2", "spot_3", "spot_4"])
}

///|
test "BANKSY: compute_agf adds H1" {
  let model = @src.banksy_compute(
    [[1.0, 2.0, 3.0, 4.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
    k_geom=2,
    compute_agf=true,
  ) catch {
    _ => abort("BANKSY AGF compute should succeed")
  }
  assert_eq(model.max_harmonic(), 1)
  assert_true(model.harmonic(0) is Some(_))
  assert_true(model.harmonic(1) is Some(_))
  assert_true(model.harmonic(2) is None)
}

///|
test "BANKSY: model owns defensive copies of input" {
  let expression = [[1.0, 2.0, 3.0]]
  let coordinates = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]
  let model = @src.banksy_compute(expression, coordinates, k_geom=1) catch {
    _ => abort("BANKSY defensive copy fixture should compute")
  }
  expression[0][0] = 99.0
  coordinates[0][0] = 99.0
  assert_eq(model.expression[0][0], 1.0)
  assert_eq(model.coordinates[0][0], 0.0)
}

///|
test "BANKSY: harmonic accessor returns a defensive copy" {
  let model = banksy_test_model()
  let first = match model.harmonic(0) {
    Some(value) => value
    None => abort("H0 should exist")
  }
  let expected = first[0][0]
  first[0][0] = 999.0
  let second = match model.harmonic(0) {
    Some(value) => value
    None => abort("H0 should exist")
  }
  assert_eq(second[0][0], expected)
}

///|
test "BANKSY: H0 preserves constant expression" {
  let model = @src.banksy_compute(
    [[7.0, 7.0, 7.0, 7.0]],
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
    k_geom=3,
  ) catch {
    _ => abort("constant BANKSY fixture should compute")
  }
  for value in model.harmonics[0][0] {
    banksy_test_close(value, 7.0, 1.0e-12)
  }
}

///|
test "BANKSY: centered H1 is zero for constant expression" {
  let model = @src.banksy_compute(
    [[4.0, 4.0, 4.0, 4.0]],
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
    k_geom=3,
    compute_agf=true,
  ) catch {
    _ => abort("constant BANKSY AGF fixture should compute")
  }
  for value in model.harmonics[1][0] {
    banksy_test_close(value, 0.0, 1.0e-12)
  }
}

///|
test "BANKSY: directional gradient produces positive H1" {
  let model = @src.banksy_compute(
    [[0.0, 1.0, 2.0, 3.0, 4.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
    k_geom=2,
    compute_agf=true,
  ) catch {
    _ => abort("gradient BANKSY fixture should compute")
  }
  assert_true(model.harmonics[1][0][2] > 0.0)
}

///|
test "BANKSY: each spatial kernel produces normalized weights" {
  let modes = [
    @src.banksy_knn_median(),
    @src.banksy_knn_inverse_distance(),
    @src.banksy_knn_inverse_power(),
    @src.banksy_knn_rank(),
    @src.banksy_knn_uniform(),
    @src.banksy_radius_gaussian(),
  ]
  let (expression, coordinates, gene_names, spot_names) = @src.banksy_example_data()
  for mode in modes {
    let model = @src.banksy_compute_with_config(
      expression,
      coordinates,
      banksy_test_config(mode),
      gene_names~,
      spot_names~,
    ) catch {
      _ => abort("all BANKSY kernels should compute")
    }
    for neighborhood in model.neighbors[0] {
      let mut total = 0.0
      for neighbor in neighborhood {
        total = total + neighbor.weight
      }
      banksy_test_close(total, 1.0, 1.0e-10)
    }
  }
}

///|
test "BANKSY: uniform kernel assigns equal weights" {
  let model = @src.banksy_compute_with_config(
    [[1.0, 2.0, 3.0, 4.0, 5.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
    @src.BanksyComputeConfig::create(
      max_harmonic=0,
      k_geom=[3],
      spatial_mode=@src.banksy_knn_uniform(),
    ) catch {
      _ => abort("uniform configuration should be valid")
    },
  ) catch {
    _ => abort("uniform kernel should compute")
  }
  for neighbor in model.neighbors[0][2] {
    banksy_test_close(neighbor.weight, 1.0 / 3.0, 1.0e-12)
  }
}

///|
test "BANKSY: inverse-distance kernel favors closer spots" {
  let model = @src.banksy_compute_with_config(
    [[1.0, 2.0, 3.0, 4.0]],
    [[0.0, 0.0], [1.0, 0.0], [3.0, 0.0], [7.0, 0.0]],
    @src.BanksyComputeConfig::create(
      k_geom=[3],
      spatial_mode=@src.banksy_knn_inverse_distance(),
    ) catch {
      _ => abort("inverse-distance configuration should be valid")
    },
  ) catch {
    _ => abort("inverse-distance kernel should compute")
  }
  assert_true(model.neighbors[0][0][0].weight > model.neighbors[0][0][1].weight)
  assert_true(model.neighbors[0][0][1].weight > model.neighbors[0][0][2].weight)
}

///|
test "BANKSY: rank kernel weights decrease with rank" {
  let model = @src.banksy_compute_with_config(
    [[1.0, 2.0, 3.0, 4.0, 5.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
    @src.BanksyComputeConfig::create(
      k_geom=[4],
      spatial_mode=@src.banksy_knn_rank(),
    ) catch {
      _ => abort("rank configuration should be valid")
    },
  ) catch {
    _ => abort("rank kernel should compute")
  }
  let neighborhood = model.neighbors[0][0]
  assert_true(neighborhood[0].weight > neighborhood[1].weight)
  assert_true(neighborhood[1].weight > neighborhood[2].weight)
  assert_true(neighborhood[2].weight > neighborhood[3].weight)
}

///|
test "BANKSY: neighbor angles use the first two spatial dimensions" {
  let model = @src.banksy_compute(
    [[1.0, 2.0, 3.0]],
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    k_geom=2,
  ) catch {
    _ => abort("angle fixture should compute")
  }
  banksy_test_close(model.neighbors[0][0][0].angle, 0.0, 1.0e-12)
  banksy_test_close(
    model.neighbors[0][0][1].angle,
    3.141592653589793 / 2.0,
    1.0e-12,
  )
}

///|
test "BANKSY: deterministic sampling selects requested neighbors" {
  let config = @src.BanksyComputeConfig::create(
    k_geom=[5],
    sample_size=2,
    sample_renormalize=true,
    seed=11,
  ) catch {
    _ => abort("sampling configuration should be valid")
  }
  let expression = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]
  let coordinates = [
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [3.0, 0.0],
    [4.0, 0.0],
    [5.0, 0.0],
  ]
  let first = @src.banksy_compute_with_config(expression, coordinates, config) catch {
    _ => abort("sampled BANKSY model should compute")
  }
  let second = @src.banksy_compute_with_config(expression, coordinates, config) catch {
    _ => abort("sampled BANKSY model should repeat")
  }
  assert_eq(first.neighbors[0][0].length(), 2)
  assert_eq(first.neighbors[0][0], second.neighbors[0][0])
  banksy_test_close(
    first.neighbors[0][0][0].weight + first.neighbors[0][0][1].weight,
    1.0,
    1.0e-12,
  )
}

///|
test "BANKSY: non-renormalized sampling retains original mass" {
  let model = @src.banksy_compute_with_config(
    [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [5.0, 0.0]],
    @src.BanksyComputeConfig::create(
      k_geom=[5],
      sample_size=2,
      sample_renormalize=false,
    ) catch {
      _ => abort("non-renormalized sampling configuration should be valid")
    },
  ) catch {
    _ => abort("non-renormalized sampling should compute")
  }
  let mut total = 0.0
  for neighbor in model.neighbors[0][0] {
    total = total + neighbor.weight
  }
  assert_true(total > 0.0)
  assert_true(total < 1.0)
}

///|
test "BANKSY: model summary reports dimensions and harmonics" {
  let summary = banksy_test_model().summary()
  assert_true(summary.contains("5 genes"))
  assert_true(summary.contains("24 spots"))
  assert_true(summary.contains("H0-H1"))
}

///|
test "BANKSY: compute rejects empty expression" {
  let failed = try {
    ignore(@src.banksy_compute([], []))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects a single spot" {
  let failed = try {
    ignore(@src.banksy_compute([[1.0]], [[0.0, 0.0]], k_geom=1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects ragged expression" {
  let failed = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 2.0, 3.0], [4.0, 5.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        k_geom=1,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects non-finite expression" {
  let failed = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 1.0e301, 3.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        k_geom=1,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects coordinate count mismatch" {
  let failed = try {
    ignore(
      @src.banksy_compute([[1.0, 2.0, 3.0]], [[0.0, 0.0], [1.0, 0.0]], k_geom=1),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects one-dimensional coordinates" {
  let failed = try {
    ignore(
      @src.banksy_compute([[1.0, 2.0, 3.0]], [[0.0], [1.0], [2.0]], k_geom=1),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects ragged coordinates" {
  let failed = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 2.0, 3.0]],
        [[0.0, 0.0], [1.0, 0.0, 1.0], [2.0, 0.0]],
        k_geom=1,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: compute rejects duplicate names" {
  let duplicate_gene = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        k_geom=1,
        gene_names=["g", "g"],
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  let duplicate_spot = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 2.0, 3.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        k_geom=1,
        spot_names=["s", "s", "t"],
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(duplicate_gene)
  assert_true(duplicate_spot)
}

///|
test "BANKSY: compute rejects excessive neighbor and sample counts" {
  let neighbors = try {
    ignore(
      @src.banksy_compute(
        [[1.0, 2.0, 3.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        k_geom=3,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  let sample = try {
    let config = @src.BanksyComputeConfig::create(k_geom=[2], sample_size=3)
    ignore(
      @src.banksy_compute_with_config(
        [[1.0, 2.0, 3.0]],
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        config,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(neighbors)
  assert_true(sample)
}

///|
test "BANKSY: M1 component weights match upstream lambda split" {
  let weights = @src.banksy_component_weights(0.2, 1) catch {
    _ => abort("valid BANKSY weights should compute")
  }
  banksy_test_close(weights[0], 0.8.sqrt(), 1.0e-12)
  banksy_test_close(weights[1], (0.4 / 3.0).sqrt(), 1.0e-12)
  banksy_test_close(weights[2], (0.2 / 3.0).sqrt(), 1.0e-12)
}

///|
test "BANKSY: lambda endpoints isolate expression and spatial components" {
  let expression = @src.banksy_component_weights(0.0, 1) catch {
    _ => abort("lambda zero should be valid")
  }
  let spatial = @src.banksy_component_weights(1.0, 1) catch {
    _ => abort("lambda one should be valid")
  }
  assert_eq(expression[0], 1.0)
  assert_eq(expression[1], 0.0)
  assert_eq(expression[2], 0.0)
  assert_eq(spatial[0], 0.0)
  banksy_test_close(
    spatial[1] * spatial[1] + spatial[2] * spatial[2],
    1.0,
    1.0e-12,
  )
}

///|
test "BANKSY: component weights reject invalid arguments" {
  let lambda = try {
    ignore(@src.banksy_component_weights(1.1, 1))
    false
  } catch {
    BanksyError(_) => true
  }
  let harmonic = try {
    ignore(@src.banksy_component_weights(0.2, -1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(lambda)
  assert_true(harmonic)
}

///|
test "BANKSY: joint matrix concatenates expression H0 and H1" {
  let matrix = @src.banksy_get_matrix(
    banksy_test_model(),
    lambda=0.2,
    max_harmonic=1,
  ) catch {
    _ => abort("BANKSY matrix should compute")
  }
  assert_eq(matrix.n_features(), 15)
  assert_eq(matrix.n_spots(), 24)
  assert_eq(matrix.feature_names[0], "left_marker")
  assert_eq(matrix.feature_names[5], "left_marker_H0")
  assert_eq(matrix.feature_names[10], "left_marker_H1")
  assert_eq(matrix.component_weights.length(), 3)
}

///|
test "BANKSY: spot matrix transposes feature orientation" {
  let matrix = @src.banksy_get_matrix(
    banksy_test_model(),
    lambda=0.2,
    max_harmonic=0,
  ) catch {
    _ => abort("BANKSY matrix should compute")
  }
  let spots = matrix.spot_matrix()
  assert_eq(spots.length(), matrix.n_spots())
  assert_eq(spots[0].length(), matrix.n_features())
  assert_eq(spots[3][2], matrix.values[2][3])
}

///|
test "BANKSY: global scaling centers every feature" {
  let matrix = banksy_test_matrix()
  for row in matrix.values {
    let mut mean = 0.0
    for value in row {
      mean = mean + value
    }
    mean = mean / row.length().to_double()
    banksy_test_close(mean, 0.0, 1.0e-10)
  }
}

///|
test "BANKSY: scaling maps constant features to zero" {
  let model = @src.banksy_compute(
    [[5.0, 5.0, 5.0, 5.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
    k_geom=2,
  ) catch {
    _ => abort("constant BANKSY model should compute")
  }
  let matrix = @src.banksy_get_matrix(model, scale=true) catch {
    _ => abort("constant BANKSY matrix should scale")
  }
  for row in matrix.values {
    for value in row {
      banksy_test_close(value, 0.0, 1.0e-12)
    }
  }
}

///|
test "BANKSY: group scaling centers each group independently" {
  let model = @src.banksy_compute(
    [[1.0, 3.0, 101.0, 103.0]],
    [[0.0, 0.0], [1.0, 0.0], [10.0, 0.0], [11.0, 0.0]],
    k_geom=1,
  ) catch {
    _ => abort("group scaling fixture should compute")
  }
  let matrix = @src.banksy_get_matrix(
    model,
    lambda=0.0,
    max_harmonic=0,
    scale=true,
    groups=["A", "A", "B", "B"],
  ) catch {
    _ => abort("group scaling should compute")
  }
  banksy_test_close(matrix.values[0][0] + matrix.values[0][1], 0.0, 1.0e-12)
  banksy_test_close(matrix.values[0][2] + matrix.values[0][3], 0.0, 1.0e-12)
}

///|
test "BANKSY: joint matrix validates groups and harmonics" {
  let model = banksy_test_model()
  let groups = try {
    ignore(@src.banksy_get_matrix(model, scale=true, groups=["A"]))
    false
  } catch {
    BanksyError(_) => true
  }
  let harmonic = try {
    ignore(@src.banksy_get_matrix(model, max_harmonic=2))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(groups)
  assert_true(harmonic)
}

///|
test "BANKSY: PCA returns requested dimensions" {
  let pca = banksy_test_pca()
  assert_eq(pca.n_components, 3)
  assert_eq(pca.n_spots(), 24)
  assert_eq(pca.scores[0].length(), 3)
  assert_eq(pca.loadings.length(), 15)
  assert_eq(pca.loadings[0].length(), 3)
  assert_eq(pca.eigenvalues.length(), 3)
}

///|
test "BANKSY: PCA scores are centered" {
  let pca = banksy_test_pca()
  for component in 0..<pca.n_components {
    let mut total = 0.0
    for spot in 0..<pca.n_spots() {
      total = total + pca.scores[spot][component]
    }
    banksy_test_close(total, 0.0, 1.0e-8)
  }
}

///|
test "BANKSY: PCA eigenvalues are ordered and variance is bounded" {
  let pca = banksy_test_pca()
  assert_true(pca.eigenvalues[0] >= pca.eigenvalues[1])
  assert_true(pca.eigenvalues[1] >= pca.eigenvalues[2])
  let mut total = 0.0
  for value in pca.variance_explained {
    assert_true(value >= 0.0)
    assert_true(value <= 1.0)
    total = total + value
  }
  assert_true(total > 0.0)
  assert_true(total <= 1.0 + 1.0e-10)
}

///|
test "BANKSY: PCA summary reports component count" {
  let summary = banksy_test_pca().summary()
  assert_true(summary.contains("3 components"))
  assert_true(summary.contains("explained variance"))
}

///|
test "BANKSY: PCA validates component count" {
  let matrix = banksy_test_matrix()
  let zero = try {
    ignore(@src.banksy_run_pca(matrix, n_components=0))
    false
  } catch {
    BanksyError(_) => true
  }
  let large = try {
    ignore(@src.banksy_run_pca(matrix, n_components=24))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(zero)
  assert_true(large)
}

///|
test "BANKSY: PCA rejects constant input" {
  let model = @src.banksy_compute(
    [[5.0, 5.0, 5.0, 5.0]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
    k_geom=2,
  ) catch {
    _ => abort("constant PCA fixture should compute")
  }
  let matrix = @src.banksy_get_matrix(model, lambda=0.0) catch {
    _ => abort("constant PCA matrix should compute")
  }
  let failed = try {
    ignore(@src.banksy_run_pca(matrix, n_components=1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: default cluster configuration exposes deterministic controls" {
  let config = @src.BanksyClusterConfig::default()
  assert_eq(config.n_clusters, 5)
  assert_eq(config.n_starts, 10)
  assert_eq(config.max_iterations, 100)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.seed, 1)
}

///|
test "BANKSY: cluster configuration preserves custom controls" {
  let config = banksy_test_cluster_config()
  assert_eq(config.n_clusters, 2)
  assert_eq(config.n_starts, 5)
  assert_eq(config.max_iterations, 100)
  assert_eq(config.tolerance, 1.0e-8)
  assert_eq(config.seed, 3)
}

///|
test "BANKSY: cluster configuration rejects invalid controls" {
  let clusters = try {
    ignore(@src.BanksyClusterConfig::create(n_clusters=0))
    false
  } catch {
    BanksyError(_) => true
  }
  let starts = try {
    ignore(@src.BanksyClusterConfig::create(n_starts=0))
    false
  } catch {
    BanksyError(_) => true
  }
  let iterations = try {
    ignore(@src.BanksyClusterConfig::create(max_iterations=0))
    false
  } catch {
    BanksyError(_) => true
  }
  let tolerance = try {
    ignore(@src.BanksyClusterConfig::create(tolerance=0.0))
    false
  } catch {
    BanksyError(_) => true
  }
  let seed = try {
    ignore(@src.BanksyClusterConfig::create(seed=-1))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(clusters)
  assert_true(starts)
  assert_true(iterations)
  assert_true(tolerance)
  assert_true(seed)
}

///|
test "BANKSY: deterministic k-means recovers two separated groups" {
  let model = @src.banksy_compute(
    [[-5.0, -4.8, -5.2, 5.0, 4.8, 5.2]],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [10.0, 0.0], [11.0, 0.0], [12.0, 0.0]],
    k_geom=1,
  ) catch {
    _ => abort("separated BANKSY model should compute")
  }
  let matrix = @src.banksy_get_matrix(model, lambda=0.0, max_harmonic=0) catch {
    _ => abort("separated BANKSY matrix should compute")
  }
  let pca = @src.banksy_run_pca(matrix, n_components=1) catch {
    _ => abort("separated BANKSY PCA should compute")
  }
  let fit = @src.banksy_cluster(pca, config=banksy_test_cluster_config()) catch {
    _ => abort("separated BANKSY clusters should fit")
  }
  assert_eq(fit.n_clusters(), 2)
  assert_eq(fit.cluster_sizes[0] + fit.cluster_sizes[1], 6)
  assert_eq(fit.labels[0], fit.labels[1])
  assert_eq(fit.labels[1], fit.labels[2])
  assert_eq(fit.labels[3], fit.labels[4])
  assert_eq(fit.labels[4], fit.labels[5])
  assert_true(fit.labels[0] != fit.labels[3])
  assert_true(fit.silhouette > 0.9)
}

///|
test "BANKSY: k-means is reproducible for fixed configuration" {
  let pca = banksy_test_pca()
  let first = @src.banksy_cluster(pca, config=banksy_test_cluster_config()) catch {
    _ => abort("first BANKSY clustering should fit")
  }
  let second = @src.banksy_cluster(pca, config=banksy_test_cluster_config()) catch {
    _ => abort("second BANKSY clustering should fit")
  }
  assert_eq(first.labels, second.labels)
  assert_eq(first.inertia, second.inertia)
}

///|
test "BANKSY: cluster result summary reports fit metrics" {
  let fit = @src.banksy_cluster(
    banksy_test_pca(),
    config=banksy_test_cluster_config(),
  ) catch {
    _ => abort("BANKSY clustering should fit")
  }
  let summary = fit.summary()
  assert_true(summary.contains("2 clusters"))
  assert_true(summary.contains("inertia"))
  assert_true(summary.contains("silhouette"))
}

///|
test "BANKSY: clustering rejects more clusters than spots" {
  let failed = try {
    ignore(
      @src.banksy_cluster(
        banksy_test_pca(),
        config=@src.BanksyClusterConfig::create(n_clusters=25),
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: silhouette is high for separated groups" {
  let value = @src.banksy_silhouette([[0.0], [0.1], [10.0], [10.1]], [
    0, 0, 1, 1,
  ]) catch {
    _ => abort("valid silhouette should compute")
  }
  assert_true(value > 0.98)
}

///|
test "BANKSY: silhouette validates labels" {
  let length = try {
    ignore(@src.banksy_silhouette([[0.0], [1.0]], [0]))
    false
  } catch {
    BanksyError(_) => true
  }
  let negative = try {
    ignore(@src.banksy_silhouette([[0.0], [1.0]], [0, -1]))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(length)
  assert_true(negative)
}

///|
test "BANKSY: label smoothing repairs an isolated label" {
  let result = @src.banksy_smooth_labels(
    [0, 0, 1, 0, 0],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]],
    k=2,
    proportion_threshold=0.5,
    max_iterations=5,
  ) catch {
    _ => abort("BANKSY smoothing should compute")
  }
  assert_eq(result.labels, [0, 0, 0, 0, 0])
  assert_eq(result.changed_labels, 1)
  assert_true(result.converged)
}

///|
test "BANKSY: strict smoothing threshold can retain labels" {
  let result = @src.banksy_smooth_labels(
    [0, 1, 0, 1],
    [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
    k=1,
    proportion_threshold=1.0,
    max_iterations=2,
  ) catch {
    _ => abort("strict BANKSY smoothing should compute")
  }
  assert_eq(result.labels, [0, 1, 0, 1])
  assert_eq(result.changed_labels, 0)
}

///|
test "BANKSY: smoothing validates controls" {
  let k = try {
    ignore(@src.banksy_smooth_labels([0, 1], [[0.0, 0.0], [1.0, 0.0]], k=2))
    false
  } catch {
    BanksyError(_) => true
  }
  let threshold = try {
    ignore(
      @src.banksy_smooth_labels(
        [0, 1],
        [[0.0, 0.0], [1.0, 0.0]],
        k=1,
        proportion_threshold=1.1,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  let iterations = try {
    ignore(
      @src.banksy_smooth_labels(
        [0, 1],
        [[0.0, 0.0], [1.0, 0.0]],
        k=1,
        max_iterations=0,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(k)
  assert_true(threshold)
  assert_true(iterations)
}

///|
test "BANKSY: spatial agreement identifies homogeneous labels" {
  let coordinates = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]
  let homogeneous = @src.banksy_spatial_agreement(
    [0, 0, 0, 0],
    coordinates,
    k=2,
  ) catch {
    _ => abort("homogeneous agreement should compute")
  }
  let alternating = @src.banksy_spatial_agreement(
    [0, 1, 0, 1],
    coordinates,
    k=1,
  ) catch {
    _ => abort("alternating agreement should compute")
  }
  assert_eq(homogeneous, 1.0)
  assert_true(alternating < homogeneous)
}

///|
test "BANKSY: ARI is one for identical and relabeled partitions" {
  let identical = @src.banksy_adjusted_rand_index([0, 0, 1, 1, 2, 2], [
    0, 0, 1, 1, 2, 2,
  ]) catch {
    _ => abort("identical ARI should compute")
  }
  let relabeled = @src.banksy_adjusted_rand_index([0, 0, 1, 1, 2, 2], [
    2, 2, 0, 0, 1, 1,
  ]) catch {
    _ => abort("relabeled ARI should compute")
  }
  assert_eq(identical, 1.0)
  assert_eq(relabeled, 1.0)
}

///|
test "BANKSY: ARI detects disagreeing partitions" {
  let value = @src.banksy_adjusted_rand_index([0, 0, 0, 1, 1, 1], [
    0, 1, 0, 1, 0, 1,
  ]) catch {
    _ => abort("disagreeing ARI should compute")
  }
  assert_true(value < 0.3)
}

///|
test "BANKSY: ARI validates label vectors" {
  let length = try {
    ignore(@src.banksy_adjusted_rand_index([0, 1], [0]))
    false
  } catch {
    BanksyError(_) => true
  }
  let negative = try {
    ignore(@src.banksy_adjusted_rand_index([0, -1], [0, 1]))
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(length)
  assert_true(negative)
}

///|
test "BANKSY: complete run exposes matrix PCA clustering and agreement" {
  let run = @src.banksy_run(
    banksy_test_model(),
    lambda=0.2,
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
    scale=true,
    smooth=false,
    smoothing_k=3,
  ) catch {
    _ => abort("complete BANKSY run should compute")
  }
  assert_eq(run.lambda, 0.2)
  assert_eq(run.matrix.n_features(), 15)
  assert_eq(run.pca.n_components, 3)
  assert_eq(run.clustering.n_clusters(), 2)
  assert_true(run.smoothing is None)
  assert_true(run.spatial_agreement >= 0.0)
  assert_true(run.spatial_agreement <= 1.0)
  assert_eq(run.final_labels(), run.clustering.labels)
}

///|
test "BANKSY: complete run can smooth labels" {
  let run = @src.banksy_run(
    banksy_test_model(),
    lambda=0.8,
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
    smooth=true,
    smoothing_k=3,
    smoothing_iterations=5,
  ) catch {
    _ => abort("smoothed BANKSY run should compute")
  }
  assert_true(run.smoothing is Some(_))
  assert_eq(run.final_labels().length(), 24)
}

///|
test "BANKSY: run summary includes lambda and spatial agreement" {
  let run = @src.banksy_run(
    banksy_test_model(),
    lambda=0.2,
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    smoothing_k=3,
  ) catch {
    _ => abort("BANKSY summary run should compute")
  }
  let summary = run.summary()
  assert_true(summary.contains("lambda=0.2"))
  assert_true(summary.contains("spatial agreement"))
}

///|
test "BANKSY: parameter sweep preserves lambda order" {
  let runs = @src.banksy_parameter_sweep(
    banksy_test_model(),
    [0.0, 0.2, 0.8],
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
  ) catch {
    _ => abort("BANKSY parameter sweep should compute")
  }
  assert_eq(runs.length(), 3)
  assert_eq(runs[0].lambda, 0.0)
  assert_eq(runs[1].lambda, 0.2)
  assert_eq(runs[2].lambda, 0.8)
}

///|
test "BANKSY: parameter sweep rejects empty lambda list" {
  let failed = try {
    ignore(
      @src.banksy_parameter_sweep(
        banksy_test_model(),
        [],
        n_components=3,
        cluster_config=banksy_test_cluster_config(),
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}

///|
test "BANKSY: SpatialExperiment wrapper writes harmonics and labels" {
  let output = @src.banksy_spatial_experiment(
    banksy_test_experiment(),
    config=banksy_test_config(@src.banksy_knn_median()),
    lambda=0.2,
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
    group_key="section",
    smooth=true,
    smoothing_k=3,
    smoothing_iterations=5,
  ) catch {
    _ => abort("BANKSY SpatialExperiment integration should compute")
  }
  assert_true(output.experiment.assay.contains("logcounts"))
  assert_true(output.experiment.assay.contains("H0"))
  assert_true(output.experiment.assay.contains("H1"))
  assert_eq(output.experiment.assay["H0"].length(), 5)
  assert_eq(output.experiment.assay["H0"][0].length(), 24)
  assert_true(output.experiment.col_data[0].contains("banksy_cluster"))
  assert_true(output.experiment.col_data[0].contains("banksy_cluster_raw"))
  assert_eq(output.experiment.metadata["banksy_assay"], "logcounts")
  assert_eq(output.experiment.metadata["banksy_lambda"], "0.2")
  assert_eq(output.experiment.metadata["banksy_max_harmonic"], "1")
  assert_eq(output.model.n_spots(), 24)
}

///|
test "BANKSY: SpatialExperiment wrapper does not mutate its input" {
  let experiment = banksy_test_experiment()
  let output = @src.banksy_spatial_experiment(
    experiment,
    config=banksy_test_config(@src.banksy_knn_median()),
    lambda=0.2,
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
    smooth=false,
    smoothing_k=3,
  ) catch {
    _ => abort("immutable BANKSY integration should compute")
  }
  assert_true(!experiment.assay.contains("H0"))
  assert_true(!experiment.col_data[0].contains("banksy_cluster"))
  assert_true(!experiment.metadata.contains("banksy_lambda"))
  output.experiment.assay["H0"][0][0] = 999.0
  assert_true(experiment.assay["logcounts"][0][0] != 999.0)
}

///|
test "BANKSY: SpatialExperiment wrapper generates missing annotations" {
  let (expression, coordinates, _, _) = @src.banksy_example_data()
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", expression))
  for coordinate in coordinates {
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(coordinate[0], coordinate[1]),
      ),
    )
  }
  let output = @src.banksy_spatial_experiment(
    experiment,
    config=banksy_test_config(@src.banksy_knn_median()),
    n_components=3,
    cluster_config=banksy_test_cluster_config(),
    max_harmonic=1,
    smooth=false,
    smoothing_k=3,
  ) catch {
    _ => abort("BANKSY should generate missing annotations")
  }
  assert_eq(output.experiment.row_data.length(), 5)
  assert_eq(output.experiment.col_data.length(), 24)
  assert_eq(output.experiment.row_data[0]["gene_name"], "gene_1")
  assert_eq(output.experiment.col_data[0]["spot_id"], "spot_1")
}

///|
test "BANKSY: SpatialExperiment wrapper validates assay and coordinates" {
  let missing_assay = try {
    ignore(
      @src.banksy_spatial_experiment(
        @src.SpatialExperiment::new(),
        assay_name="missing",
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", [[1.0, 2.0]]))
  let missing_coordinates = try {
    ignore(
      @src.banksy_spatial_experiment(
        experiment,
        config=@src.BanksyComputeConfig::create(k_geom=[1]),
        n_components=1,
        cluster_config=@src.BanksyClusterConfig::create(n_clusters=1),
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(missing_assay)
  assert_true(missing_coordinates)
}

///|
test "BANKSY: SpatialExperiment wrapper validates grouping field" {
  let failed = try {
    ignore(
      @src.banksy_spatial_experiment(
        banksy_test_experiment(),
        config=banksy_test_config(@src.banksy_knn_median()),
        n_components=3,
        cluster_config=banksy_test_cluster_config(),
        max_harmonic=1,
        group_key="missing",
        smooth=false,
        smoothing_k=3,
      ),
    )
    false
  } catch {
    BanksyError(_) => true
  }
  assert_true(failed)
}
