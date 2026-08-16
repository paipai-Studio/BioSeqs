///|
fn sling_adv_test_coordinates() -> Array[Array[Double]] {
  [
    [-0.1, 0.0],
    [0.0, 0.1],
    [0.1, -0.1],
    [0.9, 0.0],
    [1.0, 0.1],
    [1.1, -0.1],
    [1.9, 0.0],
    [2.0, 0.1],
    [2.1, -0.1],
    [2.9, 0.9],
    [3.0, 1.0],
    [3.1, 1.1],
    [2.9, -0.9],
    [3.0, -1.0],
    [3.1, -1.1],
  ]
}

///|
fn sling_adv_test_labels() -> Array[String] {
  ["A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E"]
}

///|
fn sling_adv_test_config() -> @src.SlingshotAdvancedConfig raise @src.SlingshotAdvancedError {
  @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    end_clusters=["D", "E"],
    distance=@src.slingshot_center_euclidean(),
    extension=@src.slingshot_extend_none(),
    max_iterations=8,
    tolerance=1.0e-4,
    smoother_span=0.4,
    curve_points=30,
  )
}

///|
fn sling_adv_test_result() -> @src.SlingshotAdvancedResult raise @src.SlingshotAdvancedError {
  @src.slingshot_advanced(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
}

///|
fn sling_adv_close(left : Double, right : Double, tolerance : Double) -> Bool {
  (left - right).abs() <= tolerance
}

///|
test "slingshot advanced distance accessors" {
  assert_true(
    @src.slingshot_center_euclidean() != @src.slingshot_scaled_diagonal(),
  )
  assert_true(@src.slingshot_scaled_diagonal() != @src.slingshot_scaled_full())
  assert_true(@src.slingshot_scaled_full() != @src.slingshot_center_euclidean())
}

///|
test "slingshot advanced extension accessors" {
  assert_true(@src.slingshot_extend_none() != @src.slingshot_extend_line())
  assert_true(@src.slingshot_extend_line() != @src.slingshot_extend_pc1())
  assert_true(@src.slingshot_extend_pc1() != @src.slingshot_extend_none())
}

///|
test "slingshot advanced config defaults" {
  let config = @src.SlingshotAdvancedConfig::create()
  assert_eq(config.distance, @src.slingshot_scaled_full())
  assert_eq(config.extension, @src.slingshot_extend_line())
  assert_eq(config.max_iterations, 15)
  assert_eq(config.curve_points, 100)
  assert_true(config.reweight)
  assert_true(config.reassign)
}

///|
test "slingshot advanced config explicit values" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["root"],
    end_clusters=["leaf"],
    use_median=true,
    omega=2.0,
    omega_scale=2.5,
    shrink=0.25,
    reweight=false,
    reassign=false,
    max_iterations=3,
    tolerance=0.01,
    smoother_span=0.5,
    curve_points=12,
    covariance_ridge=0.001,
  )
  assert_eq(config.start_clusters, ["root"])
  assert_eq(config.end_clusters, ["leaf"])
  assert_true(config.use_median)
  assert_eq(config.omega, 2.0)
  assert_eq(config.shrink, 0.25)
  assert_false(config.reweight)
}

///|
test "slingshot advanced rejects invalid omega" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(omega=-2.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid omega scale" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(omega_scale=0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid shrink" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(shrink=1.1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid iterations" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(max_iterations=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid tolerance" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(tolerance=0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid smoother span" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(smoother_span=1.1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid curve points" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(curve_points=1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects invalid covariance ridge" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(covariance_ridge=0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced rejects empty constraint names" {
  let mut failures = 0
  ignore(@src.SlingshotAdvancedConfig::create(start_clusters=[""])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.SlingshotAdvancedConfig::create(end_clusters=[""])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "slingshot advanced validates coordinate shape" {
  let mut failures = 0
  ignore(@src.slingshot_infer_lineages([], [], sling_adv_test_config())) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.slingshot_infer_lineages(
      [[0.0], [1.0, 2.0]],
      ["A", "B"],
      sling_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "slingshot advanced validates label count" {
  let mut failures = 0
  ignore(
    @src.slingshot_infer_lineages(
      [[0.0], [1.0]],
      ["A"],
      sling_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced validates finite coordinates" {
  let mut failures = 0
  ignore(
    @src.slingshot_infer_lineages(
      [[0.0], [0.0 / 0.0]],
      ["A", "B"],
      sling_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced validates start and end clusters" {
  let mut failures = 0
  let start_config = @src.SlingshotAdvancedConfig::create(start_clusters=[
    "missing",
  ])
  ignore(
    @src.slingshot_infer_lineages([[0.0], [1.0]], ["A", "B"], start_config),
  ) catch {
    _ => failures = failures + 1
  }
  let end_config = @src.SlingshotAdvancedConfig::create(end_clusters=["missing"])
  ignore(@src.slingshot_infer_lineages([[0.0], [1.0]], ["A", "B"], end_config)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "slingshot advanced hard labels preserve first occurrence order" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_eq(model.cluster_names, ["A", "B", "C", "D", "E"])
}

///|
test "slingshot advanced cluster centers are correct" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_true(sling_adv_close(model.clusters[0].center[0], 0.0, 1.0e-12))
  assert_true(sling_adv_close(model.clusters[2].center[0], 2.0, 1.0e-12))
  assert_true(sling_adv_close(model.clusters[3].center[1], 1.0, 1.0e-12))
  assert_eq(model.clusters[0].size, 3.0)
}

///|
test "slingshot advanced covariance is finite and regularized" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_true(model.clusters[0].covariance[0][0] > 0.0)
  assert_true(model.clusters[0].covariance[1][1] > 0.0)
  assert_true(
    sling_adv_close(
      model.clusters[0].covariance[0][1],
      model.clusters[0].covariance[1][0],
      1.0e-12,
    ),
  )
}

///|
test "slingshot advanced euclidean cluster distances" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_true(sling_adv_close(model.distance_matrix[0][1], 1.0, 1.0e-12))
  assert_true(sling_adv_close(model.distance_matrix[2][3], 2.0.sqrt(), 1.0e-12))
  assert_eq(model.distance_matrix[1][0], model.distance_matrix[0][1])
}

///|
test "slingshot advanced scaled diagonal distance" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_scaled_diagonal(),
  )
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    config,
  )
  assert_true(model.distance_matrix[0][1] > 1.0)
}

///|
test "slingshot advanced scaled full distance" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_scaled_full(),
  )
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    config,
  )
  assert_true(model.distance_matrix[0][1] > 0.0)
  assert_eq(model.distance_matrix[0][0], 0.0)
}

///|
test "slingshot advanced constrained mst has expected size" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_eq(model.edges.length(), 4)
  assert_eq(model.roots.length(), 1)
}

///|
test "slingshot advanced forced endpoints are leaves" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  let degrees = Array::make(model.cluster_names.length(), 0)
  for edge in model.edges {
    degrees[edge.from] = degrees[edge.from] + 1
    degrees[edge.to] = degrees[edge.to] + 1
  }
  assert_eq(degrees[3], 1)
  assert_eq(degrees[4], 1)
}

///|
test "slingshot advanced root honors start cluster" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_eq(model.roots, [0])
}

///|
test "slingshot advanced enumerates root to leaf lineages" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_eq(model.lineages.length(), 2)
  assert_eq(model.lineages[0][0], 0)
  assert_eq(model.lineages[1][0], 0)
  let ends = [model.lineages[0][3], model.lineages[1][3]]
  ends.sort()
  assert_eq(ends, [3, 4])
}

///|
test "slingshot advanced initial weights share common trunk" {
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    sling_adv_test_config(),
  )
  assert_eq(model.initial_weights[0], [1.0, 1.0])
  assert_eq(model.initial_weights[9], [1.0, 0.0])
  assert_eq(model.initial_weights[12], [0.0, 1.0])
}

///|
test "slingshot advanced treats minus one as unclustered" {
  let coordinates = [[0.0], [1.0], [2.0]]
  let labels = ["A", "-1", "B"]
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    curve_points=10,
    max_iterations=2,
  )
  let model = @src.slingshot_infer_lineages(coordinates, labels, config)
  assert_eq(model.cluster_names, ["A", "B"])
  assert_eq(model.cluster_weights[1], [0.0, 0.0])
  assert_eq(model.initial_weights[1], [0.0])
}

///|
test "slingshot advanced weighted memberships normalize rows" {
  let model = @src.slingshot_infer_weighted_lineages(
    [[0.0], [1.0], [2.0]],
    [[2.0, 0.0], [1.0, 1.0], [0.0, 3.0]],
    ["A", "B"],
    @src.SlingshotAdvancedConfig::create(
      start_clusters=["A"],
      distance=@src.slingshot_center_euclidean(),
    ),
  )
  assert_eq(model.cluster_weights[0], [1.0, 0.0])
  assert_eq(model.cluster_weights[1], [0.5, 0.5])
  assert_eq(model.cluster_weights[2], [0.0, 1.0])
}

///|
test "slingshot advanced rejects invalid weighted memberships" {
  let mut failures = 0
  let config = @src.SlingshotAdvancedConfig::create()
  ignore(
    @src.slingshot_infer_weighted_lineages(
      [[0.0], [1.0]],
      [[1.0], [-1.0]],
      ["A"],
      config,
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.slingshot_infer_weighted_lineages(
      [[0.0], [1.0]],
      [[1.0, 0.0], [1.0]],
      ["A", "B"],
      config,
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "slingshot advanced rejects empty weighted clusters" {
  let mut failures = 0
  ignore(
    @src.slingshot_infer_weighted_lineages(
      [[0.0], [1.0]],
      [[1.0, 0.0], [1.0, 0.0]],
      ["A", "B"],
      @src.SlingshotAdvancedConfig::create(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced fixed omega creates forest" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    omega=0.5,
  )
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    config,
  )
  assert_eq(model.edges.length(), 0)
  assert_eq(model.roots.length(), 5)
  assert_eq(model.lineages.length(), 5)
}

///|
test "slingshot advanced automatic omega uses mst median" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    automatic_omega=true,
    omega_scale=1.0,
  )
  let model = @src.slingshot_infer_lineages(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    config,
  )
  assert_true(model.omega_threshold > 0.0)
  assert_true(model.omega_threshold < 2.0)
}

///|
test "slingshot advanced single cluster lineage" {
  let config = @src.SlingshotAdvancedConfig::create(
    distance=@src.slingshot_center_euclidean(),
    curve_points=12,
    max_iterations=2,
  )
  let result = @src.slingshot_advanced(
    [[0.0, 0.0], [1.0, 0.1], [2.0, -0.1]],
    ["A", "A", "A"],
    config,
  )
  assert_eq(result.curves.length(), 1)
  assert_eq(result.curves[0].points.length(), 12)
  assert_eq(result.lineage_model.lineages[0], [0])
}

///|
test "slingshot advanced fits one curve per lineage" {
  let result = sling_adv_test_result()
  assert_eq(result.curves.length(), 2)
  assert_eq(result.curves[0].points.length(), 30)
  assert_eq(result.curves[1].points.length(), 30)
}

///|
test "slingshot advanced pseudotime matrix dimensions" {
  let result = sling_adv_test_result()
  assert_eq(result.pseudotime.length(), 15)
  assert_eq(result.pseudotime[0].length(), 2)
  assert_eq(result.weights.length(), 15)
  assert_eq(result.weights[0].length(), 2)
}

///|
test "slingshot advanced pseudotime begins near root" {
  let result = sling_adv_test_result()
  let root_time = result.average_pseudotime[0]
  assert_true(root_time < result.average_pseudotime[7])
  assert_true(root_time < result.average_pseudotime[10])
  assert_true(root_time < result.average_pseudotime[13])
}

///|
test "slingshot advanced optional pseudotime follows weights" {
  let result = sling_adv_test_result()
  for cell in 0..<result.weights.length() {
    for lineage in 0..<result.curves.length() {
      match result.pseudotime[cell][lineage] {
        Some(_) => assert_true(result.weights[cell][lineage] > 0.0)
        None => assert_eq(result.weights[cell][lineage], 0.0)
      }
    }
  }
}

///|
test "slingshot advanced reports bounded iterations" {
  let result = sling_adv_test_result()
  assert_true(result.iterations >= 1)
  assert_true(result.iterations <= result.config.max_iterations)
  assert_true(result.total_distance >= 0.0)
}

///|
test "slingshot advanced no reweight preserves memberships" {
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    end_clusters=["D", "E"],
    distance=@src.slingshot_center_euclidean(),
    extension=@src.slingshot_extend_none(),
    reweight=false,
    reassign=false,
    max_iterations=2,
    curve_points=20,
  )
  let result = @src.slingshot_advanced(
    sling_adv_test_coordinates(),
    sling_adv_test_labels(),
    config,
  )
  assert_eq(result.weights, result.lineage_model.initial_weights)
}

///|
test "slingshot advanced shrink aligns shared origin" {
  let result = sling_adv_test_result()
  let first = result.curves[0].points[0]
  let second = result.curves[1].points[0]
  assert_true(sling_adv_close(first[0], second[0], 1.0e-8))
  assert_true(sling_adv_close(first[1], second[1], 1.0e-8))
}

///|
test "slingshot advanced supports all extension modes" {
  let lengths : Array[Double] = []
  for
    extension in [
      @src.slingshot_extend_none(),
      @src.slingshot_extend_line(),
      @src.slingshot_extend_pc1(),
    ] {
    let config = @src.SlingshotAdvancedConfig::create(
      start_clusters=["A"],
      end_clusters=["D", "E"],
      distance=@src.slingshot_center_euclidean(),
      extension~,
      reweight=false,
      reassign=false,
      max_iterations=1,
      curve_points=15,
    )
    let result = @src.slingshot_advanced(
      sling_adv_test_coordinates(),
      sling_adv_test_labels(),
      config,
    )
    lengths.push(result.curves[0].length)
  }
  assert_eq(lengths.length(), 3)
  assert_true(lengths[0] > 0.0)
  assert_true(lengths[1] > 0.0)
  assert_true(lengths[2] > 0.0)
}

///|
test "slingshot advanced probability weights sum to one" {
  let probabilities = @src.slingshot_curve_weight_probabilities(
    sling_adv_test_result(),
  )
  for row in probabilities {
    let mut total = 0.0
    for value in row {
      total = total + value
    }
    assert_true(sling_adv_close(total, 1.0, 1.0e-10))
  }
}

///|
test "slingshot advanced branch ids use one based labels" {
  let result = sling_adv_test_result()
  let ids = @src.slingshot_branch_ids(result)
  assert_true(ids[0].contains("1"))
  assert_true(ids[0].contains("2"))
  assert_true(ids[10].contains("1"))
  assert_true(ids[13].contains("2"))
}

///|
test "slingshot advanced branch ids validate threshold" {
  let mut failures = 0
  ignore(@src.slingshot_branch_ids(sling_adv_test_result(), threshold=1.1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced prediction dimensions" {
  let prediction = @src.slingshot_predict(sling_adv_test_result(), [
    [0.0, 0.0],
    [3.0, 1.0],
    [3.0, -1.0],
  ])
  assert_eq(prediction.pseudotime.length(), 3)
  assert_eq(prediction.weights[0].length(), 2)
  assert_eq(prediction.distances[0].length(), 2)
  assert_eq(prediction.average_pseudotime.length(), 3)
}

///|
test "slingshot advanced prediction maps branch endpoints" {
  let prediction = @src.slingshot_predict(sling_adv_test_result(), [
    [3.0, 1.0],
    [3.0, -1.0],
  ])
  assert_true(prediction.branch_ids[0].contains("1"))
  assert_true(prediction.branch_ids[1].contains("2"))
}

///|
test "slingshot advanced prediction assigns distant cells" {
  let prediction = @src.slingshot_predict(sling_adv_test_result(), [
    [100.0, 100.0],
  ])
  let mut total = 0.0
  for weight in prediction.weights[0] {
    total = total + weight
  }
  assert_true(total > 0.0)
  assert_true(prediction.branch_ids[0].length() > 0)
}

///|
test "slingshot advanced prediction validates dimensions" {
  let mut failures = 0
  ignore(@src.slingshot_predict(sling_adv_test_result(), [[0.0, 1.0, 2.0]])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 1)
}

///|
test "slingshot advanced weighted full pipeline" {
  let coordinates = [[0.0], [0.5], [1.0], [1.5], [2.0]]
  let weights = [[1.0, 0.0], [0.75, 0.25], [0.5, 0.5], [0.25, 0.75], [0.0, 1.0]]
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    extension=@src.slingshot_extend_none(),
    curve_points=15,
    max_iterations=3,
  )
  let result = @src.slingshot_advanced_weighted(
    coordinates,
    weights,
    ["A", "B"],
    config,
  )
  assert_eq(result.curves.length(), 1)
  assert_true(result.average_pseudotime[0] < result.average_pseudotime[4])
}

///|
test "slingshot advanced summary reports model dimensions" {
  let summary = @src.slingshot_advanced_summary(sling_adv_test_result())
  assert_true(summary.contains("Slingshot 2.21.0"))
  assert_true(summary.contains("cells: 15"))
  assert_true(summary.contains("clusters: 5"))
  assert_true(summary.contains("lineages: 2"))
}

///|
test "slingshot advanced sce writes trajectory outputs" {
  let experiment = @src.SingleCellExperiment::new(
    [
      [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
        15.0,
      ],
    ],
    ["gene"],
    [
      "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12",
      "c13", "c14", "c15",
    ],
  )
  experiment.reduced_dims["PCA"] = sling_adv_test_coordinates()
  experiment.col_data["cluster"] = sling_adv_test_labels()
  let output = @src.slingshot_advanced_sce(experiment, sling_adv_test_config())
  assert_eq(output.result.curves.length(), 2)
  assert_eq(output.experiment.col_data["slingshot.branch"].length(), 15)
  assert_eq(
    output.experiment.reduced_dims["slingshot.pseudotime"][0].length(),
    2,
  )
  assert_eq(output.experiment.reduced_dims["slingshot.weights"][0].length(), 2)
}

///|
test "slingshot advanced sce preserves input immutability" {
  let experiment = @src.SingleCellExperiment::new(
    [
      [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
        15.0,
      ],
    ],
    ["gene"],
    [
      "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12",
      "c13", "c14", "c15",
    ],
  )
  experiment.reduced_dims["PCA"] = sling_adv_test_coordinates()
  experiment.col_data["cluster"] = sling_adv_test_labels()
  experiment.metadata["source"] = "original"
  let output = @src.slingshot_advanced_sce(experiment, sling_adv_test_config())
  output.experiment.assays["counts"][0][0] = 99.0
  output.experiment.reduced_dims["PCA"][0][0] = 99.0
  output.experiment.col_data["cluster"][0] = "changed"
  assert_eq(experiment.assays["counts"][0][0], 1.0)
  assert_eq(experiment.reduced_dims["PCA"][0][0], -0.1)
  assert_eq(experiment.col_data["cluster"][0], "A")
  assert_eq(experiment.metadata["source"], "original")
}

///|
test "slingshot advanced sce recursively copies alternatives" {
  let experiment = @src.SingleCellExperiment::new([[1.0, 2.0]], ["gene"], [
    "c1", "c2",
  ])
  experiment.reduced_dims["PCA"] = [[0.0], [1.0]]
  experiment.col_data["cluster"] = ["A", "B"]
  experiment.alternative_experiments["alt"] = @src.SingleCellExperiment::new(
    [[3.0, 4.0]],
    ["feature"],
    ["c1", "c2"],
  )
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    curve_points=10,
    max_iterations=2,
  )
  let output = @src.slingshot_advanced_sce(experiment, config)
  output.experiment.alternative_experiments["alt"].assays["counts"][0][0] = 8.0
  assert_eq(
    experiment.alternative_experiments["alt"].assays["counts"][0][0],
    3.0,
  )
}

///|
test "slingshot advanced sce supports custom output prefix" {
  let experiment = @src.SingleCellExperiment::new([[1.0, 2.0]], ["gene"], [
    "c1", "c2",
  ])
  experiment.reduced_dims["UMAP"] = [[0.0], [1.0]]
  experiment.col_data["group"] = ["A", "B"]
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    distance=@src.slingshot_center_euclidean(),
    curve_points=10,
    max_iterations=2,
  )
  let output = @src.slingshot_advanced_sce(
    experiment,
    config,
    reduced_dim_name="UMAP",
    cluster_column="group",
    output_prefix="trajectory",
  )
  assert_true(output.experiment.col_data.contains("trajectory.branch"))
  assert_true(output.experiment.reduced_dims.contains("trajectory.pseudotime"))
  assert_eq(output.experiment.metadata["trajectory.reduced_dim"], "UMAP")
}

///|
test "slingshot advanced sce validates names and inputs" {
  let experiment = @src.SingleCellExperiment::new([[1.0, 2.0]], ["gene"], [
    "c1", "c2",
  ])
  let config = @src.SlingshotAdvancedConfig::create()
  let mut failures = 0
  ignore(@src.slingshot_advanced_sce(experiment, config, reduced_dim_name="")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.slingshot_advanced_sce(experiment, config)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}
