///|
fn flowsom_test_close(
  left : Double,
  right : Double,
  tolerance : Double,
) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn flowsom_test_some(value : Double?) -> Double {
  match value {
    Some(actual) => actual
    None => abort("expected a numeric value")
  }
}

///|
fn flowsom_test_unique_count(values : Array[Int]) -> Int {
  let unique : Array[Int] = []
  for value in values {
    if !unique.contains(value) {
      unique.push(value)
    }
  }
  unique.length()
}

///|
fn flowsom_test_sum(values : Array[Int]) -> Int {
  let mut total = 0
  for value in values {
    total = total + value
  }
  total
}

///|
fn flowsom_test_data() -> Array[Array[Double]] {
  [
    [0.0, 0.0],
    [0.1, 0.2],
    [0.2, 0.1],
    [0.3, 0.2],
    [9.8, 10.0],
    [10.0, 9.9],
    [10.1, 10.2],
    [10.3, 10.1],
  ]
}

///|
fn flowsom_test_config(
  metaclusters? : Int = 0,
  seed? : Int = 31,
) -> @src.FlowSomConfig {
  @src.FlowSomConfig::create(
    xdim=2,
    ydim=2,
    rlen=3,
    mst_runs=2,
    alpha_start=0.1,
    alpha_end=0.02,
    radius_start=1.0,
    radius_end=0.0,
    metaclusters~,
    meta_starts=3,
    seed~,
  ) catch {
    _ => abort("valid FlowSOM test configuration should build")
  }
}

///|
fn flowsom_test_model(
  metaclusters? : Int = 0,
  seed? : Int = 31,
) -> @src.FlowSomModel {
  @src.flowsom_train(
    flowsom_test_data(),
    flowsom_test_config(metaclusters~, seed~),
    marker_names=["CD3", "CD19"],
  ) catch {
    _ => abort("valid FlowSOM test model should train")
  }
}

///|
fn flowsom_test_one_node_model(metaclusters? : Int = 0) -> @src.FlowSomModel {
  let config = @src.FlowSomConfig::create(
    xdim=1,
    ydim=1,
    rlen=2,
    alpha_start=0.1,
    alpha_end=0.02,
    radius_start=0.0,
    radius_end=0.0,
    metaclusters~,
    seed=41,
  ) catch {
    _ => abort("valid one-node FlowSOM configuration should build")
  }
  @src.flowsom_train([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], config, marker_names=[
    "A", "B",
  ]) catch {
    _ => abort("valid one-node FlowSOM model should train")
  }
}

///|
fn flowsom_test_sce() -> @src.SingleCellExperiment {
  let experiment = @src.SingleCellExperiment::new(
    [
      [0.0, 0.1, 0.2, 10.0, 10.1, 10.2],
      [1.0, 1.1, 1.2, 11.0, 11.1, 11.2],
      [2.0, 2.1, 2.2, 12.0, 12.1, 12.2],
    ],
    ["CD3", "CD19", "CD45"],
    ["C1", "C2", "C3", "C4", "C5", "C6"],
  )
  experiment.row_data["symbol"] = ["T", "B", "pan"]
  experiment.col_data["batch"] = ["A", "A", "A", "B", "B", "B"]
  experiment.reduced_dims["PCA"] = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.1],
    [10.0, 10.0],
    [10.1, 10.0],
    [10.2, 10.1],
  ]
  experiment.metadata["project"] = "flowsom-test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[5.0, 6.0, 7.0, 8.0, 9.0, 10.0]],
    ["Spike1"],
    ["C1", "C2", "C3", "C4", "C5", "C6"],
  )
  experiment
}

///|
test "FlowSOM: distance and initialization accessors expose enum values" {
  assert_true(
    @src.flowsom_manhattan_distance() != @src.flowsom_euclidean_distance(),
  )
  assert_true(
    @src.flowsom_euclidean_distance() != @src.flowsom_chebyshev_distance(),
  )
  assert_true(
    @src.flowsom_chebyshev_distance() != @src.flowsom_cosine_distance(),
  )
  assert_true(
    @src.flowsom_random_initialization() != @src.flowsom_kwsp_initialization(),
  )
  assert_true(
    @src.flowsom_kwsp_initialization() != @src.flowsom_pca_initialization(),
  )
}

///|
test "FlowSOM: configuration defaults match the upstream workflow" {
  let config = @src.FlowSomConfig::create() catch {
    _ => abort("default FlowSOM configuration should build")
  }
  assert_eq(config.xdim, 10)
  assert_eq(config.ydim, 10)
  assert_eq(config.rlen, 10)
  assert_eq(config.mst_runs, 1)
  assert_eq(config.alpha_start, 0.05)
  assert_eq(config.alpha_end, 0.01)
  assert_eq(config.radius_start, -1.0)
  assert_eq(config.radius_end, 0.0)
  assert_true(config.distance == @src.flowsom_euclidean_distance())
  assert_true(config.initialization == @src.flowsom_kwsp_initialization())
  assert_eq(config.metaclusters, 0)
  assert_eq(config.seed, 1)
}

///|
test "FlowSOM: configuration preserves explicit controls" {
  let config = @src.FlowSomConfig::create(
    xdim=3,
    ydim=2,
    rlen=7,
    mst_runs=3,
    alpha_start=0.2,
    alpha_end=0.03,
    radius_start=2.0,
    radius_end=0.5,
    distance=@src.flowsom_cosine_distance(),
    initialization=@src.flowsom_pca_initialization(),
    importance=[2.0, 0.5],
    metaclusters=3,
    meta_starts=4,
    outlier_mad=2.5,
    seed=17,
  ) catch {
    _ => abort("explicit FlowSOM configuration should build")
  }
  assert_eq(config.xdim, 3)
  assert_eq(config.ydim, 2)
  assert_eq(config.rlen, 7)
  assert_eq(config.mst_runs, 3)
  assert_eq(config.importance, [2.0, 0.5])
  assert_eq(config.metaclusters, 3)
  assert_eq(config.meta_starts, 4)
  assert_eq(config.outlier_mad, 2.5)
  assert_eq(config.seed, 17)
}

///|
test "FlowSOM: configuration rejects invalid grid and run controls" {
  let mut failures = 0
  ignore(@src.FlowSomConfig::create(xdim=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(ydim=-1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(rlen=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(mst_runs=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "FlowSOM: configuration rejects invalid learning schedules" {
  let mut failures = 0
  ignore(@src.FlowSomConfig::create(alpha_start=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(alpha_start=0.01, alpha_end=0.02)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(alpha_end=-0.01)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(alpha_start=1.0e301)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "FlowSOM: configuration rejects invalid radius schedules" {
  let mut failures = 0
  ignore(@src.FlowSomConfig::create(radius_start=-2.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(radius_end=-0.1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(radius_start=1.0, radius_end=2.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: configuration rejects invalid meta, outlier and seed controls" {
  let mut failures = 0
  ignore(@src.FlowSomConfig::create(xdim=2, ydim=2, metaclusters=5)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(meta_starts=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(outlier_mad=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(seed=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.FlowSomConfig::create(importance=[1.0, 0.0])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 5)
}

///|
test "FlowSOM: four distances follow their definitions" {
  let left = [1.0, -2.0, 3.0]
  let right = [4.0, 2.0, 3.0]
  assert_eq(
    @src.flowsom_distance(left, right, @src.flowsom_manhattan_distance()),
    7.0,
  )
  assert_eq(
    @src.flowsom_distance(left, right, @src.flowsom_euclidean_distance()),
    5.0,
  )
  assert_eq(
    @src.flowsom_distance(left, right, @src.flowsom_chebyshev_distance()),
    4.0,
  )
  let cosine = @src.flowsom_distance(
    [1.0, 0.0],
    [0.0, 1.0],
    @src.flowsom_cosine_distance(),
  )
  assert_true(flowsom_test_close(cosine, 1.0, 1.0e-12))
}

///|
test "FlowSOM: cosine distance handles zero vectors" {
  let cosine = @src.flowsom_cosine_distance()
  assert_eq(@src.flowsom_distance([0.0, 0.0], [0.0, 0.0], cosine), 0.0)
  assert_eq(@src.flowsom_distance([0.0, 0.0], [1.0, 0.0], cosine), 1.0)
  assert_true(
    flowsom_test_close(
      @src.flowsom_distance([1.0, 1.0], [2.0, 2.0], cosine),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "FlowSOM: distance validates vector dimensions and values" {
  let mut failures = 0
  ignore(@src.flowsom_distance([], [], @src.flowsom_euclidean_distance())) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.flowsom_distance([1.0], [1.0, 2.0], @src.flowsom_euclidean_distance()),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.flowsom_distance([1.0e301], [0.0], @src.flowsom_euclidean_distance()),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: KWSP initialization repeatedly selects the farthest cell" {
  let codes = @src.flowsom_initialize_kwsp([[0.0], [1.0], [10.0]], 3, seed=1) catch {
    _ => abort("KWSP initialization should work")
  }
  assert_eq(codes, [[0.0], [10.0], [1.0]])
}

///|
test "FlowSOM: KWSP initialization is reproducible and copies rows" {
  let data = [[0.0, 0.0], [1.0, 1.0], [5.0, 5.0], [10.0, 10.0]]
  let first = @src.flowsom_initialize_kwsp(data, 3, seed=23) catch {
    _ => abort("first KWSP initialization should work")
  }
  let second = @src.flowsom_initialize_kwsp(data, 3, seed=23) catch {
    _ => abort("second KWSP initialization should work")
  }
  assert_eq(first, second)
  first[0][0] = 999.0
  assert_true(data[0][0] != 999.0)
  assert_true(data[1][0] != 999.0)
  assert_true(data[2][0] != 999.0)
  assert_true(data[3][0] != 999.0)
}

///|
test "FlowSOM: KWSP initialization validates nodes and seed" {
  let data = [[0.0], [1.0]]
  let mut failures = 0
  ignore(@src.flowsom_initialize_kwsp(data, 0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_initialize_kwsp(data, 3)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_initialize_kwsp(data, 1, seed=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: one-node PCA initialization returns the column means" {
  let codes = @src.flowsom_initialize_pca(
    [[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]],
    1,
    1,
  ) catch {
    _ => abort("one-node PCA initialization should work")
  }
  assert_eq(codes.length(), 1)
  assert_true(flowsom_test_close(codes[0][0], 3.0, 1.0e-12))
  assert_true(flowsom_test_close(codes[0][1], 6.0, 1.0e-12))
}

///|
test "FlowSOM: PCA grid is symmetric around the data mean" {
  let codes = @src.flowsom_initialize_pca(
    [[0.0, 0.0], [2.0, 1.0], [4.0, 2.0]],
    2,
    1,
  ) catch {
    _ => abort("PCA grid initialization should work")
  }
  assert_eq(codes.length(), 2)
  assert_eq(codes[0].length(), 2)
  assert_true(
    flowsom_test_close((codes[0][0] + codes[1][0]) / 2.0, 2.0, 1.0e-9),
  )
  assert_true(
    flowsom_test_close((codes[0][1] + codes[1][1]) / 2.0, 1.0, 1.0e-9),
  )
}

///|
test "FlowSOM: PCA initialization validates grid dimensions" {
  let mut failures = 0
  ignore(@src.flowsom_initialize_pca([[1.0]], 0, 1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_initialize_pca([[1.0]], 1, 0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: MST connects a one-dimensional chain with minimum weight" {
  let edges = @src.flowsom_build_mst([[0.0], [1.0], [3.0], [7.0]]) catch {
    _ => abort("MST construction should work")
  }
  assert_eq(edges.length(), 3)
  assert_eq(edges[0].from, 0)
  assert_eq(edges[0].to, 1)
  assert_eq(edges[0].weight, 1.0)
  assert_eq(edges[1].from, 1)
  assert_eq(edges[1].to, 2)
  assert_eq(edges[1].weight, 2.0)
  assert_eq(edges[2].from, 2)
  assert_eq(edges[2].to, 3)
  assert_eq(edges[2].weight, 4.0)
}

///|
test "FlowSOM: MST tie breaking is stable by node index" {
  let edges = @src.flowsom_build_mst([[0.0], [1.0], [-1.0]]) catch {
    _ => abort("tied MST construction should work")
  }
  assert_eq(edges.length(), 2)
  assert_eq(edges[0].from, 0)
  assert_eq(edges[0].to, 1)
  assert_eq(edges[1].from, 0)
  assert_eq(edges[1].to, 2)
}

///|
test "FlowSOM: one-node MST has no edges" {
  let edges = @src.flowsom_build_mst([[2.0, 3.0]]) catch {
    _ => abort("one-node MST construction should work")
  }
  assert_eq(edges.length(), 0)
}

///|
test "FlowSOM: MST topology reports unweighted path lengths" {
  let distances = @src.flowsom_mst_distances([[0.0], [1.0], [3.0], [7.0]]) catch {
    _ => abort("MST topology should work")
  }
  assert_eq(distances.length(), 4)
  assert_eq(distances[0], [0, 1, 2, 3])
  assert_eq(distances[3], [3, 2, 1, 0])
}

///|
test "FlowSOM: MST validates rectangular finite codebooks" {
  let mut failures = 0
  ignore(@src.flowsom_build_mst([])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_build_mst([[1.0], [1.0, 2.0]])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_build_mst([[1.0e301]])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: data mapping returns nearest nodes and distances" {
  let mapping = @src.flowsom_map_data([[0.0, 0.0], [10.0, 10.0]], [
    [1.0, 0.0],
    [9.0, 10.0],
  ]) catch {
    _ => abort("FlowSOM mapping should work")
  }
  assert_eq(mapping.clusters, [0, 1])
  assert_eq(mapping.distances, [1.0, 1.0])
}

///|
test "FlowSOM: data mapping resolves nearest-node ties by index" {
  let mapping = @src.flowsom_map_data([[0.0], [2.0]], [[1.0]]) catch {
    _ => abort("tied FlowSOM mapping should work")
  }
  assert_eq(mapping.clusters, [0])
  assert_eq(mapping.distances, [1.0])
}

///|
test "FlowSOM: data mapping supports non-Euclidean distance" {
  let mapping = @src.flowsom_map_data(
    [[0.0, 0.0], [4.0, 4.0]],
    [[3.0, 0.0]],
    distance=@src.flowsom_manhattan_distance(),
  ) catch {
    _ => abort("Manhattan FlowSOM mapping should work")
  }
  assert_eq(mapping.clusters, [0])
  assert_eq(mapping.distances, [3.0])
}

///|
test "FlowSOM: data mapping validates marker dimensions" {
  let mut raised = false
  ignore(@src.flowsom_map_data([[0.0, 1.0]], [[0.0]])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "FlowSOM: training returns a complete grid, MST and mapping" {
  let model = flowsom_test_model()
  assert_eq(model.data.length(), 8)
  assert_eq(model.marker_names, ["CD3", "CD19"])
  assert_eq(model.codes.length(), 4)
  assert_eq(model.grid.length(), 4)
  assert_eq(model.grid[0], (0, 0))
  assert_eq(model.grid[3], (1, 1))
  assert_eq(model.mapping.clusters.length(), 8)
  assert_eq(model.mapping.distances.length(), 8)
  assert_eq(model.mst_edges.length(), 3)
  assert_eq(model.topology_distances.length(), 4)
  assert_eq(flowsom_test_sum(model.node_counts), 8)
}

///|
test "FlowSOM: fixed seeds make training reproducible" {
  let first = flowsom_test_model(seed=37)
  let second = flowsom_test_model(seed=37)
  assert_eq(first.codes, second.codes)
  assert_eq(first.mapping.clusters, second.mapping.clusters)
  assert_eq(first.mapping.distances, second.mapping.distances)
  assert_eq(first.iterations, second.iterations)
}

///|
test "FlowSOM: random initialization training path is available" {
  let config = @src.FlowSomConfig::create(
    xdim=2,
    ydim=2,
    rlen=2,
    initialization=@src.flowsom_random_initialization(),
    seed=43,
  ) catch {
    _ => abort("random FlowSOM configuration should build")
  }
  let model = @src.flowsom_train(flowsom_test_data(), config) catch {
    _ => abort("random-initialized FlowSOM should train")
  }
  assert_eq(model.codes.length(), 4)
  assert_eq(model.mapping.clusters.length(), 8)
}

///|
test "FlowSOM: PCA training permits more nodes than cells" {
  let config = @src.FlowSomConfig::create(
    xdim=2,
    ydim=2,
    rlen=1,
    initialization=@src.flowsom_pca_initialization(),
    seed=47,
  ) catch {
    _ => abort("PCA FlowSOM configuration should build")
  }
  let model = @src.flowsom_train([[0.0, 0.0], [1.0, 1.0]], config) catch {
    _ => abort("PCA FlowSOM with extra nodes should train")
  }
  assert_eq(model.codes.length(), 4)
  assert_eq(model.mapping.clusters.length(), 2)
}

///|
test "FlowSOM: multi-stage training tracks iterations and final topology" {
  let model = flowsom_test_model()
  assert_true(model.iterations > 0)
  assert_true(model.iterations <= 3 * 8 * 2)
  for node in 0..<model.codes.length() {
    assert_eq(model.topology_distances[node][node], 0)
  }
}

///|
test "FlowSOM: training generates marker names when omitted" {
  let model = @src.flowsom_train(flowsom_test_data(), flowsom_test_config()) catch {
    _ => abort("FlowSOM with generated marker names should train")
  }
  assert_eq(model.marker_names, ["marker.1", "marker.2"])
}

///|
test "FlowSOM: training validates marker names" {
  let config = flowsom_test_config()
  let mut failures = 0
  ignore(@src.flowsom_train(flowsom_test_data(), config, marker_names=["only"])) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.flowsom_train(flowsom_test_data(), config, marker_names=["A", "A"]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.flowsom_train(flowsom_test_data(), config, marker_names=["A", ""]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: training validates importance dimensions" {
  let config = @src.FlowSomConfig::create(xdim=1, ydim=1, importance=[1.0]) catch {
    _ => abort("importance configuration should build before data validation")
  }
  let mut raised = false
  ignore(@src.flowsom_train([[1.0, 2.0]], config)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "FlowSOM: training rejects empty, ragged and non-finite data" {
  let config = flowsom_test_config()
  let mut failures = 0
  ignore(@src.flowsom_train([], config)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_train([[1.0], [1.0, 2.0]], config)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_train([[1.0e301, 0.0]], config)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "FlowSOM: non-PCA initialization requires one cell per node" {
  let config = @src.FlowSomConfig::create(xdim=2, ydim=2) catch {
    _ => abort("valid FlowSOM configuration should build")
  }
  let mut raised = false
  ignore(@src.flowsom_train([[0.0], [1.0]], config)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "FlowSOM: node summaries use medians, sample SD, CV and R-style MAD" {
  let model = flowsom_test_one_node_model()
  assert_eq(model.node_counts, [3])
  assert_eq(model.node_percentages, [1.0])
  assert_eq(flowsom_test_some(model.node_medians[0][0]), 3.0)
  assert_eq(flowsom_test_some(model.node_medians[0][1]), 4.0)
  assert_true(
    flowsom_test_close(flowsom_test_some(model.node_sds[0][0]), 2.0, 1.0e-12),
  )
  assert_true(
    flowsom_test_close(
      flowsom_test_some(model.node_cvs[0][0]),
      2.0 / 3.0,
      1.0e-12,
    ),
  )
  assert_true(
    flowsom_test_close(flowsom_test_some(model.node_cvs[0][1]), 0.5, 1.0e-12),
  )
  assert_true(
    flowsom_test_close(
      flowsom_test_some(model.node_mads[0][0]),
      2.9652,
      1.0e-12,
    ),
  )
}

///|
test "FlowSOM: empty nodes retain missing optional statistics" {
  let config = @src.FlowSomConfig::create(
    xdim=3,
    ydim=1,
    rlen=1,
    initialization=@src.flowsom_pca_initialization(),
    seed=53,
  ) catch {
    _ => abort("PCA FlowSOM configuration should build")
  }
  let model = @src.flowsom_train([[0.0], [10.0]], config) catch {
    _ => abort("PCA FlowSOM should train")
  }
  let mut empty = -1
  for node in 0..<model.node_counts.length() {
    if model.node_counts[node] == 0 {
      empty = node
      break
    }
  }
  assert_true(empty >= 0)
  assert_true(model.node_medians[empty][0] is None)
  assert_true(model.node_sds[empty][0] is None)
  assert_true(model.node_mads[empty][0] is None)
}

///|
test "FlowSOM: outlier report is aligned to nodes and cells" {
  let model = flowsom_test_model()
  assert_eq(model.outliers.median_distances.length(), 4)
  assert_eq(model.outliers.mad_distances.length(), 4)
  assert_eq(model.outliers.thresholds.length(), 4)
  assert_eq(model.outliers.counts.length(), 4)
  assert_eq(model.outliers.maximum_distances.length(), 4)
  assert_eq(model.outliers.per_cell.length(), 8)
  let mut flagged = 0
  for value in model.outliers.per_cell {
    if value {
      flagged = flagged + 1
    }
  }
  assert_eq(flowsom_test_sum(model.outliers.counts), flagged)
}

///|
test "FlowSOM: quantization and topographic errors are bounded" {
  let model = flowsom_test_model()
  assert_true(model.quantization_error >= 0.0)
  assert_true(model.topographic_error >= 0.0)
  assert_true(model.topographic_error <= 1.0)
}

///|
test "FlowSOM: node positive percentages use strict cutoffs" {
  let model = flowsom_test_one_node_model()
  let percentages = @src.flowsom_node_positive_percentages(model, [2.0, 4.0]) catch {
    _ => abort("node positive percentages should work")
  }
  assert_true(
    flowsom_test_close(flowsom_test_some(percentages[0][0]), 2.0 / 3.0, 1.0e-12),
  )
  assert_true(
    flowsom_test_close(flowsom_test_some(percentages[0][1]), 1.0 / 3.0, 1.0e-12),
  )
}

///|
test "FlowSOM: node positive percentages validate cutoffs" {
  let model = flowsom_test_one_node_model()
  let mut failures = 0
  ignore(@src.flowsom_node_positive_percentages(model, [1.0])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_node_positive_percentages(model, [1.0, 1.0e301])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: one meta-cluster assigns every code to zero" {
  let labels = @src.flowsom_meta_cluster_codes([[0.0], [1.0], [10.0]], 1) catch {
    _ => abort("one FlowSOM meta-cluster should work")
  }
  assert_eq(labels, [0, 0, 0])
}

///|
test "FlowSOM: k-means meta-clustering separates distant code groups" {
  let labels = @src.flowsom_meta_cluster_codes(
    [[0.0], [0.1], [0.2], [10.0], [10.1], [10.2]],
    2,
    starts=4,
    seed=59,
  ) catch {
    _ => abort("FlowSOM meta-clustering should work")
  }
  assert_eq(flowsom_test_unique_count(labels), 2)
  assert_true(labels[0] == labels[1] && labels[1] == labels[2])
  assert_true(labels[3] == labels[4] && labels[4] == labels[5])
  assert_true(labels[0] != labels[3])
}

///|
test "FlowSOM: meta-clustering validates cluster controls" {
  let codes = [[0.0], [1.0]]
  let mut failures = 0
  ignore(@src.flowsom_meta_cluster_codes(codes, 0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_meta_cluster_codes(codes, 3)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_meta_cluster_codes(codes, 2, starts=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_meta_cluster_codes(codes, 2, seed=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "FlowSOM: automatic meta-clustering returns small maxima directly" {
  assert_eq(@src.flowsom_determine_meta_clusters([[0.0]], 1), 1)
  assert_eq(@src.flowsom_determine_meta_clusters([[0.0], [1.0]], 2), 2)
  let mut failures = 0
  ignore(@src.flowsom_determine_meta_clusters([[0.0], [1.0]], 2, starts=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_determine_meta_clusters([[0.0], [1.0]], 2, seed=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: automatic meta-clustering selects an internal elbow" {
  let selected = @src.flowsom_determine_meta_clusters(
    [[0.0], [0.1], [0.2], [10.0], [10.1], [10.2]],
    5,
    starts=3,
    seed=61,
  ) catch {
    _ => abort("automatic FlowSOM meta-clustering should work")
  }
  assert_true(selected >= 2)
  assert_true(selected < 5)
}

///|
test "FlowSOM: automatic meta-clustering validates its maximum" {
  let mut failures = 0
  ignore(@src.flowsom_determine_meta_clusters([[0.0]], 0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_determine_meta_clusters([[0.0]], 2)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: configured meta-clusters propagate to cells" {
  let model = flowsom_test_model(metaclusters=2)
  assert_eq(model.meta_clusters.length(), 4)
  assert_eq(model.cell_meta_clusters.length(), 8)
  assert_eq(flowsom_test_unique_count(model.meta_clusters), 2)
  assert_eq(
    @src.flowsom_meta_medians(model).length(),
    @src.flowsom_meta_counts(model).length(),
  )
  for cell in 0..<model.mapping.clusters.length() {
    assert_eq(
      model.cell_meta_clusters[cell],
      model.meta_clusters[model.mapping.clusters[cell]],
    )
  }
}

///|
test "FlowSOM: meta counts and medians aggregate original cells" {
  let model = flowsom_test_one_node_model(metaclusters=1)
  assert_eq(@src.flowsom_meta_counts(model), [3])
  let medians = @src.flowsom_meta_medians(model)
  assert_eq(medians.length(), 1)
  assert_eq(flowsom_test_some(medians[0][0]), 3.0)
  assert_eq(flowsom_test_some(medians[0][1]), 4.0)
}

///|
test "FlowSOM: meta summaries reject models without meta-clusters" {
  let model = flowsom_test_one_node_model()
  let mut failures = 0
  ignore(@src.flowsom_meta_counts(model)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_meta_medians(model)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: post-hoc meta-clustering leaves the original model unchanged" {
  let model = flowsom_test_model()
  let updated = @src.flowsom_with_metaclusters(model, 2, starts=3, seed=67) catch {
    _ => abort("post-hoc FlowSOM meta-clustering should work")
  }
  assert_eq(model.meta_clusters.length(), 0)
  assert_eq(updated.meta_clusters.length(), 4)
  assert_eq(updated.cell_meta_clusters.length(), 8)
  assert_eq(updated.config.metaclusters, 2)
  updated.codes[0][0] = 999.0
  assert_true(model.codes[0][0] != 999.0)
}

///|
test "FlowSOM: automatic post-hoc meta-clustering updates the model" {
  let model = flowsom_test_model()
  let updated = @src.flowsom_with_automatic_metaclusters(
    model,
    3,
    starts=3,
    seed=71,
  ) catch {
    _ => abort("automatic post-hoc meta-clustering should work")
  }
  assert_true(updated.config.metaclusters >= 2)
  assert_true(updated.config.metaclusters <= 3)
  assert_eq(updated.meta_clusters.length(), 4)
  assert_eq(updated.cell_meta_clusters.length(), 8)
}

///|
test "FlowSOM: new data mapping applies node and meta-cluster labels" {
  let model = flowsom_test_model(metaclusters=2)
  let mapped = @src.flowsom_map_new(model, [[0.05, 0.1], [10.2, 10.0]]) catch {
    _ => abort("new FlowSOM data should map")
  }
  assert_eq(mapped.mapping.clusters.length(), 2)
  assert_eq(mapped.mapping.distances.length(), 2)
  assert_eq(mapped.meta_clusters.length(), 2)
  assert_eq(mapped.outliers.length(), 2)
  for cell in 0..<2 {
    assert_eq(
      mapped.meta_clusters[cell],
      model.meta_clusters[mapped.mapping.clusters[cell]],
    )
  }
}

///|
test "FlowSOM: new data mapping validates marker dimensions" {
  let model = flowsom_test_model()
  let mut raised = false
  ignore(@src.flowsom_map_new(model, [[1.0]])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "FlowSOM: marker outliers classify low and high observations" {
  let model = flowsom_test_one_node_model()
  let outliers = @src.flowsom_marker_outliers(
    model,
    [[0.0, 4.0], [3.0, 4.0], [7.0, 8.0]],
    mad_allowed=0.0,
  ) catch {
    _ => abort("marker-level FlowSOM outliers should work")
  }
  assert_eq(outliers, [[-1, 0], [0, 0], [1, 1]])
}

///|
test "FlowSOM: marker outliers validate dimensions and MAD control" {
  let model = flowsom_test_one_node_model()
  let mut failures = 0
  ignore(@src.flowsom_marker_outliers(model, [[1.0]])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_marker_outliers(model, [[1.0, 2.0]], mad_allowed=-1.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: weighted and unweighted purity follow contingency counts" {
  let real = [0, 0, 1, 1]
  let predicted = [0, 0, 0, 1]
  let weighted = @src.flowsom_purity(real, predicted) catch {
    _ => abort("weighted FlowSOM purity should work")
  }
  let unweighted = @src.flowsom_purity(real, predicted, weighted=false) catch {
    _ => abort("unweighted FlowSOM purity should work")
  }
  assert_true(flowsom_test_close(weighted.mean, 0.75, 1.0e-12))
  assert_true(flowsom_test_close(weighted.worst, 2.0 / 3.0, 1.0e-12))
  assert_eq(weighted.below_075, 1)
  assert_true(flowsom_test_close(unweighted.mean, 5.0 / 6.0, 1.0e-12))
}

///|
test "FlowSOM: F-measure weights best matches by real-cluster size" {
  let score = @src.flowsom_f_measure([0, 0, 1, 1], [0, 0, 0, 1]) catch {
    _ => abort("FlowSOM F-measure should work")
  }
  assert_true(flowsom_test_close(score, 11.0 / 15.0, 1.0e-12))
}

///|
test "FlowSOM: purity and F-measure validate label vectors" {
  let mut failures = 0
  ignore(@src.flowsom_purity([], [])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_purity([0], [0, 1])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_f_measure([], [])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_f_measure([0], [0, 1])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "FlowSOM: summary reports model dimensions and quality" {
  let summary = @src.flowsom_summary(flowsom_test_model(metaclusters=2))
  assert_true(summary.contains("FlowSOM model"))
  assert_true(summary.contains("cells: 8"))
  assert_true(summary.contains("markers: 2"))
  assert_true(summary.contains("nodes: 4"))
  assert_true(summary.contains("metaclusters: 2"))
  assert_true(summary.contains("quantization error:"))
  assert_true(summary.contains("topographic error:"))
}

///|
test "FlowSOM: FlowFrame adapter selects markers in requested order" {
  let frame = @src.FlowFrame::new(
    [[0.0, 1.0, 2.0], [0.1, 1.1, 2.1], [10.0, 11.0, 12.0], [10.1, 11.1, 12.1]],
    [
      @src.ParameterDescription::new("CD3", 20.0, 0.0, 20.0),
      @src.ParameterDescription::new("CD19", 20.0, 0.0, 20.0),
      @src.ParameterDescription::new("CD45", 20.0, 0.0, 20.0),
    ],
  )
  let config = @src.FlowSomConfig::create(xdim=2, ydim=1, rlen=2, seed=73) catch {
    _ => abort("FlowFrame FlowSOM configuration should build")
  }
  let model = @src.flowsom_train_flow_frame(frame, config, marker_indices=[2, 0]) catch {
    _ => abort("FlowFrame FlowSOM training should work")
  }
  assert_eq(model.marker_names, ["CD45", "CD3"])
  assert_eq(model.data[0], [2.0, 0.0])
  assert_eq(model.data[3], [12.1, 10.1])
}

///|
test "FlowSOM: FlowFrame adapter validates marker selection" {
  let frame = @src.FlowFrame::new([[0.0, 1.0], [2.0, 3.0]], [
    @src.ParameterDescription::new("A", 10.0, 0.0, 10.0),
    @src.ParameterDescription::new("B", 10.0, 0.0, 10.0),
  ])
  let config = @src.FlowSomConfig::create(xdim=1, ydim=1) catch {
    _ => abort("FlowFrame FlowSOM configuration should build")
  }
  let mut failures = 0
  ignore(@src.flowsom_train_flow_frame(frame, config, marker_indices=[0, 0])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_train_flow_frame(frame, config, marker_indices=[2])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "FlowSOM: SCE adapter writes one-based labels without mutating input" {
  let experiment = flowsom_test_sce()
  let config = @src.FlowSomConfig::create(
    xdim=2,
    ydim=1,
    rlen=2,
    metaclusters=2,
    meta_starts=3,
    seed=79,
  ) catch {
    _ => abort("SCE FlowSOM configuration should build")
  }
  let result = @src.flowsom_cluster_sce(experiment, config) catch {
    _ => abort("SCE FlowSOM clustering should work")
  }
  assert_false(experiment.col_data.contains("FlowSOM.cluster"))
  assert_false(experiment.col_data.contains("FlowSOM.metacluster"))
  assert_true(result.experiment.col_data.contains("FlowSOM.cluster"))
  assert_true(result.experiment.col_data.contains("FlowSOM.metacluster"))
  assert_eq(result.experiment.col_data["FlowSOM.cluster"].length(), 6)
  assert_eq(result.experiment.col_data["FlowSOM.metacluster"].length(), 6)
  for label in result.experiment.col_data["FlowSOM.cluster"] {
    assert_true(label == "1" || label == "2")
  }
  assert_eq(result.experiment.metadata["FlowSOM.assay"], "counts")
  assert_eq(result.experiment.metadata["FlowSOM.grid"], "2x1")
}

///|
test "FlowSOM: SCE adapter transposes assays and honors marker subsets" {
  let experiment = flowsom_test_sce()
  let config = @src.FlowSomConfig::create(xdim=2, ydim=1, rlen=2, seed=83) catch {
    _ => abort("SCE FlowSOM configuration should build")
  }
  let result = @src.flowsom_cluster_sce(
    experiment,
    config,
    marker_indices=[2, 0],
    cluster_column="som",
    meta_column="meta",
  ) catch {
    _ => abort("subset SCE FlowSOM clustering should work")
  }
  assert_eq(result.model.marker_names, ["CD45", "CD3"])
  assert_eq(result.model.data[0], [2.0, 0.0])
  assert_eq(result.model.data[5], [12.2, 10.2])
  assert_true(result.experiment.col_data.contains("som"))
  assert_false(result.experiment.col_data.contains("meta"))
  assert_eq(result.experiment.assays["counts"].length(), 3)
  assert_eq(result.experiment.assays["counts"][0].length(), 6)
}

///|
test "FlowSOM: SCE adapter recursively deep-copies nested containers" {
  let experiment = flowsom_test_sce()
  let config = @src.FlowSomConfig::create(xdim=2, ydim=1, rlen=2, seed=89) catch {
    _ => abort("SCE FlowSOM configuration should build")
  }
  let result = @src.flowsom_cluster_sce(experiment, config) catch {
    _ => abort("SCE FlowSOM clustering should work")
  }
  result.experiment.assays["counts"][0][0] = 999.0
  result.experiment.row_data["symbol"][0] = "changed"
  result.experiment.col_data["batch"][0] = "changed"
  result.experiment.reduced_dims["PCA"][0][0] = 999.0
  result.experiment.metadata["project"] = "changed"
  result.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 999.0
  assert_eq(experiment.assays["counts"][0][0], 0.0)
  assert_eq(experiment.row_data["symbol"][0], "T")
  assert_eq(experiment.col_data["batch"][0], "A")
  assert_eq(experiment.reduced_dims["PCA"][0][0], 0.0)
  assert_eq(experiment.metadata["project"], "flowsom-test")
  assert_eq(
    experiment.alternative_experiments["spike"].assays["counts"][0][0],
    5.0,
  )
}

///|
test "FlowSOM: SCE adapter validates assay, marker and output names" {
  let experiment = flowsom_test_sce()
  let config = @src.FlowSomConfig::create(xdim=2, ydim=1) catch {
    _ => abort("SCE FlowSOM configuration should build")
  }
  let mut failures = 0
  ignore(@src.flowsom_cluster_sce(experiment, config, assay_name="missing")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_cluster_sce(experiment, config, marker_indices=[0, 0])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.flowsom_cluster_sce(experiment, config, cluster_column="")) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.flowsom_cluster_sce(
      experiment,
      config,
      cluster_column="same",
      meta_column="same",
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}
