///|
/// Tests for the Bioconductor miloR-inspired neighbourhood DA workflow.

///|
fn milo_test_assert_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort("miloR values differ beyond tolerance")
  }
}

///|
fn milo_test_small_graph() -> @src.MiloExperiment {
  @src.milo_build_graph(
    [[0.0], [1.0], [3.0], [10.0]],
    cell_names=["A", "B", "C", "D"],
    k=1,
  ) catch {
    _ => abort("valid small Milo graph should build")
  }
}

///|
fn milo_test_counted() -> (@src.MiloExperiment, @src.MiloDesign) {
  let data = @src.milo_sample_data()
  let graph = @src.milo_build_graph(
    data.coordinates,
    cell_names=data.cell_names,
    k=8,
    dimensions=2,
  ) catch {
    _ => abort("sample Milo graph should build")
  }
  let neighborhoods = graph.make_neighborhoods(proportion=0.5, seed=19) catch {
    _ => abort("sample Milo neighborhoods should build")
  }
  let counted = neighborhoods.count_cells(data.sample_ids) catch {
    _ => abort("sample Milo cells should count")
  }
  let design = @src.milo_binary_design(data.sample_conditions, reference="A") catch {
    _ => abort("sample Milo design should build")
  }
  (counted, design)
}

///|
test "miloR: sample data contains balanced replicated samples" {
  let data = @src.milo_sample_data()
  assert_eq(data.coordinates.length(), 32)
  assert_eq(data.cell_names.length(), 32)
  assert_eq(data.sample_ids.length(), 32)
  assert_eq(data.sample_names, ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"])
  assert_eq(data.sample_conditions, ["A", "A", "A", "A", "B", "B", "B", "B"])
}

///|
test "miloR: builds exact KNN graph and generates cell names" {
  let graph = @src.milo_build_graph([[0.0], [1.0], [3.0]], k=1) catch {
    _ => abort("valid graph should build")
  }
  assert_eq(graph.cell_names, ["Cell1", "Cell2", "Cell3"])
  assert_eq(graph.dimensions, 1)
  assert_eq(graph.k, 1)
  assert_eq(graph.knn_indices, [[1], [0], [1]])
}

///|
test "miloR: exact KNN excludes each query cell" {
  let graph = milo_test_small_graph()
  for cell in 0..<graph.knn_indices.length() {
    assert_false(graph.knn_indices[cell].contains(cell))
  }
}

///|
test "miloR: KNN distances are sorted nearest first" {
  let graph = @src.milo_build_graph([[0.0], [1.0], [3.0], [10.0]], k=3) catch {
    _ => abort("valid graph should build")
  }
  assert_eq(graph.knn_indices[0], [1, 2, 3])
  assert_eq(graph.knn_distances[0], [1.0, 3.0, 10.0])
}

///|
test "miloR: KNN ties are resolved by cell index" {
  let graph = @src.milo_build_graph([[-1.0], [0.0], [1.0]], k=2) catch {
    _ => abort("valid tied graph should build")
  }
  assert_eq(graph.knn_indices[1], [0, 2])
}

///|
test "miloR: directed KNN edges are symmetrized" {
  let graph = milo_test_small_graph()
  assert_eq(graph.graph_neighbors, [[1], [0, 2], [1, 3], [2]])
}

///|
test "miloR: selected dimensions control graph distances" {
  let graph = @src.milo_build_graph(
    [[0.0, 100.0], [1.0, 0.0], [2.0, 100.0]],
    k=1,
    dimensions=1,
  ) catch {
    _ => abort("valid dimension subset should build")
  }
  assert_eq(graph.knn_indices[0], [1])
  assert_eq(graph.knn_indices[2], [1])
}

///|
test "miloR: graph construction copies coordinate input" {
  let coordinates = [[0.0], [1.0], [2.0]]
  let graph = @src.milo_build_graph(coordinates, k=1) catch {
    _ => abort("valid graph should build")
  }
  coordinates[0][0] = 99.0
  assert_eq(graph.coordinates[0][0], 0.0)
}

///|
test "miloR: graph rejects malformed coordinates" {
  let too_small = try {
    ignore(@src.milo_build_graph([[0.0]], k=1))
    false
  } catch {
    MiloError(_) => true
  }
  let ragged = try {
    ignore(@src.milo_build_graph([[0.0], [1.0, 2.0]], k=1))
    false
  } catch {
    MiloError(_) => true
  }
  let non_finite = try {
    ignore(@src.milo_build_graph([[0.0], [@double.not_a_number]], k=1))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(too_small)
  assert_true(ragged)
  assert_true(non_finite)
}

///|
test "miloR: graph validates k and dimensions" {
  let bad_k = try {
    ignore(@src.milo_build_graph([[0.0], [1.0]], k=2))
    false
  } catch {
    MiloError(_) => true
  }
  let bad_dimensions = try {
    ignore(@src.milo_build_graph([[0.0], [1.0]], k=1, dimensions=2))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(bad_k)
  assert_true(bad_dimensions)
}

///|
test "miloR: graph validates cell names" {
  let wrong_length = try {
    ignore(@src.milo_build_graph([[0.0], [1.0], [2.0]], cell_names=["A"], k=1))
    false
  } catch {
    MiloError(_) => true
  }
  let duplicates = try {
    ignore(
      @src.milo_build_graph(
        [[0.0], [1.0], [2.0]],
        cell_names=["A", "A", "B"],
        k=1,
      ),
    )
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(wrong_length)
  assert_true(duplicates)
}

///|
test "miloR: unrefined sampling is deterministic" {
  let graph = @src.milo_build_graph(
    [[0.0], [1.0], [2.0], [3.0], [4.0], [5.0]],
    k=2,
  ) catch {
    _ => abort("valid graph should build")
  }
  let first = graph.make_neighborhoods(proportion=0.5, refined=false, seed=7) catch {
    _ => abort("valid neighborhoods should build")
  }
  let second = graph.make_neighborhoods(proportion=0.5, refined=false, seed=7) catch {
    _ => abort("valid neighborhoods should build")
  }
  assert_eq(first.neighborhood_indices, second.neighborhood_indices)
  assert_eq(first.neighborhood_indices.length(), 3)
}

///|
test "miloR: refined sampling collapses duplicate representatives" {
  let data = @src.milo_sample_data()
  let graph = @src.milo_build_graph(data.coordinates, k=4) catch {
    _ => abort("valid graph should build")
  }
  let refined = graph.make_neighborhoods(proportion=0.5, refined=true, seed=11) catch {
    _ => abort("valid refined neighborhoods should build")
  }
  assert_true(refined.neighborhoods.length() > 0)
  assert_true(refined.neighborhoods.length() <= 16)
  for index in refined.neighborhood_indices {
    assert_true(index >= 0 && index < data.coordinates.length())
  }
}

///|
test "miloR: neighborhoods include index and undirected graph neighbors" {
  let graph = milo_test_small_graph()
  let experiment = graph.make_neighborhoods(
    proportion=0.75,
    refined=false,
    seed=3,
  ) catch {
    _ => abort("valid neighborhoods should build")
  }
  for neighborhood in 0..<experiment.neighborhoods.length() {
    let index = experiment.neighborhood_indices[neighborhood]
    assert_true(experiment.neighborhoods[neighborhood].contains(index))
    for neighbor in experiment.graph_neighbors[index] {
      assert_true(experiment.neighborhoods[neighborhood].contains(neighbor))
    }
  }
}

///|
test "miloR: neighborhood sampling validates proportion" {
  let graph = milo_test_small_graph()
  let zero = try {
    ignore(graph.make_neighborhoods(proportion=0.0))
    false
  } catch {
    MiloError(_) => true
  }
  let one = try {
    ignore(graph.make_neighborhoods(proportion=1.0))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(zero)
  assert_true(one)
}

///|
test "miloR: counts neighborhood cells by first-seen sample order" {
  let graph = milo_test_small_graph()
  let neighborhoods = graph.make_neighborhoods(
    proportion=0.75,
    refined=false,
    seed=5,
  ) catch {
    _ => abort("valid neighborhoods should build")
  }
  let counted = neighborhoods.count_cells(["S2", "S1", "S1", "S2"]) catch {
    _ => abort("valid sample IDs should count")
  }
  assert_eq(counted.sample_names, ["S2", "S1"])
  for index in 0..<counted.neighborhoods.length() {
    let sum = counted.neighborhood_counts[index][0] +
      counted.neighborhood_counts[index][1]
    assert_eq(sum, counted.neighborhoods[index].length())
  }
}

///|
test "miloR: counting returns a new experiment" {
  let graph = milo_test_small_graph()
  let neighborhoods = graph.make_neighborhoods(proportion=0.5, refined=false) catch {
    _ => abort("valid neighborhoods should build")
  }
  let counted = neighborhoods.count_cells(["S1", "S1", "S2", "S2"]) catch {
    _ => abort("valid sample IDs should count")
  }
  assert_eq(neighborhoods.sample_names.length(), 0)
  assert_eq(counted.sample_names, ["S1", "S2"])
}

///|
test "miloR: counting validates workflow and sample IDs" {
  let graph = milo_test_small_graph()
  let before_neighborhoods = try {
    ignore(graph.count_cells(["S1", "S1", "S2", "S2"]))
    false
  } catch {
    MiloError(_) => true
  }
  let neighborhoods = graph.make_neighborhoods(proportion=0.5, refined=false) catch {
    _ => abort("valid neighborhoods should build")
  }
  let wrong_length = try {
    ignore(neighborhoods.count_cells(["S1"]))
    false
  } catch {
    MiloError(_) => true
  }
  let empty = try {
    ignore(neighborhoods.count_cells(["S1", "", "S2", "S2"]))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(before_neighborhoods)
  assert_true(wrong_length)
  assert_true(empty)
}

///|
test "miloR: computes mean feature expression per neighborhood" {
  let graph = milo_test_small_graph()
  let neighborhoods = graph.make_neighborhoods(
    proportion=0.5,
    refined=false,
    seed=2,
  ) catch {
    _ => abort("valid neighborhoods should build")
  }
  let expression = neighborhoods.neighborhood_expression([[0.0, 2.0, 4.0, 8.0]]) catch {
    _ => abort("valid expression should aggregate")
  }
  assert_eq(expression.length(), 1)
  assert_eq(expression[0].length(), neighborhoods.neighborhoods.length())
  for index in 0..<neighborhoods.neighborhoods.length() {
    let mut total = 0.0
    for cell in neighborhoods.neighborhoods[index] {
      total = total + [0.0, 2.0, 4.0, 8.0][cell]
    }
    milo_test_assert_close(
      expression[0][index],
      total / neighborhoods.neighborhoods[index].length().to_double(),
      1.0e-12,
    )
  }
}

///|
test "miloR: expression aggregation validates columns and values" {
  let graph = milo_test_small_graph()
  let neighborhoods = graph.make_neighborhoods(proportion=0.5, refined=false) catch {
    _ => abort("valid neighborhoods should build")
  }
  let wrong_length = try {
    ignore(neighborhoods.neighborhood_expression([[1.0]]))
    false
  } catch {
    MiloError(_) => true
  }
  let non_finite = try {
    ignore(
      neighborhoods.neighborhood_expression([
        [@double.not_a_number, 1.0, 2.0, 3.0],
      ]),
    )
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(wrong_length)
  assert_true(non_finite)
}

///|
test "miloR: neighborhood overlap matrix is symmetric" {
  let (experiment, _) = milo_test_counted()
  let overlaps = experiment.neighborhood_overlap_matrix()
  assert_eq(overlaps.length(), experiment.neighborhoods.length())
  for left in 0..<overlaps.length() {
    assert_eq(overlaps[left][left], experiment.neighborhoods[left].length())
    for right in 0..<overlaps.length() {
      assert_eq(overlaps[left][right], overlaps[right][left])
    }
  }
}

///|
test "miloR: binary design uses explicit reference condition" {
  let design = @src.milo_binary_design(
    ["treated", "control", "treated", "control"],
    reference="control",
  ) catch {
    _ => abort("valid binary design should build")
  }
  assert_eq(design.levels, ["control", "treated"])
  assert_eq(design.reference, "control")
  assert_eq(design.coefficient_name, "treated_vs_control")
  assert_eq(design.matrix, [[1.0, 1.0], [1.0, 0.0], [1.0, 1.0], [1.0, 0.0]])
}

///|
test "miloR: binary design validates conditions and reference" {
  let one_level = try {
    ignore(@src.milo_binary_design(["A", "A"]))
    false
  } catch {
    MiloError(_) => true
  }
  let three_levels = try {
    ignore(@src.milo_binary_design(["A", "B", "C"]))
    false
  } catch {
    MiloError(_) => true
  }
  let absent_reference = try {
    ignore(@src.milo_binary_design(["A", "B"], reference="C"))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(one_level)
  assert_true(three_levels)
  assert_true(absent_reference)
}

///|
test "miloR: graph spatial FDR supports k-distance weighting" {
  let (experiment, _) = milo_test_counted()
  let p_values = Array::make(experiment.neighborhoods.length(), 0.05)
  let adjusted = @src.milo_graph_spatial_fdr(experiment, p_values) catch {
    _ => abort("valid spatial FDR should compute")
  }
  assert_eq(adjusted.length(), p_values.length())
  for value in adjusted {
    assert_true(value >= 0.0 && value <= 1.0)
  }
}

///|
test "miloR: spatial FDR supports all graph weighting schemes" {
  let (experiment, _) = milo_test_counted()
  let p_values : Array[Double] = []
  for index in 0..<experiment.neighborhoods.length() {
    p_values.push(
      (index + 1).to_double() /
      (experiment.neighborhoods.length() + 1).to_double(),
    )
  }
  let neighbor = @src.milo_graph_spatial_fdr(
    experiment,
    p_values,
    weighting=@src.milo_neighbour_distance_weighting(),
  ) catch {
    _ => abort("neighbor-distance FDR should compute")
  }
  let maximum = @src.milo_graph_spatial_fdr(
    experiment,
    p_values,
    weighting=@src.milo_max_distance_weighting(),
  ) catch {
    _ => abort("maximum-distance FDR should compute")
  }
  let overlap = @src.milo_graph_spatial_fdr(
    experiment,
    p_values,
    weighting=@src.milo_graph_overlap_weighting(),
  ) catch {
    _ => abort("graph-overlap FDR should compute")
  }
  assert_eq(neighbor.length(), p_values.length())
  assert_eq(maximum.length(), p_values.length())
  assert_eq(overlap.length(), p_values.length())
}

///|
test "miloR: disabled spatial weighting returns missing values" {
  let (experiment, _) = milo_test_counted()
  let adjusted = @src.milo_graph_spatial_fdr(
    experiment,
    Array::make(experiment.neighborhoods.length(), 0.5),
    weighting=@src.milo_no_spatial_weighting(),
  ) catch {
    _ => abort("disabled spatial FDR should return missing values")
  }
  for value in adjusted {
    assert_true(value.is_nan())
  }
}

///|
test "miloR: spatial FDR validates p-values" {
  let (experiment, _) = milo_test_counted()
  let wrong_length = try {
    ignore(@src.milo_graph_spatial_fdr(experiment, []))
    false
  } catch {
    MiloError(_) => true
  }
  let invalid = Array::make(experiment.neighborhoods.length(), 0.5)
  invalid[0] = 1.5
  let out_of_range = try {
    ignore(@src.milo_graph_spatial_fdr(experiment, invalid))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(wrong_length)
  assert_true(out_of_range)
}

///|
test "miloR: complete sample pipeline reports expected dimensions" {
  let (experiment, _) = milo_test_counted()
  assert_eq(experiment.coordinates.length(), 32)
  assert_eq(experiment.sample_names.length(), 8)
  assert_true(experiment.neighborhoods.length() > 1)
  assert_true(experiment.summary().contains("cells=32"))
  assert_true(experiment.summary().contains("samples=8"))
}

///|
test "miloR: NB GLM returns one finite result per neighborhood" {
  let (experiment, design) = milo_test_counted()
  let results = experiment.test_neighborhoods(design.matrix, 1) catch {
    _ => abort("valid Milo DA test should fit")
  }
  assert_eq(results.length(), experiment.neighborhoods.length())
  for result in results {
    assert_true(result.tested)
    assert_true(!result.log_fc.is_nan())
    assert_true(!result.statistic.is_nan())
    assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
    assert_true(result.fdr >= 0.0 && result.fdr <= 1.0)
    assert_true(result.spatial_fdr >= 0.0 && result.spatial_fdr <= 1.0)
    assert_true(result.dispersion > 0.0)
  }
}

///|
test "miloR: DA direction follows left and right state enrichment" {
  let (experiment, design) = milo_test_counted()
  let results = experiment.test_neighborhoods(design.matrix, 1, cell_sizes=[
    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
  ]) catch {
    _ => abort("valid Milo DA test should fit")
  }
  let mut left = 0
  let mut right = 0
  for result in results {
    let x = experiment.coordinates[result.index_cell][0]
    let counts = experiment.neighborhood_counts[result.neighborhood]
    let condition_a = counts[0] + counts[1] + counts[2] + counts[3]
    let condition_b = counts[4] + counts[5] + counts[6] + counts[7]
    if x < 2.5 {
      assert_true(condition_a >= condition_b)
      if condition_a > condition_b {
        assert_true(result.log_fc < 0.0)
        left = left + 1
      }
    } else {
      assert_true(condition_b >= condition_a)
      if condition_b > condition_a {
        assert_true(result.log_fc > 0.0)
        right = right + 1
      }
    }
  }
  assert_true(left > 0)
  assert_true(right > 0)
}

///|
test "miloR: mean threshold marks low-count neighborhoods untested" {
  let (experiment, design) = milo_test_counted()
  let results = experiment.test_neighborhoods(design.matrix, 1, min_mean=100.0) catch {
    _ => abort("valid threshold should apply")
  }
  for result in results {
    assert_false(result.tested)
    assert_eq(result.p_value, 1.0)
    assert_eq(result.fdr, 1.0)
  }
}

///|
test "miloR: DA supports explicit cell-size normalization" {
  let (experiment, design) = milo_test_counted()
  let results = experiment.test_neighborhoods(design.matrix, 1, cell_sizes=[
    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
  ]) catch {
    _ => abort("valid cell sizes should normalize")
  }
  assert_eq(results.length(), experiment.neighborhoods.length())
}

///|
test "miloR: DA supports disabling spatial correction" {
  let (experiment, design) = milo_test_counted()
  let results = experiment.test_neighborhoods(
    design.matrix,
    1,
    weighting=@src.milo_no_spatial_weighting(),
  ) catch {
    _ => abort("valid DA test should fit without spatial correction")
  }
  for result in results {
    assert_true(result.spatial_fdr.is_nan())
  }
}

///|
test "miloR: DA validates design dimensions rank and contrast" {
  let (experiment, _) = milo_test_counted()
  let wrong_rows = try {
    ignore(experiment.test_neighborhoods([[1.0, 0.0]], 1))
    false
  } catch {
    MiloError(_) => true
  }
  let singular : Array[Array[Double]] = []
  for _ in 0..<8 {
    singular.push([1.0, 1.0])
  }
  let bad_rank = try {
    ignore(experiment.test_neighborhoods(singular, 1))
    false
  } catch {
    MiloError(_) => true
  }
  let design = @src.milo_binary_design(["A", "A", "A", "A", "B", "B", "B", "B"]) catch {
    _ => abort("valid design should build")
  }
  let bad_contrast = try {
    ignore(experiment.test_neighborhoods(design.matrix, 2))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(wrong_rows)
  assert_true(bad_rank)
  assert_true(bad_contrast)
}

///|
test "miloR: DA validates minimum mean and cell sizes" {
  let (experiment, design) = milo_test_counted()
  let bad_mean = try {
    ignore(experiment.test_neighborhoods(design.matrix, 1, min_mean=-1.0))
    false
  } catch {
    MiloError(_) => true
  }
  let bad_sizes = try {
    ignore(experiment.test_neighborhoods(design.matrix, 1, cell_sizes=[1.0]))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(bad_mean)
  assert_true(bad_sizes)
}

///|
test "miloR: SCE constructor uses named reduced dimensions and cells" {
  let sce = @src.SingleCellExperiment::new([[1.0, 2.0, 3.0, 4.0]], ["Gene1"], [
    "A", "B", "C", "D",
  ])
  let with_pca = @src.sce_set_reduced_dim(sce, "PCA", [
    [0.0],
    [1.0],
    [2.0],
    [3.0],
  ])
  let graph = @src.milo_from_single_cell_experiment(with_pca, k=1) catch {
    _ => abort("valid SCE reduced dimensions should build a graph")
  }
  assert_eq(graph.cell_names, ["A", "B", "C", "D"])
  assert_eq(graph.knn_indices.length(), 4)
}

///|
test "miloR: SCE constructor rejects missing reduced dimensions" {
  let sce = @src.SingleCellExperiment::new([[1.0, 2.0]], ["Gene1"], ["A", "B"])
  let missing = try {
    ignore(@src.milo_from_single_cell_experiment(sce, k=1))
    false
  } catch {
    MiloError(_) => true
  }
  assert_true(missing)
}
