// Tests for the Bioconductor nnSVG-inspired NNGP spatial model.

///|
fn nnsvg_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "nnSVG value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn nnsvg_test_config() -> @src.NnsvgConfig {
  @src.NnsvgConfig::create(
    n_neighbors=3,
    length_scale_min=0.03,
    length_scale_max=1.5,
    length_scale_grid=5,
    proportion_grid=5,
    refinement_steps=1,
    minimum_variance=1.0e-9,
    fdr_threshold=0.1,
  ) catch {
    _ => abort("nnSVG test configuration should be valid")
  }
}

///|
fn nnsvg_test_result() -> @src.NnsvgResult {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  @src.nnsvg_with_config(
    expression,
    coordinates,
    nnsvg_test_config(),
    gene_names~,
  ) catch {
    _ => abort("nnSVG example fit should succeed")
  }
}

///|
fn nnsvg_test_spatial_experiment() -> @src.SpatialExperiment {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", expression))
  let counts : Array[Array[Double]] = []
  for row in expression {
    let count_row : Array[Double] = []
    for value in row {
      count_row.push((value * 4.0).round())
    }
    counts.push(count_row)
  }
  ignore(@src.se_add_assay(experiment, "counts", counts))
  for gene_name in gene_names {
    ignore(@src.se_add_row(experiment, Map([("gene_name", gene_name)])))
  }
  for spot in 0..<coordinates.length() {
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", "spot_" + (spot + 1).to_string())]),
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
test "nnSVG: default configuration exposes upstream-style controls" {
  let config = @src.NnsvgConfig::default()
  assert_eq(config.n_neighbors, 10)
  assert_true(config.ordering is @src.NnsvgOrdering::Ammd)
  assert_eq(config.length_scale_min, 0.02)
  assert_eq(config.length_scale_max, 2.0)
  assert_eq(config.length_scale_grid, 10)
  assert_eq(config.proportion_grid, 10)
  assert_eq(config.refinement_steps, 2)
  assert_eq(config.fdr_threshold, 0.05)
}

///|
test "nnSVG: custom configuration preserves parameters" {
  let ordering = @src.nnsvg_sum_coordinates_ordering()
  let config = @src.NnsvgConfig::create(
    n_neighbors=4,
    ordering~,
    length_scale_min=0.05,
    length_scale_max=3.0,
    length_scale_grid=7,
    proportion_grid=8,
    refinement_steps=3,
    minimum_variance=1.0e-8,
    fdr_threshold=0.2,
  ) catch {
    _ => abort("custom nnSVG configuration should be valid")
  }
  assert_eq(config.n_neighbors, 4)
  assert_true(config.ordering is @src.NnsvgOrdering::SumCoordinates)
  assert_eq(config.length_scale_grid, 7)
  assert_eq(config.proportion_grid, 8)
  assert_eq(config.refinement_steps, 3)
  assert_eq(config.minimum_variance, 1.0e-8)
}

///|
test "nnSVG: ordering helpers expose both ordering modes" {
  assert_true(@src.nnsvg_ammd_ordering() is @src.NnsvgOrdering::Ammd)
  assert_true(
    @src.nnsvg_sum_coordinates_ordering() is @src.NnsvgOrdering::SumCoordinates,
  )
}

///|
test "nnSVG: configuration rejects non-positive neighbor count" {
  let failed = try {
    ignore(@src.NnsvgConfig::create(n_neighbors=0))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: configuration rejects non-positive minimum length scale" {
  let failed = try {
    ignore(@src.NnsvgConfig::create(length_scale_min=0.0))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: configuration rejects reversed length-scale bounds" {
  let failed = try {
    ignore(@src.NnsvgConfig::create(length_scale_min=2.0, length_scale_max=1.0))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: configuration rejects undersized parameter grids" {
  let length_grid = try {
    ignore(@src.NnsvgConfig::create(length_scale_grid=1))
    false
  } catch {
    NnsvgError(_) => true
  }
  let proportion_grid = try {
    ignore(@src.NnsvgConfig::create(proportion_grid=1))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(length_grid)
  assert_true(proportion_grid)
}

///|
test "nnSVG: configuration rejects negative refinement count" {
  let failed = try {
    ignore(@src.NnsvgConfig::create(refinement_steps=-1))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: configuration rejects invalid variance floor" {
  let failed = try {
    ignore(@src.NnsvgConfig::create(minimum_variance=0.0))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: configuration rejects invalid FDR threshold" {
  let zero = try {
    ignore(@src.NnsvgConfig::create(fdr_threshold=0.0))
    false
  } catch {
    NnsvgError(_) => true
  }
  let large = try {
    ignore(@src.NnsvgConfig::create(fdr_threshold=1.1))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(zero)
  assert_true(large)
}

///|
test "nnSVG: example data uses gene by spot orientation" {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  assert_eq(expression.length(), 4)
  assert_eq(expression[0].length(), 12)
  assert_eq(coordinates.length(), 12)
  assert_eq(coordinates[0].length(), 2)
  assert_eq(gene_names.length(), 4)
}

///|
test "nnSVG: neighbor graph scales coordinates by global range" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=3) catch {
    _ => abort("neighbor graph should build")
  }
  nnsvg_test_close(graph.range_scale, 3.0, 1.0e-12)
  nnsvg_test_close(graph.coordinates[3][0], 1.0, 1.0e-12)
  nnsvg_test_close(graph.coordinates[4][1], 1.0 / 3.0, 1.0e-12)
}

///|
test "nnSVG: sum-coordinate ordering is deterministic" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(
    coordinates,
    n_neighbors=2,
    ordering=@src.nnsvg_sum_coordinates_ordering(),
  ) catch {
    _ => abort("sum-coordinate graph should build")
  }
  assert_eq(graph.order[0], 0)
  assert_eq(graph.order[1], 1)
  assert_eq(graph.order[2], 4)
  assert_eq(graph.order.length(), coordinates.length())
}

///|
test "nnSVG: AMMD ordering is a full permutation" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=3) catch {
    _ => abort("AMMD graph should build")
  }
  let seen = Array::make(coordinates.length(), false)
  for index in graph.order {
    assert_true(index >= 0 && index < coordinates.length())
    assert_false(seen[index])
    seen[index] = true
  }
  for value in seen {
    assert_true(value)
  }
}

///|
test "nnSVG: AMMD ordering is reproducible" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let first = @src.nnsvg_build_neighbors(coordinates, n_neighbors=3) catch {
    _ => abort("first graph should build")
  }
  let second = @src.nnsvg_build_neighbors(coordinates, n_neighbors=3) catch {
    _ => abort("second graph should build")
  }
  assert_eq(first.order, second.order)
  assert_eq(first.neighbors, second.neighbors)
}

///|
test "nnSVG: inverse ordering maps originals to processing positions" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates) catch {
    _ => abort("graph should build")
  }
  for position in 0..<graph.order.length() {
    assert_eq(graph.inverse_order[graph.order[position]], position)
  }
}

///|
test "nnSVG: all neighbors are predecessors" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=4) catch {
    _ => abort("graph should build")
  }
  for original in 0..<coordinates.length() {
    for neighbor in graph.neighbors[original] {
      assert_true(graph.inverse_order[neighbor] < graph.inverse_order[original])
    }
  }
}

///|
test "nnSVG: predecessor neighbor count is capped" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=100) catch {
    _ => abort("graph should build")
  }
  assert_eq(graph.n_neighbors, coordinates.length() - 1)
  for original in 0..<coordinates.length() {
    assert_true(
      graph.neighbors[original].length() <= graph.inverse_order[original],
    )
  }
}

///|
test "nnSVG: stored neighbor distances are non-decreasing" {
  let (_, coordinates, _) = @src.nnsvg_example_data()
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=5) catch {
    _ => abort("graph should build")
  }
  for original in 0..<coordinates.length() {
    for index in 1..<graph.distances[original].length() {
      assert_true(
        graph.distances[original][index] >= graph.distances[original][index - 1],
      )
    }
  }
}

///|
test "nnSVG: duplicate coordinates remain valid with a nugget" {
  let coordinates = [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
  let graph = @src.nnsvg_build_neighbors(coordinates, n_neighbors=2) catch {
    _ => abort("duplicate-coordinate graph should build")
  }
  assert_eq(graph.order.length(), 4)
  assert_eq(graph.neighbors.length(), 4)
}

///|
test "nnSVG: neighbor graph rejects too few spots" {
  let failed = try {
    ignore(@src.nnsvg_build_neighbors([[0.0, 0.0], [1.0, 1.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: neighbor graph rejects one-dimensional coordinates" {
  let failed = try {
    ignore(@src.nnsvg_build_neighbors([[0.0], [1.0], [2.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: neighbor graph rejects ragged coordinates" {
  let failed = try {
    ignore(@src.nnsvg_build_neighbors([[0.0, 0.0], [1.0], [2.0, 0.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: neighbor graph rejects non-finite coordinates" {
  let failed = try {
    ignore(@src.nnsvg_build_neighbors([[0.0, 0.0], [1.0e301, 1.0], [2.0, 0.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: neighbor graph rejects invariant coordinates" {
  let failed = try {
    ignore(@src.nnsvg_build_neighbors([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit dimensions and names match input" {
  let result = nnsvg_test_result()
  assert_eq(result.n_genes(), 4)
  assert_eq(result.n_spots(), 12)
  assert_eq(result.gene_names[0], "x_y_gradient")
  assert_eq(result.genes[3].gene_name, "constant")
  assert_eq(result.design_columns, 1)
}

///|
test "nnSVG: default function generates gene names" {
  let (expression, coordinates, _) = @src.nnsvg_example_data()
  let result = @src.nnsvg(expression, coordinates, n_neighbors=3) catch {
    _ => abort("default nnSVG fit should succeed")
  }
  assert_eq(result.gene_names[0], "gene_1")
  assert_eq(result.gene_names[3], "gene_4")
}

///|
test "nnSVG: custom design coefficients are retained" {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  let design : Array[Array[Double]] = []
  for coordinate in coordinates {
    design.push([1.0, coordinate[0]])
  }
  let result = @src.nnsvg_with_config(
    expression,
    coordinates,
    nnsvg_test_config(),
    design~,
    gene_names~,
  ) catch {
    _ => abort("covariate-adjusted nnSVG fit should succeed")
  }
  assert_eq(result.design_columns, 2)
  for gene in result.genes {
    assert_eq(gene.beta.length(), 2)
  }
}

///|
test "nnSVG: likelihood-ratio identity holds" {
  let result = nnsvg_test_result()
  for gene in result.genes {
    nnsvg_test_close(
      gene.likelihood_ratio,
      2.0 * (gene.log_likelihood - gene.linear_log_likelihood),
      1.0e-8,
    )
  }
}

///|
test "nnSVG: p-values use the two-degree chi-square tail" {
  let result = nnsvg_test_result()
  for gene in result.genes {
    nnsvg_test_close(
      gene.p_value,
      @math.exp(-0.5 * gene.likelihood_ratio),
      1.0e-12,
    )
  }
}

///|
test "nnSVG: variance components reconstruct total fitted variance" {
  let result = nnsvg_test_result()
  for gene in result.genes {
    let total = gene.sigma_sq + gene.tau_sq
    assert_true(total > 0.0)
    nnsvg_test_close(
      gene.proportion_spatial_variance,
      gene.sigma_sq / total,
      1.0e-10,
    )
  }
}

///|
test "nnSVG: positive spatial fits report reciprocal phi" {
  let result = nnsvg_test_result()
  for gene in result.genes {
    if gene.length_scale > 0.0 {
      nnsvg_test_close(gene.phi, 1.0 / gene.length_scale, 1.0e-10)
    } else {
      assert_eq(gene.phi, 0.0)
    }
  }
}

///|
test "nnSVG: adjusted p-values are bounded and not below raw values" {
  let result = nnsvg_test_result()
  for gene in result.genes {
    assert_true(gene.p_value >= 0.0 && gene.p_value <= 1.0)
    assert_true(
      gene.adjusted_p_value >= gene.p_value - 1.0e-12 &&
      gene.adjusted_p_value <= 1.0,
    )
  }
}

///|
test "nnSVG: constant gene falls back to non-spatial model" {
  let result = nnsvg_test_result()
  let constant = result.genes[3]
  nnsvg_test_close(constant.likelihood_ratio, 0.0, 1.0e-12)
  nnsvg_test_close(constant.p_value, 1.0, 1.0e-12)
  nnsvg_test_close(constant.sigma_sq, 0.0, 1.0e-12)
  assert_eq(constant.length_scale, 0.0)
}

///|
test "nnSVG: spatial gradients outrank constant expression" {
  let result = nnsvg_test_result()
  assert_true(result.genes[0].rank < result.genes[3].rank)
  assert_true(result.genes[1].rank < result.genes[3].rank)
  assert_true(
    result.genes[0].likelihood_ratio > result.genes[3].likelihood_ratio,
  )
}

///|
test "nnSVG: gene lookup finds names and rejects unknown names" {
  let result = nnsvg_test_result()
  let found = match result.gene("checkerboard") {
    Some(value) => value
    None => abort("known gene should be found")
  }
  assert_eq(found.gene_name, "checkerboard")
  assert_true(result.gene("missing") is None)
}

///|
test "nnSVG: top results are ordered by rank" {
  let top = nnsvg_test_result().top(3)
  assert_eq(top.length(), 3)
  assert_true(top[0].rank <= top[1].rank)
  assert_true(top[1].rank <= top[2].rank)
}

///|
test "nnSVG: top handles zero and oversized requests" {
  let result = nnsvg_test_result()
  assert_eq(result.top(0).length(), 0)
  assert_eq(result.top(100).length(), result.n_genes())
}

///|
test "nnSVG: significant supports an explicit threshold" {
  let result = nnsvg_test_result()
  let all = result.significant(threshold=1.0)
  let none = result.significant(threshold=0.0)
  assert_eq(all.length(), result.n_genes())
  assert_eq(none.length(), 0)
}

///|
test "nnSVG: summary reports dimensions and neighbor count" {
  let summary = nnsvg_test_result().summary()
  assert_true(summary.contains("4 genes x 12 spots"))
  assert_true(summary.contains("neighbors=3"))
  assert_true(summary.contains("FDR=0.1"))
}

///|
test "nnSVG: fit rejects empty expression" {
  let failed = try {
    ignore(@src.nnsvg([], [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects too few spots" {
  let failed = try {
    ignore(@src.nnsvg([[1.0, 2.0]], [[0.0, 0.0], [1.0, 0.0]]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects ragged expression" {
  let failed = try {
    ignore(
      @src.nnsvg([[1.0, 2.0, 3.0], [1.0, 2.0]], [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
      ]),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects non-finite expression" {
  let failed = try {
    ignore(
      @src.nnsvg([[1.0, 1.0e301, 3.0]], [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects coordinate count mismatch" {
  let (expression, coordinates, _) = @src.nnsvg_example_data()
  let shortened : Array[Array[Double]] = []
  for index in 0..<(coordinates.length() - 1) {
    shortened.push(coordinates[index])
  }
  let failed = try {
    ignore(@src.nnsvg(expression, shortened))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects gene name count mismatch" {
  let (expression, coordinates, _) = @src.nnsvg_example_data()
  let failed = try {
    ignore(@src.nnsvg(expression, coordinates, gene_names=["only_one"]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects design row mismatch" {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  let failed = try {
    ignore(
      @src.nnsvg(expression, coordinates, design=[[1.0], [1.0]], gene_names~),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: fit rejects rank-deficient design" {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  let design : Array[Array[Double]] = []
  for _ in coordinates {
    design.push([1.0, 1.0])
  }
  let failed = try {
    ignore(@src.nnsvg(expression, coordinates, design~, gene_names~))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: gene filtering applies expression threshold" {
  let counts = [
    [4.0, 4.0, 0.0, 0.0],
    [2.0, 2.0, 2.0, 2.0],
    [5.0, 5.0, 5.0, 5.0],
  ]
  let filter = @src.nnsvg_filter_genes(
    counts,
    ["A", "B", "C"],
    minimum_count=3.0,
    minimum_spot_percentage=50.0,
    filter_mitochondrial=false,
  ) catch {
    _ => abort("gene filtering should succeed")
  }
  assert_eq(filter.kept_indices, [0, 2])
  assert_eq(filter.removed_low_expression, [1])
  assert_eq(filter.counts.length(), 2)
}

///|
test "nnSVG: gene filtering removes mitochondrial prefixes" {
  let filter = @src.nnsvg_filter_genes(
    [[4.0, 4.0], [4.0, 4.0], [4.0, 4.0]],
    ["MT-ND1", "mt-Co1", "ACTB"],
    minimum_spot_percentage=100.0,
  ) catch {
    _ => abort("mitochondrial filtering should succeed")
  }
  assert_eq(filter.kept_indices, [2])
  assert_eq(filter.removed_mitochondrial, [0, 1])
}

///|
test "nnSVG: mitochondrial filtering can be disabled" {
  let filter = @src.nnsvg_filter_genes(
    [[4.0, 4.0], [4.0, 4.0]],
    ["MT-ND1", "ACTB"],
    minimum_spot_percentage=100.0,
    filter_mitochondrial=false,
  ) catch {
    _ => abort("disabled mitochondrial filtering should succeed")
  }
  assert_eq(filter.kept_indices, [0, 1])
  assert_eq(filter.removed_mitochondrial.length(), 0)
}

///|
test "nnSVG: default filter percentage requires one detected spot" {
  let filter = @src.nnsvg_filter_genes(
    [[0.0, 0.0, 0.0, 3.0], [0.0, 0.0, 0.0, 0.0]],
    ["detected", "absent"],
    filter_mitochondrial=false,
  ) catch {
    _ => abort("default filtering should succeed")
  }
  assert_eq(filter.kept_indices, [0])
  assert_eq(filter.removed_low_expression, [1])
}

///|
test "nnSVG: gene filtering rejects name mismatch" {
  let failed = try {
    ignore(@src.nnsvg_filter_genes([[1.0, 2.0]], []))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: gene filtering rejects ragged counts" {
  let failed = try {
    ignore(@src.nnsvg_filter_genes([[1.0, 2.0], [1.0]], ["A", "B"]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: gene filtering rejects negative counts" {
  let failed = try {
    ignore(@src.nnsvg_filter_genes([[-1.0, 2.0]], ["A"]))
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: gene filtering rejects invalid thresholds" {
  let count = try {
    ignore(@src.nnsvg_filter_genes([[1.0, 2.0]], ["A"], minimum_count=-1.0))
    false
  } catch {
    NnsvgError(_) => true
  }
  let percentage = try {
    ignore(
      @src.nnsvg_filter_genes(
        [[1.0, 2.0]],
        ["A"],
        minimum_spot_percentage=101.0,
      ),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(count)
  assert_true(percentage)
}

///|
test "nnSVG: SpatialExperiment integration preserves assays" {
  let experiment = nnsvg_test_spatial_experiment()
  let output = @src.nnsvg_spatial_experiment(
    experiment,
    config=nnsvg_test_config(),
  ) catch {
    _ => abort("SpatialExperiment integration should succeed")
  }
  assert_eq(output.result.n_genes(), 4)
  assert_true(output.experiment.assay.contains("logcounts"))
  assert_true(output.experiment.assay.contains("counts"))
}

///|
test "nnSVG: SpatialExperiment integration writes rowData statistics" {
  let output = @src.nnsvg_spatial_experiment(
    nnsvg_test_spatial_experiment(),
    config=nnsvg_test_config(),
  ) catch {
    _ => abort("SpatialExperiment integration should succeed")
  }
  let row = output.experiment.row_data[0]
  assert_true(row.contains("nnsvg_sigma_sq"))
  assert_true(row.contains("nnsvg_tau_sq"))
  assert_true(row.contains("nnsvg_prop_sv"))
  assert_true(row.contains("nnsvg_LR_stat"))
  assert_true(row.contains("nnsvg_pval"))
  assert_true(row.contains("nnsvg_padj"))
}

///|
test "nnSVG: SpatialExperiment integration writes metadata" {
  let output = @src.nnsvg_spatial_experiment(
    nnsvg_test_spatial_experiment(),
    config=nnsvg_test_config(),
  ) catch {
    _ => abort("SpatialExperiment integration should succeed")
  }
  assert_eq(output.experiment.metadata["nnsvg_assay"], "logcounts")
  assert_eq(output.experiment.metadata["nnsvg_neighbors"], "3")
  assert_eq(output.experiment.metadata["nnsvg_genes"], "4")
}

///|
test "nnSVG: SpatialExperiment integration leaves input rowData unchanged" {
  let experiment = nnsvg_test_spatial_experiment()
  ignore(
    @src.nnsvg_spatial_experiment(experiment, config=nnsvg_test_config()) catch {
      _ => abort("SpatialExperiment integration should succeed")
    },
  )
  assert_false(experiment.row_data[0].contains("nnsvg_pval"))
}

///|
test "nnSVG: SpatialExperiment integration rejects missing assay" {
  let failed = try {
    ignore(
      @src.nnsvg_spatial_experiment(
        nnsvg_test_spatial_experiment(),
        assay_name="missing",
        config=nnsvg_test_config(),
      ),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: SpatialExperiment integration rejects coordinate mismatch" {
  let experiment = nnsvg_test_spatial_experiment()
  ignore(experiment.spatial_coords.pop())
  let failed = try {
    ignore(
      @src.nnsvg_spatial_experiment(experiment, config=nnsvg_test_config()),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: SpatialExperiment filter subsets every assay" {
  let output = @src.nnsvg_filter_spatial_experiment(
    nnsvg_test_spatial_experiment(),
    minimum_count=5.0,
    minimum_spot_percentage=50.0,
    filter_mitochondrial=false,
  ) catch {
    _ => abort("SpatialExperiment filtering should succeed")
  }
  assert_eq(
    output.experiment.assay["counts"].length(),
    output.filter.kept_indices.length(),
  )
  assert_eq(
    output.experiment.assay["logcounts"].length(),
    output.filter.kept_indices.length(),
  )
  assert_eq(
    output.experiment.row_data.length(),
    output.filter.kept_indices.length(),
  )
}

///|
test "nnSVG: SpatialExperiment filter preserves spot data" {
  let experiment = nnsvg_test_spatial_experiment()
  let output = @src.nnsvg_filter_spatial_experiment(
    experiment,
    minimum_count=5.0,
    minimum_spot_percentage=50.0,
    filter_mitochondrial=false,
  ) catch {
    _ => abort("SpatialExperiment filtering should succeed")
  }
  assert_eq(
    output.experiment.spatial_coords.length(),
    experiment.spatial_coords.length(),
  )
  assert_eq(output.experiment.col_data.length(), experiment.col_data.length())
}

///|
test "nnSVG: SpatialExperiment filter records metadata" {
  let output = @src.nnsvg_filter_spatial_experiment(
    nnsvg_test_spatial_experiment(),
    minimum_count=5.0,
    minimum_spot_percentage=50.0,
    filter_mitochondrial=false,
  ) catch {
    _ => abort("SpatialExperiment filtering should succeed")
  }
  assert_true(output.experiment.metadata.contains("nnsvg_filter_kept"))
  assert_true(output.experiment.metadata.contains("nnsvg_filter_removed"))
}

///|
test "nnSVG: SpatialExperiment filter rejects missing count assay" {
  let failed = try {
    ignore(
      @src.nnsvg_filter_spatial_experiment(
        nnsvg_test_spatial_experiment(),
        count_assay="missing",
      ),
    )
    false
  } catch {
    NnsvgError(_) => true
  }
  assert_true(failed)
}

///|
test "nnSVG: three-dimensional SpatialExperiment coordinates are supported" {
  let experiment = nnsvg_test_spatial_experiment()
  experiment.spatial_coords[1] = @src.SpatialCoord::new(1.0, 0.0, 1.0)
  let output = @src.nnsvg_spatial_experiment(
    experiment,
    config=nnsvg_test_config(),
  ) catch {
    _ => abort("3D SpatialExperiment integration should succeed")
  }
  assert_eq(output.result.graph.coordinates[0].length(), 3)
}
