// Tests for Bioconductor decontX-inspired ambient RNA decontamination.

///|
fn decontx_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "decontX value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn decontx_test_result() -> @src.DecontXResult {
  let (counts, genes, cells, clusters, _) = @src.decontx_example_data()
  @src.decontx(counts, genes, cells, clusters) catch {
    _ => abort("example decontX analysis should succeed")
  }
}

///|
fn decontx_test_sce() -> @src.SingleCellExperiment {
  let (counts, _, _, clusters, _) = @src.decontx_example_data()
  @src.SCEBuilder::new()
  |> @src.SCEBuilder::add_assay("counts", counts)
  |> @src.SCEBuilder::add_row_data("symbol", [
    "TR1", "TR2", "BR1", "BR2", "HK", "LOW",
  ])
  |> @src.SCEBuilder::add_col_data("cluster", clusters)
  |> @src.SCEBuilder::build
}

///|
test "decontX: default configuration is scientifically conservative" {
  let config = @src.DecontXConfig::default()
  assert_eq(config.max_iterations, 100)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.contamination_prior_alpha, 1.0)
  assert_eq(config.contamination_prior_beta, 9.0)
  assert_eq(config.initial_contamination, 0.1)
}

///|
test "decontX: custom configuration preserves all parameters" {
  let config = @src.DecontXConfig::create(
    max_iterations=25,
    tolerance=1.0e-5,
    contamination_prior_alpha=2.0,
    contamination_prior_beta=8.0,
    profile_pseudocount=0.25,
    initial_contamination=0.2,
    minimum_probability=1.0e-10,
  ) catch {
    _ => abort("valid custom configuration should build")
  }
  assert_eq(config.max_iterations, 25)
  assert_eq(config.profile_pseudocount, 0.25)
  assert_eq(config.minimum_probability, 1.0e-10)
}

///|
test "decontX: configuration rejects invalid iteration and tolerance" {
  let iterations = try {
    ignore(@src.DecontXConfig::create(max_iterations=0))
    false
  } catch {
    DecontXError(_) => true
  }
  let tolerance = try {
    ignore(@src.DecontXConfig::create(tolerance=0.0))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(iterations)
  assert_true(tolerance)
}

///|
test "decontX: configuration rejects invalid priors" {
  let beta = try {
    ignore(@src.DecontXConfig::create(contamination_prior_beta=0.0))
    false
  } catch {
    DecontXError(_) => true
  }
  let profile = try {
    ignore(@src.DecontXConfig::create(profile_pseudocount=-1.0))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(beta)
  assert_true(profile)
}

///|
test "decontX: configuration rejects invalid probabilities" {
  let initial = try {
    ignore(@src.DecontXConfig::create(initial_contamination=1.0))
    false
  } catch {
    DecontXError(_) => true
  }
  let minimum = try {
    ignore(@src.DecontXConfig::create(minimum_probability=0.0))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(initial)
  assert_true(minimum)
}

///|
test "decontX: example data uses gene by cell orientation" {
  let (counts, genes, cells, clusters, background) = @src.decontx_example_data()
  assert_eq(counts.length(), 6)
  assert_eq(counts[0].length(), 8)
  assert_eq(genes.length(), 6)
  assert_eq(cells.length(), 8)
  assert_eq(clusters.length(), 8)
  assert_eq(background.length(), 6)
}

///|
test "decontX: result dimensions match the source matrix" {
  let result = decontx_test_result()
  assert_eq(result.n_genes(), 6)
  assert_eq(result.n_cells(), 8)
  assert_eq(result.n_clusters(), 2)
  assert_eq(result.contamination_counts.length(), 6)
  assert_eq(result.contamination_counts[0].length(), 8)
}

///|
test "decontX: omitted identifiers are generated deterministically" {
  let result = @src.decontx([[10.0, 1.0], [1.0, 10.0]], [], [], ["A", "B"]) catch {
    _ => abort("generated names should be accepted")
  }
  assert_eq(result.gene_names, ["Gene1", "Gene2"])
  assert_eq(result.cell_names, ["Cell1", "Cell2"])
}

///|
test "decontX: cluster order follows first appearance" {
  let result = @src.decontx([[10.0, 1.0, 5.0], [1.0, 10.0, 4.0]], [], [], [
    "B", "A", "B",
  ]) catch {
    _ => abort("valid labels should be accepted")
  }
  assert_eq(result.cluster_names, ["B", "A"])
  assert_eq(result.cluster_indices, [0, 1, 0])
}

///|
test "decontX: native and contaminant counts reconstruct observations" {
  let (counts, _, _, _, _) = @src.decontx_example_data()
  let result = decontx_test_result()
  for gene in 0..<counts.length() {
    for cell in 0..<counts[gene].length() {
      decontx_test_close(
        result.corrected_counts[gene][cell] +
        result.contamination_counts[gene][cell],
        counts[gene][cell],
        1.0e-10,
      )
    }
  }
}

///|
test "decontX: per-cell contamination estimates are bounded" {
  let result = decontx_test_result()
  for value in result.contamination {
    assert_true(value > 0.0 && value < 1.0)
  }
}

///|
test "decontX: native profiles are normalized" {
  let result = decontx_test_result()
  for profile in result.native_profiles {
    let mut total = 0.0
    for value in profile {
      total = total + value
    }
    decontx_test_close(total, 1.0, 1.0e-10)
  }
}

///|
test "decontX: contaminant profiles are normalized" {
  let result = decontx_test_result()
  for profile in result.contaminant_profiles {
    let mut total = 0.0
    for value in profile {
      total = total + value
    }
    decontx_test_close(total, 1.0, 1.0e-10)
  }
}

///|
test "decontX: B-cell marker contamination is removed from T cells" {
  let result = decontx_test_result()
  assert_true(result.corrected_counts[2][0] < 8.0)
  assert_true(result.contamination_counts[2][0] > 0.0)
}

///|
test "decontX: T-cell marker contamination is removed from B cells" {
  let result = decontx_test_result()
  assert_true(result.corrected_counts[0][4] < 8.0)
  assert_true(result.contamination_counts[0][4] > 0.0)
}

///|
test "decontX: diagnostics record bounded inference work" {
  let result = decontx_test_result()
  assert_true(result.diagnostics.iterations >= 1)
  assert_true(result.diagnostics.iterations <= 100)
  assert_eq(
    result.diagnostics.log_likelihoods.length(),
    result.diagnostics.iterations,
  )
}

///|
test "decontX: diagnostic likelihoods are finite" {
  let result = decontx_test_result()
  assert_true(result.diagnostics.initial_log_likelihood.abs() < 1.0e300)
  assert_true(result.diagnostics.final_log_likelihood.abs() < 1.0e300)
  for value in result.diagnostics.log_likelihoods {
    assert_true(value == value)
    assert_true(value.abs() < 1.0e300)
  }
}

///|
test "decontX: cell estimate reports decomposed library size" {
  let result = decontx_test_result()
  match result.cell_estimate(0) {
    Some(estimate) => {
      assert_eq(estimate.cell_name, "T1")
      assert_eq(estimate.cluster, "T")
      decontx_test_close(
        estimate.native_counts + estimate.contaminant_counts,
        estimate.library_size,
        1.0e-10,
      )
      decontx_test_close(estimate.library_size, 214.0, 1.0e-10)
    }
    None => abort("valid cell estimate should exist")
  }
}

///|
test "decontX: cell estimate validates indices" {
  let result = decontx_test_result()
  assert_true(result.cell_estimate(-1) is None)
  assert_true(result.cell_estimate(8) is None)
}

///|
test "decontX: gene and cell lookup use preserved identifiers" {
  let result = decontx_test_result()
  assert_eq(result.gene_index("B_marker_1"), Some(2))
  assert_eq(result.cell_index("B3"), Some(6))
  assert_true(result.gene_index("missing") is None)
  assert_true(result.cell_index("missing") is None)
}

///|
test "decontX: cluster estimates partition all cells" {
  let result = decontx_test_result()
  let estimates = result.cluster_estimates()
  assert_eq(estimates.length(), 2)
  assert_eq(estimates[0].cluster, "T")
  assert_eq(estimates[0].n_cells, 4)
  assert_eq(estimates[1].cluster, "B")
  assert_eq(estimates[1].n_cells, 4)
  decontx_test_close(
    estimates[0].native_counts + estimates[0].contaminant_counts,
    estimates[0].total_counts,
    1.0e-10,
  )
}

///|
test "decontX: most contaminated cells are sorted descending" {
  let result = decontx_test_result()
  let top = result.most_contaminated_cells(4)
  assert_eq(top.length(), 4)
  for index in 1..<top.length() {
    assert_true(top[index - 1].contamination >= top[index].contamination)
  }
}

///|
test "decontX: top-cell query handles zero and oversized limits" {
  let result = decontx_test_result()
  assert_eq(result.most_contaminated_cells(0).length(), 0)
  assert_eq(result.most_contaminated_cells(100).length(), 8)
}

///|
test "decontX: summary exposes model dimensions and convergence" {
  let summary = decontx_test_result().summary()
  assert_true(summary.contains("6 genes x 8 cells"))
  assert_true(summary.contains("2 clusters"))
  assert_true(summary.contains("mean contamination="))
  assert_true(summary.contains("converged="))
}

///|
test "decontX: explicit background fixes a common ambient profile" {
  let (counts, genes, cells, clusters, background) = @src.decontx_example_data()
  let result = @src.decontx(counts, genes, cells, clusters, background~) catch {
    _ => abort("valid background should be accepted")
  }
  assert_true(result.used_background)
  assert_eq(result.contaminant_profiles[0], result.contaminant_profiles[1])
}

///|
test "decontX: background changes the inferred contaminant distribution" {
  let (counts, genes, cells, clusters, background) = @src.decontx_example_data()
  let without = @src.decontx(counts, genes, cells, clusters) catch {
    _ => abort("analysis without background should succeed")
  }
  let with_background = @src.decontx(
    counts,
    genes,
    cells,
    clusters,
    background~,
  ) catch {
    _ => abort("analysis with background should succeed")
  }
  assert_true(
    (without.contaminant_profiles[0][0] -
    with_background.contaminant_profiles[0][0]).abs() >
    1.0e-6,
  )
}

///|
test "decontX: one cluster is identifiable with external background" {
  let result = @src.decontx(
    [[10.0, 12.0], [2.0, 3.0]],
    ["A", "B"],
    ["C1", "C2"],
    ["one", "one"],
    background=[[1.0], [9.0]],
  ) catch {
    _ => abort("background should identify one-cluster contamination")
  }
  assert_eq(result.n_clusters(), 1)
  assert_true(result.used_background)
}

///|
test "decontX: one cluster without background is rejected" {
  let raised = try {
    ignore(@src.decontx([[10.0, 12.0], [2.0, 3.0]], [], [], ["one", "one"]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(raised)
}

///|
test "decontX: empty and ragged matrices are rejected" {
  let empty = try {
    ignore(@src.decontx([], [], [], []))
    false
  } catch {
    DecontXError(_) => true
  }
  let ragged = try {
    ignore(@src.decontx([[1.0], [2.0, 3.0]], [], [], ["A"]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(empty)
  assert_true(ragged)
}

///|
test "decontX: negative and non-finite counts are rejected" {
  let negative = try {
    ignore(@src.decontx([[1.0, -1.0]], [], [], ["A", "B"]))
    false
  } catch {
    DecontXError(_) => true
  }
  let non_finite = try {
    ignore(@src.decontx([[1.0, @double.not_a_number]], [], [], ["A", "B"]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(negative)
  assert_true(non_finite)
}

///|
test "decontX: identifier dimensions are validated" {
  let genes = try {
    ignore(@src.decontx([[1.0, 2.0], [2.0, 1.0]], ["only_one"], [], ["A", "B"]))
    false
  } catch {
    DecontXError(_) => true
  }
  let cells = try {
    ignore(@src.decontx([[1.0, 2.0], [2.0, 1.0]], [], ["only_one"], ["A", "B"]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(genes)
  assert_true(cells)
}

///|
test "decontX: duplicate and blank identifiers are rejected" {
  let duplicate = try {
    ignore(
      @src.decontx([[1.0, 2.0], [2.0, 1.0]], ["same", "same"], [], ["A", "B"]),
    )
    false
  } catch {
    DecontXError(_) => true
  }
  let blank = try {
    ignore(@src.decontx([[1.0, 2.0], [2.0, 1.0]], [], ["C1", " "], ["A", "B"]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(duplicate)
  assert_true(blank)
}

///|
test "decontX: cluster labels are validated" {
  let length = try {
    ignore(@src.decontx([[1.0, 2.0], [2.0, 1.0]], [], [], ["A"]))
    false
  } catch {
    DecontXError(_) => true
  }
  let blank = try {
    ignore(@src.decontx([[1.0, 2.0], [2.0, 1.0]], [], [], ["A", ""]))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(length)
  assert_true(blank)
}

///|
test "decontX: background dimensions and values are validated" {
  let dimensions = try {
    ignore(
      @src.decontx([[1.0, 2.0], [2.0, 1.0]], [], [], ["A", "B"], background=[
        [1.0],
      ]),
    )
    false
  } catch {
    DecontXError(_) => true
  }
  let negative = try {
    ignore(
      @src.decontx([[1.0, 2.0], [2.0, 1.0]], [], [], ["A", "B"], background=[
        [1.0],
        [-1.0],
      ]),
    )
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(dimensions)
  assert_true(negative)
}

///|
test "decontX: zero-library cells remain zero" {
  let result = @src.decontx([[0.0, 10.0, 1.0], [0.0, 1.0, 10.0]], [], [], [
    "A", "A", "B",
  ]) catch {
    _ => abort("zero-library cells should be supported")
  }
  assert_eq(result.corrected_counts[0][0], 0.0)
  assert_eq(result.corrected_counts[1][0], 0.0)
  assert_eq(result.contamination_counts[0][0], 0.0)
  assert_eq(result.contamination_counts[1][0], 0.0)
}

///|
test "decontX: configured maximum iteration is respected" {
  let (counts, genes, cells, clusters, _) = @src.decontx_example_data()
  let config = @src.DecontXConfig::create(max_iterations=1, tolerance=1.0e-30) catch {
    _ => abort("valid iteration limit should build")
  }
  let result = @src.decontx(counts, genes, cells, clusters, config~) catch {
    _ => abort("bounded analysis should succeed")
  }
  assert_eq(result.diagnostics.iterations, 1)
  assert_true(!result.diagnostics.converged)
}

///|
test "decontX: stronger low-contamination prior shrinks estimates" {
  let (counts, genes, cells, clusters, _) = @src.decontx_example_data()
  let weak = @src.DecontXConfig::create(
    max_iterations=1,
    contamination_prior_alpha=1.0,
    contamination_prior_beta=1.0,
  ) catch {
    _ => abort("weak prior should build")
  }
  let strong = @src.DecontXConfig::create(
    max_iterations=1,
    contamination_prior_alpha=1.0,
    contamination_prior_beta=99.0,
  ) catch {
    _ => abort("strong prior should build")
  }
  let weak_result = @src.decontx(counts, genes, cells, clusters, config=weak) catch {
    _ => abort("weak-prior analysis should succeed")
  }
  let strong_result = @src.decontx(
    counts,
    genes,
    cells,
    clusters,
    config=strong,
  ) catch {
    _ => abort("strong-prior analysis should succeed")
  }
  assert_true(
    strong_result.mean_contamination() < weak_result.mean_contamination(),
  )
}

///|
test "decontX: automatic clustering finds both example populations" {
  let (counts, genes, cells, _, _) = @src.decontx_example_data()
  let result = @src.decontx_auto(counts, genes, cells, 2) catch {
    _ => abort("automatic clustering should succeed")
  }
  assert_eq(result.n_cells(), 8)
  assert_eq(result.n_clusters(), 2)
  assert_true(result.cluster_names[0].has_prefix("cluster_"))
}

///|
test "decontX: automatic clustering validates requested cluster count" {
  let (counts, genes, cells, _, _) = @src.decontx_example_data()
  let one = try {
    ignore(@src.decontx_auto(counts, genes, cells, 1))
    false
  } catch {
    DecontXError(_) => true
  }
  let too_many = try {
    ignore(@src.decontx_auto(counts, genes, cells, 9))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(one)
  assert_true(too_many)
}

///|
test "decontX: SCE integration adds assay and diagnostics" {
  let output = @src.decontx_sce(decontx_test_sce()) catch {
    _ => abort("valid SCE integration should succeed")
  }
  let corrected = @src.sce_get_assay(output.experiment, "decontXcounts")
  assert_eq(corrected.length(), 6)
  assert_eq(corrected[0].length(), 8)
  assert_eq(
    @src.sce_get_col_data(output.experiment, "decontX_contamination").length(),
    8,
  )
  assert_eq(@src.sce_get_col_data(output.experiment, "decontX_cluster"), [
    "T", "T", "T", "T", "B", "B", "B", "B",
  ])
}

///|
test "decontX: SCE integration does not mutate the source container" {
  let source = decontx_test_sce()
  let _ = @src.decontx_sce(source) catch {
    _ => abort("valid SCE integration should succeed")
  }
  assert_eq(@src.sce_get_assay(source, "decontXcounts").length(), 0)
  assert_eq(@src.sce_get_col_data(source, "decontX_contamination").length(), 0)
}

///|
test "decontX: SCE integration validates assay and cluster metadata" {
  let missing_cluster = @src.SingleCellExperiment::new(
    [[1.0, 2.0], [2.0, 1.0]],
    ["G1", "G2"],
    ["C1", "C2"],
  )
  let cluster_error = try {
    ignore(@src.decontx_sce(missing_cluster))
    false
  } catch {
    DecontXError(_) => true
  }
  let missing_assay = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("other", [[1.0, 2.0], [2.0, 1.0]])
    |> @src.SCEBuilder::add_col_data("cluster", ["A", "B"])
    |> @src.SCEBuilder::build
  let assay_error = try {
    ignore(@src.decontx_sce(missing_assay))
    false
  } catch {
    DecontXError(_) => true
  }
  assert_true(cluster_error)
  assert_true(assay_error)
}

///|
test "decontX: SCE integration supports a custom output assay" {
  let output = @src.decontx_sce(
    decontx_test_sce(),
    output_assay="ambient_corrected",
  ) catch {
    _ => abort("custom output assay should be accepted")
  }
  assert_eq(
    @src.sce_get_assay(output.experiment, "ambient_corrected").length(),
    6,
  )
  assert_eq(@src.sce_get_assay(output.experiment, "decontXcounts").length(), 0)
}
