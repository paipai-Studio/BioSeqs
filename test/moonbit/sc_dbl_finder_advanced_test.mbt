///|
fn scdfa_test_small_data() -> @src.SingleCellData {
  @src.SingleCellData::create(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [2.0, 4.0, 6.0]],
    ["c1", "c2", "c3", "c4"],
    ["g1", "g2", "g3"],
  ) catch {
    _ => abort("valid small scDblFinder data should build")
  }
}

///|
fn scdfa_test_config() -> @src.ScDblFinderConfig {
  @src.ScDblFinderConfig::create(
    expected_doublet_rate=0.15,
    artificial_doublets=48,
    selected_features=18,
    dimensions=5,
    neighbors=8,
    iterations=2,
    classifier_steps=120,
    seed=17,
  ) catch {
    _ => abort("valid scDblFinder configuration should build")
  }
}

///|
fn scdfa_test_result() -> @src.ScDblFinderResult {
  let (data, clusters, samples, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("valid scDblFinder example should build")
  }
  @src.sc_dbl_finder(
    data,
    clusters~,
    samples~,
    known_doublets=truth,
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid scDblFinder workflow should fit")
  }
}

///|
fn scdfa_test_experiment() -> @src.SingleCellExperiment {
  let (data, clusters, samples, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("valid scDblFinder example should build")
  }
  let assay : Array[Array[Double]] = []
  for gene in 0..<data.gene_names.length() {
    let row : Array[Double] = []
    for cell in 0..<data.cell_names.length() {
      row.push(data.counts[cell][gene])
    }
    assay.push(row)
  }
  let experiment = @src.SingleCellExperiment::new(
    assay,
    data.gene_names,
    data.cell_names,
  )
  experiment.col_data["cluster"] = clusters
  experiment.col_data["capture"] = samples
  let known : Array[String] = []
  for value in truth {
    known.push(if value { "doublet" } else { "singlet" })
  }
  experiment.col_data["known"] = known
  experiment
}

///|
test "scDblFinder advanced example has typed truth and captures" {
  let (data, clusters, samples, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  assert_eq(data.counts.length(), 42)
  assert_eq(data.gene_names.length(), 24)
  assert_eq(clusters.length(), 42)
  assert_eq(samples.length(), 42)
  assert_eq(truth.length(), 42)
  let mut doublets = 0
  for value in truth {
    if value {
      doublets = doublets + 1
    }
  }
  assert_eq(doublets, 6)
}

///|
test "scDblFinder SingleCellData constructor copies inputs" {
  let counts = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let cells = ["a", "b", "c"]
  let genes = ["x", "y"]
  let data = @src.SingleCellData::create(counts, cells, genes) catch {
    _ => abort("valid data should build")
  }
  counts[0][0] = 99.0
  cells[0] = "changed"
  genes[0] = "changed"
  assert_eq(data.counts[0][0], 1.0)
  assert_eq(data.cell_names[0], "a")
  assert_eq(data.gene_names[0], "x")
}

///|
test "scDblFinder rejects fewer than three cells" {
  let raised = try {
    ignore(
      @src.SingleCellData::create([[1.0, 2.0], [3.0, 4.0]], ["a", "b"], [
        "x", "y",
      ]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder rejects fewer than two genes" {
  let raised = try {
    ignore(
      @src.SingleCellData::create([[1.0], [2.0], [3.0]], ["a", "b", "c"], ["x"]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder rejects a ragged count matrix" {
  let raised = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, 2.0], [3.0], [4.0, 5.0]],
        ["a", "b", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder validates cell and gene identifiers" {
  let wrong_length = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  let empty_name = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        ["a", "", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  let duplicate_gene = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b", "c"],
        ["x", "x"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(wrong_length && empty_name && duplicate_gene)
}

///|
test "scDblFinder rejects non-finite and negative counts" {
  let negative = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, -1.0], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  let nan = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, Double::nan()], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  let infinite = try {
    ignore(
      @src.SingleCellData::create(
        [[1.0, 1.0e301], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(negative && nan && infinite)
}

///|
test "scDblFinder rejects zero-library cells" {
  let raised = try {
    ignore(
      @src.SingleCellData::create(
        [[0.0, 0.0], [3.0, 4.0], [5.0, 6.0]],
        ["a", "b", "c"],
        ["x", "y"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder default configuration matches portable workflow" {
  let config = @src.ScDblFinderConfig::default()
  assert_eq(config.expected_doublet_rate, -1.0)
  assert_eq(config.doublet_rate_per_1000, 0.008)
  assert_eq(config.selected_features, 1000)
  assert_eq(config.dimensions, 20)
  assert_eq(config.neighbors, 15)
  assert_eq(config.iterations, 3)
  assert_eq(config.classifier_steps, 300)
  assert_eq(config.cluster_count, 0)
}

///|
test "scDblFinder configuration preserves custom values" {
  let config = @src.ScDblFinderConfig::create(
    expected_doublet_rate=0.12,
    rate_uncertainty=0.03,
    artificial_doublets=40,
    selected_features=12,
    dimensions=4,
    neighbors=7,
    iterations=2,
    classifier_steps=50,
    random_fraction=0.2,
    half_size_fraction=0.4,
    cluster_count=3,
    seed=9,
  ) catch {
    _ => abort("valid custom configuration should build")
  }
  assert_eq(config.expected_doublet_rate, 0.12)
  assert_eq(config.artificial_doublets, 40)
  assert_eq(config.dimensions, 4)
  assert_eq(config.cluster_count, 3)
  assert_eq(config.seed, 9)
}

///|
test "scDblFinder configuration rejects invalid rates" {
  let expected = try {
    ignore(@src.ScDblFinderConfig::create(expected_doublet_rate=1.1))
    false
  } catch {
    _ => true
  }
  let per_thousand = try {
    ignore(@src.ScDblFinderConfig::create(doublet_rate_per_1000=-0.1))
    false
  } catch {
    _ => true
  }
  let uncertainty = try {
    ignore(@src.ScDblFinderConfig::create(rate_uncertainty=-0.2))
    false
  } catch {
    _ => true
  }
  assert_true(expected && per_thousand && uncertainty)
}

///|
test "scDblFinder configuration rejects invalid dimensions" {
  let features = try {
    ignore(@src.ScDblFinderConfig::create(selected_features=1))
    false
  } catch {
    _ => true
  }
  let dimensions = try {
    ignore(@src.ScDblFinderConfig::create(dimensions=0))
    false
  } catch {
    _ => true
  }
  let neighbors = try {
    ignore(@src.ScDblFinderConfig::create(neighbors=0))
    false
  } catch {
    _ => true
  }
  assert_true(features && dimensions && neighbors)
}

///|
test "scDblFinder configuration validates classifier controls" {
  let iterations = try {
    ignore(@src.ScDblFinderConfig::create(iterations=0))
    false
  } catch {
    _ => true
  }
  let steps = try {
    ignore(@src.ScDblFinderConfig::create(classifier_steps=0))
    false
  } catch {
    _ => true
  }
  let learning_rate = try {
    ignore(@src.ScDblFinderConfig::create(learning_rate=1.1))
    false
  } catch {
    _ => true
  }
  let regularization = try {
    ignore(@src.ScDblFinderConfig::create(regularization=-0.01))
    false
  } catch {
    _ => true
  }
  assert_true(iterations && steps && learning_rate && regularization)
}

///|
test "scDblFinder configuration validates threshold and mixing fractions" {
  let stringency = try {
    ignore(@src.ScDblFinderConfig::create(stringency=1.0))
    false
  } catch {
    _ => true
  }
  let random_fraction = try {
    ignore(@src.ScDblFinderConfig::create(random_fraction=-0.1))
    false
  } catch {
    _ => true
  }
  let half_fraction = try {
    ignore(@src.ScDblFinderConfig::create(half_size_fraction=1.1))
    false
  } catch {
    _ => true
  }
  let unidentifiable = try {
    ignore(@src.ScDblFinderConfig::create(unidentifiable_threshold=1.1))
    false
  } catch {
    _ => true
  }
  assert_true(stringency && random_fraction && half_fraction && unidentifiable)
}

///|
test "scDblFinder configuration rejects a single requested cluster" {
  let raised = try {
    ignore(@src.ScDblFinderConfig::create(cluster_count=1))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder expected rate follows cells per thousand" {
  let rate = @src.scdf_expected_doublet_rate(5000) catch {
    _ => abort("valid expected rate should compute")
  }
  assert_true((rate - 0.04).abs() < 1.0e-12)
}

///|
test "scDblFinder expected rate is capped and validated" {
  let capped = @src.scdf_expected_doublet_rate(2000, rate_per_1000=1.0) catch {
    _ => abort("valid capped rate should compute")
  }
  assert_eq(capped, 1.0)
  let raised = try {
    ignore(@src.scdf_expected_doublet_rate(0))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder homotypic proportion uses squared cluster frequencies" {
  let proportion = @src.scdf_homotypic_proportion(["A", "A", "B", "B"]) catch {
    _ => abort("valid homotypic proportion should compute")
  }
  assert_true((proportion - 0.5).abs() < 1.0e-12)
}

///|
test "scDblFinder homotypic proportion validates labels" {
  let empty = try {
    ignore(@src.scdf_homotypic_proportion([]))
    false
  } catch {
    _ => true
  }
  let blank = try {
    ignore(@src.scdf_homotypic_proportion(["A", ""]))
    false
  } catch {
    _ => true
  }
  assert_true(empty && blank)
}

///|
test "scDblFinder random artificial doublets have requested dimensions" {
  let artificial = @src.scdf_generate_artificial_doublets(
    scdfa_test_small_data(),
    7,
    random_fraction=1.0,
    half_size_fraction=0.0,
    seed=3,
  ) catch {
    _ => abort("valid artificial doublets should build")
  }
  assert_eq(artificial.counts.length(), 7)
  assert_eq(artificial.names.length(), 7)
  assert_eq(artificial.parent_one.length(), 7)
  assert_eq(artificial.origins, ["", "", "", "", "", "", ""])
  for row in artificial.counts {
    assert_eq(row.length(), 3)
  }
}

///|
test "scDblFinder artificial counts sum parent profiles" {
  let data = scdfa_test_small_data()
  let artificial = @src.scdf_generate_artificial_doublets(
    data,
    5,
    half_size_fraction=0.0,
    seed=4,
  ) catch {
    _ => abort("valid artificial doublets should build")
  }
  for index in 0..<artificial.counts.length() {
    for gene in 0..<data.gene_names.length() {
      assert_eq(
        artificial.counts[index][gene],
        data.counts[artificial.parent_one[index]][gene] +
        data.counts[artificial.parent_two[index]][gene],
      )
    }
  }
}

///|
test "scDblFinder half-size artificial doublets are scaled" {
  let data = scdfa_test_small_data()
  let artificial = @src.scdf_generate_artificial_doublets(
    data,
    4,
    half_size_fraction=0.5,
    seed=2,
  ) catch {
    _ => abort("valid half-size artificial doublets should build")
  }
  for index in 0..<4 {
    let expected_scale = if index < 2 { 0.5 } else { 1.0 }
    let expected = (
        data.counts[artificial.parent_one[index]][0] +
        data.counts[artificial.parent_two[index]][0]
      ) *
      expected_scale
    assert_eq(artificial.counts[index][0], expected)
  }
}

///|
test "scDblFinder artificial parents are distinct and in bounds" {
  let artificial = @src.scdf_generate_artificial_doublets(
    scdfa_test_small_data(),
    20,
    seed=8,
  ) catch {
    _ => abort("valid artificial doublets should build")
  }
  for index in 0..<artificial.counts.length() {
    assert_true(artificial.parent_one[index] >= 0)
    assert_true(artificial.parent_two[index] < 4)
    assert_true(artificial.parent_one[index] != artificial.parent_two[index])
  }
}

///|
test "scDblFinder cluster-aware artificial doublets use cross-cluster pairs" {
  let data = scdfa_test_small_data()
  let clusters = ["A", "A", "B", "B"]
  let artificial = @src.scdf_generate_artificial_doublets(
    data,
    12,
    clusters~,
    random_fraction=0.0,
    half_size_fraction=0.0,
  ) catch {
    _ => abort("valid cluster-aware artificial doublets should build")
  }
  for index in 0..<artificial.counts.length() {
    assert_true(
      clusters[artificial.parent_one[index]] !=
      clusters[artificial.parent_two[index]],
    )
    assert_eq(artificial.origins[index], "A+B")
  }
}

///|
test "scDblFinder artificial generation is deterministic for a seed" {
  let data = scdfa_test_small_data()
  let first = @src.scdf_generate_artificial_doublets(data, 10, seed=5) catch {
    _ => abort("valid artificial doublets should build")
  }
  let second = @src.scdf_generate_artificial_doublets(data, 10, seed=5) catch {
    _ => abort("valid artificial doublets should build")
  }
  assert_eq(first.counts, second.counts)
  assert_eq(first.parent_one, second.parent_one)
  assert_eq(first.parent_two, second.parent_two)
}

///|
test "scDblFinder artificial generation responds to its seed" {
  let data = scdfa_test_small_data()
  let first = @src.scdf_generate_artificial_doublets(data, 10, seed=1) catch {
    _ => abort("valid artificial doublets should build")
  }
  let second = @src.scdf_generate_artificial_doublets(data, 10, seed=2) catch {
    _ => abort("valid artificial doublets should build")
  }
  assert_true(
    first.parent_one != second.parent_one ||
    first.parent_two != second.parent_two,
  )
}

///|
test "scDblFinder artificial generation validates arguments" {
  let data = scdfa_test_small_data()
  let number = try {
    ignore(@src.scdf_generate_artificial_doublets(data, 0))
    false
  } catch {
    _ => true
  }
  let clusters = try {
    ignore(@src.scdf_generate_artificial_doublets(data, 2, clusters=["A", "B"]))
    false
  } catch {
    _ => true
  }
  let fraction = try {
    ignore(@src.scdf_generate_artificial_doublets(data, 2, random_fraction=1.1))
    false
  } catch {
    _ => true
  }
  assert_true(number && clusters && fraction)
}

///|
test "scDblFinder threshold separates ideal score classes" {
  let threshold = @src.scdf_optimize_threshold(
    [0.05, 0.10, 0.15, 0.20],
    [0.80, 0.90, 0.95],
    0.0,
    uncertainty=0.0,
  ) catch {
    _ => abort("valid threshold should optimize")
  }
  assert_true(threshold.threshold > 0.20)
  assert_true(threshold.threshold < 0.80)
  assert_eq(threshold.false_positive_rate, 0.0)
  assert_eq(threshold.false_negative_rate, 0.0)
}

///|
test "scDblFinder threshold diagnostics remain bounded with overlap" {
  let threshold = @src.scdf_optimize_threshold(
    [0.1, 0.3, 0.6, 0.8],
    [0.4, 0.5, 0.7, 0.9],
    0.25,
  ) catch {
    _ => abort("valid threshold should optimize")
  }
  assert_true(threshold.threshold >= 0.0 && threshold.threshold <= 1.0)
  assert_true(threshold.called_rate >= 0.0 && threshold.called_rate <= 1.0)
  assert_true(
    threshold.false_negative_rate >= 0.0 && threshold.false_negative_rate <= 1.0,
  )
  assert_true(threshold.cost >= 0.0)
}

///|
test "scDblFinder threshold optimization is deterministic" {
  let first = @src.scdf_optimize_threshold(
    [0.1, 0.2, 0.7],
    [0.4, 0.8, 0.9],
    0.2,
  ) catch {
    _ => abort("valid threshold should optimize")
  }
  let second = @src.scdf_optimize_threshold(
    [0.1, 0.2, 0.7],
    [0.4, 0.8, 0.9],
    0.2,
  ) catch {
    _ => abort("valid threshold should optimize")
  }
  assert_eq(first.threshold, second.threshold)
  assert_eq(first.cost, second.cost)
}

///|
test "scDblFinder threshold requires both classes" {
  let raised = try {
    ignore(@src.scdf_optimize_threshold([], [0.8], 0.1))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder threshold rejects invalid scores" {
  let nan = try {
    ignore(@src.scdf_optimize_threshold([Double::nan()], [0.8], 0.1))
    false
  } catch {
    _ => true
  }
  let outside = try {
    ignore(@src.scdf_optimize_threshold([0.2], [1.1], 0.1))
    false
  } catch {
    _ => true
  }
  assert_true(nan && outside)
}

///|
test "scDblFinder threshold validates rate and stringency" {
  let rate = try {
    ignore(@src.scdf_optimize_threshold([0.2], [0.8], 1.1))
    false
  } catch {
    _ => true
  }
  let stringency = try {
    ignore(@src.scdf_optimize_threshold([0.2], [0.8], 0.1, stringency=0.0))
    false
  } catch {
    _ => true
  }
  assert_true(rate && stringency)
}

///|
test "scDblFinder fast clustering labels every cell" {
  let (data, _, _, _) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let clusters = @src.scdf_fast_cluster(
    data,
    3,
    dimensions=5,
    selected_features=18,
  ) catch {
    _ => abort("valid fast clustering should run")
  }
  assert_eq(clusters.length(), data.cell_names.length())
  for label in clusters {
    assert_true(label.has_prefix("cluster_"))
  }
}

///|
test "scDblFinder fast clustering is deterministic" {
  let (data, _, _, _) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let first = @src.scdf_fast_cluster(data, 3) catch {
    _ => abort("valid fast clustering should run")
  }
  let second = @src.scdf_fast_cluster(data, 3) catch {
    _ => abort("valid fast clustering should run")
  }
  assert_eq(first, second)
}

///|
test "scDblFinder fast clustering validates cluster count" {
  let data = scdfa_test_small_data()
  let one = try {
    ignore(@src.scdf_fast_cluster(data, 1))
    false
  } catch {
    _ => true
  }
  let too_many = try {
    ignore(@src.scdf_fast_cluster(data, 5))
    false
  } catch {
    _ => true
  }
  assert_true(one && too_many)
}

///|
test "scDblFinder workflow returns aligned cell-level arrays" {
  let result = scdfa_test_result()
  assert_eq(result.n_cells(), 42)
  assert_eq(result.scores.length(), 42)
  assert_eq(result.calls.length(), 42)
  assert_eq(result.weighted_ratios.length(), 42)
  assert_eq(result.neighbor_ratios.length(), 42)
  assert_eq(result.cxds_scores.length(), 42)
  assert_eq(result.origins.length(), 42)
  assert_eq(result.library_sizes.length(), 42)
}

///|
test "scDblFinder classifier scores are finite probabilities" {
  let result = scdfa_test_result()
  for score in result.scores {
    assert_true(score == score)
    assert_true(score >= 0.0 && score <= 1.0)
  }
}

///|
test "scDblFinder neighborhood and cxds features are bounded" {
  let result = scdfa_test_result()
  for cell in 0..<result.n_cells() {
    assert_true(
      result.weighted_ratios[cell] >= 0.0 && result.weighted_ratios[cell] <= 1.0,
    )
    assert_true(
      result.neighbor_ratios[cell] >= 0.0 && result.neighbor_ratios[cell] <= 1.0,
    )
    assert_true(
      result.cxds_scores[cell] >= 0.0 && result.cxds_scores[cell] <= 1.0,
    )
    assert_true(
      result.origin_ambiguities[cell] >= 0.0 &&
      result.origin_ambiguities[cell] <= 1.0,
    )
  }
}

///|
test "scDblFinder kNN distances and classes are populated" {
  let result = scdfa_test_result()
  for cell in 0..<result.n_cells() {
    assert_true(result.distances_to_nearest[cell] >= 0.0)
    assert_true(result.distances_to_doublet[cell] >= 0.0)
    assert_true(result.distances_to_real[cell] >= 0.0)
    assert_true(
      result.nearest_classes[cell] == "cell" ||
      result.nearest_classes[cell] == "artificialDoublet",
    )
  }
}

///|
test "scDblFinder PCA has requested portable dimensions" {
  let result = scdfa_test_result()
  assert_eq(result.pca.length(), 42)
  for coordinates in result.pca {
    assert_eq(coordinates.length(), 5)
    for value in coordinates {
      assert_true(value == value && value.abs() <= 1.0e300)
    }
  }
}

///|
test "scDblFinder feature selection records selected genes" {
  let result = scdfa_test_result()
  assert_true(result.selected_gene_names.length() > 0)
  assert_true(result.selected_gene_names.length() <= result.gene_names.length())
  for gene in result.selected_gene_names {
    assert_true(result.gene_names.contains(gene))
  }
}

///|
test "scDblFinder performs independent capture fits" {
  let result = scdfa_test_result()
  assert_eq(result.sample_ids, ["capture_A", "capture_B"])
  assert_eq(result.thresholds.length(), 2)
  assert_eq(result.expected_doublet_rates.length(), 2)
  assert_eq(result.artificial_doublets, 96)
  for threshold in result.thresholds {
    assert_true(threshold.threshold >= 0.0 && threshold.threshold <= 1.0)
  }
}

///|
test "scDblFinder preserves known labels and workflow configuration" {
  let result = scdfa_test_result()
  let mut known = 0
  for value in result.known_doublets {
    if value {
      known = known + 1
    }
  }
  assert_eq(known, 6)
  assert_eq(result.classifier_iterations, 2)
  assert_eq(result.config.seed, 17)
}

///|
test "scDblFinder ranks top cells by descending score" {
  let result = scdfa_test_result()
  let top = result.top_doublets(limit=7) catch {
    _ => abort("valid top-doublet query should run")
  }
  assert_eq(top.length(), 7)
  for index in 1..<top.length() {
    assert_true(top[index - 1].score >= top[index].score)
  }
}

///|
test "scDblFinder cell lookup returns typed diagnostics" {
  let result = scdfa_test_result()
  match result.cell("known_doublet_1") {
    Some(cell) => {
      assert_eq(cell.cell_name, "known_doublet_1")
      assert_true(cell.score >= 0.0 && cell.score <= 1.0)
      assert_true(cell.origin != "")
    }
    None => abort("known cell should be found")
  }
  assert_true(result.cell("missing") is None)
}

///|
test "scDblFinder summary reports workflow metadata" {
  let summary = scdfa_test_result().summary()
  assert_true(summary.contains("scDblFinder Result"))
  assert_true(summary.contains("cells=42"))
  assert_true(summary.contains("samples=2"))
  assert_true(summary.contains("iterations=2"))
}

///|
test "scDblFinder singlet filtering follows calls" {
  let (data, _, _, _) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let result = scdfa_test_result()
  let filtered = result.filter_singlets(data) catch {
    _ => abort("valid singlet filtering should run")
  }
  assert_eq(
    filtered.cell_names.length(),
    result.n_cells() - result.n_doublets(),
  )
  assert_eq(filtered.gene_names, data.gene_names)
}

///|
test "scDblFinder synthetic classification exceeds baseline accuracy" {
  let (_, _, _, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let result = scdfa_test_result()
  let mut correct = 0
  for cell in 0..<truth.length() {
    if truth[cell] == result.calls[cell] {
      correct = correct + 1
    }
  }
  assert_true(correct.to_double() / truth.length().to_double() >= 0.75)
}

///|
test "scDblFinder workflow is deterministic" {
  let first = scdfa_test_result()
  let second = scdfa_test_result()
  assert_eq(first.scores, second.scores)
  assert_eq(first.calls, second.calls)
  assert_eq(first.origins, second.origins)
  assert_eq(first.pca, second.pca)
}

///|
test "scDblFinder supports a single combined capture" {
  let (data, clusters, _, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let result = @src.sc_dbl_finder(
    data,
    clusters~,
    known_doublets=truth,
    config=scdfa_test_config(),
  ) catch {
    _ => abort("single-capture workflow should fit")
  }
  assert_eq(result.sample_ids, ["all"])
  assert_eq(result.thresholds.length(), 1)
  assert_eq(result.artificial_doublets, 48)
}

///|
test "scDblFinder can infer clusters before fitting" {
  let (data, _, _, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let config = @src.ScDblFinderConfig::create(
    expected_doublet_rate=0.15,
    artificial_doublets=36,
    selected_features=18,
    dimensions=5,
    neighbors=8,
    iterations=1,
    classifier_steps=80,
    cluster_count=3,
  ) catch {
    _ => abort("valid auto-cluster configuration should build")
  }
  let result = @src.sc_dbl_finder(data, known_doublets=truth, config~) catch {
    _ => abort("auto-cluster workflow should fit")
  }
  assert_eq(result.clusters.length(), 42)
  for cluster in result.clusters {
    assert_true(cluster != "")
  }
}

///|
test "scDblFinder workflow validates metadata lengths" {
  let data = scdfa_test_small_data()
  let clusters = try {
    ignore(@src.sc_dbl_finder(data, clusters=["A", "B"]))
    false
  } catch {
    _ => true
  }
  let samples = try {
    ignore(@src.sc_dbl_finder(data, samples=["one"]))
    false
  } catch {
    _ => true
  }
  let known = try {
    ignore(@src.sc_dbl_finder(data, known_doublets=[true]))
    false
  } catch {
    _ => true
  }
  assert_true(clusters && samples && known)
}

///|
test "scDblFinder rejects captures with fewer than three cells" {
  let data = scdfa_test_small_data()
  let raised = try {
    ignore(
      @src.sc_dbl_finder(data, samples=["small", "small", "other", "other"]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder cluster-aware fitting assigns origins" {
  let result = scdfa_test_result()
  for origin in result.origins {
    assert_true(origin != "")
    assert_true(origin.contains("+"))
  }
}

///|
test "scDblFinder pairwise enrichment covers cluster pairs" {
  let enrichment = @src.scdf_pairwise_enrichment(scdfa_test_result()) catch {
    _ => abort("valid enrichment should compute")
  }
  assert_eq(enrichment.length(), 3)
  for index in 0..<enrichment.length() {
    let row = enrichment[index]
    assert_true(row.combination.contains("+"))
    assert_true(row.observed >= 0)
    assert_true(row.expected >= 0.0)
    assert_true(row.p_value >= 0.0 && row.p_value <= 1.0)
    assert_true(
      row.adjusted_p_value >= row.p_value && row.adjusted_p_value <= 1.0,
    )
    if index > 0 {
      assert_true(enrichment[index - 1].p_value <= row.p_value)
    }
  }
}

///|
test "scDblFinder enrichment requires cluster-aware fitting" {
  let (data, _, _, truth) = @src.scdf_create_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let result = @src.sc_dbl_finder(
    data,
    known_doublets=truth,
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid cluster-free workflow should fit")
  }
  let raised = try {
    ignore(@src.scdf_pairwise_enrichment(result))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder SingleCellExperiment writes cell diagnostics" {
  let output = @src.sc_dbl_finder_single_cell_experiment(
    scdfa_test_experiment(),
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known",
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid SingleCellExperiment workflow should fit")
  }
  assert_eq(output.result.n_cells(), 42)
  assert_eq(output.experiment.col_data["scDblFinder.score"].length(), 42)
  assert_eq(output.experiment.col_data["scDblFinder.class"].length(), 42)
  assert_eq(
    output.experiment.col_data["scDblFinder.mostLikelyOrigin"].length(),
    42,
  )
}

///|
test "scDblFinder SingleCellExperiment writes feature and PCA diagnostics" {
  let output = @src.sc_dbl_finder_single_cell_experiment(
    scdfa_test_experiment(),
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known",
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid SingleCellExperiment workflow should fit")
  }
  assert_eq(output.experiment.row_data["scDblFinder.selected"].length(), 24)
  assert_eq(output.experiment.reduced_dims["scDblFinder.PCA"].length(), 42)
  assert_true(output.experiment.metadata.contains("scDblFinder.doublets"))
  assert_true(output.experiment.metadata.contains("scDblFinder.rate"))
  assert_true(output.experiment.metadata.contains("scDblFinder.artificial"))
}

///|
test "scDblFinder SingleCellExperiment leaves the input unchanged" {
  let experiment = scdfa_test_experiment()
  let output = @src.sc_dbl_finder_single_cell_experiment(
    experiment,
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known",
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid SingleCellExperiment workflow should fit")
  }
  assert_true(!experiment.col_data.contains("scDblFinder.score"))
  output.experiment.col_data["capture"][0] = "changed"
  assert_eq(experiment.col_data["capture"][0], "capture_A")
}

///|
test "scDblFinder SingleCellExperiment supports a custom output prefix" {
  let output = @src.sc_dbl_finder_single_cell_experiment(
    scdfa_test_experiment(),
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known",
    output_prefix="doublet",
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid custom-prefix workflow should fit")
  }
  assert_true(output.experiment.col_data.contains("doublet.score"))
  assert_true(output.experiment.row_data.contains("doublet.selected"))
  assert_true(output.experiment.reduced_dims.contains("doublet.PCA"))
  assert_true(!output.experiment.col_data.contains("scDblFinder.score"))
}

///|
test "scDblFinder SingleCellExperiment transposes gene-by-cell assays" {
  let experiment = scdfa_test_experiment()
  let output = @src.sc_dbl_finder_single_cell_experiment(
    experiment,
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known",
    config=scdfa_test_config(),
  ) catch {
    _ => abort("valid SingleCellExperiment workflow should fit")
  }
  let mut expected_library = 0.0
  for gene in 0..<experiment.row_names.length() {
    expected_library = expected_library + experiment.assays["counts"][gene][0]
  }
  assert_true(
    (output.result.library_sizes[0] - expected_library).abs() < 1.0e-12,
  )
}

///|
test "scDblFinder SingleCellExperiment validates assay and columns" {
  let experiment = scdfa_test_experiment()
  let assay = try {
    ignore(
      @src.sc_dbl_finder_single_cell_experiment(
        experiment,
        assay_name="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  let column = try {
    ignore(
      @src.sc_dbl_finder_single_cell_experiment(
        experiment,
        cluster_column="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(assay && column)
}

///|
test "scDblFinder SingleCellExperiment rejects invalid known labels" {
  let experiment = scdfa_test_experiment()
  experiment.col_data["known"][0] = "unknown"
  let raised = try {
    ignore(
      @src.sc_dbl_finder_single_cell_experiment(
        experiment,
        known_doublet_column="known",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "scDblFinder SingleCellExperiment validates assay dimensions" {
  let experiment = @src.SingleCellExperiment::new(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    ["only_one_name"],
    ["c1", "c2", "c3"],
  )
  let raised = try {
    ignore(@src.sc_dbl_finder_single_cell_experiment(experiment))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
