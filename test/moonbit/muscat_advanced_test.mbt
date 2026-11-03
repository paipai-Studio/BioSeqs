///|
fn muscat_adv_test_data() -> @src.MuscatAdvancedData {
  @src.muscat_advanced_example() catch {
    _ => abort("valid muscat advanced example should build")
  }
}

///|
fn muscat_adv_test_config() -> @src.MuscatAdvancedConfig {
  @src.MuscatAdvancedConfig::create(
    min_cells=10,
    min_count=1.0,
    min_samples=2,
    max_iterations=100,
    tolerance=1.0e-8,
    minimum_dispersion=1.0e-8,
    maximum_dispersion=50.0,
    dispersion_prior_df=10.0,
    ridge=1.0e-6,
    fdr_threshold=0.1,
    detection_filter=0.9,
  ) catch {
    _ => abort("valid muscat advanced configuration should build")
  }
}

///|
fn muscat_adv_test_sum() -> @src.MuscatAdvancedPseudoBulk {
  @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_sum_counts(),
  ) catch {
    _ => abort("sum-count pseudobulk should build")
  }
}

///|
fn muscat_adv_test_detection() -> @src.MuscatAdvancedPseudoBulk {
  @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_number_detected(),
  ) catch {
    _ => abort("number-detected pseudobulk should build")
  }
}

///|
fn muscat_adv_test_design() -> @src.MuscatAdvancedDesign {
  @src.muscat_group_design_advanced(muscat_adv_test_sum(), reference="ctrl") catch {
    _ => abort("valid muscat group design should build")
  }
}

///|
fn muscat_adv_test_contrasts() -> Array[@src.MuscatAdvancedContrast] {
  @src.muscat_default_contrasts_advanced(muscat_adv_test_design()) catch {
    _ => abort("valid muscat contrast should build")
  }
}

///|
fn muscat_adv_test_ds() -> @src.MuscatAdvancedResults {
  @src.muscat_pbds_advanced(
    muscat_adv_test_sum(),
    muscat_adv_test_design(),
    muscat_adv_test_contrasts(),
    config=muscat_adv_test_config(),
  ) catch {
    _ => abort("valid muscat DS model should fit")
  }
}

///|
fn muscat_adv_test_dd() -> @src.MuscatAdvancedResults {
  @src.muscat_pbdd_advanced(
    muscat_adv_test_detection(),
    muscat_adv_test_design(),
    muscat_adv_test_contrasts(),
    config=muscat_adv_test_config(),
  ) catch {
    _ => abort("valid muscat DD model should fit")
  }
}

///|
fn muscat_adv_test_stagewise() -> @src.MuscatStagewiseResults {
  @src.muscat_stagewise_ds_dd(
    muscat_adv_test_ds(),
    muscat_adv_test_dd(),
    alpha=0.1,
  ) catch {
    _ => abort("valid muscat stagewise analysis should run")
  }
}

///|
fn muscat_adv_test_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e299
}

///|
fn muscat_adv_test_sce() -> @src.SingleCellExperiment {
  let data = muscat_adv_test_data()
  let experiment = @src.SingleCellExperiment::new(
    data.counts,
    data.gene_names,
    data.cell_names,
  )
  experiment.col_data["sample"] = data.sample_ids.copy()
  experiment.col_data["cluster"] = data.cluster_ids.copy()
  experiment.col_data["group"] = data.group_ids.copy()
  experiment.col_data["batch"] = Array::make(data.cell_names.length(), "one")
  experiment.row_data["symbol"] = data.gene_names.copy()
  experiment.reduced_dims["PCA"] = [[1.0, 2.0], [3.0, 4.0]]
  experiment.metadata["source"] = "test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[1.0, 2.0]],
    ["spike1"],
    ["s1", "s2"],
  )
  experiment
}

///|
test "muscat advanced default configuration follows pseudobulk defaults" {
  let config = @src.MuscatAdvancedConfig::default()
  assert_eq(config.min_cells, 10)
  assert_eq(config.min_count, 1.0)
  assert_eq(config.min_samples, 2)
  assert_eq(config.max_iterations, 100)
  assert_eq(config.tolerance, 1.0e-8)
  assert_eq(config.minimum_dispersion, 1.0e-8)
  assert_eq(config.maximum_dispersion, 100.0)
  assert_eq(config.dispersion_prior_df, 10.0)
  assert_eq(config.ridge, 1.0e-8)
  assert_eq(config.fdr_threshold, 0.05)
  assert_eq(config.lfc_threshold, 0.0)
  assert_eq(config.detection_filter, 0.9)
}

///|
test "muscat advanced configuration preserves explicit controls" {
  let config = @src.MuscatAdvancedConfig::create(
    min_cells=3,
    min_count=2.0,
    min_samples=3,
    max_iterations=40,
    tolerance=1.0e-5,
    minimum_dispersion=1.0e-4,
    maximum_dispersion=20.0,
    dispersion_prior_df=5.0,
    ridge=1.0e-5,
    fdr_threshold=0.2,
    lfc_threshold=1.0,
    detection_filter=0.8,
  ) catch {
    _ => abort("explicit muscat controls should be valid")
  }
  assert_eq(config.min_cells, 3)
  assert_eq(config.min_count, 2.0)
  assert_eq(config.min_samples, 3)
  assert_eq(config.max_iterations, 40)
  assert_eq(config.maximum_dispersion, 20.0)
  assert_eq(config.dispersion_prior_df, 5.0)
  assert_eq(config.fdr_threshold, 0.2)
  assert_eq(config.lfc_threshold, 1.0)
  assert_eq(config.detection_filter, 0.8)
}

///|
test "muscat advanced rejects invalid filtering controls" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedConfig::create(min_cells=-1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(min_count=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(min_count=0.0 / 0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(min_samples=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "muscat advanced rejects invalid iteration controls" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedConfig::create(max_iterations=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(tolerance=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(tolerance=0.0 / 0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced rejects invalid dispersion controls" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedConfig::create(minimum_dispersion=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedConfig::create(
      minimum_dispersion=2.0,
      maximum_dispersion=1.0,
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(dispersion_prior_df=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(ridge=-1.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "muscat advanced rejects invalid testing controls" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedConfig::create(fdr_threshold=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(fdr_threshold=1.1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(lfc_threshold=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(detection_filter=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedConfig::create(detection_filter=1.1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 5)
}

///|
test "muscat advanced aggregation helpers expose all upstream statistics" {
  assert_eq(@src.muscat_sum_counts(), @src.muscat_sum_counts())
  assert_eq(@src.muscat_mean_expression(), @src.muscat_mean_expression())
  assert_eq(@src.muscat_median_expression(), @src.muscat_median_expression())
  assert_eq(
    @src.muscat_proportion_detected(),
    @src.muscat_proportion_detected(),
  )
  assert_eq(@src.muscat_number_detected(), @src.muscat_number_detected())
  assert_true(@src.muscat_sum_counts() != @src.muscat_mean_expression())
}

///|
test "muscat advanced example exposes gene by cell metadata contract" {
  let data = muscat_adv_test_data()
  assert_eq(data.counts.length(), 10)
  assert_eq(data.counts[0].length(), 160)
  assert_eq(data.gene_names.length(), 10)
  assert_eq(data.cell_names.length(), 160)
  assert_eq(data.sample_names, ["C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"])
  assert_eq(data.cluster_names, ["A", "B"])
  assert_eq(data.group_names, ["ctrl", "stim"])
  assert_eq(data.sample_groups, [
    "ctrl", "ctrl", "ctrl", "ctrl", "stim", "stim", "stim", "stim",
  ])
}

///|
test "muscat advanced validates count matrix dimensions" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedData::create([], [], [], [])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedData::create([[]], [], [], [])) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create(
      [[1.0, 2.0], [1.0]],
      ["s1", "s2"],
      ["A", "A"],
      ["c", "t"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced validates count values" {
  let mut failures = 0
  ignore(
    @src.MuscatAdvancedData::create([[-1.0, 2.0]], ["s1", "s2"], ["A", "A"], [
      "c", "t",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create([[1.5, 2.0]], ["s1", "s2"], ["A", "A"], [
      "c", "t",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create(
      [[0.0 / 0.0, 2.0]],
      ["s1", "s2"],
      ["A", "A"],
      ["c", "t"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced validates cell metadata dimensions and values" {
  let mut failures = 0
  ignore(
    @src.MuscatAdvancedData::create([[1.0, 2.0]], ["s1"], ["A", "A"], ["c", "t"]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create([[1.0, 2.0]], ["s1", "s2"], ["A", " "], [
      "c", "t",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "muscat advanced requires replicated sample and group identities" {
  let mut failures = 0
  ignore(
    @src.MuscatAdvancedData::create([[1.0, 2.0]], ["s1", "s1"], ["A", "A"], [
      "c", "t",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create([[1.0, 2.0]], ["s1", "s2"], ["A", "A"], [
      "c", "c",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "muscat advanced requires each sample to have one group" {
  let failed = try {
    ignore(
      @src.MuscatAdvancedData::create(
        [[1.0, 2.0, 3.0, 4.0]],
        ["s1", "s1", "s2", "s2"],
        ["A", "A", "A", "A"],
        ["c", "t", "t", "t"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced validates gene and cell names" {
  let mut failures = 0
  ignore(
    @src.MuscatAdvancedData::create(
      [[1.0, 2.0], [2.0, 3.0]],
      ["s1", "s2"],
      ["A", "A"],
      ["c", "t"],
      gene_names=["one"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create(
      [[1.0, 2.0]],
      ["s1", "s2"],
      ["A", "A"],
      ["c", "t"],
      gene_names=[" "],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedData::create(
      [[1.0, 2.0]],
      ["s1", "s2"],
      ["A", "A"],
      ["c", "t"],
      cell_names=["cell", "cell"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced data owns copies of caller arrays" {
  let counts = [[1.0, 2.0]]
  let samples = ["s1", "s2"]
  let clusters = ["A", "A"]
  let groups = ["c", "t"]
  let data = @src.MuscatAdvancedData::create(counts, samples, clusters, groups) catch {
    _ => abort("valid copied data should build")
  }
  data.counts[0][0] = 99.0
  data.sample_ids[0] = "changed"
  assert_eq(counts[0][0], 1.0)
  assert_eq(samples[0], "s1")
}

///|
test "muscat advanced sum pseudobulk has cluster gene sample orientation" {
  let pseudobulk = muscat_adv_test_sum()
  assert_eq(pseudobulk.values.length(), 2)
  assert_eq(pseudobulk.values[0].length(), 10)
  assert_eq(pseudobulk.values[0][0].length(), 8)
  assert_eq(pseudobulk.cluster_names, ["A", "B"])
  assert_eq(pseudobulk.cell_counts[0], [10, 10, 10, 10, 10, 10, 10, 10])
  assert_eq(pseudobulk.aggregation, @src.muscat_sum_counts())
  assert_false(pseudobulk.scaled_cpm)
}

///|
test "muscat advanced sum aggregation preserves designed signals" {
  let pseudobulk = muscat_adv_test_sum()
  assert_eq(pseudobulk.values[0][0][0], 1009.0)
  assert_eq(pseudobulk.values[0][1][0], 39.0)
  assert_eq(pseudobulk.values[0][1][4], 189.0)
  assert_eq(pseudobulk.values[1][2][0], 100.0)
  assert_eq(pseudobulk.values[1][2][4], 100.0)
  assert_eq(pseudobulk.values[0][3][0], 30.0)
  assert_eq(pseudobulk.values[0][3][4], 100.0)
}

///|
test "muscat advanced mean aggregation divides by cluster sample cells" {
  let pseudobulk = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_mean_expression(),
  )
  assert_true((pseudobulk.values[0][0][0] - 100.9).abs() < 1.0e-12)
  assert_eq(pseudobulk.values[1][2][0], 10.0)
  assert_eq(pseudobulk.values[1][2][4], 10.0)
}

///|
test "muscat advanced median aggregation handles even cell counts" {
  let pseudobulk = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_median_expression(),
  )
  assert_eq(pseudobulk.values[0][0][0], 101.0)
  assert_eq(pseudobulk.values[1][2][0], 0.0)
  assert_eq(pseudobulk.values[1][2][4], 10.0)
}

///|
test "muscat advanced proportion detected stays on unit interval" {
  let pseudobulk = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_proportion_detected(),
  )
  assert_eq(pseudobulk.values[0][0][0], 1.0)
  assert_eq(pseudobulk.values[1][2][0], 0.1)
  assert_eq(pseudobulk.values[1][2][4], 1.0)
  for cluster in pseudobulk.values {
    for gene in cluster {
      for value in gene {
        assert_true(value >= 0.0 && value <= 1.0)
      }
    }
  }
}

///|
test "muscat advanced number detected counts positive cells" {
  let pseudobulk = muscat_adv_test_detection()
  assert_eq(pseudobulk.values[0][0][0], 10.0)
  assert_eq(pseudobulk.values[1][2][0], 1.0)
  assert_eq(pseudobulk.values[1][2][4], 10.0)
  assert_eq(pseudobulk.values[0][3][0], 1.0)
  assert_eq(pseudobulk.values[0][3][4], 10.0)
}

///|
test "muscat advanced aggregation fills missing cluster sample combinations" {
  let data = @src.MuscatAdvancedData::create(
    [[1.0, 2.0]],
    ["s1", "s2"],
    ["A", "B"],
    ["ctrl", "stim"],
    gene_names=["g"],
  )
  let pseudobulk = @src.muscat_aggregate_advanced(data)
  assert_eq(pseudobulk.values[0][0], [1.0, 0.0])
  assert_eq(pseudobulk.values[1][0], [0.0, 2.0])
  assert_eq(pseudobulk.cell_counts[0], [1, 0])
  assert_eq(pseudobulk.cell_counts[1], [0, 1])
}

///|
test "muscat advanced library sizes retain raw sum counts" {
  let sums = muscat_adv_test_sum()
  let means = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_mean_expression(),
  )
  assert_eq(sums.library_sizes[0][0], 1107.0)
  assert_eq(sums.library_sizes[1][0], 1177.0)
  assert_eq(means.library_sizes, sums.library_sizes)
}

///|
test "muscat advanced CPM scaling uses raw cluster sample library size" {
  let pseudobulk = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_sum_counts(),
    scale_cpm=true,
  )
  assert_true(pseudobulk.scaled_cpm)
  assert_true(
    (pseudobulk.values[0][0][0] - 1009.0 / 1107.0 * 1.0e6).abs() < 1.0e-8,
  )
}

///|
test "muscat advanced rejects CPM scaling of detection counts" {
  let failed = try {
    ignore(
      @src.muscat_aggregate_advanced(
        muscat_adv_test_data(),
        aggregation=@src.muscat_number_detected(),
        scale_cpm=true,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced group design uses reference coding" {
  let design = muscat_adv_test_design()
  assert_eq(design.matrix.length(), 8)
  assert_eq(design.matrix[0], [1.0, 0.0])
  assert_eq(design.matrix[3], [1.0, 0.0])
  assert_eq(design.matrix[4], [1.0, 1.0])
  assert_eq(design.coefficient_names, ["(Intercept)", "stim-ctrl"])
  assert_eq(design.sample_names, muscat_adv_test_sum().sample_names)
}

///|
test "muscat advanced group design supports alternate reference" {
  let design = @src.muscat_group_design_advanced(
    muscat_adv_test_sum(),
    reference="stim",
  )
  assert_eq(design.matrix[0], [1.0, 1.0])
  assert_eq(design.matrix[4], [1.0, 0.0])
  assert_eq(design.coefficient_names, ["(Intercept)", "ctrl-stim"])
}

///|
test "muscat advanced group design rejects absent reference" {
  let failed = try {
    ignore(
      @src.muscat_group_design_advanced(
        muscat_adv_test_sum(),
        reference="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced custom design preserves arbitrary covariates" {
  let design = @src.MuscatAdvancedDesign::create(
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
    ["s1", "s2", "s3", "s4"],
    ["intercept", "treatment"],
  )
  assert_eq(design.matrix[2], [1.0, 1.0])
  assert_eq(design.coefficient_names, ["intercept", "treatment"])
}

///|
test "muscat advanced validates design dimensions and names" {
  let mut failures = 0
  ignore(@src.MuscatAdvancedDesign::create([], [], [])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedDesign::create([[1.0], [1.0]], ["s1"], ["i"])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.MuscatAdvancedDesign::create([[1.0], [1.0]], ["s1", "s1"], ["i"])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced requires finite full-rank residual design" {
  let mut failures = 0
  ignore(
    @src.MuscatAdvancedDesign::create([[1.0, 0.0], [0.0, 1.0]], ["s1", "s2"], [
      "a", "b",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedDesign::create(
      [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
      ["s1", "s2", "s3"],
      ["a", "b"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.MuscatAdvancedDesign::create([[1.0], [0.0 / 0.0]], ["s1", "s2"], ["a"]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced contrast preserves coefficient vector" {
  let contrast = @src.muscat_contrast_advanced(
    muscat_adv_test_design(),
    [0.0, 1.0],
    "stimulus",
  )
  assert_eq(contrast.values, [0.0, 1.0])
  assert_eq(contrast.name, "stimulus")
}

///|
test "muscat advanced validates contrast dimensions values and name" {
  let design = muscat_adv_test_design()
  let mut failures = 0
  ignore(@src.muscat_contrast_advanced(design, [1.0], "short")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.muscat_contrast_advanced(design, [0.0, 0.0], "zero")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.muscat_contrast_advanced(design, [0.0, 0.0 / 0.0], "nan")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.muscat_contrast_advanced(design, [0.0, 1.0], " ")) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "muscat advanced default contrasts cover all non-intercept coefficients" {
  let design = @src.MuscatAdvancedDesign::create(
    [
      [1.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 0.0, 1.0],
      [1.0, 0.0, 1.0],
    ],
    ["s1", "s2", "s3", "s4", "s5", "s6"],
    ["intercept", "B-A", "C-A"],
  )
  let contrasts = @src.muscat_default_contrasts_advanced(design)
  assert_eq(contrasts.length(), 2)
  assert_eq(contrasts[0].name, "B-A")
  assert_eq(contrasts[0].values, [0.0, 1.0, 0.0])
  assert_eq(contrasts[1].values, [0.0, 0.0, 1.0])
}

///|
test "muscat advanced testing requires exact pseudobulk sample order" {
  let design = muscat_adv_test_design()
  let reversed = @src.MuscatAdvancedDesign::create(
    design.matrix,
    ["T4", "T3", "T2", "T1", "C4", "C3", "C2", "C1"],
    design.coefficient_names,
  )
  let failed = try {
    ignore(
      @src.muscat_pbds_advanced(
        muscat_adv_test_sum(),
        reversed,
        muscat_adv_test_contrasts(),
        config=muscat_adv_test_config(),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced testing validates contrast collection" {
  let design = muscat_adv_test_design()
  let contrast = muscat_adv_test_contrasts()[0]
  let mut failures = 0
  ignore(
    @src.muscat_pbds_advanced(
      muscat_adv_test_sum(),
      design,
      [],
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_pbds_advanced(
      muscat_adv_test_sum(),
      design,
      [contrast, contrast],
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "muscat advanced DS requires unscaled sum counts" {
  let mean = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    aggregation=@src.muscat_mean_expression(),
  )
  let scaled = @src.muscat_aggregate_advanced(
    muscat_adv_test_data(),
    scale_cpm=true,
  )
  let mut failures = 0
  ignore(
    @src.muscat_pbds_advanced(
      mean,
      muscat_adv_test_design(),
      muscat_adv_test_contrasts(),
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_pbds_advanced(
      scaled,
      muscat_adv_test_design(),
      muscat_adv_test_contrasts(),
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "muscat advanced DD requires number-detected pseudobulk" {
  let failed = try {
    ignore(
      @src.muscat_pbdd_advanced(
        muscat_adv_test_sum(),
        muscat_adv_test_design(),
        muscat_adv_test_contrasts(),
        config=muscat_adv_test_config(),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced reports when filtering removes every model" {
  let failed = try {
    ignore(
      @src.muscat_pbds_advanced(
        muscat_adv_test_sum(),
        muscat_adv_test_design(),
        muscat_adv_test_contrasts(),
        config=@src.MuscatAdvancedConfig::create(min_cells=11),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "muscat advanced DS returns cluster gene contrast results" {
  let result = muscat_adv_test_ds()
  assert_eq(result.mode, "DS")
  assert_eq(result.results.length(), 20)
  assert_eq(result.gene_names, muscat_adv_test_data().gene_names)
  assert_eq(result.cluster_names, ["A", "B"])
  assert_eq(result.contrast_names, ["stim-ctrl"])
  assert_eq(result.n_tested(), 18)
  assert_true(result.summary().contains("10 genes x 2 clusters x 1 contrasts"))
}

///|
test "muscat advanced DS detects cluster-specific abundance signal" {
  let result = muscat_adv_test_ds()
  let signal = match result.get("ds_a", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("DS signal result should exist")
  }
  let inactive = match result.get("ds_a", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("inactive DS result should exist")
  }
  assert_true(signal.tested)
  assert_true(signal.log_fc > 1.5)
  assert_true(signal.statistic > inactive.statistic)
  assert_true(signal.p_value < inactive.p_value)
  assert_true(signal.local_fdr <= 0.1)
}

///|
test "muscat advanced DS leaves equal-sum DD signal unchanged" {
  let result = muscat_adv_test_ds()
  let pure_dd = match result.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("pure DD result should exist")
  }
  assert_true(pure_dd.tested)
  assert_true(pure_dd.log_fc.abs() < 1.0e-5)
  assert_true(pure_dd.p_value > 0.9)
}

///|
test "muscat advanced DS statistics and dispersion are finite" {
  let result = muscat_adv_test_ds()
  let config = muscat_adv_test_config()
  for value in result.results {
    assert_true(muscat_adv_test_finite(value.log_fc))
    assert_true(muscat_adv_test_finite(value.log_cpm))
    assert_true(muscat_adv_test_finite(value.statistic))
    assert_true(value.p_value >= 0.0 && value.p_value <= 1.0)
    assert_true(value.local_fdr >= value.p_value - 1.0e-12)
    assert_true(value.global_fdr >= value.p_value - 1.0e-12)
    assert_true(value.dispersion >= config.minimum_dispersion)
    assert_true(value.dispersion <= config.maximum_dispersion)
    assert_eq(value.n_samples, 8)
  }
}

///|
test "muscat advanced filters all-zero genes without nonfinite values" {
  let result = muscat_adv_test_ds()
  let zero = match result.get("zero", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("zero-gene result should exist")
  }
  assert_false(zero.tested)
  assert_eq(zero.statistic, 0.0)
  assert_eq(zero.p_value, 1.0)
  assert_eq(zero.local_fdr, 1.0)
  assert_eq(zero.global_fdr, 1.0)
}

///|
test "muscat advanced fold-change threshold cannot increase Wald evidence" {
  let baseline = muscat_adv_test_ds()
  let thresholded = @src.muscat_pbds_advanced(
    muscat_adv_test_sum(),
    muscat_adv_test_design(),
    muscat_adv_test_contrasts(),
    config=@src.MuscatAdvancedConfig::create(
      min_cells=10,
      ridge=1.0e-6,
      fdr_threshold=0.1,
      lfc_threshold=1.0,
    ),
  )
  for index in 0..<baseline.results.length() {
    assert_true(
      thresholded.results[index].statistic <=
      baseline.results[index].statistic + 1.0e-10,
    )
    assert_true(
      thresholded.results[index].p_value + 1.0e-10 >=
      baseline.results[index].p_value,
    )
  }
}

///|
test "muscat advanced DD filters ubiquitously detected genes" {
  let result = muscat_adv_test_dd()
  let stable = match result.get("stable", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("stable DD result should exist")
  }
  let ds = match result.get("ds_a", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("ubiquitous DS gene result should exist")
  }
  assert_false(stable.tested)
  assert_false(ds.tested)
  assert_eq(stable.p_value, 1.0)
  assert_eq(ds.p_value, 1.0)
}

///|
test "muscat advanced DD detects equal-sum detection signal" {
  let result = muscat_adv_test_dd()
  let signal = match result.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("DD signal result should exist")
  }
  let inactive = match result.get("dd_b", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("inactive DD result should exist")
  }
  assert_true(signal.tested)
  assert_false(inactive.tested)
  assert_true(signal.log_fc > 1.5)
  assert_true(signal.p_value < 0.01)
  assert_true(signal.local_fdr < 0.1)
  assert_true(signal.global_fdr < 0.1)
}

///|
test "muscat advanced DD detects combined abundance and detection signal" {
  let result = muscat_adv_test_dd()
  let signal = match result.get("both_a", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("combined DD signal result should exist")
  }
  assert_true(signal.tested)
  assert_true(signal.log_fc > 1.5)
  assert_true(signal.statistic > 5.0)
  assert_true(signal.local_fdr < 0.1)
}

///|
test "muscat advanced DD statistics remain finite after CDR normalization" {
  let result = muscat_adv_test_dd()
  assert_eq(result.mode, "DD")
  assert_eq(result.results.length(), 20)
  assert_eq(result.n_tested(), 12)
  for value in result.results {
    assert_true(muscat_adv_test_finite(value.log_fc))
    assert_true(muscat_adv_test_finite(value.statistic))
    assert_true(muscat_adv_test_finite(value.dispersion))
    assert_true(value.p_value >= 0.0 && value.p_value <= 1.0)
  }
}

///|
test "muscat advanced local and global BH adjustments are bounded" {
  for result in [muscat_adv_test_ds(), muscat_adv_test_dd()] {
    for value in result.results {
      assert_true(value.local_fdr >= 0.0 && value.local_fdr <= 1.0)
      assert_true(value.global_fdr >= 0.0 && value.global_fdr <= 1.0)
      if value.tested {
        assert_true(value.local_fdr + 1.0e-12 >= value.p_value)
        assert_true(value.global_fdr + 1.0e-12 >= value.p_value)
      }
    }
  }
}

///|
test "muscat advanced stagewise output aligns DS and DD hypotheses" {
  let result = muscat_adv_test_stagewise()
  assert_eq(result.results.length(), 20)
  assert_eq(result.gene_names, muscat_adv_test_data().gene_names)
  assert_eq(result.cluster_names, ["A", "B"])
  assert_eq(result.contrast_names, ["stim-ctrl"])
  assert_eq(result.alpha, 0.1)
  for value in result.results {
    assert_true(value.screen_p_value >= 0.0 && value.screen_p_value <= 1.0)
    assert_true(value.screen_fdr >= 0.0 && value.screen_fdr <= 1.0)
    assert_true(value.ds_confirmation >= 0.0 && value.ds_confirmation <= 1.0)
    assert_true(value.dd_confirmation >= 0.0 && value.dd_confirmation <= 1.0)
  }
}

///|
test "muscat advanced stagewise classifies DS DD and combined signals" {
  let result = muscat_adv_test_stagewise()
  let ds = match result.get("ds_a", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("stagewise DS result should exist")
  }
  let dd = match result.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("stagewise DD result should exist")
  }
  let both = match result.get("both_a", "A", "stim-ctrl") {
    Some(value) => value
    None => abort("stagewise combined result should exist")
  }
  assert_eq(ds.classification, "DS")
  assert_eq(dd.classification, "DD")
  assert_eq(both.classification, "both")
  assert_true(result.n_confirmed() >= 3)
}

///|
test "muscat advanced stagewise screening uses harmonic mean p-value" {
  let ds = muscat_adv_test_ds()
  let dd = muscat_adv_test_dd()
  let result = @src.muscat_stagewise_ds_dd(ds, dd, alpha=0.1)
  let ds_gene = match ds.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("DS component should exist")
  }
  let dd_gene = match dd.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("DD component should exist")
  }
  let combined = match result.get("dd_b", "B", "stim-ctrl") {
    Some(value) => value
    None => abort("stagewise component should exist")
  }
  let expected = 2.0 / (1.0 / ds_gene.p_value + 1.0 / dd_gene.p_value)
  assert_true((combined.screen_p_value - expected).abs() < 1.0e-12)
}

///|
test "muscat advanced stagewise validates alpha and result modes" {
  let ds = muscat_adv_test_ds()
  let dd = muscat_adv_test_dd()
  let mut failures = 0
  ignore(@src.muscat_stagewise_ds_dd(ds, dd, alpha=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.muscat_stagewise_ds_dd(ds, dd, alpha=1.1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.muscat_stagewise_ds_dd(ds, ds)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "muscat advanced result lookup returns none for absent keys" {
  assert_true(muscat_adv_test_ds().get("missing", "A", "stim-ctrl") is None)
  assert_true(
    muscat_adv_test_stagewise().get("missing", "A", "stim-ctrl") is None,
  )
}

///|
test "muscat advanced result summaries report mode and significance" {
  let ds = muscat_adv_test_ds()
  let dd = muscat_adv_test_dd()
  assert_true(ds.summary().contains("muscat advanced DS"))
  assert_true(dd.summary().contains("muscat advanced DD"))
  assert_true(ds.n_significant() >= 2)
  assert_true(ds.n_significant(global=true) >= 2)
  assert_eq(dd.n_significant(), 2)
  assert_eq(dd.n_significant(global=true), 2)
}

///|
test "muscat advanced SCE writes DS DD and stagewise annotations" {
  let output = @src.muscat_advanced_sce(
    muscat_adv_test_sce(),
    "sample",
    "cluster",
    "group",
    reference="ctrl",
    config=muscat_adv_test_config(),
  ) catch {
    _ => abort("valid muscat SCE should analyze")
  }
  assert_eq(output.sum_pseudobulk.values.length(), 2)
  assert_eq(output.detection_pseudobulk.values.length(), 2)
  assert_eq(output.ds_result.mode, "DS")
  assert_eq(output.dd_result.mode, "DD")
  assert_eq(output.experiment.row_data["muscat.A.dsLogFC"].length(), 10)
  assert_eq(output.experiment.row_data["muscat.A.dsFdr"].length(), 10)
  assert_eq(output.experiment.row_data["muscat.B.ddLogFC"].length(), 10)
  assert_eq(output.experiment.row_data["muscat.B.stageClass"].length(), 10)
  assert_eq(output.experiment.metadata["muscat.version"], "muscat 1.27.4")
  assert_eq(output.experiment.metadata["muscat.reference"], "ctrl")
  assert_eq(output.experiment.metadata["muscat.contrast"], "stim-ctrl")
}

///|
test "muscat advanced SCE preserves input immutability" {
  let experiment = muscat_adv_test_sce()
  let output = @src.muscat_advanced_sce(
    experiment,
    "sample",
    "cluster",
    "group",
    config=muscat_adv_test_config(),
  )
  output.experiment.assays["counts"][0][0] = 999.0
  output.experiment.row_data["symbol"][0] = "changed"
  output.experiment.col_data["batch"][0] = "changed"
  output.experiment.reduced_dims["PCA"][0][0] = 999.0
  output.experiment.metadata["source"] = "changed"
  assert_eq(experiment.assays["counts"][0][0], 100.0)
  assert_eq(experiment.row_data["symbol"][0], "stable")
  assert_eq(experiment.col_data["batch"][0], "one")
  assert_eq(experiment.reduced_dims["PCA"][0][0], 1.0)
  assert_eq(experiment.metadata["source"], "test")
  assert_false(experiment.row_data.contains("muscat.A.dsLogFC"))
}

///|
test "muscat advanced SCE recursively copies alternative experiments" {
  let experiment = muscat_adv_test_sce()
  let output = @src.muscat_advanced_sce(
    experiment,
    "sample",
    "cluster",
    "group",
    config=muscat_adv_test_config(),
  )
  output.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 8.0
  assert_eq(
    experiment.alternative_experiments["spike"].assays["counts"][0][0],
    1.0,
  )
}

///|
test "muscat advanced SCE supports custom assay and output prefix" {
  let experiment = muscat_adv_test_sce()
  experiment.assays["raw"] = experiment.assays["counts"].map(fn(row) {
    row.copy()
  })
  let output = @src.muscat_advanced_sce(
    experiment,
    "sample",
    "cluster",
    "group",
    assay_name="raw",
    reference="ctrl",
    output_prefix="mx",
    config=muscat_adv_test_config(),
  )
  assert_true(output.experiment.row_data.contains("mx.A.dsLogFC"))
  assert_true(output.experiment.row_data.contains("mx.B.stageClass"))
  assert_true(output.experiment.col_data.contains("mx.sample_id"))
  assert_eq(output.experiment.metadata["mx.assay"], "raw")
  assert_eq(output.experiment.metadata["mx.clusters"], "A,B")
}

///|
test "muscat advanced SCE validates assay columns and prefix" {
  let experiment = muscat_adv_test_sce()
  let mut failures = 0
  ignore(
    @src.muscat_advanced_sce(
      experiment,
      "sample",
      "cluster",
      "group",
      assay_name="missing",
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_advanced_sce(
      experiment,
      "missing",
      "cluster",
      "group",
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_advanced_sce(
      experiment,
      "sample",
      "missing",
      "group",
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_advanced_sce(
      experiment,
      "sample",
      "cluster",
      "missing",
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.muscat_advanced_sce(
      experiment,
      "sample",
      "cluster",
      "group",
      output_prefix=" ",
      config=muscat_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 5)
}

///|
test "muscat advanced SCE validates reference group" {
  let failed = try {
    ignore(
      @src.muscat_advanced_sce(
        muscat_adv_test_sce(),
        "sample",
        "cluster",
        "group",
        reference="missing",
        config=muscat_adv_test_config(),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}
