///|
fn mast_adv_test_group_data() -> @src.MastAdvancedData {
  @src.mast_advanced_from_groups(
    [
      [0.0, 1.0, 0.0, 1.1, 0.0, 2.0, 0.0, 2.1],
      [0.0, 0.0, 1.0, 1.1, 2.0, 2.1, 2.2, 2.3],
      [1.0, 1.1, 0.9, 1.2, 1.0, 1.1, 0.9, 1.2],
    ],
    ["A", "A", "A", "A", "B", "B", "B", "B"],
    reference="A",
    gene_names=["continuous", "detection", "stable"],
    cell_names=["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"],
    include_cdr=false,
  ) catch {
    _ => abort("valid grouped MAST data should build")
  }
}

///|
fn mast_adv_test_example_result() -> @src.MastAdvancedResult {
  let (data, tested) = @src.mast_advanced_example() catch {
    _ => abort("MAST advanced example should build")
  }
  @src.mast_advanced_zlm(data, tested) catch {
    _ => abort("MAST advanced example should fit")
  }
}

///|
fn mast_adv_test_sce() -> @src.SingleCellExperiment {
  let (data, _) = @src.mast_advanced_example() catch {
    _ => abort("MAST advanced example should build")
  }
  let experiment = @src.SingleCellExperiment::new(
    data.expression,
    data.gene_names,
    data.cell_names,
  )
  experiment.assays["logcounts"] = data.expression.map(fn(row) { row.copy() })
  let groups : Array[String] = []
  for cell in 0..<data.n_cells {
    groups.push(if cell < 24 { "control" } else { "stimulated" })
  }
  experiment.col_data["condition"] = groups
  experiment.row_data["symbol"] = data.gene_names.copy()
  experiment.reduced_dims["PCA"] = [[0.0, 1.0], [1.0, 0.0]]
  experiment.metadata["source"] = "test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[1.0, 2.0]],
    ["spike1"],
    ["s1", "s2"],
  )
  experiment
}

///|
fn mast_adv_test_na_result() -> @src.MastAdvancedResult {
  let data = @src.mast_advanced_from_groups(
    [
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      [1.0, 1.1, 0.9, 1.2, 1.0, 1.1, 2.0, 2.1, 1.9, 2.2, 2.0, 2.1],
      [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 1.1, 0.0, 1.2, 2.0, 0.0, 2.1, 0.0, 2.2, 0.0],
    ],
    ["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"],
    reference="A",
    gene_names=["zero", "all_positive", "sparse", "mixed"],
    include_cdr=false,
  ) catch {
    _ => abort("valid MAST NA test data should build")
  }
  let config = @src.MastAdvancedConfig::create(empirical_bayes=false) catch {
    _ => abort("valid MAST configuration should build")
  }
  @src.mast_advanced_zlm(data, [1], config~) catch {
    _ => abort("MAST NA test data should fit")
  }
}

///|
test "MAST advanced default configuration matches zlm workflow" {
  let config = @src.MastAdvancedConfig::default()
  assert_eq(config.max_iterations, 100)
  assert_eq(config.cauchy_prior_scale, 2.5)
  assert_true(config.empirical_bayes)
  assert_eq(config.ebayes_use_full_model, false)
  assert_eq(config.minimum_detected, 2)
  assert_eq(config.fdr_threshold, 0.05)
}

///|
test "MAST advanced configuration preserves explicit controls" {
  let config = @src.MastAdvancedConfig::create(
    max_iterations=25,
    tolerance=1.0e-6,
    cauchy_prior_scale=1.5,
    ridge=1.0e-5,
    empirical_bayes=false,
    ebayes_use_full_model=true,
    maximum_prior_df=50.0,
    minimum_detected=3,
    minimum_probability=1.0e-7,
    fdr_threshold=0.1,
  ) catch {
    _ => abort("valid explicit MAST configuration should build")
  }
  assert_eq(config.max_iterations, 25)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.cauchy_prior_scale, 1.5)
  assert_eq(config.ridge, 1.0e-5)
  assert_eq(config.empirical_bayes, false)
  assert_true(config.ebayes_use_full_model)
  assert_eq(config.maximum_prior_df, 50.0)
  assert_eq(config.minimum_detected, 3)
  assert_eq(config.minimum_probability, 1.0e-7)
  assert_eq(config.fdr_threshold, 0.1)
}

///|
test "MAST advanced rejects invalid iteration and tolerance controls" {
  let iterations = try {
    ignore(@src.MastAdvancedConfig::create(max_iterations=0))
    false
  } catch {
    _ => true
  }
  let tolerance = try {
    ignore(@src.MastAdvancedConfig::create(tolerance=0.0))
    false
  } catch {
    _ => true
  }
  assert_true(iterations)
  assert_true(tolerance)
}

///|
test "MAST advanced rejects invalid prior and ridge controls" {
  let scale = try {
    ignore(@src.MastAdvancedConfig::create(cauchy_prior_scale=0.0))
    false
  } catch {
    _ => true
  }
  let ridge = try {
    ignore(@src.MastAdvancedConfig::create(ridge=-1.0))
    false
  } catch {
    _ => true
  }
  let prior_df = try {
    ignore(@src.MastAdvancedConfig::create(maximum_prior_df=0.0))
    false
  } catch {
    _ => true
  }
  assert_true(scale)
  assert_true(ridge)
  assert_true(prior_df)
}

///|
test "MAST advanced rejects invalid probability controls" {
  let minimum = try {
    ignore(@src.MastAdvancedConfig::create(minimum_probability=0.5))
    false
  } catch {
    _ => true
  }
  let fdr = try {
    ignore(@src.MastAdvancedConfig::create(fdr_threshold=1.1))
    false
  } catch {
    _ => true
  }
  let detected = try {
    ignore(@src.MastAdvancedConfig::create(minimum_detected=0))
    false
  } catch {
    _ => true
  }
  assert_true(minimum)
  assert_true(fdr)
  assert_true(detected)
}

///|
test "MAST advanced grouped constructor uses feature by cell orientation" {
  let data = mast_adv_test_group_data()
  assert_eq(data.n_genes, 3)
  assert_eq(data.n_cells, 8)
  assert_eq(data.expression[0].length(), 8)
  assert_eq(data.gene_names, ["continuous", "detection", "stable"])
  assert_eq(data.cell_names[7], "c8")
}

///|
test "MAST advanced grouped constructor uses treatment coding" {
  let data = mast_adv_test_group_data()
  assert_eq(data.coefficient_names, ["(Intercept)", "group:B"])
  assert_eq(data.design[0], [1.0, 0.0])
  assert_eq(data.design[4], [1.0, 1.0])
}

///|
test "MAST advanced grouped constructor honors reference level" {
  let data = @src.mast_advanced_from_groups(
    [[1.0, 2.0, 3.0, 4.0]],
    ["A", "A", "B", "B"],
    reference="B",
    include_cdr=false,
  ) catch {
    _ => abort("valid reference level should build")
  }
  assert_eq(data.coefficient_names, ["(Intercept)", "group:A"])
  assert_eq(data.design[0], [1.0, 1.0])
  assert_eq(data.design[2], [1.0, 0.0])
}

///|
test "MAST advanced constructor appends cell detection rate" {
  let data = @src.MastAdvancedData::create(
    [[1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]],
    [[1.0], [1.0], [1.0], [1.0]],
    ["(Intercept)"],
  ) catch {
    _ => abort("variable CDR should produce a full-rank design")
  }
  assert_eq(data.coefficient_names, ["(Intercept)", "cngeneson"])
  assert_eq(data.cdr, [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0])
  assert_eq(data.design[3], [1.0, 1.0])
}

///|
test "MAST advanced constructor supports arbitrary design matrices" {
  let data = @src.MastAdvancedData::create(
    [[1.0, 1.2, 1.4, 2.0, 2.2, 2.4]],
    [
      [1.0, 0.0, -1.0],
      [1.0, 0.0, 0.0],
      [1.0, 0.0, 1.0],
      [1.0, 1.0, -1.0],
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 1.0],
    ],
    ["(Intercept)", "condition", "batch_score"],
    include_cdr=false,
  ) catch {
    _ => abort("full-rank arbitrary design should build")
  }
  assert_eq(data.n_coefficients, 3)
  assert_eq(data.design[5][2], 1.0)
}

///|
test "MAST advanced constructor defensively copies inputs" {
  let expression = [[1.0, 0.0, 2.0, 0.0]]
  let design = [[1.0], [1.0], [1.0], [1.0]]
  let genes = ["g1"]
  let cells = ["c1", "c2", "c3", "c4"]
  let data = @src.MastAdvancedData::create(
    expression,
    design,
    ["(Intercept)"],
    gene_names=genes,
    cell_names=cells,
    include_cdr=false,
  ) catch {
    _ => abort("valid MAST data should build")
  }
  expression[0][0] = 99.0
  design[0][0] = 99.0
  genes[0] = "changed"
  cells[0] = "changed"
  assert_eq(data.expression[0][0], 1.0)
  assert_eq(data.design[0][0], 1.0)
  assert_eq(data.gene_names[0], "g1")
  assert_eq(data.cell_names[0], "c1")
}

///|
test "MAST advanced constructor rejects empty and ragged expression" {
  let empty = try {
    ignore(@src.MastAdvancedData::create([], [], [], include_cdr=false))
    false
  } catch {
    _ => true
  }
  let ragged = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0], [3.0]],
        [[1.0], [1.0]],
        ["intercept"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(empty)
  assert_true(ragged)
}

///|
test "MAST advanced constructor rejects invalid expression values" {
  let negative = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, -1.0]],
        [[1.0], [1.0]],
        ["intercept"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let non_finite = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, @double.not_a_number]],
        [[1.0], [1.0]],
        ["intercept"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(negative)
  assert_true(non_finite)
}

///|
test "MAST advanced constructor validates design dimensions" {
  let rows = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0]],
        [[1.0]],
        ["intercept"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let columns = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0]],
        [[1.0], [1.0, 0.0]],
        ["intercept"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let names = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0]],
        [[1.0], [1.0]],
        [],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(rows)
  assert_true(columns)
  assert_true(names)
}

///|
test "MAST advanced constructor rejects rank deficient designs" {
  let raised = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0, 3.0]],
        [[1.0, 2.0], [1.0, 2.0], [1.0, 2.0]],
        ["intercept", "constant"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "MAST advanced constructor validates identifiers" {
  let length = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0]],
        [[1.0], [1.0]],
        ["intercept"],
        gene_names=["g1", "g2"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let duplicate = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0], [2.0, 3.0]],
        [[1.0], [1.0]],
        ["intercept"],
        gene_names=["g1", "g1"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let empty = try {
    ignore(
      @src.MastAdvancedData::create(
        [[1.0, 2.0]],
        [[1.0], [1.0]],
        [" "],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(length)
  assert_true(duplicate)
  assert_true(empty)
}

///|
test "MAST advanced grouped constructor validates groups and reference" {
  let length = try {
    ignore(
      @src.mast_advanced_from_groups([[1.0, 2.0]], ["A"], include_cdr=false),
    )
    false
  } catch {
    _ => true
  }
  let one_level = try {
    ignore(
      @src.mast_advanced_from_groups(
        [[1.0, 2.0]],
        ["A", "A"],
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  let reference = try {
    ignore(
      @src.mast_advanced_from_groups(
        [[1.0, 2.0]],
        ["A", "B"],
        reference="C",
        include_cdr=false,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(length)
  assert_true(one_level)
  assert_true(reference)
}

///|
test "MAST advanced example exposes intended dimensions and contrast" {
  let (data, tested) = @src.mast_advanced_example() catch {
    _ => abort("MAST advanced example should build")
  }
  assert_eq(data.n_genes, 8)
  assert_eq(data.n_cells, 48)
  assert_eq(data.coefficient_names[1], "group:stimulated")
  assert_eq(tested, [1])
}

///|
test "MAST advanced zlm returns aligned result matrices" {
  let result = mast_adv_test_example_result()
  assert_eq(result.gene_names.length(), 8)
  assert_eq(result.discrete_coefficients.length(), 8)
  assert_eq(result.continuous_coefficients.length(), 8)
  assert_eq(
    result.discrete_coefficients[0].length(),
    result.coefficient_names.length(),
  )
  assert_eq(result.tested_names, ["group:stimulated"])
}

///|
test "MAST advanced recovers discrete and continuous effect directions" {
  let result = mast_adv_test_example_result()
  assert_true(result.discrete_coefficients[0][1] > 0.0)
  assert_true(result.continuous_coefficients[0][1] > 0.0)
  assert_true(result.continuous_coefficients[1][1] > 0.0)
}

///|
test "MAST advanced computes valid component and hurdle probabilities" {
  let result = mast_adv_test_example_result()
  for gene in 0..<result.gene_names.length() {
    match result.hurdle_p_values[gene] {
      Some(value) => assert_true(value >= 0.0 && value <= 1.0)
      None => ()
    }
    match result.hurdle_fdr[gene] {
      Some(value) => assert_true(value >= 0.0 && value <= 1.0)
      None => ()
    }
  }
  assert_true(result.n_tested() > 0)
}

///|
test "MAST advanced hurdle statistic sums testable components" {
  let result = mast_adv_test_example_result()
  let discrete = result.discrete_statistics[0].unwrap_or(-1.0)
  let continuous = result.continuous_statistics[0].unwrap_or(-1.0)
  let hurdle = result.hurdle_statistics[0].unwrap_or(-1.0)
  assert_true(discrete >= 0.0)
  assert_true(continuous >= 0.0)
  assert_true((hurdle - discrete - continuous).abs() < 1.0e-8)
}

///|
test "MAST advanced empirical Bayes estimates a finite variance prior" {
  let result = mast_adv_test_example_result()
  assert_true(result.empirical_bayes)
  assert_true(result.prior_variance > 0.0)
  assert_true(result.prior_variance < 1.0e300)
  assert_true(result.prior_df > 0.0)
  assert_true(result.prior_df <= 1000.0)
}

///|
test "MAST advanced empirical Bayes supports full-model residuals" {
  let (data, tested) = @src.mast_advanced_example() catch {
    _ => abort("MAST advanced example should build")
  }
  let config = @src.MastAdvancedConfig::create(
    ebayes_use_full_model=true,
    maximum_prior_df=80.0,
  ) catch {
    _ => abort("valid H1 eBayes configuration should build")
  }
  let result = @src.mast_advanced_zlm(data, tested, config~) catch {
    _ => abort("full-model eBayes workflow should fit")
  }
  assert_true(result.empirical_bayes)
  assert_true(result.prior_variance > 0.0)
  assert_true(result.prior_df > 0.0 && result.prior_df <= 80.0)
}

///|
test "MAST advanced can disable empirical Bayes moderation" {
  let (data, tested) = @src.mast_advanced_example() catch {
    _ => abort("MAST advanced example should build")
  }
  let config = @src.MastAdvancedConfig::create(empirical_bayes=false) catch {
    _ => abort("valid unmoderated configuration should build")
  }
  let result = @src.mast_advanced_zlm(data, tested, config~) catch {
    _ => abort("unmoderated MAST workflow should fit")
  }
  assert_eq(result.empirical_bayes, false)
  assert_eq(result.prior_variance, 0.0)
  assert_eq(result.prior_df, 0.0)
  for gene in 0..<result.gene_names.length() {
    match (result.continuous_variance[gene], result.moderated_variance[gene]) {
      (Some(raw), Some(moderated)) =>
        assert_true((raw - moderated).abs() < 1.0e-10)
      _ => ()
    }
  }
}

///|
test "MAST advanced moderated variances are finite and positive" {
  let result = mast_adv_test_example_result()
  for value in result.moderated_variance {
    match value {
      Some(variance) => {
        assert_true(variance > 0.0)
        assert_true(variance < 1.0e300)
      }
      None => ()
    }
  }
}

///|
test "MAST advanced preserves NA semantics for untestable genes" {
  let result = mast_adv_test_na_result()
  assert_true(result.discrete_p_values[0] is None)
  assert_true(result.continuous_p_values[0] is None)
  assert_true(result.hurdle_p_values[0] is None)
  assert_true(result.discrete_p_values[1] is None)
  assert_true(result.continuous_p_values[1] is Some(_))
  assert_true(result.continuous_p_values[2] is None)
}

///|
test "MAST advanced combines whichever hurdle components are testable" {
  let result = mast_adv_test_na_result()
  assert_true(result.hurdle_p_values[1] is Some(_))
  assert_true(result.hurdle_p_values[2] is Some(_))
  assert_true(result.hurdle_p_values[3] is Some(_))
}

///|
test "MAST advanced BH correction preserves missing values" {
  let result = mast_adv_test_na_result()
  assert_true(result.hurdle_fdr[0] is None)
  for gene in 1..<result.gene_names.length() {
    match (result.hurdle_p_values[gene], result.hurdle_fdr[gene]) {
      (Some(p), Some(fdr)) => assert_true(fdr + 1.0e-12 >= p)
      _ => ()
    }
  }
}

///|
test "MAST advanced marginal effects and delta variances are finite" {
  let result = mast_adv_test_example_result()
  assert_true(result.marginal_log_fc[0].unwrap_or(0.0) > 0.0)
  assert_true(result.marginal_log_fc[1] is Some(_))
  assert_true(result.marginal_log_fc[1].unwrap_or(1.0e300).abs() < 1.0e300)
  assert_true(result.marginal_log_fc_variance[0].unwrap_or(-1.0) >= 0.0)
}

///|
test "MAST advanced fitting is deterministic" {
  let first = mast_adv_test_example_result()
  let second = mast_adv_test_example_result()
  assert_eq(first.discrete_coefficients, second.discrete_coefficients)
  assert_eq(first.continuous_coefficients, second.continuous_coefficients)
  assert_eq(first.hurdle_p_values, second.hurdle_p_values)
  assert_eq(first.hurdle_fdr, second.hurdle_fdr)
}

///|
test "MAST advanced supports joint deletion of design columns" {
  let data = @src.MastAdvancedData::create(
    [[1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3, 3.0, 3.1, 3.2, 3.3]],
    [
      [1.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 0.0, 1.0],
      [1.0, 0.0, 1.0],
      [1.0, 0.0, 1.0],
      [1.0, 0.0, 1.0],
    ],
    ["(Intercept)", "group:B", "group:C"],
    include_cdr=false,
  ) catch {
    _ => abort("valid multi-level design should build")
  }
  let config = @src.MastAdvancedConfig::create(empirical_bayes=false) catch {
    _ => abort("valid MAST configuration should build")
  }
  let result = @src.mast_advanced_zlm(data, [1, 2], config~) catch {
    _ => abort("joint reduced-model test should fit")
  }
  assert_eq(result.tested_names, ["group:B", "group:C"])
  assert_true(result.continuous_p_values[0] is Some(_))
  assert_true(result.hurdle_statistics[0].unwrap_or(-1.0) >= 0.0)
}

///|
test "MAST advanced validates tested coefficient indices" {
  let data = mast_adv_test_group_data()
  let empty = try {
    ignore(@src.mast_advanced_zlm(data, []))
    false
  } catch {
    _ => true
  }
  let intercept = try {
    ignore(@src.mast_advanced_zlm(data, [0]))
    false
  } catch {
    _ => true
  }
  let range = try {
    ignore(@src.mast_advanced_zlm(data, [2]))
    false
  } catch {
    _ => true
  }
  let duplicate = try {
    ignore(@src.mast_advanced_zlm(data, [1, 1]))
    false
  } catch {
    _ => true
  }
  assert_true(empty)
  assert_true(intercept)
  assert_true(range)
  assert_true(duplicate)
}

///|
test "MAST advanced result query helpers are consistent" {
  let result = mast_adv_test_example_result()
  assert_eq(result.index_of("both_components"), Some(0))
  assert_true(result.index_of("missing") is None)
  assert_true(result.n_tested() <= result.gene_names.length())
  assert_eq(result.significant_genes().length(), result.n_significant())
  assert_true(result.top_genes(3).length() <= 3)
  assert_eq(result.top_genes(0), [])
}

///|
test "MAST advanced result summary reports contrast and tested count" {
  let result = mast_adv_test_example_result()
  let summary = result.summary()
  assert_true(summary.contains("MAST advanced zlm"))
  assert_true(summary.contains("group:stimulated"))
  assert_true(summary.contains(result.n_tested().to_string()))
}

///|
test "MAST advanced SingleCellExperiment writes row and cell results" {
  let output = @src.mast_advanced_zlm_sce(
    mast_adv_test_sce(),
    "condition",
    reference="control",
    output_prefix="mast39",
  ) catch {
    _ => abort("valid MAST SCE workflow should fit")
  }
  assert_eq(output.experiment.row_data["mast39.detectionRate"].length(), 8)
  assert_eq(output.experiment.row_data["mast39.nDetected"].length(), 8)
  assert_eq(output.experiment.row_data["mast39.hurdleP"].length(), 8)
  assert_eq(output.experiment.row_data["mast39.hurdleFdr"].length(), 8)
  assert_eq(output.experiment.row_data["mast39.class"].length(), 8)
  assert_eq(output.experiment.col_data["mast39.cdr"].length(), 48)
}

///|
test "MAST advanced SingleCellExperiment records provenance metadata" {
  let output = @src.mast_advanced_zlm_sce(
    mast_adv_test_sce(),
    "condition",
    reference="control",
    output_prefix="mast39",
  ) catch {
    _ => abort("valid MAST SCE workflow should fit")
  }
  assert_eq(output.experiment.metadata["mast39.assay"], "logcounts")
  assert_eq(output.experiment.metadata["mast39.group"], "condition")
  assert_eq(output.experiment.metadata["mast39.reference"], "control")
  assert_eq(output.experiment.metadata["mast39.contrast"], "group:stimulated")
  assert_eq(output.experiment.metadata["source"], "test")
}

///|
test "MAST advanced SingleCellExperiment does not mutate input" {
  let input = mast_adv_test_sce()
  let output = @src.mast_advanced_zlm_sce(
    input,
    "condition",
    reference="control",
    output_prefix="mast39",
  ) catch {
    _ => abort("valid MAST SCE workflow should fit")
  }
  assert_eq(input.row_data.contains("mast39.hurdleP"), false)
  assert_eq(input.col_data.contains("mast39.cdr"), false)
  output.experiment.assays["logcounts"][0][0] = 999.0
  output.experiment.row_data["symbol"][0] = "changed"
  output.experiment.reduced_dims["PCA"][0][0] = 999.0
  output.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 999.0
  assert_true(input.assays["logcounts"][0][0] != 999.0)
  assert_eq(input.row_data["symbol"][0], "both_components")
  assert_eq(input.reduced_dims["PCA"][0][0], 0.0)
  assert_eq(input.alternative_experiments["spike"].assays["counts"][0][0], 1.0)
}

///|
test "MAST advanced SingleCellExperiment supports custom assay and prefix" {
  let input = mast_adv_test_sce()
  let output = @src.mast_advanced_zlm_sce(
    input,
    "condition",
    assay_name="counts",
    reference="control",
    output_prefix="custom",
    include_cdr=false,
  ) catch {
    _ => abort("custom MAST SCE workflow should fit")
  }
  assert_true(output.experiment.row_data.contains("custom.logFC"))
  assert_eq(output.experiment.metadata["custom.assay"], "counts")
  assert_eq(output.result.coefficient_names, ["(Intercept)", "group:stimulated"])
}

///|
test "MAST advanced SingleCellExperiment writes NA for untestable values" {
  let result = mast_adv_test_na_result()
  let experiment = @src.SingleCellExperiment::new(
    [
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      [1.0, 1.1, 0.9, 1.2, 1.0, 1.1, 2.0, 2.1, 1.9, 2.2, 2.0, 2.1],
      [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 1.1, 0.0, 1.2, 2.0, 0.0, 2.1, 0.0, 2.2, 0.0],
    ],
    result.gene_names,
    ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12"],
  )
  experiment.assays["logcounts"] = experiment.assays["counts"]
  experiment.col_data["condition"] = [
    "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
  ]
  let config = @src.MastAdvancedConfig::create(empirical_bayes=false) catch {
    _ => abort("valid MAST configuration should build")
  }
  let output = @src.mast_advanced_zlm_sce(
    experiment,
    "condition",
    reference="A",
    include_cdr=false,
    config~,
  ) catch {
    _ => abort("MAST NA SCE workflow should fit")
  }
  assert_eq(output.experiment.row_data["mast.hurdleP"][0], "NA")
  assert_eq(output.experiment.row_data["mast.class"][0], "not_tested")
}

///|
test "MAST advanced SingleCellExperiment rejects missing inputs" {
  let assay = try {
    ignore(
      @src.mast_advanced_zlm_sce(
        mast_adv_test_sce(),
        "condition",
        assay_name="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  let group = try {
    ignore(@src.mast_advanced_zlm_sce(mast_adv_test_sce(), "missing"))
    false
  } catch {
    _ => true
  }
  assert_true(assay)
  assert_true(group)
}

///|
test "MAST advanced SingleCellExperiment validates prefix and reference" {
  let prefix = try {
    ignore(
      @src.mast_advanced_zlm_sce(
        mast_adv_test_sce(),
        "condition",
        output_prefix=" ",
      ),
    )
    false
  } catch {
    _ => true
  }
  let reference = try {
    ignore(
      @src.mast_advanced_zlm_sce(
        mast_adv_test_sce(),
        "condition",
        reference="absent",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(prefix)
  assert_true(reference)
}
