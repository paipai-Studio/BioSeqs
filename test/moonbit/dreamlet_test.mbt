///|
/// Tests for the Bioconductor dreamlet-inspired pseudobulk mixed-model
/// workflow.

///|
fn dreamlet_test_close(
  actual : Double,
  expected : Double,
  tolerance? : Double = 1.0e-8,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "dreamlet value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn dreamlet_test_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[String],
  Array[String],
  Map[String, Array[String]],
) {
  let genes = [
    "response_up", "stable", "cell_type", "rare", "donor", "response_down",
  ]
  let counts : Array[Array[Double]] = []
  for _ in genes {
    counts.push([])
  }
  let cells : Array[String] = []
  let samples : Array[String] = []
  let clusters : Array[String] = []
  let treatments : Array[String] = []
  let donors : Array[String] = []
  let ages : Array[String] = []
  for donor in 0..<4 {
    for treatment in 0..<2 {
      let sample = "D" +
        (donor + 1).to_string() +
        (if treatment == 0 { "_C" } else { "_T" })
      for cluster in 0..<2 {
        let cluster_name = if cluster == 0 { "T_cell" } else { "Monocyte" }
        for replicate in 0..<3 {
          cells.push(
            sample + "_" + cluster_name + "_" + (replicate + 1).to_string(),
          )
          samples.push(sample)
          clusters.push(cluster_name)
          treatments.push(if treatment == 0 { "Control" } else { "Treated" })
          donors.push("D" + (donor + 1).to_string())
          ages.push((30 + donor * 5).to_string())
          counts[0].push(
            (10 + treatment * 10 + donor + (if replicate == 2 { 1 } else { 0 })).to_double(),
          )
          counts[1].push((8 + replicate % 2).to_double())
          counts[2].push(
            (if cluster == 0 { 18 + treatment } else { 4 + treatment }).to_double(),
          )
          counts[3].push(
            if donor == 0 && treatment == 0 && cluster == 0 && replicate == 0 {
              1.0
            } else {
              0.0
            },
          )
          counts[4].push((6 + donor * 2 + replicate).to_double())
          counts[5].push(
            (18 - treatment * 8 + cluster + (if replicate == 1 { 1 } else { 0 })).to_double(),
          )
        }
      }
    }
  }
  let metadata : Map[String, Array[String]] = Map([
    ("Treatment", treatments),
    ("Donor", donors),
    ("Age", ages),
  ])
  (counts, genes, cells, samples, clusters, metadata)
}

///|
fn dreamlet_test_pseudobulk() -> @src.DreamletPseudobulk {
  let (counts, genes, cells, samples, clusters, metadata) = dreamlet_test_data()
  @src.dreamlet_aggregate_to_pseudobulk(
    counts,
    genes,
    cells,
    samples,
    clusters,
    cell_metadata=metadata,
  ) catch {
    _ => abort("valid dreamlet pseudobulk aggregation should succeed")
  }
}

///|
fn dreamlet_test_config() -> @src.DreamletProcessConfig {
  @src.DreamletProcessConfig::create(
    min_cells=2,
    min_count=1.0,
    min_samples=6,
    min_prop=0.5,
    min_total_count=10.0,
    span=0.7,
    max_iterations=30,
    tolerance=0.02,
  ) catch {
    _ => abort("valid dreamlet test configuration should build")
  }
}

///|
fn dreamlet_test_fixed_model() -> @src.DreamletModelSpec {
  @src.dreamlet_model([
    @src.dreamlet_categorical_effect("Treatment", reference="Control"),
  ]) catch {
    _ => abort("valid dreamlet fixed model should build")
  }
}

///|
fn dreamlet_test_mixed_model() -> @src.DreamletModelSpec {
  @src.dreamlet_model([
    @src.dreamlet_categorical_effect("Treatment", reference="Control"),
    @src.dreamlet_random_effect("Donor"),
  ]) catch {
    _ => abort("valid dreamlet mixed model should build")
  }
}

///|
fn dreamlet_test_processed(mixed? : Bool = false) -> @src.DreamletProcessedData {
  let model = if mixed {
    dreamlet_test_mixed_model()
  } else {
    dreamlet_test_fixed_model()
  }
  @src.dreamlet_process_assays(
    dreamlet_test_pseudobulk(),
    model,
    config=dreamlet_test_config(),
  ) catch {
    _ => abort("valid dreamlet processing should succeed")
  }
}

///|
fn dreamlet_test_result() -> @src.DreamletResult {
  @src.dreamlet(
    dreamlet_test_processed(),
    "Treatment:Treated",
    contrast_name="Treated-Control",
    max_iterations=30,
    tolerance=0.02,
  ) catch {
    _ => abort("valid dreamlet differential expression should succeed")
  }
}

///|
test "dreamlet: effect factories preserve kinds and references" {
  let numeric = @src.dreamlet_numeric_effect("Age")
  let categorical = @src.dreamlet_categorical_effect(
    "Treatment",
    reference="Control",
  )
  let random = @src.dreamlet_random_effect("Donor")
  assert_true(numeric.kind is Numeric)
  assert_true(categorical.kind is Categorical)
  assert_true(random.kind is Random)
  assert_eq(categorical.reference, "Control")
}

///|
test "dreamlet: typed model preserves effect order" {
  let model = dreamlet_test_mixed_model()
  assert_eq(model.effects.length(), 2)
  assert_eq(model.effects[0].name, "Treatment")
  assert_eq(model.effects[1].name, "Donor")
  assert_true(model.intercept)
}

///|
test "dreamlet: model rejects duplicate effect names" {
  let raised = try {
    ignore(
      @src.dreamlet_model([
        @src.dreamlet_numeric_effect("Age"),
        @src.dreamlet_random_effect("Age"),
      ]),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: model rejects empty effect names" {
  let raised = try {
    ignore(@src.dreamlet_model([@src.dreamlet_numeric_effect("")]))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: model rejects a no-intercept design without fixed effects" {
  let raised = try {
    ignore(
      @src.dreamlet_model(
        [@src.dreamlet_random_effect("Donor")],
        intercept=false,
      ),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: default process configuration matches upstream defaults" {
  let config = @src.DreamletProcessConfig::default()
  assert_eq(config.min_cells, 5)
  assert_eq(config.min_count, 5.0)
  assert_eq(config.min_samples, 4)
  assert_eq(config.min_prop, 0.4)
  assert_eq(config.prior_count, 0.5)
  assert_eq(config.logratio_trim, 0.3)
  assert_eq(config.sum_trim, 0.05)
}

///|
test "dreamlet: custom process configuration preserves controls" {
  let config = dreamlet_test_config()
  assert_eq(config.min_cells, 2)
  assert_eq(config.min_samples, 6)
  assert_eq(config.span, 0.7)
  assert_true(config.use_poisson_initial_weights)
  assert_false(config.rescale_initial_weights)
}

///|
test "dreamlet: process configuration rejects invalid sample filters" {
  let cells = try {
    ignore(@src.DreamletProcessConfig::create(min_cells=0))
    false
  } catch {
    DreamletError(_) => true
  }
  let samples = try {
    ignore(@src.DreamletProcessConfig::create(min_samples=2))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(cells)
  assert_true(samples)
}

///|
test "dreamlet: process configuration rejects invalid proportions" {
  let low = try {
    ignore(@src.DreamletProcessConfig::create(min_prop=0.0))
    false
  } catch {
    DreamletError(_) => true
  }
  let high = try {
    ignore(@src.DreamletProcessConfig::create(min_prop=1.1))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "dreamlet: process configuration rejects invalid TMM trims" {
  let ratio = try {
    ignore(@src.DreamletProcessConfig::create(logratio_trim=0.5))
    false
  } catch {
    DreamletError(_) => true
  }
  let abundance = try {
    ignore(@src.DreamletProcessConfig::create(sum_trim=-0.1))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(ratio)
  assert_true(abundance)
}

///|
test "dreamlet: pseudobulk uses gene by sample orientation" {
  let pseudobulk = dreamlet_test_pseudobulk()
  assert_eq(pseudobulk.gene_names.length(), 6)
  assert_eq(pseudobulk.sample_names.length(), 8)
  assert_eq(pseudobulk.assays.length(), 2)
  assert_eq(pseudobulk.assays[0].counts.length(), 6)
  assert_eq(pseudobulk.assays[0].counts[0].length(), 8)
}

///|
test "dreamlet: pseudobulk preserves first-observed cell type order" {
  let pseudobulk = dreamlet_test_pseudobulk()
  assert_eq(pseudobulk.assays[0].cluster_id, "T_cell")
  assert_eq(pseudobulk.assays[1].cluster_id, "Monocyte")
}

///|
test "dreamlet: pseudobulk preserves first-observed sample order" {
  let pseudobulk = dreamlet_test_pseudobulk()
  assert_eq(pseudobulk.sample_names[0], "D1_C")
  assert_eq(pseudobulk.sample_names[1], "D1_T")
  assert_eq(pseudobulk.sample_names[7], "D4_T")
}

///|
test "dreamlet: pseudobulk sums raw counts within sample and cell type" {
  let pseudobulk = dreamlet_test_pseudobulk()
  let assay = match pseudobulk.assay("T_cell") {
    Some(value) => value
    None => abort("T_cell assay should exist")
  }
  dreamlet_test_close(assay.counts[0][0], 31.0)
  dreamlet_test_close(assay.counts[0][1], 61.0)
}

///|
test "dreamlet: pseudobulk records observed cell counts" {
  let pseudobulk = dreamlet_test_pseudobulk()
  for assay in pseudobulk.assays {
    for count in assay.cell_counts {
      assert_eq(count, 3)
    }
  }
}

///|
test "dreamlet: pseudobulk library sizes equal column sums" {
  let assay = dreamlet_test_pseudobulk().assays[0]
  for sample in 0..<assay.sample_names.length() {
    let mut total = 0.0
    for gene in assay.counts {
      total = total + gene[sample]
    }
    dreamlet_test_close(total, assay.library_sizes[sample])
  }
}

///|
test "dreamlet: pseudobulk collapses sample metadata" {
  let assay = dreamlet_test_pseudobulk().assays[0]
  assert_eq(assay.metadata["Treatment"][0], "Control")
  assert_eq(assay.metadata["Treatment"][1], "Treated")
  assert_eq(assay.metadata["Donor"][0], "D1")
  assert_eq(assay.metadata["Donor"][1], "D1")
}

///|
test "dreamlet: pseudobulk fills missing sample-cell type combinations" {
  let metadata : Map[String, Array[String]] = Map([])
  let result = @src.dreamlet_aggregate_to_pseudobulk(
    [[3.0, 4.0, 5.0]],
    ["g"],
    ["c1", "c2", "c3"],
    ["s1", "s2", "s2"],
    ["A", "A", "B"],
    cell_metadata=metadata,
  ) catch {
    _ => abort("sparse combinations should aggregate")
  }
  let b = match result.assay("B") {
    Some(value) => value
    None => abort("B assay should exist")
  }
  assert_eq(b.cell_counts, [0, 1])
  assert_eq(b.counts[0], [0.0, 5.0])
}

///|
test "dreamlet: pseudobulk summary reports dimensions" {
  let summary = dreamlet_test_pseudobulk().summary()
  assert_true(summary.contains("Genes: 6"))
  assert_true(summary.contains("Samples: 8"))
  assert_true(summary.contains("Cell types: 2"))
}

///|
test "dreamlet: aggregation rejects negative counts" {
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_to_pseudobulk([[-1.0]], ["g"], ["c"], ["s"], ["A"]),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: aggregation rejects non-integer raw counts" {
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_to_pseudobulk([[1.5]], ["g"], ["c"], ["s"], ["A"]),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: aggregation rejects duplicate gene names" {
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_to_pseudobulk(
        [[1.0], [2.0]],
        ["g", "g"],
        ["c"],
        ["s"],
        ["A"],
      ),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: aggregation rejects metadata varying within sample" {
  let metadata : Map[String, Array[String]] = Map([
    ("Treatment", ["Control", "Treated"]),
  ])
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_to_pseudobulk(
        [[1.0, 2.0]],
        ["g"],
        ["c1", "c2"],
        ["s", "s"],
        ["A", "A"],
        cell_metadata=metadata,
      ),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: aggregation rejects metadata dimension mismatch" {
  let metadata : Map[String, Array[String]] = Map([("Group", ["A"])])
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_to_pseudobulk(
        [[1.0, 2.0]],
        ["g"],
        ["c1", "c2"],
        ["s1", "s2"],
        ["A", "A"],
        cell_metadata=metadata,
      ),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: TMM factors have geometric mean one" {
  let counts = [
    [100.0, 200.0, 90.0],
    [50.0, 100.0, 55.0],
    [30.0, 60.0, 25.0],
    [10.0, 20.0, 12.0],
  ]
  let result = @src.dreamlet_tmm(counts) catch {
    _ => abort("valid TMM normalization should succeed")
  }
  let mut product = 1.0
  for factor in result.factors {
    product = product * factor
  }
  dreamlet_test_close(product, 1.0, tolerance=1.0e-10)
}

///|
test "dreamlet: TMM identical compositions produce unit factors" {
  let result = @src.dreamlet_tmm([
    [10.0, 20.0, 30.0],
    [5.0, 10.0, 15.0],
    [2.0, 4.0, 6.0],
  ]) catch {
    _ => abort("proportional libraries should normalize")
  }
  for factor in result.factors {
    dreamlet_test_close(factor, 1.0, tolerance=1.0e-10)
  }
}

///|
test "dreamlet: TMM reference index is valid" {
  let result = @src.dreamlet_tmm([
    [10.0, 12.0, 9.0],
    [3.0, 7.0, 5.0],
    [8.0, 4.0, 6.0],
  ]) catch {
    _ => abort("valid TMM normalization should succeed")
  }
  assert_true(result.reference_sample >= 0)
  assert_true(result.reference_sample < 3)
}

///|
test "dreamlet: TMM is invariant to a common count multiplier" {
  let first = @src.dreamlet_tmm([
    [10.0, 12.0, 9.0],
    [3.0, 7.0, 5.0],
    [8.0, 4.0, 6.0],
  ]) catch {
    _ => abort("first TMM normalization should succeed")
  }
  let second = @src.dreamlet_tmm([
    [20.0, 24.0, 18.0],
    [6.0, 14.0, 10.0],
    [16.0, 8.0, 12.0],
  ]) catch {
    _ => abort("second TMM normalization should succeed")
  }
  for sample in 0..<3 {
    dreamlet_test_close(
      first.factors[sample],
      second.factors[sample],
      tolerance=1.0e-10,
    )
  }
}

///|
test "dreamlet: TMM rejects zero-size libraries" {
  let raised = try {
    ignore(@src.dreamlet_tmm([[1.0, 0.0], [2.0, 0.0]]))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: normalized counts are counts per million" {
  let normalized = @src.dreamlet_compute_norm_counts(
    [[10.0, 20.0], [90.0, 180.0]],
    [100.0, 200.0],
  ) catch {
    _ => abort("valid normalized counts should succeed")
  }
  assert_eq(normalized[0], [100000.0, 100000.0])
  assert_eq(normalized[1], [900000.0, 900000.0])
}

///|
test "dreamlet: log CPM uses prior counts for zeros" {
  let transformed = @src.dreamlet_compute_log_cpm(
    [[0.0, 0.0]],
    [100.0, 200.0],
    prior_count=0.5,
  ) catch {
    _ => abort("valid log CPM should succeed")
  }
  assert_true(transformed[0][0] == transformed[0][0])
  assert_true(transformed[0][1] == transformed[0][1])
  assert_true(transformed[0][0] > transformed[0][1])
}

///|
test "dreamlet: count transformations reject invalid library sizes" {
  let normalized = try {
    ignore(@src.dreamlet_compute_norm_counts([[1.0]], [0.0]))
    false
  } catch {
    DreamletError(_) => true
  }
  let logged = try {
    ignore(@src.dreamlet_compute_log_cpm([[1.0]], [-1.0]))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(normalized)
  assert_true(logged)
}

///|
test "dreamlet: processAssays retains both well-powered cell types" {
  let processed = dreamlet_test_processed()
  assert_eq(processed.assays.length(), 2)
  assert_eq(processed.exclusions.length(), 0)
}

///|
test "dreamlet: processAssays filters rare genes" {
  let assay = dreamlet_test_processed().assays[0]
  assert_false(assay.gene_names.contains("rare"))
  assert_true(assay.dropped_gene_names.contains("rare"))
  assert_eq(assay.gene_names.length(), 5)
}

///|
test "dreamlet: processAssays builds treatment design" {
  let assay = dreamlet_test_processed().assays[0]
  assert_eq(assay.design.coefficient_names, ["(Intercept)", "Treatment:Treated"])
  assert_eq(assay.design.sample_names, assay.sample_names)
}

///|
test "dreamlet: mixed processAssays reuses random-effect solver" {
  let assay = dreamlet_test_processed(mixed=true).assays[0]
  assert_eq(assay.design.random_effects.length(), 1)
  assert_eq(assay.design.random_effects[0].name, "Donor")
  assert_eq(assay.design.random_effects[0].levels.length(), 4)
}

///|
test "dreamlet: processAssays stores TMM effective libraries" {
  let assay = dreamlet_test_processed().assays[0]
  assert_eq(assay.norm_factors.length(), 8)
  assert_eq(assay.effective_library_sizes.length(), 8)
  for sample in 0..<8 {
    dreamlet_test_close(
      assay.effective_library_sizes[sample],
      assay.library_sizes[sample] * assay.norm_factors[sample],
    )
  }
}

///|
test "dreamlet: processAssays stores normalized expression dimensions" {
  let assay = dreamlet_test_processed().assays[0]
  assert_eq(assay.normalized_cpm.length(), assay.gene_names.length())
  assert_eq(assay.log_cpm.length(), assay.gene_names.length())
  assert_eq(assay.log_cpm[0].length(), assay.sample_names.length())
}

///|
test "dreamlet: Poisson initial weights are row-centered" {
  let assay = dreamlet_test_processed().assays[0]
  for row in assay.initial_weights {
    let mut total = 0.0
    for weight in row {
      assert_true(weight > 0.0)
      total = total + weight
    }
    dreamlet_test_close(total / row.length().to_double(), 1.0)
  }
}

///|
test "dreamlet: voom precision weights are finite and positive" {
  let assay = dreamlet_test_processed().assays[0]
  assert_eq(assay.precision_weights.length(), assay.gene_names.length())
  for row in assay.precision_weights {
    assert_eq(row.length(), assay.sample_names.length())
    for weight in row {
      assert_true(weight > 0.0)
      assert_true(weight <= 1.0e8)
      assert_true(weight == weight)
    }
  }
}

///|
test "dreamlet: voom trend coordinates are sorted" {
  let assay = dreamlet_test_processed().assays[0]
  assert_eq(assay.voom_x.length(), assay.gene_names.length())
  assert_eq(assay.voom_y.length(), assay.gene_names.length())
  assert_eq(assay.voom_fitted.length(), assay.gene_names.length())
  for index in 1..<assay.voom_x.length() {
    assert_true(assay.voom_x[index - 1] <= assay.voom_x[index])
  }
}

///|
test "dreamlet: processed assay lookup uses cell type name" {
  let processed = dreamlet_test_processed()
  assert_true(processed.assay("T_cell") is Some(_))
  assert_true(processed.assay("missing") is None)
}

///|
test "dreamlet: processAssays drops constant numeric effects" {
  let model = @src.dreamlet_model([
    @src.dreamlet_numeric_effect("Constant"),
    @src.dreamlet_categorical_effect("Treatment", reference="Control"),
  ]) catch {
    _ => abort("constant-effect model should build")
  }
  let pseudobulk = dreamlet_test_pseudobulk()
  for assay in pseudobulk.assays {
    assay.metadata["Constant"] = Array::make(8, "1")
  }
  let processed = @src.dreamlet_process_assays(
    pseudobulk,
    model,
    config=dreamlet_test_config(),
  ) catch {
    _ => abort("constant effects should be dropped")
  }
  assert_true(processed.assays[0].dropped_effects[0].contains("Constant"))
}

///|
test "dreamlet: processAssays reports missing model metadata" {
  let model = @src.dreamlet_model([@src.dreamlet_numeric_effect("Missing")]) catch {
    _ => abort("model itself should build")
  }
  let raised = try {
    ignore(
      @src.dreamlet_process_assays(
        dreamlet_test_pseudobulk(),
        model,
        config=dreamlet_test_config(),
      ),
    )
    false
  } catch {
    DreamletError(message) => message.contains("retained no assays")
  }
  assert_true(raised)
}

///|
test "dreamlet: processed summary reports model count" {
  let summary = dreamlet_test_processed().summary()
  assert_true(summary.contains("Retained cell types: 2"))
  assert_true(summary.contains("Excluded cell types: 0"))
  assert_true(summary.contains("Gene-cell type models: 10"))
}

///|
test "dreamlet: differential expression fits every retained cell type" {
  let result = dreamlet_test_result()
  assert_eq(result.assays.length(), 2)
  assert_eq(result.genes.length(), 10)
  assert_eq(result.contrast_name, "Treated-Control")
}

///|
test "dreamlet: response-up effect is positive in each cell type" {
  let result = dreamlet_test_result()
  let mut found = 0
  for gene in result.genes {
    if gene.gene_name == "response_up" {
      assert_true(gene.estimate > 0.5)
      found = found + 1
    }
  }
  assert_eq(found, 2)
}

///|
test "dreamlet: response-down effect is negative in each cell type" {
  let result = dreamlet_test_result()
  let mut found = 0
  for gene in result.genes {
    if gene.gene_name == "response_down" {
      assert_true(gene.estimate < -0.3)
      found = found + 1
    }
  }
  assert_eq(found, 2)
}

///|
test "dreamlet: result stores within-cluster and study-wide FDR" {
  let result = dreamlet_test_result()
  for gene in result.genes {
    assert_true(gene.within_cluster_fdr >= gene.p_value)
    assert_true(gene.study_wide_fdr >= gene.p_value)
    assert_true(gene.study_wide_fdr <= 1.0)
  }
}

///|
test "dreamlet: top table sorts all cell types by p-value" {
  let table = dreamlet_test_result().top_table() catch {
    _ => abort("valid dreamlet top table should succeed")
  }
  assert_eq(table.length(), 10)
  for index in 1..<table.length() {
    assert_true(table[index - 1].p_value <= table[index].p_value)
  }
}

///|
test "dreamlet: top table filters by cell type" {
  let table = dreamlet_test_result().top_table(cluster_id="Monocyte") catch {
    _ => abort("cell-type dreamlet top table should succeed")
  }
  assert_eq(table.length(), 5)
  for gene in table {
    assert_eq(gene.cluster_id, "Monocyte")
  }
}

///|
test "dreamlet: result assay lookup exposes underlying dream result" {
  let result = dreamlet_test_result()
  match result.assay("T_cell") {
    Some(value) => {
      assert_eq(value.contrast_name, "Treated-Control")
      assert_eq(value.genes.length(), 5)
    }
    None => abort("T_cell dream result should exist")
  }
  assert_true(result.assay("missing") is None)
}

///|
test "dreamlet: result summary reports study-wide analysis" {
  let summary = dreamlet_test_result().summary()
  assert_true(summary.contains("dreamlet Differential Expression Summary"))
  assert_true(summary.contains("Contrast: Treated-Control"))
  assert_true(summary.contains("Cell types: 2"))
  assert_true(summary.contains("Tests: 10"))
}

///|
test "dreamlet: differential expression rejects absent coefficient" {
  let raised = try {
    ignore(
      @src.dreamlet(
        dreamlet_test_processed(),
        "Missing",
        max_iterations=10,
        tolerance=0.05,
      ),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: top table validates FDR threshold" {
  let raised = try {
    ignore(dreamlet_test_result().top_table(maximum_fdr=1.1))
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: SCE integration preserves selected metadata" {
  let (counts, genes, cells, samples, clusters, metadata) = dreamlet_test_data()
  let sce = @src.SingleCellExperiment::new(counts, genes, cells)
  sce.col_data["sample_id"] = samples
  sce.col_data["cluster_id"] = clusters
  sce.col_data["Treatment"] = metadata["Treatment"]
  sce.col_data["Donor"] = metadata["Donor"]
  let pseudobulk = @src.dreamlet_aggregate_sce(sce, "sample_id", "cluster_id", metadata_fields=[
    "Treatment", "Donor",
  ]) catch {
    _ => abort("valid SCE dreamlet aggregation should succeed")
  }
  assert_eq(pseudobulk.assays.length(), 2)
  assert_eq(pseudobulk.assays[0].metadata["Treatment"][1], "Treated")
}

///|
test "dreamlet: SCE integration rejects absent assay" {
  let sce = @src.SingleCellExperiment::new([[1.0]], ["g"], ["c"])
  sce.col_data["sample"] = ["s"]
  sce.col_data["cluster"] = ["A"]
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_sce(sce, "sample", "cluster", assay="missing"),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}

///|
test "dreamlet: SCE integration rejects absent metadata field" {
  let sce = @src.SingleCellExperiment::new([[1.0]], ["g"], ["c"])
  sce.col_data["sample"] = ["s"]
  sce.col_data["cluster"] = ["A"]
  let raised = try {
    ignore(
      @src.dreamlet_aggregate_sce(sce, "sample", "cluster", metadata_fields=[
        "Missing",
      ]),
    )
    false
  } catch {
    DreamletError(_) => true
  }
  assert_true(raised)
}
