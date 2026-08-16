///|
/// Tests for Bioconductor variancePartition-inspired mixed models.

///|
fn vp_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "variancePartition value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn vp_test_samples() -> Array[String] {
  let samples : Array[String] = []
  for subject in 0..<4 {
    for replicate in 0..<4 {
      samples.push(
        "S" + (subject + 1).to_string() + "_R" + (replicate + 1).to_string(),
      )
    }
  }
  samples
}

///|
fn vp_test_subjects() -> Array[String] {
  let subjects : Array[String] = []
  for subject in 0..<4 {
    for _ in 0..<4 {
      subjects.push("S" + (subject + 1).to_string())
    }
  }
  subjects
}

///|
fn vp_test_treatments() -> Array[String] {
  let treatments : Array[String] = []
  for _ in 0..<4 {
    treatments.push("Control")
    treatments.push("Treated")
    treatments.push("Control")
    treatments.push("Treated")
  }
  treatments
}

///|
fn vp_test_batches() -> Array[String] {
  let batches : Array[String] = []
  for _ in 0..<4 {
    batches.push("B1")
    batches.push("B1")
    batches.push("B2")
    batches.push("B2")
  }
  batches
}

///|
fn vp_test_design(
  include_batch? : Bool = false,
) -> @src.VariancePartitionDesign {
  let treatment = @src.vp_categorical_effect(
    "Treatment",
    vp_test_treatments(),
    reference="Control",
  ) catch {
    _ => abort("valid treatment effect should build")
  }
  let subject = @src.vp_random_effect("Subject", vp_test_subjects()) catch {
    _ => abort("valid subject effect should build")
  }
  let random_effects = [subject]
  if include_batch {
    let batch = @src.vp_random_effect("Batch", vp_test_batches()) catch {
      _ => abort("valid batch effect should build")
    }
    random_effects.push(batch)
  }
  @src.vp_design(vp_test_samples(), [treatment], random_effects) catch {
    _ => abort("valid variancePartition design should build")
  }
}

///|
fn vp_test_subject_gene() -> Array[Double] {
  let values : Array[Double] = []
  let noise = [0.08, -0.08, 0.04, -0.04]
  for subject in 0..<4 {
    let subject_effect = (subject.to_double() - 1.5) * 2.0
    for replicate in 0..<4 {
      let treatment_effect = if replicate % 2 == 1 { 0.1 } else { 0.0 }
      values.push(10.0 + subject_effect + treatment_effect + noise[replicate])
    }
  }
  values
}

///|
fn vp_test_treatment_gene() -> Array[Double] {
  let values : Array[Double] = []
  let noise = [0.05, -0.03, -0.05, 0.03]
  for subject in 0..<4 {
    let subject_effect = (subject.to_double() - 1.5) * 0.08
    for replicate in 0..<4 {
      let treatment_effect = if replicate % 2 == 1 { 4.0 } else { 0.0 }
      values.push(6.0 + subject_effect + treatment_effect + noise[replicate])
    }
  }
  values
}

///|
fn vp_test_batch_gene() -> Array[Double] {
  let values : Array[Double] = []
  let noise = [0.03, -0.03, 0.02, -0.02]
  for subject in 0..<4 {
    let subject_effect = (subject.to_double() - 1.5) * 0.1
    for replicate in 0..<4 {
      let batch_effect = if replicate >= 2 { 3.0 } else { 0.0 }
      values.push(8.0 + subject_effect + batch_effect + noise[replicate])
    }
  }
  values
}

///|
fn vp_test_residual_gene() -> Array[Double] {
  let pattern = [1.0, -1.0, -1.0, 1.0]
  let values : Array[Double] = []
  for subject in 0..<4 {
    let scale = 1.0 + subject.to_double() * 0.1
    for replicate in 0..<4 {
      values.push(5.0 + scale * pattern[replicate])
    }
  }
  values
}

///|
fn vp_test_expression() -> Array[Array[Double]] {
  [vp_test_subject_gene(), vp_test_treatment_gene(), vp_test_residual_gene()]
}

///|
fn vp_test_assay() -> @src.SummarizedExperiment {
  let expression = vp_test_expression()
  let weights : Array[Array[Double]] = []
  for _ in 0..<expression.length() {
    weights.push(Array::make(vp_test_samples().length(), 1.0))
  }
  let assays : Map[String, Array[Array[Double]]] = Map([
    ("logcounts", expression),
    ("weights", weights),
  ])
  @src.summarized_experiment(assays, [], [], Map([]))
}

///|
test "variancePartition: numeric fixed effect preserves values" {
  let effect = @src.vp_numeric_effect("Age", [20.0, 30.0, 40.0]) catch {
    _ => abort("valid numeric effect should build")
  }
  assert_eq(effect.name, "Age")
  assert_eq(effect.column_names, ["Age"])
  assert_eq(effect.columns, [[20.0], [30.0], [40.0]])
}

///|
test "variancePartition: categorical effect uses requested reference" {
  let effect = @src.vp_categorical_effect(
    "Disease",
    ["Case", "Control", "Case", "Other"],
    reference="Control",
  ) catch {
    _ => abort("valid categorical effect should build")
  }
  assert_eq(effect.column_names, ["Disease:Case", "Disease:Other"])
  assert_eq(effect.columns[0], [1.0, 0.0])
  assert_eq(effect.columns[1], [0.0, 0.0])
  assert_eq(effect.columns[3], [0.0, 1.0])
}

///|
test "variancePartition: random effect encoding is stable" {
  let effect = @src.vp_random_effect("Subject", ["B", "A", "B", "A"]) catch {
    _ => abort("valid random effect should build")
  }
  assert_eq(effect.levels, ["B", "A"])
  assert_eq(effect.level_indices, [0, 1, 0, 1])
}

///|
test "variancePartition: design expands intercept and treatment" {
  let design = vp_test_design()
  assert_eq(design.sample_names.length(), 16)
  assert_eq(design.coefficient_names, ["(Intercept)", "Treatment:Treated"])
  assert_eq(design.matrix[0], [1.0, 0.0])
  assert_eq(design.matrix[1], [1.0, 1.0])
  assert_eq(design.random_effects[0].levels.length(), 4)
}

///|
test "variancePartition: constructors copy caller arrays" {
  let ages = [20.0, 30.0, 40.0]
  let effect = @src.vp_numeric_effect("Age", ages) catch {
    _ => abort("valid numeric effect should build")
  }
  ages[0] = 99.0
  assert_eq(effect.columns[0][0], 20.0)
  let samples = vp_test_samples()
  let design = vp_test_design()
  samples[0] = "changed"
  assert_eq(design.sample_names[0], "S1_R1")
}

///|
test "variancePartition: subject random variance dominates subject gene" {
  let fit = @src.vp_fit_gene(
    vp_test_subject_gene(),
    vp_test_design(),
    gene_name="subject_gene",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("subject mixed model should fit")
  }
  assert_true(fit.converged)
  assert_true(fit.random_variances[0] > fit.residual_variance)
  let fraction = fit.fraction("Subject")
  match fraction {
    Some(value) => assert_true(value > 0.8)
    None => abort("Subject fraction should be present")
  }
}

///|
test "variancePartition: fixed variance dominates treatment gene" {
  let fit = @src.vp_fit_gene(
    vp_test_treatment_gene(),
    vp_test_design(),
    gene_name="treatment_gene",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("treatment mixed model should fit")
  }
  vp_test_close(fit.coefficients[1], 4.0, 0.15)
  match fit.fraction("Treatment") {
    Some(value) => assert_true(value > 0.8)
    None => abort("Treatment fraction should be present")
  }
}

///|
test "variancePartition: variance fractions sum to one" {
  let fit = @src.vp_fit_gene(
    vp_test_subject_gene(),
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("mixed model should fit")
  }
  let mut total = 0.0
  for value in fit.variance_fractions {
    assert_true(value >= 0.0)
    total = total + value
  }
  vp_test_close(total, 1.0, 1.0e-10)
}

///|
test "variancePartition: BLUPs recover subject ordering" {
  let fit = @src.vp_fit_gene(
    vp_test_subject_gene(),
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("mixed model should fit")
  }
  let first = match fit.blup("Subject", "S1") {
    Some(value) => value
    None => abort("S1 BLUP should be present")
  }
  let last = match fit.blup("Subject", "S4") {
    Some(value) => value
    None => abort("S4 BLUP should be present")
  }
  assert_true(first < 0.0)
  assert_true(last > 0.0)
  assert_true(last > first)
}

///|
test "variancePartition: fitted values and residuals preserve sample order" {
  let expression = vp_test_treatment_gene()
  let fit = @src.vp_fit_gene(
    expression,
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("mixed model should fit")
  }
  assert_eq(fit.fitted_values.length(), expression.length())
  assert_eq(fit.residuals.length(), expression.length())
  for index in 0..<expression.length() {
    vp_test_close(
      fit.fitted_values[index] + fit.residuals[index],
      expression[index],
      1.0e-10,
    )
  }
}

///|
test "variancePartition: uniform precision weights match unweighted fit" {
  let expression = vp_test_treatment_gene()
  let unweighted = @src.vp_fit_gene(
    expression,
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("unweighted model should fit")
  }
  let weighted = @src.vp_fit_gene(
    expression,
    vp_test_design(),
    weights=Array::make(expression.length(), 5.0),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("weighted model should fit")
  }
  vp_test_close(unweighted.coefficients[1], weighted.coefficients[1], 1.0e-8)
  vp_test_close(
    unweighted.residual_variance,
    weighted.residual_variance,
    1.0e-8,
  )
}

///|
test "variancePartition: heterogeneous precision weights are accepted" {
  let expression = vp_test_treatment_gene()
  let weights = Array::make(expression.length(), 1.0)
  weights[0] = 0.2
  let fit = @src.vp_fit_gene(
    expression,
    vp_test_design(),
    weights~,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("heterogeneous weighted model should fit")
  }
  assert_true(fit.residual_variance > 0.0)
  assert_true(fit.log_likelihood < 1.0e300)
}

///|
test "variancePartition: ML and REML fits are identified" {
  let reml_fit = @src.vp_fit_gene(
    vp_test_subject_gene(),
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("REML model should fit")
  }
  let ml_fit = @src.vp_fit_gene(
    vp_test_subject_gene(),
    vp_test_design(),
    reml=false,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("ML model should fit")
  }
  assert_true(reml_fit.reml)
  assert_false(ml_fit.reml)
  assert_true(reml_fit.log_likelihood != ml_fit.log_likelihood)
}

///|
test "variancePartition: multiple random intercepts are fitted" {
  let fit = @src.vp_fit_gene(
    vp_test_batch_gene(),
    vp_test_design(include_batch=true),
    gene_name="batch_gene",
    max_iterations=100,
    tolerance=0.01,
  ) catch {
    _ => abort("multi-random-effect model should fit")
  }
  assert_eq(fit.random_effect_names, ["Subject", "Batch"])
  assert_eq(fit.random_variances.length(), 2)
  assert_true(fit.random_variances[1] > fit.random_variances[0])
  match fit.fraction("Batch") {
    Some(value) => assert_true(value > 0.5)
    None => abort("Batch fraction should be present")
  }
}

///|
test "variancePartition: fixed-only design reduces to linear model" {
  let treatment = @src.vp_categorical_effect(
    "Treatment",
    vp_test_treatments(),
    reference="Control",
  ) catch {
    _ => abort("valid fixed effect should build")
  }
  let design = @src.vp_design(vp_test_samples(), [treatment], []) catch {
    _ => abort("fixed-only design should build")
  }
  let fit = @src.vp_fit_gene(
    vp_test_treatment_gene(),
    design,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("fixed-only model should fit")
  }
  vp_test_close(fit.coefficients[1], 4.0, 0.15)
  assert_eq(fit.random_variances.length(), 0)
}

///|
test "variancePartition: matrix fit generates gene names" {
  let result = @src.fit_extract_variance_partition(
    vp_test_expression(),
    vp_test_design(),
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("variance partition matrix should fit")
  }
  assert_eq(result.gene_names, ["Gene1", "Gene2", "Gene3"])
  assert_eq(result.fractions.length(), 3)
  assert_eq(result.fits.length(), 3)
  assert_false(result.fits[0].reml)
}

///|
test "variancePartition: matrix result supports named lookup" {
  let result = @src.fit_extract_variance_partition(
    vp_test_expression(),
    vp_test_design(),
    gene_names=["subject", "treatment", "residual"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("variance partition matrix should fit")
  }
  match result.fraction("subject", "Subject") {
    Some(value) => assert_true(value > 0.8)
    None => abort("named variance fraction should be present")
  }
  assert_true(result.fraction("missing", "Subject") is None)
}

///|
test "variancePartition: matrix result reports component medians" {
  let result = @src.fit_extract_variance_partition(
    vp_test_expression(),
    vp_test_design(),
    gene_names=["subject", "treatment", "residual"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("variance partition matrix should fit")
  }
  assert_eq(result.component_names, ["Treatment", "Subject", "Residuals"])
  assert_eq(result.median_fractions.length(), 3)
  for value in result.median_fractions {
    assert_true(value >= 0.0 && value <= 1.0)
  }
}

///|
test "variancePartition: result summary contains component names" {
  let result = @src.fit_extract_variance_partition(
    [vp_test_subject_gene()],
    vp_test_design(),
    gene_names=["subject"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("variance partition matrix should fit")
  }
  let summary = result.summary()
  assert_true(summary.contains("variancePartition Summary"))
  assert_true(summary.contains("Subject median fraction"))
  assert_true(summary.contains("Genes: 1"))
}

///|
test "variancePartition: dream estimates repeated-measures contrast" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    vp_test_expression(),
    design,
    contrast,
    gene_names=["subject", "treatment", "residual"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  vp_test_close(result.genes[1].estimate, 4.0, 0.15)
  assert_true(result.genes[1].standard_error > 0.0)
  assert_true(result.genes[1].p_value < 0.001)
}

///|
test "variancePartition: dream uses bounded small-sample degrees of freedom" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    [vp_test_treatment_gene()],
    design,
    contrast,
    gene_names=["treatment"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  let degrees = result.genes[0].degrees_of_freedom
  assert_true(degrees >= 1.0)
  assert_true(degrees <= 14.0)
}

///|
test "variancePartition: dream applies monotone BH correction" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    vp_test_expression(),
    design,
    contrast,
    gene_names=["subject", "treatment", "residual"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  for gene in result.genes {
    assert_true(gene.adjusted_p_value >= gene.p_value)
    assert_true(gene.adjusted_p_value <= 1.0)
  }
}

///|
test "variancePartition: dream top table sorts by p-value" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    vp_test_expression(),
    design,
    contrast,
    gene_names=["subject", "treatment", "residual"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  let table = result.top_table() catch {
    _ => abort("valid top table should build")
  }
  assert_eq(table.length(), 3)
  for index in 1..<table.length() {
    assert_true(table[index - 1].p_value <= table[index].p_value)
  }
}

///|
test "variancePartition: custom dream contrast is supported" {
  let contrast = @src.vp_contrast("Treatment", [0.0, 1.0]) catch {
    _ => abort("valid custom contrast should build")
  }
  assert_eq(contrast.name, "Treatment")
  assert_eq(contrast.coefficients, [0.0, 1.0])
}

///|
test "variancePartition: dream summary reports contrast and significance" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    [vp_test_treatment_gene()],
    design,
    contrast,
    gene_names=["treatment"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  let summary = result.summary()
  assert_true(summary.contains("dream Differential Expression Summary"))
  assert_true(summary.contains("Contrast: Treated-Control"))
  assert_true(summary.contains("FDR <= 0.05: 1"))
}

///|
test "variancePartition: SummarizedExperiment variance entry uses assays" {
  let result = @src.fit_extract_variance_partition_se(
    vp_test_assay(),
    "logcounts",
    vp_test_design(),
    gene_names=["subject", "treatment", "residual"],
    weights_assay="weights",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("SummarizedExperiment variance analysis should fit")
  }
  assert_eq(result.gene_names.length(), 3)
  match result.fraction("subject", "Subject") {
    Some(value) => assert_true(value > 0.8)
    None => abort("Subject fraction should be present")
  }
}

///|
test "variancePartition: SummarizedExperiment dream entry uses assays" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream_se(
    vp_test_assay(),
    "logcounts",
    design,
    contrast,
    gene_names=["subject", "treatment", "residual"],
    weights_assay="weights",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("SummarizedExperiment dream analysis should fit")
  }
  assert_true(result.genes[1].p_value < 0.001)
}

///|
test "variancePartition: rejects categorical reference that is absent" {
  let failed = try {
    ignore(
      @src.vp_categorical_effect(
        "Disease",
        ["Case", "Control"],
        reference="Other",
      ),
    )
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects unreplicated random effects" {
  let failed = try {
    ignore(@src.vp_random_effect("Subject", ["A", "B", "C"]))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects duplicate sample names" {
  let age = @src.vp_numeric_effect("Age", [1.0, 2.0, 3.0]) catch {
    _ => abort("numeric effect should build")
  }
  let failed = try {
    ignore(@src.vp_design(["A", "A", "B"], [age], []))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects rank-deficient fixed design" {
  let first = @src.vp_numeric_effect("First", [0.0, 1.0, 2.0, 3.0]) catch {
    _ => abort("numeric effect should build")
  }
  let second = @src.vp_numeric_effect("Second", [0.0, 2.0, 4.0, 6.0]) catch {
    _ => abort("numeric effect should build")
  }
  let failed = try {
    ignore(@src.vp_design(["A", "B", "C", "D"], [first, second], []))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects expression dimension mismatch" {
  let failed = try {
    ignore(@src.vp_fit_gene([1.0, 2.0], vp_test_design()))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects nonpositive precision weights" {
  let weights = Array::make(16, 1.0)
  weights[3] = 0.0
  let failed = try {
    ignore(@src.vp_fit_gene(vp_test_subject_gene(), vp_test_design(), weights~))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects duplicate gene names" {
  let failed = try {
    ignore(
      @src.fit_extract_variance_partition(
        [vp_test_subject_gene(), vp_test_treatment_gene()],
        vp_test_design(),
        gene_names=["gene", "gene"],
      ),
    )
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects zero contrast" {
  let failed = try {
    ignore(@src.vp_contrast("zero", [0.0, 0.0]))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects contrast dimension mismatch" {
  let failed = try {
    ignore(
      @src.dream(
        [vp_test_treatment_gene()],
        vp_test_design(),
        @src.vp_contrast("bad", [1.0]) catch {
          _ => abort("nonzero contrast should build")
        },
      ),
    )
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: rejects invalid top-table FDR" {
  let design = vp_test_design()
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("valid coefficient contrast should build")
  }
  let result = @src.dream(
    [vp_test_treatment_gene()],
    design,
    contrast,
    gene_names=["treatment"],
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("dream analysis should fit")
  }
  let failed = try {
    ignore(result.top_table(maximum_fdr=1.5))
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}

///|
test "variancePartition: reports missing SummarizedExperiment assay" {
  let failed = try {
    ignore(
      @src.fit_extract_variance_partition_se(
        vp_test_assay(),
        "missing",
        vp_test_design(),
      ),
    )
    false
  } catch {
    VariancePartitionError(_) => true
  }
  assert_true(failed)
}
