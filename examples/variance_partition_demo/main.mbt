///|
fn demo_samples() -> Array[String] {
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
fn demo_subjects() -> Array[String] {
  let subjects : Array[String] = []
  for subject in 0..<4 {
    for _ in 0..<4 {
      subjects.push("S" + (subject + 1).to_string())
    }
  }
  subjects
}

///|
fn demo_treatments() -> Array[String] {
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
fn demo_subject_gene() -> Array[Double] {
  let values : Array[Double] = []
  let noise = [0.08, 0.04, -0.08, -0.04]
  for subject in 0..<4 {
    let subject_effect = (subject.to_double() - 1.5) * 2.0
    for replicate in 0..<4 {
      values.push(10.0 + subject_effect + noise[replicate])
    }
  }
  values
}

///|
fn demo_treatment_gene() -> Array[Double] {
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
fn main {
  println("=== Bioconductor variancePartition Demo ===")
  let treatment = @src.vp_categorical_effect(
    "Treatment",
    demo_treatments(),
    reference="Control",
  ) catch {
    _ => abort("failed to encode treatment")
  }
  let subject = @src.vp_random_effect("Subject", demo_subjects()) catch {
    _ => abort("failed to encode subject")
  }
  let design = @src.vp_design(demo_samples(), [treatment], [subject]) catch {
    _ => abort("failed to build mixed-model design")
  }
  let expression = [demo_subject_gene(), demo_treatment_gene()]
  let weights = [
    Array::make(16, 1.0),
    [
      0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0,
    ],
  ]

  println("\n1. Typed fixed and random effects")
  println("   coefficients: " + design.coefficient_names.to_string())
  println("   random effect: " + design.random_effects[0].name)
  println("   levels: " + design.random_effects[0].levels.to_string())

  println("\n2. ML variance fractions with observation-level weights")
  let partition = @src.fit_extract_variance_partition(
    expression,
    design,
    gene_names=["subject_gene", "treatment_gene"],
    weights~,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("failed to partition variance")
  }
  println(partition.summary())
  for fit in partition.fits {
    println("   " + fit.gene_name + ":")
    for component in 0..<fit.component_names.length() {
      println(
        "     " +
        fit.component_names[component] +
        "=" +
        fit.variance_fractions[component].to_string(),
      )
    }
  }

  println("\n3. Random-intercept BLUPs")
  let subject_fit = partition.fits[0]
  for index in 0..<subject_fit.blups[0].levels.length() {
    println(
      "   " +
      subject_fit.blups[0].levels[index] +
      ": " +
      subject_fit.blups[0].estimates[index].to_string(),
    )
  }

  println("\n4. dream repeated-measures contrast")
  let contrast = @src.vp_coefficient_contrast(
    design, "Treated-Control", "Treatment:Treated",
  ) catch {
    _ => abort("failed to build treatment contrast")
  }
  let differential = @src.dream(
    expression,
    design,
    contrast,
    gene_names=["subject_gene", "treatment_gene"],
    weights~,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("failed to run dream")
  }
  for gene in differential.genes {
    println(
      "   " +
      gene.gene_name +
      ": estimate=" +
      gene.estimate.to_string() +
      ", t=" +
      gene.statistic.to_string() +
      ", df=" +
      gene.degrees_of_freedom.to_string() +
      ", FDR=" +
      gene.adjusted_p_value.to_string(),
    )
  }

  println("\n5. SummarizedExperiment assay integration")
  let experiment = @src.summarized_experiment(
    Map([("logcounts", expression), ("weights", weights)]),
    [],
    [],
    Map([("study", "repeated-measures-demo")]),
  )
  let assay_result = @src.dream_se(
    experiment,
    "logcounts",
    design,
    contrast,
    gene_names=["subject_gene", "treatment_gene"],
    weights_assay="weights",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("failed to run assay-backed dream")
  }
  println("   " + assay_result.summary())

  println(
    "\nScope: random intercepts and numerical Satterthwaite tests; random slopes, Kenward-Roger, voom/eBayes and plotting are not included.",
  )
  println("\n=== Demo Complete ===")
}
