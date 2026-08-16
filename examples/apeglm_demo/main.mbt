///|
fn print_gene(result : @src.ApeglmResult, name : String) -> Unit {
  match result.find_gene(name) {
    Some(gene) =>
      println(
        "   " +
        name +
        ": MLE=" +
        (gene.mle[result.coefficient] / @math.ln(2.0)).to_string() +
        ", MAP=" +
        gene.effect(result.coefficient, log2_scale=true).to_string() +
        ", posterior SD=" +
        gene.standard_error(result.coefficient, log2_scale=true).to_string() +
        ", FSR=" +
        gene.fsr.to_string() +
        ", s-value=" +
        gene.s_value.to_string() +
        ", FSOS=" +
        gene.threshold_probability.to_string(),
      )
    None => abort("example gene not found: " + name)
  }
}

///|
fn main {
  println("=== Bioconductor apeglm Demo ===")
  let (counts, design, dispersions, genes, coefficients) = @src.apeglm_example_data()
  let config = @src.ApeglmConfig::create(coefficient=2, threshold=@math.ln(2.0)) catch {
    _ => abort("failed to create apeglm configuration")
  }

  println("\n1. Fit a batch-adjusted negative-binomial GLM")
  let result = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("failed to fit apeglm model")
  }
  let summary = result.summary(maximum_s_value=0.1)
  println(
    "   coefficient=" +
    summary.coefficient_name +
    ", genes=" +
    summary.gene_count.to_string() +
    ", converged=" +
    summary.converged_count.to_string() +
    ", selected=" +
    summary.selected_count.to_string(),
  )
  println(
    "   adaptive Cauchy scale=" +
    result.prior_control.prior_scale.to_string() +
    " (natural-log scale)",
  )

  println("\n2. Compare unshrunk MLEs with posterior MAP estimates")
  print_gene(result, "strong_up")
  print_gene(result, "strong_down")
  print_gene(result, "low_noisy")
  print_gene(result, "stable")

  println("\n3. Rank and select effects by directional error control")
  let ranked = result.ranked()
  println(
    "   top-ranked gene=" +
    ranked[0].gene_name +
    ", absolute MAP=" +
    ranked[0].map[result.coefficient].abs().to_string(),
  )
  println(
    "   genes with s-value <= 0.10: " +
    result.select(maximum_s_value=0.1).length().to_string(),
  )

  println("\n4. Add posterior summaries to an immutable container copy")
  let assays : Map[String, Array[Array[Double]]] = Map([])
  assays["counts"] = counts
  let metadata : Map[String, String] = Map([])
  metadata["source"] = "apeglm_demo"
  let experiment = @src.summarized_experiment(assays, [], [], metadata)
  let output = @src.apeglm_summarized_experiment(
    experiment,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("failed to integrate apeglm with SummarizedExperiment")
  }
  let posterior_map = match @src.se_assay(output.experiment, "apeglm_map") {
    Some(value) => value
    None => abort("posterior MAP assay missing")
  }
  println(
    "   posterior MAP rows=" +
    posterior_map.length().to_string() +
    ", source object unchanged=" +
    (@src.se_assay(experiment, "apeglm_map") is None).to_string(),
  )

  println("\n5. Export a DESeq2-style log2 table")
  println(result.to_tsv(log2_scale=true))

  println("6. Invalid model inputs produce explicit diagnostics")
  let rejected = try {
    ignore(
      @src.apeglm_fit([[1.0, -1.0], [2.0, 3.0]], [[1.0, 0.0], [1.0, 1.0]], [
        0.1, 0.1,
      ]),
    )
    false
  } catch {
    ApeglmError(message) => {
      println("   " + message)
      true
    }
  }
  println("   malformed counts rejected=" + rejected.to_string())
  println("\n=== Demo Complete ===")
}
