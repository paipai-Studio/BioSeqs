///|
fn mast_demo_option(value : Double?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "NA"
  }
}

///|
fn main {
  println("=== Bioconductor MAST 1.39.0 Advanced Demo ===")

  let (data, tested) = @src.mast_advanced_example() catch {
    _ => abort("failed to construct the MAST example")
  }
  println("\n1. Feature x cell model and treatment-coded design")
  println(
    "  genes=" +
    data.n_genes.to_string() +
    ", cells=" +
    data.n_cells.to_string() +
    ", coefficients=" +
    data.coefficient_names.join(", "),
  )
  println(
    "  tested coefficient=" +
    data.coefficient_names[tested[0]] +
    ", first-cell CDR=" +
    data.cdr[0].to_string(),
  )

  let config = @src.MastAdvancedConfig::create(
    empirical_bayes=true,
    ebayes_use_full_model=false,
    fdr_threshold=0.05,
  ) catch {
    _ => abort("failed to construct the MAST configuration")
  }
  let result = @src.mast_advanced_zlm(data, tested, config~) catch {
    _ => abort("advanced MAST zlm fit failed")
  }
  println("\n2. Bayesian logistic + positive Gaussian hurdle fit")
  println("  " + result.summary())
  println(
    "  eBayes prior variance=" +
    result.prior_variance.to_string() +
    ", prior df=" +
    result.prior_df.to_string(),
  )

  println("\n3. Component LRT, hurdle LRT, FDR and marginal logFC")
  for gene in 0..<result.gene_names.length() {
    println(
      "  " +
      result.gene_names[gene] +
      ": detected=" +
      result.n_detected[gene].to_string() +
      "/" +
      data.n_cells.to_string() +
      ", pD=" +
      mast_demo_option(result.discrete_p_values[gene]) +
      ", pC=" +
      mast_demo_option(result.continuous_p_values[gene]) +
      ", pH=" +
      mast_demo_option(result.hurdle_p_values[gene]) +
      ", FDR=" +
      mast_demo_option(result.hurdle_fdr[gene]) +
      ", logFC=" +
      mast_demo_option(result.marginal_log_fc[gene]),
    )
  }

  println("\n4. Top genes by hurdle p-value")
  let top = result.top_genes(3)
  for index in 0..<top.length() {
    println("  " + (index + 1).to_string() + ". " + top[index])
  }

  println("\n5. Immutable SingleCellExperiment write-back")
  let experiment = @src.SingleCellExperiment::new(
    data.expression,
    data.gene_names,
    data.cell_names,
  )
  experiment.assays["logcounts"] = data.expression.map(fn(row) { row.copy() })
  experiment.col_data["condition"] = Array::make(data.n_cells, "control")
  for cell in 24..<data.n_cells {
    experiment.col_data["condition"][cell] = "stimulated"
  }
  let output = @src.mast_advanced_zlm_sce(
    experiment,
    "condition",
    reference="control",
    output_prefix="mast39",
    config~,
  ) catch {
    _ => abort("advanced MAST SingleCellExperiment fit failed")
  }
  println(
    "  original has mast39.hurdleFdr: " +
    experiment.row_data.contains("mast39.hurdleFdr").to_string(),
  )
  println(
    "  enriched row columns include hurdle FDR: " +
    output.experiment.row_data.contains("mast39.hurdleFdr").to_string(),
  )
  println(
    "  contrast metadata: " +
    output.experiment.metadata["mast39.contrast"] +
    ", tested genes=" +
    output.experiment.metadata["mast39.tested"],
  )

  println("\nMAST advanced demo completed.")
}
