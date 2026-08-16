///|
fn main {
  println("=== Bioconductor zinbwave Demo ===")
  let (counts, genes, cells, batch) = @src.zinbwave_example_data()
  let config = @src.ZinbWaveConfig::create(
    n_factors=2,
    max_iterations=35,
    inner_iterations=2,
    tolerance=1.0e-5,
    ridge=0.01,
  ) catch {
    _ => abort("failed to create zinbwave configuration")
  }

  println("\n1. Fit a ZINB latent-factor model with a batch covariate")
  let model = @src.zinbwave_fit(
    counts,
    config~,
    cell_covariates=batch,
    gene_names=genes,
    cell_names=cells,
  ) catch {
    _ => abort("failed to fit zinbwave model")
  }
  println("   " + model.summary())
  println(
    "   log-likelihood: " +
    model.diagnostics.initial_log_likelihood.to_string() +
    " -> " +
    model.diagnostics.final_log_likelihood.to_string(),
  )

  println("\n2. Inspect the inferred cell representation")
  for cell in 0..<model.n_cells() {
    println(
      "   " + model.cell_names[cell] + ": " + model.factors[cell].to_string(),
    )
  }

  println("\n3. Inspect zero-inflation diagnostics")
  match model.gene_summary(counts, 5) {
    Some(summary) =>
      println(
        "   " +
        summary.gene_name +
        ": observed zero fraction=" +
        summary.observed_zero_fraction.to_string() +
        ", fitted=" +
        summary.fitted_zero_fraction.to_string() +
        ", dispersion=" +
        summary.dispersion.to_string(),
      )
    None => abort("dropout gene summary should exist")
  }
  println(
    "   structural-zero posterior weight at dropout_gene/A1: " +
    (1.0 - model.observational_weights[5][0]).to_string(),
  )

  println("\n4. Generate normalized, imputed, and residual matrices")
  let normalized = model.normalized_values()
  let imputed = model.impute_zeros(counts) catch {
    _ => abort("failed to impute model zeros")
  }
  println("   normalized A_marker_1/A1=" + normalized[0][0].to_string())
  println("   imputed dropout_gene/A1=" + imputed[5][0].to_string())
  println(
    "   deviance residual A_marker_1/B1=" +
    model.deviance_residuals[0][4].to_string(),
  )

  println("\n5. Add outputs to an immutable SingleCellExperiment copy")
  let sce = @src.SingleCellExperiment::new(counts, genes, cells)
  let output = @src.zinbwave_sce(sce, config~, cell_covariates=batch) catch {
    _ => abort("failed to integrate zinbwave with SingleCellExperiment")
  }
  println(
    "   reduced dimensions: " +
    @src.sce_get_reduced_dim(output.experiment, "zinbwave").length().to_string() +
    " cells",
  )
  println(
    "   observation-weight assay: " +
    @src.sce_get_assay(output.experiment, "zinbwave_weights")
    .length()
    .to_string() +
    " genes",
  )
  println(
    "   source object unchanged: " +
    (@src.sce_get_assay(sce, "zinbwave_weights").length() == 0).to_string(),
  )

  println("\n=== Demo Complete ===")
}
