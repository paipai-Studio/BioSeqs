///|
/// limma Demo - Differential Expression Analysis with Linear Models.

///|
fn main {
  println("=== limma Differential Expression Analysis Demo ===\n")

  println("1. Creating simulated gene expression data...")
  let counts_data = Array::make(100 * 6, 0.0)
  let mut i = 0
  while i < 100 {
    let mut j = 0
    while j < 6 {
      if i < 10 {
        if j < 3 {
          counts_data[i * 6 + j] = 100.0 +
            i.to_double() * 5.0 +
            (j - 1).to_double() * 10.0
        } else {
          counts_data[i * 6 + j] = 500.0 +
            i.to_double() * 20.0 +
            (j - 4).to_double() * 20.0
        }
      } else {
        counts_data[i * 6 + j] = 100.0 +
          (i * i).to_double() % 50.0 +
          j.to_double() * 2.0
      }
      j = j + 1
    }
    i = i + 1
  }
  let counts = @bio.Matrix::new(100, 6, counts_data)

  println("   Created 100 genes x 6 samples matrix")
  println("   First 10 genes are differentially expressed (treatment effect)\n")

  println("2. Applying voom transformation...")
  let design = @bio.Matrix::from_2d(6, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
  ])
  let log_cpm = @bio.voom(counts, design)
  println("   voom transformation complete\n")

  // Batch effect correction example
  println("3. Removing batch effects...")
  let batch = ["B1", "B1", "B1", "B2", "B2", "B2"]
  let group = ["Control", "Control", "Control", "Treatment", "Treatment", "Treatment"]
  let batch_removed = @bio.limma_remove_batch_effect(log_cpm, batch, group = group)
  println("   Batch effect correction applied")
  println(
    "   Corrected matrix: " +
    batch_removed.rows.to_string() +
    " x " +
    batch_removed.cols.to_string() +
    "\n",
  )

  // MDS plot distances
  println("4. Computing MDS distances for visualization...")
  let mds = @bio.plot_mds_distances(log_cpm, 50)
  println("   MDS distances computed (top 50 variable genes)")
  println("   Distance matrix: " + mds.rows.to_string() + " x " + mds.cols.to_string() + "\n")

  println("5. Fitting linear model (lmFit)...")
  let fit = @bio.lm_fit(log_cpm, design)
  println(
    "   Coefficients matrix: " +
    fit.coefficients.rows.to_string() +
    " x " +
    fit.coefficients.cols.to_string(),
  )
  println(
    "   Residual degrees of freedom: " + fit.df_residual.to_string() + "\n",
  )

  println("6. Applying empirical Bayes moderation (eBayes)...")
  let eb_result = @bio.ebayes(fit, 0.01)
  println("   Prior variance: " + eb_result.prior_var.to_string())
  println("   Prior degrees of freedom: " + eb_result.prior_df.to_string())
  println(
    "   Total degrees of freedom: " + eb_result.df_total.to_string() + "\n",
  )

  // TREAT test example
  println("7. Running TREAT test...")
  let treat_result = @bio.treat(eb_result, 0.05, "both")
  println("   TREAT test complete (lfc=0.05, direction=both)")
  println(
    "   TREAT p-values: " +
    treat_result.p_value.rows.to_string() +
    " rows\n",
  )

  println("8. Extracting top differentially expressed genes (topTable)...")
  let top = @bio.top_table(eb_result, 1, 10, "BH")
  println("   Top 10 differentially expressed genes:\n")
  println("   Rank | logFC | AveExpr | t-value | P.Value | adj.P.Val | B")
  println("   -----|-------|---------|---------|---------|-----------|-----")

  i = 0
  while i < top.length() {
    let row = top[i]
    println(
      "   " +
      (i + 1).to_string().pad_start(4, ' ') +
      " | " +
      row.log_fc.to_string().pad_start(5, ' ') +
      " | " +
      row.ave_exp.to_string().pad_start(7, ' ') +
      " | " +
      row.t.to_string().pad_start(7, ' ') +
      " | " +
      row.p_value.to_string().pad_start(7, ' ') +
      " | " +
      row.adj_p_value.to_string().pad_start(9, ' ') +
      " | " +
      row.b.to_string().pad_start(4, ' '),
    )
    i = i + 1
  }
  println("\n")

  println("9. Creating contrast matrix (makeContrasts)...")
  let contrast_list = [[0.0, 1.0]]
  let contrasts = @bio.make_contrasts(contrast_list, design)
  println(
    "   Contrast matrix: " +
    contrasts.rows.to_string() +
    " x " +
    contrasts.cols.to_string() +
    "\n",
  )

  println("10. Fitting contrasts (contrasts.fit)...")
  let contrast_result = @bio.contrasts_fit(eb_result, contrasts)
  println("   Contrast fit complete")
  println(
    "   Contrast coefficients: " +
    contrast_result.coefficients.rows.to_string() +
    " x " +
    contrast_result.coefficients.cols.to_string() +
    "\n",
  )

  println("11. Summary statistics:")
  let mut sig_count = 0
  i = 0
  while i < eb_result.p_value.rows {
    if eb_result.p_value.get(i, 1) < 0.05 {
      sig_count = sig_count + 1
    }
    i = i + 1
  }
  println("   Total genes: " + eb_result.p_value.rows.to_string())
  println("   Significant genes (p < 0.05): " + sig_count.to_string())
  println("   Expected DE genes: 10 (first 10 genes have treatment effect)\n")

  // RPKM normalization
  println("12. RPKM normalization...")
  let gene_lengths = Array::make(100, 1000.0)
  let rpkm_result = @bio.limma_rpkm(counts, gene_lengths)
  println("   RPKM matrix: " + rpkm_result.rows.to_string() + " x " + rpkm_result.cols.to_string())

  // CPM normalization
  println("13. CPM normalization...")
  let cpm_result = @bio.limma_cpm(counts)
  println("   log2 CPM matrix: " + cpm_result.rows.to_string() + " x " + cpm_result.cols.to_string())

  // Quantile normalization
  println("14. Quantile normalization...")
  let quant_norm = @bio.limma_quantile_normalize(log_cpm)
  println("   Quantile normalized matrix: " + quant_norm.rows.to_string() + " x " + quant_norm.cols.to_string())

  // ComBat batch correction
  println("15. ComBat batch correction...")
  let combat_result = @bio.limma_combat(log_cpm, batch, group = group)
  println("   ComBat corrected matrix: " + combat_result.rows.to_string() + " x " + combat_result.cols.to_string())

  println("\n=== limma Analysis Complete ===")
}
