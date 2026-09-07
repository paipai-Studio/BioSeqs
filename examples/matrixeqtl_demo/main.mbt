///|
/// MatrixEQTL Demo - Fast eQTL Analysis for Large-Scale Data.

///|
fn main {
  println("=== MatrixEQTL: Fast eQTL Analysis Demo ===\n")

  // 1. Using sample data
  println("1. Loading sample eQTL dataset...")
  let input = @src.sample_matrixeqtl_input()
  println("   SNPs: " + input.snp_names.length().to_string())
  println("   Genes: " + input.gene_names.length().to_string())
  println("   Samples: " + input.sample_names.length().to_string())
  println(
    "   SNP names: " + input.snp_names.iter().map(fn (s) { s }).collect().join(", "),
  )
  println(
    "   Gene names: " +
    input.gene_names.iter().map(fn (s) { s }).collect().join(", ") +
    "\n",
  )

  // 2. Run eQTL analysis with default threshold
  println("2. Running eQTL analysis (p-value threshold = 0.05)...")
  let output = @src.matrixeqtl_run(input)
  println("   Total tests performed: " + output.n_tests.to_string())
  println("   Significant associations: " + output.results.length().to_string() + "\n")

  // 3. Print top associations
  println("3. Top eQTL associations (sorted by p-value):")
  println("   SNP      | Gene   | Beta   | t-stat | P-value | FDR")
  println("   ---------|--------|--------|--------|---------|--------")
  let mut i = 0
  while i < output.results.length() {
    let r = output.results[i]
    println(
      "   " +
      r.snp_name.pad_start(8, ' ') +
      " | " +
      r.gene_name.pad_start(6, ' ') +
      " | " +
      r.beta.to_string().pad_start(6, ' ') +
      " | " +
      r.t_stat.to_string().pad_start(6, ' ') +
      " | " +
      r.p_value.to_string().pad_start(7, ' ') +
      " | " +
      r.fdr.to_string().pad_start(6, ' '),
    )
    i = i + 1
  }
  println("")

  // 4. Create custom SNP and gene expression matrices
  println("4. Creating custom SNP and gene expression matrices...")
  let snp_names = ["rs001", "rs002", "rs003"]
  let gene_names = ["TP53", "EGFR", "MYC", "BRCA1"]
  let sample_names = [
    "Ind1", "Ind2", "Ind3", "Ind4", "Ind5", "Ind6", "Ind7", "Ind8", "Ind9", "Ind10",
  ]

  // SNP genotype data (dosage: 0, 1, 2)
  let snps : Array[Array[Double]] = [
    // rs001: correlated with TP53
    [0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0],
    // rs002: correlated with EGFR
    [2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0],
    // rs003: no strong correlation
    [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
  ]

  // Gene expression data
  let gene_expression : Array[Array[Double]] = [
    // TP53: correlated with rs001 (expression increases with dosage)
    [4.0, 4.2, 6.1, 5.8, 6.3, 8.0, 7.9, 8.1, 6.0, 4.1],
    // EGFR: correlated with rs002 (expression increases with dosage)
    [12.0, 11.8, 9.5, 9.8, 7.0, 7.2, 9.6, 9.4, 11.9, 12.1],
    // MYC: no strong correlation with any SNP
    [5.0, 5.2, 4.8, 5.1, 5.3, 4.9, 5.0, 5.1, 4.8, 5.2],
    // BRCA1: weakly correlated with rs001
    [10.0, 10.1, 11.0, 10.8, 11.2, 12.0, 11.9, 12.1, 11.0, 10.2],
  ]

  let custom_input = @src.MatrixEQTLInput::new(
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  )
  println("   Custom dataset: 3 SNPs x 4 Genes x 10 Samples\n")

  // 5. Run analysis with stricter threshold
  println("5. Running eQTL analysis on custom data (p-value threshold = 0.01)...")
  let custom_output = @src.matrixeqtl_run(custom_input, pv_threshold=0.01)
  println("   Total tests performed: " + custom_output.n_tests.to_string())
  println("   Significant associations (p < 0.01): " + custom_output.results.length().to_string() + "\n")

  // 6. Print custom results
  println("6. Significant eQTL associations from custom data:")
  if custom_output.results.length() == 0 {
    println("   No associations passed the threshold.")
  } else {
    println("   SNP    | Gene  | Beta   | t-stat | P-value | FDR")
    println("   -------|-------|--------|--------|---------|--------")
    i = 0
    while i < custom_output.results.length() {
      let r = custom_output.results[i]
      println(
        "   " +
        r.snp_name.pad_start(6, ' ') +
        " | " +
        r.gene_name.pad_start(5, ' ') +
        " | " +
        r.beta.to_string().pad_start(6, ' ') +
        " | " +
        r.t_stat.to_string().pad_start(6, ' ') +
        " | " +
        r.p_value.to_string().pad_start(7, ' ') +
        " | " +
        r.fdr.to_string().pad_start(6, ' '),
      )
      i = i + 1
    }
  }
  println("")

  // 7. Demonstrate BH FDR correction directly
  println("7. Benjamini-Hochberg FDR correction example...")
  let raw_p_values = [0.001, 0.01, 0.03, 0.05, 0.12, 0.45, 0.67, 0.89]
  let fdr_adjusted = @src.matrixeqtl_bh_fdr(raw_p_values)
  println("   Raw p-values:  " + raw_p_values.iter().map(fn (v) { v.to_string() }).collect().join(", "))
  println(
    "   FDR adjusted:  " +
    fdr_adjusted.iter().map(fn (v) { v.to_string() }).collect().join(", ") +
    "\n",
  )

  // 8. Demonstrate linear regression directly
  println("8. Direct linear regression example...")
  let snp_dosage = [0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1.0]
  let expr_values = [5.0, 6.5, 6.2, 8.0, 4.8, 6.3, 7.8, 6.4]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(
    snp_dosage,
    expr_values,
  )
  println("   SNP dosage:  " + snp_dosage.iter().map(fn (v) { v.to_string() }).collect().join(", "))
  println("   Expression:  " + expr_values.iter().map(fn (v) { v.to_string() }).collect().join(", "))
  println("   Beta (effect size): " + beta.to_string())
  println("   t-statistic: " + t_stat.to_string())
  println("   p-value: " + p_value.to_string())

  println("\n=== MatrixEQTL Analysis Complete ===")
}
