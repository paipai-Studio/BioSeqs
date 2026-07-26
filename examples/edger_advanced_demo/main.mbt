///|
/// edgeR Advanced Demo - Quasi-Likelihood F-test and Gene Set Tests
///
/// This demo shows how to:
/// 1. Create a DGEList from count data
/// 2. Normalize and fit quasi-likelihood GLM
/// 3. Perform QL F-test for differential expression
/// 4. Identify significant genes
/// 5. Run gene set tests (camera and roast)

fn main {
  println("=== edgeR Advanced: QL F-test and Gene Set Tests Demo ===\n")

  // Step 1: Create a DGEList with RNA-seq count data
  println("Step 1: Creating DGEList...")
  let genes = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE",
              "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ"]
  let groups = ["Ctrl", "Ctrl", "Ctrl", "Treat", "Treat", "Treat"]

  let counts = [
    [100, 105, 98, 200, 210, 195],   // GeneA: upregulated
    [50, 52, 48, 25, 28, 22],        // GeneB: downregulated
    [80, 82, 78, 82, 80, 79],        // GeneC: no change
    [300, 310, 295, 600, 590, 610],  // GeneD: upregulated
    [20, 22, 18, 21, 19, 23],        // GeneE: no change
    [150, 145, 148, 140, 142, 138],  // GeneF: no change
    [40, 42, 38, 80, 82, 78],        // GeneG: upregulated
    [90, 88, 92, 45, 48, 42],        // GeneH: downregulated
    [500, 510, 495, 505, 498, 512],  // GeneI: no change
    [60, 58, 62, 120, 115, 118],     // GeneJ: upregulated
  ]

  let dge = @bio.dge_list(counts, groups, genes)
  println("  Genes: \{dge.counts.length()}, Samples: \{dge.group.length()}")

  // Step 2: Normalize
  println("\nStep 2: Normalizing with TMM...")
  let dge_norm = @bio.calc_norm_factors(dge)
  println("  Norm factors: \{dge_norm.norm_factors}")

  // Step 3: Fit quasi-likelihood GLM
  println("\nStep 3: Fitting quasi-likelihood GLM...")
  let qlf_data = @bio.glm_qlf_fit(dge_norm)
  println("  Trend dispersion: \{qlf_data.trend_dispersion}")
  println("  QL dispersions: \{qlf_data.ql_dispersion}")

  // Step 4: Perform QL F-test
  println("\nStep 4: Performing QL F-test for differential expression...")
  let qlf_result = @bio.glm_qlf_test(qlf_data)
  println("  Degrees of freedom: numerator=\{qlf_result.df_num}, denominator=\{qlf_result.df_denom}")

  println("\n  Results:")
  println("  Gene    LogFC    F-stat   P-value    FDR")
  for i = 0; i < qlf_result.genes.length(); i = i + 1 {
    println("  \{qlf_result.genes[i]}  \{qlf_result.log_fc[i].to_string()}  \{qlf_result.f_stat[i].to_string()}  \{qlf_result.p_value[i].to_string()}  \{qlf_result.fdr[i].to_string()}")
  }

  // Step 5: Get significant genes
  println("\nStep 5: Significant genes (FDR < 0.05):")
  let sig_genes = @bio.qlf_significant_genes(qlf_result, alpha=0.05)
  for g in sig_genes {
    println("  - \{g}")
  }

  // Step 6: Top genes
  println("\nStep 6: Top 5 genes by FDR:")
  let top = @bio.qlf_top_genes(qlf_result, n=5)
  for entry in top {
    let (gene, lfc, fstat, fdr) = entry
    println("  \{gene}: LogFC=\{lfc.to_string()}, F=\{fstat.to_string()}, FDR=\{fdr.to_string()}")
  }

  // Step 7: Gene set tests
  println("\nStep 7: Gene set tests with camera and roast...")
  let gene_set = ["GeneA", "GeneD", "GeneG", "GeneJ"]  // Upregulated genes

  let (cam_stat, cam_p, cam_dir) = @bio.camera(qlf_result, gene_set)
  println("  Camera test for upregulated gene set:")
  println("    Statistic: \{cam_stat.to_string()}")
  println("    P-value: \{cam_p.to_string()}")
  println("    Direction: \{cam_dir}")

  let (roast_stat, roast_p, roast_dir) = @bio.roast(qlf_result, gene_set)
  println("\n  Roast test for upregulated gene set:")
  println("    Statistic: \{roast_stat.to_string()}")
  println("    P-value: \{roast_p.to_string()}")
  println("    Direction: \{roast_dir}")

  // Test with downregulated gene set
  let down_set = ["GeneB", "GeneH"]
  let (cam2_stat, cam2_p, cam2_dir) = @bio.camera(qlf_result, down_set)
  println("\n  Camera test for downregulated gene set:")
  println("    Statistic: \{cam2_stat.to_string()}")
  println("    P-value: \{cam2_p.to_string()}")
  println("    Direction: \{cam2_dir}")

  println("\n=== edgeR QL F-test Demo Complete ===")
}
