///|
/// DESeq2 Advanced Demo - Variance Stabilizing Transformation and PCA
///
/// This demo shows how to:
/// 1. Create a DESeqDataSet from count data
/// 2. Estimate size factors and dispersions
/// 3. Apply variance stabilizing transformation (VST)
/// 4. Perform PCA on transformed data for sample visualization
/// 5. Identify highly variable genes

fn main {
  println("=== DESeq2 Advanced: VST and PCA Demo ===\n")

  // Step 1: Create a DESeqDataSet with RNA-seq count data
  println("Step 1: Creating DESeqDataSet...")
  let genes = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE",
              "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ"]
  let samples = ["Ctrl1", "Ctrl2", "Treat1", "Treat2"]

  let counts = [
    [100, 105, 200, 210],   // GeneA: upregulated
    [50, 55, 25, 30],       // GeneB: downregulated
    [80, 85, 82, 78],       // GeneC: no change
    [300, 310, 600, 590],   // GeneD: upregulated
    [20, 22, 18, 21],       // GeneE: no change
    [150, 145, 140, 148],   // GeneF: no change
    [40, 42, 80, 85],       // GeneG: upregulated
    [90, 95, 45, 50],       // GeneH: downregulated
    [500, 510, 505, 498],   // GeneI: no change
    [60, 58, 120, 115],     // GeneJ: upregulated
  ]

  let design = [
    [1.0, 0.0],  // Ctrl1
    [1.0, 0.0],  // Ctrl2
    [1.0, 1.0],  // Treat1
    [1.0, 1.0],  // Treat2
  ]

  let dds = @bio.DESeqDataSet::new(counts, genes, samples, design)
  println("  Genes: \{dds.counts.length()}, Samples: \{dds.col_names.length()}")

  // Step 2: Estimate size factors and dispersions
  println("\nStep 2: Estimating size factors...")
  let dds_with_sf = @bio.estimate_size_factors(dds)
  println("  Size factors: \{dds_with_sf.size_factors}")

  println("\nStep 3: Estimating dispersions...")
  let dds_with_disp = @bio.estimate_dispersions(dds_with_sf)
  let mut disp_output = "  Dispersions (first 5): ["
  for i = 0; i < 5 && i < dds_with_disp.dispersions.length(); i = i + 1 {
    if i > 0 {
      disp_output = disp_output + ", "
    }
    disp_output = disp_output + dds_with_disp.dispersions[i].to_string()
  }
  disp_output = disp_output + "]"
  println(disp_output)

  // Step 3: Apply variance stabilizing transformation
  println("\nStep 4: Applying variance stabilizing transformation (VST)...")
  let vst_result = @bio.variance_stabilizing_transformation(dds_with_disp, blind=false)
  println("  Dispersion used: \{vst_result.dispersion}")
  println("  Transformed matrix dimensions: \{vst_result.transformed.length()} x \{vst_result.transformed[0].length()}")

  // Show some transformed values
  println("\n  Transformed values (first 3 genes, all samples):")
  for i = 0; i < 3; i = i + 1 {
    let row = vst_result.transformed[i]
    println("  \{vst_result.row_names[i]}: \{row[0].to_string()} \{row[1].to_string()} \{row[2].to_string()} \{row[3].to_string()}")
  }

  // Step 4: Fast VST with subset
  println("\nStep 5: Applying fast VST (subset=5)...")
  let fast_vst = @bio.vst(dds_with_sf, nsub=5)
  println("  Estimated dispersion: \{fast_vst.dispersion}")

  // Step 5: PCA on transformed data
  println("\nStep 6: Performing PCA on VST-transformed data...")
  let pca = @bio.plot_pca(vst_result, ntop=10)
  println("  Variance explained - PC1: \{pca.var_explained_pc1.to_string()}%, PC2: \{pca.var_explained_pc2.to_string()}%")

  println("\n  PCA coordinates:")
  for i = 0; i < pca.sample_names.length(); i = i + 1 {
    println("  \{pca.sample_names[i]}: PC1=\{pca.pc1[i].to_string()}, PC2=\{pca.pc2[i].to_string()}")
  }

  // Step 6: Identify highly variable genes
  println("\nStep 7: Gene variability (standard deviation) after VST:")
  let sds = @bio.gene_sd(vst_result)
  // Sort by SD descending (simple selection of top 3)
  let top_genes : Array[(String, Double)] = Array::new()
  for g in sds {
    top_genes.push(g)
  }
  // Simple sort by SD
  let n = top_genes.length()
  for i = 0; i < n - 1; i = i + 1 {
    for j = 0; j < n - 1 - i; j = j + 1 {
      if top_genes[j].1 < top_genes[j + 1].1 {
        let tmp = top_genes[j]
        top_genes[j] = top_genes[j + 1]
        top_genes[j + 1] = tmp
      }
    }
  }

  println("  Top 5 most variable genes:")
  for i = 0; i < 5 && i < top_genes.length(); i = i + 1 {
    println("  \{top_genes[i].0}: SD=\{top_genes[i].1.to_string()}")
  }

  println("\n=== VST Demo Complete ===")
}
