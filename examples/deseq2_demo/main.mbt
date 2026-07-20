///|
/// DESeq2 Demo - Differential Expression Analysis for RNA-seq data.
/// Inspired by Bioconductor DESeq2 package.

///|
fn main {
  println("=== DESeq2 Differential Expression Analysis Demo ===")
  println("Inspired by Bioconductor DESeq2 package")
  println("")
  
  println("1. Creating sample dataset...")
  let dds = @bio.sample_deseq_dataset()
  println("   Genes: \{dds.row_names.length().to_string()}")
  println("   Samples: \{dds.col_names.length().to_string()}")
  println("   Group: 3 Control, 3 Treatment")
  println("")
  
  println("2. Estimating size factors (median ratio method)...")
  let dds_sf = @bio.estimate_size_factors(dds)
  let mut i = 0
  while i < dds_sf.size_factors.length() {
    println("   \{dds.col_names[i]}: \{dds_sf.size_factors[i].to_string()}")
    i = i + 1
  }
  println("")
  
  println("3. Estimating dispersions (parametric fit)...")
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let mut mean_disp = 0.0
  for disp in dds_disp.dispersions {
    mean_disp = mean_disp + disp
  }
  mean_disp = mean_disp / dds_disp.dispersions.length().to_double()
  println("   Mean dispersion: \{mean_disp.to_string()}")
  println("")
  
  println("4. Fitting Negative Binomial GLM and Wald test...")
  let res = @bio.deseq(dds_disp)
  println("   Differential expression results computed")
  println("")
  
  println("5. Differential expression testing results:")
  let mut sig_count = 0
  let mut up_count = 0
  let mut down_count = 0
  for i = 0; i < res.padj.length(); i = i + 1 {
    if res.padj[i] < 0.05 {
      sig_count = sig_count + 1
      if res.log2_fold_change[i] > 0.0 {
        up_count = up_count + 1
      } else {
        down_count = down_count + 1
      }
    }
  }
  println("   Total significant genes (padj < 0.05): \{sig_count.to_string()}")
  println("   Up-regulated: \{up_count.to_string()}")
  println("   Down-regulated: \{down_count.to_string()}")
  println("")
  
  println("6. LFC shrinkage (apeglm-like method)...")
  let res_shrunk = @bio.lfc_shrink(dds_disp, res)
  let mut sig_count_shrunk = 0
  for i = 0; i < res_shrunk.padj.length(); i = i + 1 {
    if res_shrunk.padj[i] < 0.05 && res_shrunk.log2_fold_change[i].abs() > 0.5 {
      sig_count_shrunk = sig_count_shrunk + 1
    }
  }
  println("   Significant after shrinkage (padj < 0.05, |lfc| > 0.5): \{sig_count_shrunk.to_string()}")
  println("")
  
  println("7. Top differentially expressed genes:")
  let top = @bio.top_genes(res, n=5)
  println("   Gene      log2FC    padj")
  println("   ------------------------")
  for t in top {
    let (gene, lfc, padj) = t
    println("   \{gene}  \{lfc.to_string()}  \{padj.to_string()}")
  }
  println("")
  
  println("8. Significant genes with log2FC > 1:")
  let sig_genes = @bio.significant_genes(res, alpha=0.05, lfc_threshold=1.0)
  println("   Significant genes with |log2FC| > 1: \{sig_genes.length().to_string()}")
  if sig_genes.length() > 0 {
    println("   Genes: \{sig_genes.join(", ")}")
  }
  println("")
  
  println("9. Results summary:")
  println("   - BaseMean: average normalized count")
  println("   - log2FC: log2 fold change (Treatment vs Control)")
  println("   - lfcSE: standard error of log2FC")
  println("   - stat: Wald test statistic")
  println("   - pvalue: raw p-value")
  println("   - padj: Benjamini-Hochberg adjusted p-value")
  println("")
  
  println("=== DESeq2 demo completed ===")
}
