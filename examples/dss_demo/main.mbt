///|
/// DSS Demo - Differential Expression and Methylation Analysis.
/// Inspired by Bioconductor DSS package.

fn main {
  println("=== DSS Demo ===")

  // 1. Load example count data
  println("")
  println("1. Load example count data:")
  let data = @bio.dss_create_example_counts()
  println("   Genes: " + data.dss_n_genes().to_string())
  println("   Samples: " + data.dss_n_samples().to_string())

  // 2. Run dispersion estimation
  println("")
  println("2. Estimate dispersion:")
  let disp = @bio.dss_estimate_dispersion(data)
  let dispersions = disp.dss_get_dispersions()
  let shrunken = disp.dss_get_shrunken_dispersions()
  let mut i = 0
  while i < dispersions.length() {
    println("   Gene " + (i + 1).to_string() +
      ": dispersion=" + dispersions[i].to_string() +
      ", shrunken=" + shrunken[i].to_string())
    i = i + 1
  }

  // 3. Run DE analysis
  println("")
  println("3. Differential expression analysis (Wald test):")
  let de_result = @bio.dss_deseq(data, disp)
  println("   Total genes: " + de_result.dss_n_genes().to_string())
  println("   Significant genes (FDR < 0.05): " + de_result.dss_n_significant().to_string())
  let mut j = 0
  while j < de_result.dss_n_genes() {
    println("   " + de_result.dss_get_gene_id(j) +
      ": wald=" + de_result.wald_stat[j].to_string() +
      ", pval=" + de_result.dss_get_pvalue(j).to_string() +
      ", adj_pval=" + de_result.dss_get_adj_pvalue(j).to_string())
    j = j + 1
  }

  // 4. Load methylation data
  println("")
  println("4. Load example methylation data:")
  let methyl = @bio.dss_create_example_methylation()
  println("   CpG sites: " + methyl.dss_n_cpgs().to_string())
  println("   Samples: " + methyl.dss_n_samples().to_string())

  // 5. Run DM analysis
  println("")
  println("5. Differential methylation analysis:")
  let dm_result = @bio.dss_dml(methyl)
  println("   Total CpG sites: " + dm_result.dss_n_sites().to_string())
  println("   Significant sites (FDR < 0.05): " + dm_result.dss_n_significant().to_string())
  let mut k = 0
  while k < dm_result.dss_n_sites() {
    println("   " + dm_result.cpg_ids[k] +
      ": delta_meth=" + dm_result.dss_get_delta_meth(k).to_string() +
      ", stat=" + dm_result.stats[k].to_string() +
      ", adj_pval=" + dm_result.adj_p_values[k].to_string())
    k = k + 1
  }

  // 6. Run DMR detection
  println("")
  println("6. Detect differentially methylated regions (min_cpgs=2, min_delta=0.1):")
  let dmr_result = @bio.dss_dmr(methyl, dm_result, 2, 0.1)
  println("   DMR count: " + dmr_result.dss_n_dmrs().to_string())
  let mut m = 0
  while m < dmr_result.dss_n_dmrs() {
    let (id, chrom, start, end, n_cpgs, mean_delta) = dmr_result.dss_get_dmr(m)
    println("   " + id +
      ": " + chrom + ":" + start.to_string() + "-" + end.to_string() +
      ", CpGs=" + n_cpgs.to_string() +
      ", mean_delta=" + mean_delta.to_string())
    m = m + 1
  }

  println("")
  println("=== DSS Demo Completed ===")
}