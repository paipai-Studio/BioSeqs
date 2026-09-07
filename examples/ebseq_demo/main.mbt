///|
/// EBSeq Demo - Empirical Bayesian approach for differential expression analysis.
/// Inspired by Bioconductor EBSeq package.

///|
fn main {
  println("=== EBSeq Differential Expression Analysis Demo ===")
  println("Empirical Bayesian approach for identifying DE genes")
  println("")

  println("1. Creating sample dataset...")
  let data = @src.ebseq_sample_dataset()
  println("   Genes: \{data.gene_names.length().to_string()}")
  println("   Samples: \{data.sample_names.length().to_string()}")
  println("   Conditions: 3 replicates per condition (0=Control, 1=Treatment)")
  println("   DE genes: GeneA-GeneE (up), GeneF-GeneJ (down), GeneK-GeneT (equal)")
  println("")

  println("2. Estimating size factors (median-of-ratios method)...")
  let data_sf = @src.ebseq_estimate_size_factors(data)
  for i = 0; i < data_sf.sample_names.length(); i = i + 1 {
    println("   \{data_sf.sample_names[i]}: \{data_sf.size_factors[i].to_string()}")
  }
  println("")

  println("3. Estimating dispersions (method of moments with shrinkage)...")
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let mut mean_disp = 0.0
  for d in data_disp.dispersions {
    mean_disp = mean_disp + d
  }
  mean_disp = mean_disp / data_disp.dispersions.length().to_double()
  println("   Mean dispersion: \{mean_disp.to_string()}")
  println("")

  println("4. Running EBSeq DE test (PPDE threshold = 0.95)...")
  let res = @src.ebseq_test(data_disp, threshold=0.95)
  println("   Results computed for \{res.gene_names.length().to_string()} genes")
  println("")

  println("5. Differential expression summary:")
  let de_count = @src.ebseq_count_de(res)
  let mut up_count = 0
  let mut down_count = 0
  for i = 0; i < res.flag_de.length(); i = i + 1 {
    if res.flag_de[i] {
      if res.log2_fc[i] > 0.0 {
        up_count = up_count + 1
      } else {
        down_count = down_count + 1
      }
    }
  }
  println("   DE genes (PPDE >= 0.95): \{de_count.to_string()}")
  println("   Up-regulated in condition 2: \{up_count.to_string()}")
  println("   Down-regulated in condition 2: \{down_count.to_string()}")
  println("")

  println("6. Top DE genes by PPDE:")
  let top = @src.ebseq_top_genes(res, n=10)
  println("   Gene      PPDE      log2FC")
  println("   --------------------------")
  for t in top {
    let (gene, ppde, lfc) = t
    println("   \{gene}    \{ppde.to_string()}  \{lfc.to_string()}")
  }
  println("")

  println("7. DE gene names:")
  let de_genes = @src.ebseq_get_de_genes(res)
  if de_genes.length() > 0 {
    println("   \{de_genes.join(", ")}")
  } else {
    println("   No genes passed the PPDE >= 0.95 threshold")
  }
  println("")

  println("8. Creating a custom two-condition dataset...")
  let counts = [
    [100, 110, 105, 300, 310, 290],
    [200, 210, 205, 200, 195, 210],
    [300, 290, 310, 100, 110, 90],
  ]
  let gene_names = ["UpGene", "EqGene", "DownGene"]
  let sample_names = ["C1_R1", "C1_R2", "C1_R3", "C2_R1", "C2_R2", "C2_R3"]
  let conditions = [0, 0, 0, 1, 1, 1]
  let custom = @src.EBSeqDataSet::new(counts, gene_names, sample_names, conditions)
  let custom_sf = @src.ebseq_estimate_size_factors(custom)
  let custom_disp = @src.ebseq_estimate_dispersion(custom_sf)
  let custom_res = @src.ebseq_test(custom_disp, threshold=0.5)
  println("   Custom dataset: 3 genes, 6 samples")
  for i = 0; i < custom_res.gene_names.length(); i = i + 1 {
    println(
      "   \{custom_res.gene_names[i]}: PPDE=\{custom_res.ppde[i].to_string()}, log2FC=\{custom_res.log2_fc[i].to_string()}",
    )
  }
  println("")

  println("=== EBSeq demo completed ===")
}
