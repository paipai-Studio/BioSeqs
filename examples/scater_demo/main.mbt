///|
fn main {
  println("=== scater Demo ===")
  
  let n_cells = 10
  let n_genes = 20
  
  let counts : Array[Array[Int]] = Array::make(n_cells, Array::make(n_genes, 0))
  let cell_ids : Array[String] = Array::make(n_cells, "")
  let gene_names : Array[String] = Array::make(n_genes, "")
  
  let mut i = 0
  while i < n_cells {
    cell_ids[i] = "Cell" + i.to_string()
    let row : Array[Int] = Array::make(n_genes, 0)
    let mut j = 0
    while j < n_genes {
      row[j] = (i * 7 + j * 13) % 100
      j = j + 1
    }
    counts[i] = row
    i = i + 1
  }
  
  let mut g = 0
  while g < n_genes {
    if g < 2 {
      gene_names[g] = "mt-Gene" + g.to_string()
    } else if g < 4 {
      gene_names[g] = "Rpl" + g.to_string()
    } else {
      gene_names[g] = "Gene" + g.to_string()
    }
    g = g + 1
  }
  
  println("\n1. Calculate QC Metrics:")
  let cell_stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  println(cell_stats.to_string())
  
  println("\n2. QC metrics for first 5 cells:")
  let mut k = 0
  while k < 5 && k < cell_stats.cell_id.length() {
    println("  " + cell_stats.cell_id[k] + ":")
    println("    Total counts: " + cell_stats.total_counts[k].to_string())
    println("    Detected genes: " + cell_stats.detected_genes[k].to_string())
    println("    Percent mt: " + cell_stats.percent_mito[k].to_string() + "%")
    println("    Percent ribo: " + cell_stats.percent_ribo[k].to_string() + "%")
    k = k + 1
  }
  
  println("\n3. Feature QC:")
  let feature_stats = @src.scater_calculate_feature_qc(counts, gene_names)
  println("  Total genes: " + feature_stats.gene_id.length().to_string())
  println("  First gene mean expression: " + feature_stats.mean_expression[0].to_string())
  println("  First gene detection rate: " + feature_stats.percent_detected[0].to_string() + "%")
  
  println("\n4. Filter cells (default params):")
  let params = @src.QCFilterParams::new()
  let (kept_indices, keep_cells) = @src.scater_filter_cells(cell_stats, params)
  println("  Cells kept: " + kept_indices.length().to_string() + " / " + n_cells.to_string())
  
  println("\n5. Filter features (min 3 cells detected):")
  let (kept_genes, keep_genes) = @src.scater_filter_features(feature_stats, 3)
  println("  Genes kept: " + kept_genes.length().to_string() + " / " + n_genes.to_string())
  
  println("\n6. CPM normalization:")
  let cpm = @src.scater_normalize_cpm(counts)
  println("  First cell CPM sum: " + (cpm[0][0] + cpm[0][1] + cpm[0][2]).to_string() + "...")
  
  println("\n7. Log-CPM normalization:")
  let log_cpm = @src.scater_log_normalize_cpm(counts)
  println("  First gene log-CPM: " + log_cpm[0][0].to_string())
  
  println("\n8. Find highly variable genes (top 10):")
  let hvgs = @src.scater_find_hvgs(log_cpm, gene_names, 10)
  println("  Top HVGs:")
  let mut h = 0
  while h < hvgs.length() {
    println("    " + hvgs[h])
    h = h + 1
  }
  
  println("\n9. PCA (3 components):")
  let pca = @src.scater_run_pca(log_cpm, 3)
  println("  PCA result shape: " + pca.length().to_string() + " cells x " + pca[0].length().to_string() + " components")
  println("  First cell PC1: " + pca[0][0].to_string())
  println("  First cell PC2: " + pca[0][1].to_string())
  
  println("\n10. QC Report:")
  let report = @src.ScaterQCReport::new(cell_stats, feature_stats, params)
  println(report.to_string())
  
  println("\nscater demo completed!")
}
