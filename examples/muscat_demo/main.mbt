///|
/// muscat module example - Differential state analysis for single-cell RNA-seq.
///
/// Run: moon run examples/muscat_demo
fn main {
  println("=== muscat Module Demo ===")
  println("Differential State Analysis for Single-Cell RNA-seq")
  println("")

  // Generate sample single-cell data
  println("=== Sample Data ===")
  let cells = @bio.muscat_sample_data()
  println("Total single cells: \{cells.length()}")
  let mut n_ctrl = 0
  let mut n_stim = 0
  let mut i = 0
  while i < cells.length() {
    if cells[i].group_id == "ctrl" {
      n_ctrl = n_ctrl + 1
    } else {
      n_stim = n_stim + 1
    }
    i = i + 1
  }
  println("  Control cells: \{n_ctrl}")
  println("  Stimulation cells: \{n_stim}")
  println("")

  // Aggregate cells to pseudo-bulk
  println("=== Pseudo-bulk Aggregation (Sum) ===")
  let pb = @bio.aggregate_cells(cells, method=@bio.aggregation_sum())
  println("Pseudo-bulk samples: \{pb.length()}")
  let mut j = 0
  while j < pb.length() {
    println("  Sample: \{pb[j].sample_id}, Cluster: \{pb[j].cluster_id}, Group: \{pb[j].group_id}, Cells: \{pb[j].n_cells}")
    j = j + 1
  }
  println("")

  // QC summary
  println("=== QC Summary ===")
  let qc = @bio.compute_qc(pb)
  let mut k = 0
  while k < qc.length() {
    println("  \{qc[k].sample_id}/\{qc[k].cluster_id}: \{qc[k].n_cells} cells, \{qc[k].n_genes} genes, \{qc[k].total_counts} total counts")
    k = k + 1
  }
  println("")

  // Run DS analysis
  println("=== Differential State Analysis ===")
  let results = @bio.run_ds_analysis(pb, method=@bio.ds_method_edger(), fdr_threshold=0.1, log2fc_threshold=0.5)
  println(results.summary())
  println("")

  // Top differentially expressed genes
  println("=== Top 5 DE Genes ===")
  let top = results.get_top_genes(5)
  let mut m = 0
  while m < top.length() {
    let g = top[m]
    println("  \{g.gene_id} [\{g.cluster_id}]: log2FC=\{g.log2fc}, p_adj=\{g.p_adj}, sig=\{g.significant}")
    m = m + 1
  }
  println("")

  // Significant results
  println("=== Significant Results ===")
  let sig = results.get_significant()
  println("Total significant genes: \{sig.length()}")
  println("")

  // Aggregation methods comparison
  println("=== Aggregation Methods ===")
  let pb_mean = @bio.aggregate_cells(cells, method=@bio.aggregation_mean())
  let pb_med = @bio.aggregate_cells(cells, method=@bio.aggregation_median())
  println("Sum aggregation: \{pb.length()} samples")
  println("Mean aggregation: \{pb_mean.length()} samples")
  println("Median aggregation: \{pb_med.length()} samples")
}
