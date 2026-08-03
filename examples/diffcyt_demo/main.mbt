///|
/// diffcyt module demo - Differential discovery in high-dimensional cytometry.
///
/// Demonstrates FlowSOM clustering, differential abundance (DA) testing,
/// and differential state (DS) testing on simulated CyTOF data.
///
/// Reference: Weber LM & Robinson MD (2016) diffcyt. Nat Methods 17:792.
/// Run: moon run examples/diffcyt_demo
fn main {
  println("=== diffcyt Module Demo ===")
  println("Differential discovery in high-dimensional cytometry")
  println("Reference: Bioconductor diffcyt")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let cells = @bio.diffcyt_sample_data()
  println("Total cells: " + cells.length().to_string())
  let samples = @bio.diffcyt_unique_samples(cells)
  let conditions = @bio.diffcyt_unique_conditions(cells)
  println("Samples: " + samples.length().to_string())
  println("Conditions: " + conditions.to_string())
  println("Markers: " + cells[0].n_markers().to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Cluster cells (FlowSOM)
  // ------------------------------------------------------------------
  println("=== 2. Cluster Cells (FlowSOM) ===")
  let n_clusters = 5
  let codebooks = @bio.diffcyt_cluster_cells(cells, n_clusters, 10)
  println("Number of clusters: " + n_clusters.to_string())
  // Count cells per cluster.
  let cluster_counts = Array::make(n_clusters, 0)
  for cell in cells {
    cluster_counts[cell.cluster_id()] = cluster_counts[cell.cluster_id()] + 1
  }
  for i in 0..<n_clusters {
    println("  Cluster " + i.to_string() + ": " + cluster_counts[i].to_string() + " cells")
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Compute cluster counts per sample
  // ------------------------------------------------------------------
  println("=== 3. Cluster Counts per Sample ===")
  let counts = @bio.diffcyt_calc_counts(cells, n_clusters, samples)
  println("Sample | " + "Cluster counts")
  for si in 0..<samples.length() {
    let mut row_str = "  " + samples[si] + " | "
    for ci in 0..<n_clusters {
      row_str = row_str + counts[si][ci].to_string() + " "
    }
    println(row_str)
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Build design matrix
  // ------------------------------------------------------------------
  println("=== 4. Design Matrix ===")
  let conditions_arr : Array[String] = Array::new()
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id(), cell.condition())
  }
  for sid in samples {
    conditions_arr.push(cond_map.get(sid).unwrap_or(""))
  }
  let design = @bio.diffcyt_create_design_matrix(samples, conditions_arr)
  println("Design matrix dimensions: " + design.length().to_string() + " x " + design[0].length().to_string())
  println("Columns: [intercept, condition_stim]")
  println("")

  // ------------------------------------------------------------------
  // 5. Run full pipeline (DA + DS testing)
  // ------------------------------------------------------------------
  println("=== 5. Run Full Pipeline ===")
  let summary = @bio.diffcyt_run(cells, n_clusters, 0.05)
  println("Cells: " + summary.n_cells.to_string())
  println("Clusters: " + summary.n_clusters.to_string())
  println("Markers: " + summary.n_markers.to_string())
  println("Samples: " + summary.n_samples.to_string())
  println("Conditions: " + summary.n_conditions.to_string())
  println("FDR threshold: " + summary.fdr_threshold.to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. DA results
  // ------------------------------------------------------------------
  println("=== 6. Differential Abundance (DA) Results ===")
  let top_da = @bio.diffcyt_top_clusters_da(summary.da_results, n_clusters)
  for r in top_da {
    let sig = if r.significant() { "*" } else { "" }
    println(
      "  Cluster " +
        r.cluster_id().to_string() +
        ": logFC=" +
        r.log_fc().to_string() +
        " p=" +
        r.p_value().to_string() +
        " FDR=" +
        r.fdr().to_string() +
        " " +
        sig,
    )
  }
  println("Significant DA clusters: " + summary.n_sig_da().to_string())
  println("")

  // ------------------------------------------------------------------
  // 7. DS results
  // ------------------------------------------------------------------
  println("=== 7. Differential State (DS) Results ===")
  let top_ds = @bio.diffcyt_top_table_ds(summary.ds_results, 5)
  for r in top_ds {
    let sig = if r.significant() { "*" } else { "" }
    println(
      "  Cluster " +
        r.cluster_id().to_string() +
        " " +
        r.marker_name() +
        ": logFC=" +
        r.log_fc().to_string() +
        " p=" +
        r.p_value().to_string() +
        " FDR=" +
        r.fdr().to_string() +
        " " +
        sig,
    )
  }
  println("Significant DS results: " + summary.n_sig_ds().to_string())
  println("")

  println("=== Demo complete ===")
}
