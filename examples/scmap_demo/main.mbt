///|
/// scmap module demo - single-cell cell type annotation by reference projection.
///
/// Run: moon run examples/scmap_demo
fn main {
  println("=== scmap Module Demo ===")
  println("Single-cell RNA-seq cell type annotation")
  println("Reference: Kiselev VY et al. (2018) Nature Methods 15:359-362")
  println("")

  // ------------------------------------------------------------------
  // 1. Load reference atlas and query cells
  // ------------------------------------------------------------------
  println("=== 1. Load reference atlas and query data ===")
  let reference = @bio.scmap_sample_reference()
  println("Reference: \{reference.gene_names().length()} genes, \{reference.cell_names().length()} cells")
  println("Cell types: \{reference.unique_cell_types()}")
  let query = @bio.scmap_sample_query()
  println("Query: \{query.gene_names().length()} genes, \{query.cell_names().length()} cells")
  println("Query cells: \{query.cell_names()}")
  println("")

  // ------------------------------------------------------------------
  // 2. scmap-cluster: centroid-based correlation
  // ------------------------------------------------------------------
  println("=== 2. scmap-cluster (centroid correlation, threshold=0.7) ===")
  let cluster_assignments = @bio.scmap_cluster(reference, query, 0.7)
  for a in cluster_assignments {
    println(
      "  \{a.cell_name()} -> \{a.assigned_type()} (best_corr=\{a.best_correlation()})",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 3. scmap-cluster with stricter threshold
  // ------------------------------------------------------------------
  println("=== 3. scmap-cluster with strict threshold (0.95) ===")
  let strict_assignments = @bio.scmap_cluster(reference, query, 0.95)
  for a in strict_assignments {
    println(
      "  \{a.cell_name()} -> \{a.assigned_type()} (best_corr=\{a.best_correlation()})",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 4. scmap-cell: nearest-neighbour voting
  // ------------------------------------------------------------------
  println("=== 4. scmap-cell (k=5 neighbours, fraction=0.5) ===")
  let cell_assignments = @bio.scmap_cell(reference, query, 5, 0.5)
  for a in cell_assignments {
    println(
      "  \{a.cell_name()} -> \{a.assigned_type()} (best_corr=\{a.best_correlation()}, method=\{a.method()})",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Summary of cluster-based assignments
  // ------------------------------------------------------------------
  println("=== 5. Summary of cluster-based assignments ===")
  let summary = @bio.scmap_summary(cluster_assignments)
  for entry in summary {
    let (ctype, count) = entry
    println("  \{ctype}: \{count} cell(s)")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Full assignments table
  // ------------------------------------------------------------------
  println("=== 6. Full assignments table (cluster method) ===")
  println(@bio.scmap_assignments_to_string(cluster_assignments))

  println("=== Demo complete ===")
}
