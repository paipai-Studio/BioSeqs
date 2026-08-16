///|
fn format_doubles(values : Array[Double]) -> String {
  let mut result = "["
  let mut index = 0
  while index < values.length() {
    if index > 0 {
      result = result + ", "
    }
    result = result + values[index].to_string()
    index = index + 1
  }
  result + "]"
}

///|
fn main {
  println("=== Bioconductor decontX Demo ===")
  let (counts, genes, cells, clusters, background) = @src.decontx_example_data()

  println("\n1. Infer cluster-aware ambient RNA contamination")
  let result = @src.decontx(counts, genes, cells, clusters) catch {
    _ => abort("failed to run cluster-aware decontX")
  }
  println("   " + result.summary())
  println("   per-cell contamination: " + format_doubles(result.contamination))

  println("\n2. Inspect cross-population marker correction")
  println(
    "   B marker in T1: observed=" +
    counts[2][0].to_string() +
    ", native=" +
    result.corrected_counts[2][0].to_string() +
    ", contaminant=" +
    result.contamination_counts[2][0].to_string(),
  )
  println(
    "   T marker in B1: observed=" +
    counts[0][4].to_string() +
    ", native=" +
    result.corrected_counts[0][4].to_string() +
    ", contaminant=" +
    result.contamination_counts[0][4].to_string(),
  )

  println("\n3. Summarize contamination by cluster")
  for estimate in result.cluster_estimates() {
    println(
      "   " +
      estimate.cluster +
      ": cells=" +
      estimate.n_cells.to_string() +
      ", mean contamination=" +
      estimate.mean_contamination.to_string() +
      ", contaminant counts=" +
      estimate.contaminant_counts.to_string(),
    )
  }
  let most_contaminated = result.most_contaminated_cells(2)
  for estimate in most_contaminated {
    println(
      "   high-contamination cell " +
      estimate.cell_name +
      ": " +
      estimate.contamination.to_string(),
    )
  }

  println("\n4. Use empty droplets as an explicit ambient profile")
  let background_result = @src.decontx(
    counts,
    genes,
    cells,
    clusters,
    background~,
  ) catch {
    _ => abort("failed to run background-aware decontX")
  }
  println("   " + background_result.summary())
  println(
    "   shared ambient profile: " +
    format_doubles(background_result.contaminant_profiles[0]),
  )

  println("\n5. Add corrected counts and diagnostics to an SCE copy")
  let sce = @src.SingleCellExperiment::new(counts, genes, cells)
  sce.col_data["cluster"] = clusters
  let output = @src.decontx_sce(sce, background~) catch {
    _ => abort("failed to integrate decontX with SingleCellExperiment")
  }
  println(
    "   corrected assay dimensions: " +
    @src.sce_get_assay(output.experiment, "decontXcounts").length().to_string() +
    " genes x " +
    output.result.n_cells().to_string() +
    " cells",
  )
  println(
    "   contamination metadata entries: " +
    @src.sce_get_col_data(output.experiment, "decontX_contamination")
    .length()
    .to_string(),
  )
  println(
    "   source object unchanged: " +
    (@src.sce_get_assay(sce, "decontXcounts").length() == 0).to_string(),
  )

  println("\n=== Demo Complete ===")
}
