///|
/// SingleCellExperiment Demo.
fn main {
  println("=== SingleCellExperiment Demo ===")

  // 1. Create example SCE
  println("\n1. Create example SCE:")
  let sce = @src.sce_create_example()
  println(sce.to_string())
  println("  n_genes = " + @src.sce_get_n_genes(sce).to_string())
  println("  n_cells = " + @src.sce_get_n_cells(sce).to_string())

  // 2. Demonstrate Builder pattern
  println("\n2. Builder pattern:")
  let counts = @src.sce_create_example_counts()
  let built_sce = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("counts", counts)
    |> @src.SCEBuilder::add_row_data("gene_type", [
      "TF", "Kinase", "Receptor", "Enzyme", "Other",
    ])
    |> @src.SCEBuilder::add_col_data("batch", ["b1", "b1", "b2", "b2"])
    |> @src.SCEBuilder::add_reduced_dim("PCA", [
      [0.1, 0.2],
      [0.3, 0.4],
      [0.5, 0.6],
      [0.7, 0.8],
    ])
    |> @src.SCEBuilder::build
  println(
    "  Built SCE n_genes = " + @src.sce_get_n_genes(built_sce).to_string(),
  )
  println(
    "  Built SCE n_cells = " + @src.sce_get_n_cells(built_sce).to_string(),
  )

  // 3. Demonstrate assay access and setting
  println("\n3. Assay access and set:")
  let counts_assay = @src.sce_get_assay(sce, "counts")
  println(
    "  counts assay: " +
    counts_assay.length().to_string() +
    " genes x " +
    (if counts_assay.length() > 0 {
      counts_assay[0].length().to_string()
    } else {
      "0".to_string()
    }) +
    " cells",
  )
  let logcounts_assay = @src.sce_get_assay(sce, "logcounts")
  println("  logcounts assay length: " + logcounts_assay.length().to_string())
  let updated_sce = @src.sce_set_assay(sce, "custom", [[1.0, 2.0], [3.0, 4.0]])
  let custom = @src.sce_get_assay(updated_sce, "custom")
  println("  custom assay length: " + custom.length().to_string())

  // 4. Demonstrate counts and logcounts conversion
  println("\n4. Counts and logcounts conversion:")
  let raw_counts = @src.sce_get_counts(sce)
  let logcounts = @src.sce_get_logcounts(sce)
  println("  raw counts genes: " + raw_counts.length().to_string())
  println("  logcounts genes: " + logcounts.length().to_string())
  let computed_logcounts = @src.sce_compute_logcounts(raw_counts, [
    1.0, 1.0, 1.0, 1.0,
  ])
  let sce_with_log = @src.sce_set_logcounts(sce, computed_logcounts)
  println(
    "  set logcounts: " +
    @src.sce_get_logcounts(sce_with_log).length().to_string() +
    " genes",
  )

  // 5. Demonstrate PCA dimensionality reduction
  // Note: sce_compute_pca uses a square-matrix multiplication path; we pass a
  // square (4x4) logcounts matrix so the source implementation runs cleanly.
  println("\n5. PCA:")
  let square_logcounts : Array[Array[Double]] = [
    [100.0, 200.0, 50.0, 150.0],
    [120.0, 180.0, 60.0, 140.0],
    [80.0, 220.0, 40.0, 160.0],
    [90.0, 210.0, 45.0, 155.0],
  ]
  let pca_result = @src.sce_compute_pca(square_logcounts, 2)
  println(pca_result.to_string())

  // 6. Demonstrate t-SNE dimensionality reduction
  println("\n6. t-SNE:")
  let tsne_result = @src.sce_compute_tsne(square_logcounts, 3)
  println(tsne_result.to_string())

  // 7. Demonstrate UMAP dimensionality reduction
  println("\n7. UMAP:")
  let umap_result = @src.sce_compute_umap(square_logcounts, 2)
  println(umap_result.to_string())

  // 8. Demonstrate size factors computation
  println("\n8. Size factors:")
  let lib_sf = @src.sce_compute_library_size_factors(raw_counts)
  println("  library size factors: " + lib_sf.length().to_string())
  let mut i = 0
  while i < lib_sf.length() {
    println("    cell " + (i + 1).to_string() + ": " + lib_sf[i].to_string())
    i = i + 1
  }
  let deconv_sf = @src.sce_compute_deconv_size_factors(raw_counts)
  println("  deconv size factors: " + deconv_sf.length().to_string())
  let sf_result = @src.sce_compute_size_factors(raw_counts, "library")
  println("  size factor method: " + sf_result.method)

  // 9. Demonstrate normalization
  println("\n9. Normalize:")
  let norm_library = @src.sce_normalize(sce, "library")
  let norm_logcounts = @src.sce_get_logcounts(norm_library)
  println(
    "  library-normalized logcounts: " +
    norm_logcounts.length().to_string() +
    " genes",
  )
  let norm_deconv = @src.sce_normalize(sce, "deconv")
  println(
    "  deconv-normalized logcounts: " +
    @src.sce_get_logcounts(norm_deconv).length().to_string() +
    " genes",
  )

  // 10. Demonstrate gene filtering
  println("\n10. Filter genes (min_cells=2):")
  let filtered_genes = @src.sce_filter_genes(sce, 2)
  println(
    "  genes before: " +
    @src.sce_get_n_genes(sce).to_string() +
    ", after: " +
    @src.sce_get_n_genes(filtered_genes).to_string(),
  )

  // 11. Demonstrate cell filtering
  println("\n11. Filter cells (min_genes=2):")
  let filtered_cells = @src.sce_filter_cells(sce, 2)
  println(
    "  cells before: " +
    @src.sce_get_n_cells(sce).to_string() +
    ", after: " +
    @src.sce_get_n_cells(filtered_cells).to_string(),
  )

  // 12. Demonstrate subset operations
  println("\n12. Subset operations:")
  let gene_subset = @src.sce_subset_genes(sce, [0, 2, 4])
  println(
    "  subset genes [0, 2, 4] -> n_genes=" +
    @src.sce_get_n_genes(gene_subset).to_string(),
  )
  let cell_subset = @src.sce_subset_cells(sce, [0, 1])
  println(
    "  subset cells [0, 1] -> n_cells=" +
    @src.sce_get_n_cells(cell_subset).to_string(),
  )
  let row_data = @src.sce_get_row_data(built_sce, "gene_type")
  println("  row_data 'gene_type' length: " + row_data.length().to_string())
  let mut g = 0
  while g < row_data.length() {
    println("    gene " + (g + 1).to_string() + ": " + row_data[g])
    g = g + 1
  }
  let col_data = @src.sce_get_col_data(built_sce, "batch")
  println("  col_data 'batch' length: " + col_data.length().to_string())
  let mut c = 0
  while c < col_data.length() {
    println("    cell " + (c + 1).to_string() + ": " + col_data[c])
    c = c + 1
  }
  let pca_dim = @src.sce_get_reduced_dim(built_sce, "PCA")
  println("  PCA reduced dim length: " + pca_dim.length().to_string())
  let sce_with_set_dim = @src.sce_set_reduced_dim(sce, "custom_dim", [
    [0.5, 0.5],
    [0.6, 0.6],
  ])
  println(
    "  set custom_dim: " +
    @src.sce_get_reduced_dim(sce_with_set_dim, "custom_dim")
    .length()
    .to_string() +
    " rows",
  )

  // 13. Demonstrate merging multiple SCE
  println("\n13. Merge multiple SCE:")
  let sce1 = @src.sce_create_example()
  let sce2 = @src.sce_create_example()
  let merged = @src.sce_merge([sce1, sce2])
  println(
    "  merged n_genes=" +
    @src.sce_get_n_genes(merged).to_string() +
    ", n_cells=" +
    @src.sce_get_n_cells(merged).to_string(),
  )

  // 14. Demonstrate setting reduced dims on an SCE
  // Note: sce_create_example_with_reduced_dims() internally calls
  // sce_compute_pca on the (5x4) example matrix which hits the source-side
  // square-matrix assumption shown above. We instead demonstrate the
  // sce_set_reduced_dim / sce_get_reduced_dim accessors directly, attaching
  // the reduced dims computed in sections 5-7 to a fresh SCE.
  println("\n14. SCE with reduced dims (manual):")
  let sce_rd = @src.sce_set_reduced_dim(sce, "PCA", pca_result.coordinates)
    |> @src.sce_set_reduced_dim("tSNE", tsne_result.coordinates)
    |> @src.sce_set_reduced_dim("UMAP", umap_result.coordinates)
  let pca_coords = @src.sce_get_reduced_dim(sce_rd, "PCA")
  let tsne_coords = @src.sce_get_reduced_dim(sce_rd, "tSNE")
  let umap_coords = @src.sce_get_reduced_dim(sce_rd, "UMAP")
  println("  PCA coords: " + pca_coords.length().to_string() + " rows")
  println("  tSNE coords: " + tsne_coords.length().to_string() + " rows")
  println("  UMAP coords: " + umap_coords.length().to_string() + " rows")

  // 15. Demonstrate data summary
  println("\n15. SCE Summary:")
  println(@src.sce_summary(sce))

  println("\nSingleCellExperiment demo completed!")
}
