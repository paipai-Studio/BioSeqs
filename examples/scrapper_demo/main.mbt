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
fn format_ints(values : Array[Int]) -> String {
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
fn format_strings(values : Array[String]) -> String {
  let mut result = "["
  let mut index = 0
  while index < values.length() {
    if index > 0 {
      result = result + ", "
    }
    result = result + values[index]
    index = index + 1
  }
  result + "]"
}

///|
fn main {
  println("=== Bioconductor scrapper Demo ===")
  let counts = @src.scrapper_sample_counts()
  let blocks = ["batch1", "batch1", "batch1", "batch2", "batch2", "batch2"]
  let mitochondrial = @src.ScrapperNamedSubset::new("mito", [2, 4])

  println("\n1. Batch-aware RNA quality control")
  let qc = @src.scrapper_quick_rna_qc(counts, subsets=[mitochondrial], blocks~) catch {
    _ => abort("failed to compute scrapper RNA QC")
  }
  println("   " + qc.summary())
  println("   library sums: " + format_doubles(qc.metrics.sums))
  println("   detected genes: " + format_ints(qc.metrics.detected))
  println(
    "   mitochondrial proportions: " +
    format_doubles(qc.metrics.subset_proportions[0]),
  )

  println("\n2. Library-size factors and log-normalization")
  let size_factors = @src.scrapper_library_size_factors(
    counts,
    blocks~,
    mode=@src.scrapper_center_per_block(),
  ) catch {
    _ => abort("failed to compute library size factors")
  }
  let logcounts = @src.scrapper_normalize_counts(counts, size_factors) catch {
    _ => abort("failed to normalize counts")
  }
  println("   size factors: " + format_doubles(size_factors))
  println("   first normalized feature: " + format_doubles(logcounts[0]))

  println("\n3. Mean-variance trend and highly variable genes")
  let model = @src.scrapper_model_gene_variances(
    logcounts,
    mean_filter=false,
    min_window_count=3,
  ) catch {
    _ => abort("failed to model gene variances")
  }
  let highly_variable = model.highly_variable_genes(top=3) catch {
    _ => abort("failed to choose highly variable genes")
  }
  println("   means: " + format_doubles(model.means))
  println("   residual variances: " + format_doubles(model.residuals))
  println("   selected 0-based genes: " + format_ints(highly_variable))

  println("\n4. Multi-factor pseudo-bulk aggregation")
  let aggregate = @src.scrapper_aggregate_across_cells(
    counts,
    [
      @src.ScrapperFactor::new("cluster", ["T", "T", "B", "T", "B", "B"]),
      @src.ScrapperFactor::new("batch", blocks),
    ],
    compute_median=true,
  ) catch {
    _ => abort("failed to aggregate cells")
  }
  let pseudo_bulk_means = aggregate.means() catch {
    _ => abort("failed to compute pseudo-bulk means")
  }
  println("   " + aggregate.summary())
  println("   groups: " + format_strings(aggregate.group_names))
  println("   cells per group: " + format_ints(aggregate.counts))
  println("   first feature means: " + format_doubles(pseudo_bulk_means[0]))

  println("\n5. Immutable SingleCellExperiment integration")
  let sce = @src.SingleCellExperiment::new(
    counts,
    ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6"],
    ["Cell1", "Cell2", "Cell3", "Cell4", "Cell5", "Cell6"],
  )
  let normalized_sce = @src.scrapper_normalize_rna_counts_sce(
    sce,
    blocks~,
    mode=@src.scrapper_center_per_block(),
  ) catch {
    _ => abort("failed to normalize SingleCellExperiment")
  }
  let (annotated_sce, sce_qc) = @src.scrapper_quick_rna_qc_sce(
    normalized_sce,
    subsets=[mitochondrial],
    blocks~,
  ) catch {
    _ => abort("failed to annotate SingleCellExperiment QC")
  }
  println(
    "   original logcounts rows: " +
    @src.sce_get_assay(sce, "logcounts").length().to_string(),
  )
  println(
    "   copied logcounts rows: " +
    @src.sce_get_assay(annotated_sce, "logcounts").length().to_string(),
  )
  println(
    "   copied QC metadata columns populated: " +
    @src.sce_get_col_data(annotated_sce, "scrapper_keep").length().to_string(),
  )
  println("   SCE QC retained cells: " + sce_qc.retained_count().to_string())

  println("\n=== Demo Complete ===")
}
