///|
/// HTSFilter module example - RNA-seq count filtering by CPM threshold.
///
/// Run: moon run examples/htsfilter_demo
fn main {
  println("=== Bioconductor HTSFilter Module Demo ===")
  println("RNA-seq Count Filtering by CPM Threshold")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Sample RNA-seq Count Data ===")
  let (counts, groups, gene_names) = @bio.hts_filter_sample_data()
  println("Input dimensions: " + counts.length().to_string() + " genes x " +
    counts[0].length().to_string() + " samples")
  let mut ctrl = 0
  let mut treat = 0
  let mut i = 0
  while i < groups.length() {
    if groups[i] == 0 { ctrl = ctrl + 1 } else { treat = treat + 1 }
    i = i + 1
  }
  println("Control samples: " + ctrl.to_string() + ", Treatment samples: " + treat.to_string())
  println("")

  // ---- Section 2: Compute CPM ----
  println("=== 2. Calculate CPM Matrix ===")
  let lib_sizes = @bio.hts_filter_library_sizes(counts)
  println("Library sizes (first 3): " +
    lib_sizes[0].to_string() + ", " +
    lib_sizes[1].to_string() + ", " +
    lib_sizes[2].to_string())
  let cpm_mat = @bio.hts_filter_cpm(counts, lib_sizes)
  println("CPM matrix computed: " + cpm_mat.length().to_string() + " genes")
  println("")

  // ---- Section 3: Filter with different thresholds ----
  println("=== 3. Filter with Different CPM Thresholds ===")
  let thresholds = [10.0, 50.0, 100.0, 500.0]
  let mut ti = 0
  while ti < thresholds.length() {
    let threshold = thresholds[ti]
    let result = @bio.hts_filter(counts, groups, threshold, 2)
    let rate = @bio.hts_filter_retention_rate(result) * 100.0
    println("  CPM >= " + threshold.to_string() +
      " | Kept: " + result.n_genes_kept.to_string() +
      "/" + result.n_genes_input.to_string() +
      " (" + rate.to_string() + "%)")
    ti = ti + 1
  }
  println("")

  // ---- Section 4: Filtered data application ----
  println("=== 4. Apply Filter to Get Filtered Counts ===")
  let result = @bio.hts_filter(counts, groups, 100.0, 2)
  let filtered_counts = @bio.hts_filter_apply(counts, result.keep)
  let filtered_names = @bio.hts_filter_apply_names(gene_names, result.keep)
  println("Filtered matrix: " + filtered_counts.length().to_string() + " genes x " +
    filtered_counts[0].length().to_string() + " samples")
  println("First 5 kept genes:")
  let mut fi = 0
  while fi < 5 && fi < filtered_names.length() {
    println("  " + filtered_names[fi] + " (mean CPM: " +
      (filtered_counts[fi][0] + filtered_counts[fi][1] + filtered_counts[fi][2]).to_string() + ")")
    fi = fi + 1
  }
  println("")

  // ---- Section 5: Summary ----
  println("=== 5. Filtering Summary ===")
  println(@bio.hts_filter_summary(result))
}
