///|
/// csaw ChIP-seq Differential Binding Analysis Demo.
fn main {
  println("=== csaw ChIP-seq Differential Binding Analysis Demo ===\n")

  // Create windows with pre-assigned counts for multiple samples
  println("1. Creating windows with counts...")
  let windows = [
    @src.CswWindow::new("win_1", "chr1", 1000, 2000, [5.0, 10.0]),
    @src.CswWindow::new("win_2", "chr1", 2000, 3000, [15.0, 20.0]),
    @src.CswWindow::new("win_3", "chr1", 3000, 4000, [25.0, 50.0]),
    @src.CswWindow::new("win_4", "chr1", 4000, 5000, [8.0, 12.0]),
    @src.CswWindow::new("win_5", "chr1", 5000, 6000, [30.0, 45.0]),
    @src.CswWindow::new("win_6", "chr1", 6000, 7000, [12.0, 18.0]),
    @src.CswWindow::new("win_7", "chr1", 7000, 8000, [20.0, 35.0]),
    @src.CswWindow::new("win_8", "chr1", 8000, 9000, [3.0, 5.0]),
    @src.CswWindow::new("win_9", "chr1", 9000, 10000, [40.0, 60.0])
  ]
  println("   Created \{windows.length()} windows with 2 samples each")

  // Create dataset
  println("\n2. Creating CswDataSet...")
  let samples = ["control", "treatment"]
  let lib_sizes = [5000.0, 6000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  println("   Windows: \{data.csw_n_windows()}, Samples: \{data.csw_n_samples()}")

  // Library size normalization
  println("\n3. Library size normalization...")
  let norm_factors = @src.csw_norm_library(lib_sizes)
  let norm_data = @src.csw_apply_normalization(data, norm_factors)
  println("   Applied normalization factors")

  // Filtering
  println("\n4. Filtering low-abundance windows...")
  let filtered = @src.csw_filter_abundance(norm_data, 1.0)
  println("   Windows after filtering: \{filtered.csw_n_windows()}")

  // Differential binding test
  println("\n5. Differential binding testing (control vs treatment)...")
  let result = @src.csw_test_differential(filtered, [0], [1], 0.1)
  println("   Tested \{result.log_fc.length()} windows")
  println("   Significant windows (adj-p < 0.05): \{result.significant.filter(fn(x) { x }).length()}")

  // Find regions
  println("\n6. Finding differential regions...")
  let regions = @src.csw_find_regions(filtered, result, 1000)
  println("   Found \{regions.length()} differential regions")

  // Print results
  println("\n7. Summary of top results:")
  let mut i = 0
  while i < result.log_fc.length() && i < 5 {
    let logfc = result.log_fc[i].to_string()
    let pval = result.p_values[i].to_string()
    let adjp = result.adj_p_values[i].to_string()
    println("   Window \{i}: logFC=\{logfc}, p-val=\{pval}, adj-p=\{adjp}")
    i = i + 1
  }

  println("\n=== csaw Demo Complete ===")
}
