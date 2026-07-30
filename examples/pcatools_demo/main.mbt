fn main {
  println("=== Bioconductor-style PCAtools Demo ===\n")

  println("--- 1. Run PCA on sample data ---")
  let data = [
    [1.0,  2.0,  3.0],
    [2.0,  4.0,  5.0],
    [3.0,  6.0,  7.0],
    [4.0,  8.0,  9.0],
    [5.0, 10.0, 11.0],
    [6.0, 12.0, 10.0],
    [7.0, 14.0,  8.0],
    [8.0, 16.0,  5.0],
    [9.0, 18.0,  3.0],
    [50.0, 50.0, 50.0],  // outlier
  ]
  let n_rows = data.length()
  let n_cols = if n_rows > 0 { data[0].length() } else { 0 }
  println("  Input: \{n_rows} samples x \{n_cols} variables")

  let res = @bio.pcatools_run_pca(data, 3, true)
  println(@bio.pcatools_summary(res))

  println("\n--- 2. Scree plot (ASCII) ---")
  let scree = @bio.scree_plot_ascii(res, height=10, width_per_pc=6)
  println(scree)

  println("\n--- 3. Biplot (ASCII) ---")
  let opts = @bio.BiplotOptions::new()
  let biplot = @bio.pca_biplot_ascii(res, opts=opts)
  println(biplot)

  println("\n--- 4. Outlier detection ---")
  let outlier_res = @bio.find_pca_outliers(res, n_pcs=2, quantile_cutoff=0.9)
  println("  Tested samples: \{outlier_res.n_tested}")
  println("  Distance cutoff: \{outlier_res.cutoff}")
  let n_outliers = outlier_res.indices.length()
  println("  Outliers found: \{n_outliers}")
  if outlier_res.indices.length() > 0 {
    let mut idxs = "  Outlier indices: [".to_string()
    let mut oi = 0
    while oi < outlier_res.indices.length() {
      if oi > 0 { idxs = idxs + ", " }
      idxs = idxs + outlier_res.indices[oi].to_string()
      oi = oi + 1
    }
    idxs = idxs + "]"
    println(idxs)
  }

  println("\n--- 5. Variable Correlations with PCs ---")
  let corr = @bio.variable_correlations(data, res, which_pcs=2)
  let cr = corr.length()
  let cc = if cr > 0 { corr[0].length() } else { 0 }
  println("  Variables x PCs shape: \{cr} x \{cc}")
  let mut vi = 0
  while vi < corr.length() {
    let row = corr[vi]
    let mut s = "    var\{vi + 1} vs PCs: [".to_string()
    let mut pj = 0
    while pj < row.length() {
      if pj > 0 { s = s + ", " }
      s = s + row[pj].to_string()
      pj = pj + 1
    }
    s = s + "]"
    println(s)
    vi = vi + 1
  }

  println("\nPCAtools demo completed successfully!")
}
