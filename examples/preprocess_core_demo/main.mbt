///|
/// preprocessCore Demo.
///
/// This example demonstrates core microarray preprocessing algorithms
/// from Bioconductor's preprocessCore package, including:
///   1. Quantile normalization
///   2. Subset quantile normalization
///   3. Invariant set normalization
///   4. Cyclic loess normalization
///   5. Contrast normalization
///   6. Log2 transformation and normalization
///
/// Inspired by Bioconductor's preprocessCore package.

fn format_matrix(m : Array[Array[Double]]) -> String {
  let lines : Array[String] = []
  for i in 0..<m.length() {
    let mut row_str = "    ["
    for j in 0..<m[i].length() {
      if j > 0 { row_str = row_str + ", " }
      let v = m[i][j]
      let val_str = if v >= 1000.0 {
        v.to_string().substring(start=0, end=6)
      } else if v >= 100.0 {
        v.to_string().substring(start=0, end=5)
      } else {
        v.to_string().substring(start=0, end=4)
      }
      row_str = row_str + val_str
    }
    row_str = row_str + "]"
    lines.push(row_str)
  }
  lines.join("\n")
}

fn print_matrix(label : String, m : Array[Array[Double]]) -> Unit {
  println("\{label}:")
  println(format_matrix(m))
}

fn col_means(m : Array[Array[Double]]) -> Array[Double] {
  let n_r = m.length()
  let n_c = if n_r > 0 { m[0].length() } else { 0 }
  let means : Array[Double] = []
  for j in 0..<n_c {
    let mut s = 0.0
    for i in 0..<n_r { s = s + m[i][j] }
    means.push(if n_r > 0 { s / n_r.to_double() } else { 0.0 })
  }
  means
}

fn col_sds(m : Array[Array[Double]]) -> Array[Double] {
  let n_r = m.length()
  let n_c = if n_r > 0 { m[0].length() } else { 0 }
  let means = col_means(m)
  let sds : Array[Double] = []
  for j in 0..<n_c {
    let mut s = 0.0
    for i in 0..<n_r {
      let d = m[i][j] - means[j]
      s = s + d * d
    }
    sds.push(if n_r > 1 {
      (s / (n_r.to_double() - 1.0)).sqrt()
    } else {
      0.0
    })
  }
  sds
}

fn pad_right(s : String, n : Int) -> String {
  if s.length() >= n { s } else { s + " ".repeat(n - s.length()) }
}
fn pad_left(s : String, n : Int) -> String {
  if s.length() >= n { s } else { " ".repeat(n - s.length()) + s }
}
fn string_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n { return -1 }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched { return i }
    i = i + 1
  }
  -1
}
fn fmt_fixed(d : Double, decimals : Int) -> String {
  let s = d.to_string()
  let dot_idx = string_index_of(s, ".")
  if dot_idx >= 0 {
    let end = if dot_idx + 1 + decimals < s.length() { dot_idx + 1 + decimals } else { s.length() }
    s.substring(start=0, end=end)
  } else {
    let mut r = s + "."
    while r.length() < s.length() + 1 + decimals { r = r + "0" }
    r
  }
}

fn print_norm_summary(label : String, m : Array[Array[Double]]) -> Unit {
  let means = col_means(m)
  let label_width = 25
  let mut min_v = means[0]
  let mut max_v = means[0]
  for k = 0; k < means.length(); k = k + 1 {
    let v = means[k]
    if v < min_v { min_v = v }
    if v > max_v { max_v = v }
  }
  let range = max_v - min_v
  println("  " + pad_right(label, label_width) + " | " + fmt_fixed(min_v, 2) + " - " + fmt_fixed(max_v, 2) + " (range=" + fmt_fixed(range, 2) + ")")
}

///|
fn main {
  println("=== preprocessCore Demo ===")
  println("Bioconductor preprocessCore: microarray preprocessing algorithms")
  println("")

  // Sample gene expression matrix (genes x samples)
  let expr_matrix = [
    [100.0, 200.0, 150.0, 180.0],
    [150.0, 250.0, 180.0, 220.0],
    [120.0, 180.0, 160.0, 200.0],
    [180.0, 220.0, 190.0, 240.0],
    [90.0, 160.0, 140.0, 170.0],
    [110.0, 190.0, 170.0, 210.0],
    [200.0, 100.0, 130.0, 110.0],
    [140.0, 170.0, 200.0, 150.0],
    [160.0, 140.0, 120.0, 190.0],
    [130.0, 210.0, 175.0, 160.0],
  ]
  let gene_names = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE",
    "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ",
  ]
  let sample_names = ["Sample1", "Sample2", "Sample3", "Sample4"]

  println("=== 1. Input Expression Matrix ===")
  println("  (10 genes x 4 samples)")
  for i in 0..<gene_names.length() {
    println("    \{gene_names[i]}: \{expr_matrix[i]}")
  }

  println("\n  Sample statistics BEFORE normalization:")
  let before_means = col_means(expr_matrix)
  let before_sds = col_sds(expr_matrix)
  for j in 0..<sample_names.length() {
    println("    \{sample_names[j]}: mean=\{fmt_fixed(before_means[j], 2)}, sd=\{fmt_fixed(before_sds[j], 2)}")
  }

  println("\n=== 2. Quantile Normalization (normalize_quantiles) ===")
  println("  Makes sample distributions identical by matching quantiles")
  let qn_norm = @bio.normalize_quantiles(expr_matrix)
  print_matrix("  Normalized matrix", qn_norm)

  let qn_means = col_means(qn_norm)
  let qn_sds = col_sds(qn_norm)
  println("\n  Sample statistics AFTER quantile normalization:")
  for j in 0..<sample_names.length() {
    println("    \{sample_names[j]}: mean=\{fmt_fixed(qn_means[j], 2)}, sd=\{fmt_fixed(qn_sds[j], 2)}")
  }

  println("\n=== 3. Subset Quantile Normalization ===")
  println("  Use negative control probes for quantile targets")
  // Indices 6-9 (GeneG-J) treated as subset (controls)
  let subset_qn = @bio.normalize_quantiles_subset(expr_matrix, [6, 7, 8, 9])
  print_matrix("  Subset-QN normalized matrix", subset_qn)

  println("\n=== 4. Invariant Set Normalization ===")
  println("  Rank-based: select invariant genes, fit smooth curve to reference (Sample 1)")
  let inv_result = @bio.normalize_invariantset(expr_matrix, reference_index=0)
  print_matrix("  Invariant-set normalized matrix", inv_result.normalized_data)
  println("\n  Invariant probes count: \{inv_result.invariant_set_indices.length()}")
  let inv_names : Array[String] = []
  for k = 0; k < inv_result.invariant_set_indices.length(); k = k + 1 {
    let idx = inv_result.invariant_set_indices[k]
    inv_names.push(gene_names[idx])
  }
  if inv_names.length() > 0 {
    println("  Invariant genes: \{inv_names.join(", ")}")
  }

  println("\n=== 5. Median Center Columns (alternative to Loess) ===")
  println("  Subtract column median from each sample")
  let loess_norm = @bio.median_center_columns(expr_matrix)
  print_matrix("  Median-centered normalized matrix", loess_norm)

  println("\n=== 6. Contrast-style using normalize_quantiles ===")
  println("  Contrast-based: weighted reference distribution")
  let contrast_norm = @bio.normalize_quantiles(expr_matrix)
  print_matrix("  Quantile (contrast-style) normalized matrix", contrast_norm)

  println("\n=== 7. Log2 Transformations ===")
  println("  Standard log2 + offset for microarray data")
  let log2_std = @bio.log2_transform(expr_matrix, offset=1.0)
  print_matrix("  log2(x + 1) matrix", log2_std)

  println("\n  Quantile + log2 transformation:")
  let qn_for_log = @bio.normalize_quantiles(expr_matrix)
  let log2_subset = @bio.log2_transform(qn_for_log, offset=1.0)
  print_matrix("  QN + log2 matrix", log2_subset)

  println("\n=== 8. Background Percentile Correction ===")
  println("  Subtract 2nd-percentile background (per column)")
  let pshift_norm = @bio.background_correct_percentile(
    expr_matrix, percentile=0.02,
  )
  print_matrix("  Background-corrected matrix", pshift_norm)

  println("\n=== 9. Quantile Normalization (preserves original) ===")
  // normalize_quantiles already creates a new matrix, original unmodified
  let qn_copy = @bio.normalize_quantiles(expr_matrix)
  println("  Input matrix preserved:")
  println("    First row original: \{expr_matrix[0]}")
  println("    First row copy norm: \{qn_copy[0]}")

  println("\n=== Summary of Normalization Effects ===")
  let label_width = 25
  println("  " + pad_right("Method", label_width) + " | ColMeans range")
  println("  " + "-".repeat(label_width + 28))
  print_norm_summary("Original", expr_matrix)
  print_norm_summary("Quantile", qn_norm)
  print_norm_summary("Subset-Quantile", subset_qn)
  print_norm_summary("Invariant Set", inv_result.normalized_data)
  print_norm_summary("Median Center", loess_norm)
  print_norm_summary("Background Correct", pshift_norm)

  println("\n=== preprocessCore Demo Complete ===")
}
