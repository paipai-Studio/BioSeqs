fn main {
  println("=== Bioconductor-style vsn Demo ===\n")

  println("--- 1. glog (Generalized Logarithm) ---")
  let xs = [0.1, 1.0, 10.0, 100.0, 1000.0]
  println("  vs glog transform with a=1.0, b=1.0:")
  let mut i = 0
  while i < xs.length() {
    let x = xs[i]
    let y = @bio.glog(x, 1.0, 1.0)
    let yinv = @bio.glog_inv(y, 1.0, 1.0)
    let diff_abs = (x - yinv).abs()
    println("    glog(\{x}) = \{y}  (inverse check: |x - glog_inv(y)| = \{diff_abs})")
    i = i + 1
  }

  println("\n--- 2. VSN2 Variance Stabilization ---")
  let data = [
    [100.0, 200.0, 150.0, 300.0],
    [150.0, 180.0, 160.0, 290.0],
    [90.0,  210.0, 140.0, 310.0],
    [110.0, 220.0, 155.0, 305.0],
    [120.0, 190.0, 145.0, 295.0],
    [5000.0, 8000.0, 6000.0, 9000.0],
    [6000.0, 7500.0, 6500.0, 8500.0],
    [5500.0, 8500.0, 5800.0, 9500.0],
  ]
  let fit = @bio.vsn_fit_and_report(data)
  println(@bio.summarize_vsn_fit(fit))
  let normalized = @bio.vsn2(data)
  println("\n  Normalized data (first 3 rows):")
  let mut r2 = 0
  while r2 < 3 && r2 < normalized.length() {
    let row = normalized[r2]
    let mut s2 = "    [".to_string()
    let mut c2 = 0
    while c2 < row.length() {
      if c2 > 0 { s2 = s2 + ", " }
      s2 = s2 + row[c2].to_string()
      c2 = c2 + 1
    }
    s2 = s2 + "]"
    println(s2)
    r2 = r2 + 1
  }

  println("\n--- 3. VSN Denoise ---")
  let denoised = @bio.vsn_denoise(data)
  let d_rows = denoised.length()
  let d_cols = if d_rows > 0 { denoised[0].length() } else { 0 }
  println("  Denoised data has shape [\{d_rows} x \{d_cols}]")
  let mut r = 0
  while r < 3 && r < denoised.length() {
    let row = denoised[r]
    let mut s = "    row[\{r}]: [".to_string()
    let mut c = 0
    while c < row.length() {
      if c > 0 { s = s + ", " }
      s = s + row[c].to_string()
      c = c + 1
    }
    s = s + "]"
    println(s)
    r = r + 1
  }

  println("\n--- 4. Mean-SD Bins ---")
  let bins = @bio.mean_sd_bins(data, 4)
  println("  Bins count: \{bins.length()}")
  let mut b = 0
  while b < bins.length() {
    let bin = bins[b]
    println("    bin[\{b}]: mean_x=\{bin.mean_x}, sd_y=\{bin.sd_y}, count=\{bin.count}")
    b = b + 1
  }

  println("\nVSN demo completed successfully!")
}
