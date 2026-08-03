///|
/// NanoString nCounter demo - Digital barcode counting analysis.
///
/// Demonstrates sample/container construction, positive/negative/housekeeping
/// and content normalization, QC (imaging, linearity, LOQ), differential
/// expression analysis with t-tests and BH-FDR, and textual summaries.
///
/// Reference: RTC (2018) NanoStringNCTools: NanoString nCounter Tools.
/// Bioconductor package.
/// Run: moon run examples/nanostring_demo
fn main {
  println("=== NanoString nCounter Demo ===")
  println("Digital molecular barcode counting analysis")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load sample data ===")
  let data = @bio.ns_sample_data()
  println("Samples: \{data.n_samples()}")
  println("Endogenous genes: \{data.n_genes()}")
  let samples = data.samples()
  for s in samples {
    println(
      "  " +
        s.sample_id() +
        ": raw=" +
        s.raw_counts().length().to_string() +
        " pos_ctrl=" +
        s.positive_controls().length().to_string() +
        " neg_ctrl=" +
        s.negative_controls().length().to_string() +
        " hk=" +
        s.housekeeping().length().to_string() +
        " endo=" +
        s.endogenous().length().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Geometric mean and log2 utility demonstration
  // ------------------------------------------------------------------
  println("=== 2. Utility functions ===")
  let gm = @bio.ns_geometric_mean([2.0, 8.0])
  println("geometric_mean([2, 8]) = " + gm.to_string() + " (expected 4.0)")
  let log2_vals = @bio.ns_log2_counts([0.0, 1.0, 7.5, 15.5])
  println("log2_counts([0, 1, 7.5, 15.5]) = " + log2_vals.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Positive control normalization
  // ------------------------------------------------------------------
  println("=== 3. Positive control normalization ===")
  let pos_results = @bio.ns_positive_control_norm(data)
  for r in pos_results {
    println(
      "  " +
        r.sample_id() +
        ": factor=" +
        r.norm_factors()[0].to_string() +
        " method=" +
        r.method() +
        " first_count=" +
        r.normalized_counts()[0].to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Negative control (background) normalization
  // ------------------------------------------------------------------
  println("=== 4. Negative control normalization ===")
  let neg_results = @bio.ns_negative_control_norm(data)
  for r in neg_results {
    println(
      "  " +
        r.sample_id() +
        ": factor=" +
        r.norm_factors()[0].to_string() +
        " method=" +
        r.method() +
        " first_count=" +
        r.normalized_counts()[0].to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Housekeeping gene normalization
  // ------------------------------------------------------------------
  println("=== 5. Housekeeping normalization ===")
  let hk_results = @bio.ns_housekeeping_norm(data)
  for r in hk_results {
    println(
      "  " +
        r.sample_id() +
        ": factor=" +
        r.norm_factors()[0].to_string() +
        " method=" +
        r.method() +
        " first_count=" +
        r.normalized_counts()[0].to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Content normalization (positive + housekeeping)
  // ------------------------------------------------------------------
  println("=== 6. Content normalization ===")
  let content_results = @bio.ns_content_norm(data)
  for r in content_results {
    let n_factors = r.norm_factors().length()
    let factors_str = if n_factors >= 2 {
      "pos=" + r.norm_factors()[0].to_string() + ", hk=" + r.norm_factors()[1].to_string()
    } else {
      "factor=" + r.norm_factors()[0].to_string()
    }
    println(
      "  " +
        r.sample_id() +
        ": " +
        factors_str +
        " method=" +
        r.method() +
        " first_count=" +
        r.normalized_counts()[0].to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Quality control
  // ------------------------------------------------------------------
  println("=== 7. Quality control ===")
  let imaging = @bio.ns_imaging_qc(data, 100)
  let (r_sq, slope) = @bio.ns_pos_linearity_qc(data)
  let loq = @bio.ns_loq_qc(data, 50)
  println("Imaging QC (threshold=100):")
  for i in 0..<samples.length() {
    println(
      "  " +
        samples[i].sample_id() +
        ": pass=" +
        imaging[i].to_string() +
        ", below_LOQ=" +
        loq[i].to_string(),
    )
  }
  println("Positive control linearity: R^2=" + r_sq.to_string() + ", slope=" + slope.to_string())
  println("")

  // ------------------------------------------------------------------
  // 8. Differential expression analysis
  // ------------------------------------------------------------------
  println("=== 8. Differential expression analysis ===")
  let group : Array[Bool] = [
    false, false, false, true, true, true,
  ]
  println("Group: 3 control (false) + 3 treatment (true)")
  let de_results = @bio.ns_differential_expression(
    content_results, group, data.gene_names(),
  )
  let mut n_sig = 0
  let mut n_up = 0
  let mut n_down = 0
  for r in de_results {
    if r.significant() {
      n_sig = n_sig + 1
      if r.log_fc() > 0.0 {
        n_up = n_up + 1
      } else {
        n_down = n_down + 1
      }
    }
  }
  println("Genes tested: \{de_results.length()}")
  println("Significant (adj p < 0.05): \{n_sig}")
  println("  Up-regulated: \{n_up}")
  println("  Down-regulated: \{n_down}")
  println("")

  // Manual bubble sort by adjusted p-value (ascending) for display.
  let sorted_de = de_results.copy()
  let m = sorted_de.length()
  let mut i = 0
  while i < m - 1 {
    let mut j = 0
    while j < m - 1 - i {
      if sorted_de[j].adj_p_value() > sorted_de[j + 1].adj_p_value() {
        let tmp = sorted_de[j]
        sorted_de[j] = sorted_de[j + 1]
        sorted_de[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  println("Top 10 DE genes by adjusted p-value:")
  let top_n = if m < 10 { m } else { 10 }
  for k in 0..<top_n {
    let r = sorted_de[k]
    let flag = if r.significant() { "*" } else { " " }
    println(
      "  " +
        flag +
        " " +
        r.gene() +
        " logFC=" +
        r.log_fc().to_string() +
        " p=" +
        r.p_value().to_string() +
        " adjP=" +
        r.adj_p_value().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 9. Textual summaries
  // ------------------------------------------------------------------
  println("=== 9. Summaries ===")
  println(@bio.ns_summary(data))
  println(@bio.ns_qc_summary(data))
  println(@bio.ns_de_summary(de_results))

  println("=== NanoString nCounter Demo complete ===")
}
