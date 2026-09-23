///|
/// DMRcate module demo - Differentially Methylated Region calling.
///
/// Run: moon run examples/dmrcate_demo
fn main raise {
  println("=== DMRcate Module Demo ===")
  println("Differentially Methylated Region calling with spatial kernel smoothing")
  println("Reference: Peters TJ et al. (2015) Epigenetics & Chromatin 8:6")
  println("")

  // ------------------------------------------------------------------
  // 1. Create methylation beta-value matrix with positions
  // ------------------------------------------------------------------
  println("=== 1. Create synthetic methylation data ===")
  // 10 CpG sites, 6 samples (3 cases + 3 controls)
  // DMR signal: cases have elevated methylation at CpGs 3-7 on chr1
  let n_cpgs = 10
  let n_samples = 6
  let beta_values : Array[Array[Double]] = Array(capacity=0)
  let positions : Array[Int] = Array(capacity=0)
  let chr : Array[String] = Array(capacity=0)

  let mut i = 0
  while i < n_cpgs {
    chr.push("chr1")
    positions.push(10000 + i * 200)
    let vals : Array[Double] = Array(capacity=0)
    let mut s = 0
    while s < n_samples {
      let mut beta = 0.4
      // Cases: elevated methylation for CpGs 3-7
      if s < 3 && i >= 3 && i <= 7 {
        beta = beta + 0.35
      }
      // Small deterministic noise
      beta = beta + ((i * 3 + s * 7) % 5).to_double() * 0.01 - 0.02
      if beta < 0.0 {
        beta = 0.0
      }
      if beta > 1.0 {
        beta = 1.0
      }
      vals.push(beta)
      s = s + 1
    }
    beta_values.push(vals)
    i = i + 1
  }

  // Design matrix: intercept + group indicator
  let intercept : Array[Double] = Array(capacity=0)
  let group : Array[Double] = Array(capacity=0)
  let mut s = 0
  while s < n_samples {
    intercept.push(1.0)
    group.push(if s < 3 { 1.0 } else { 0.0 })
    s = s + 1
  }
  let design : Array[Array[Double]] = [intercept, group]
  let sample_names = ["case1", "case2", "case3", "ctrl1", "ctrl2", "ctrl3"]

  let input = @src.DMRcateInput::new(
    beta_values, positions, chr, design, sample_names,
  )
  println("Created input: \{input.n_cpgs()} CpGs, \{input.n_samples()} samples")
  println("")

  // ------------------------------------------------------------------
  // 2. Run DMRcate analysis
  // ------------------------------------------------------------------
  println("=== 2. Run DMRcate analysis ===")
  let result = @src.dmrcate_test(input, lambda=1000, c_scale=2.0, fdr_threshold=0.05)
  println("CpG sites tested: \{result.n_cpgs()}")
  println("DMRs called: \{result.n_dmrs()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Print significant CpGs
  // ------------------------------------------------------------------
  println("=== 3. Significant CpGs (FDR < 0.05) ===")
  println("CpG\t\tchr\tpos\tbeta\tt_stat\tp_value")
  let cpg_results = result.cpg_results()
  let mut n_sig = 0
  i = 0
  while i < cpg_results.length() {
    let cr = cpg_results[i]
    if cr.p_value() < 0.05 {
      n_sig = n_sig + 1
      println(
        "\{cr.cpg_name()}\t\{cr.chr()}\t\{cr.pos()}\t\{cr.beta()}\t\{cr.t_stat()}\t\{cr.p_value()}",
      )
    }
    i = i + 1
  }
  println("Total significant CpGs: \{n_sig}")
  println("")

  // ------------------------------------------------------------------
  // 4. Print DMRs
  // ------------------------------------------------------------------
  println("=== 4. Differentially Methylated Regions ===")
  let dmrs = result.dmrs()
  if dmrs.length() == 0 {
    println("No DMRs found at FDR < 0.05")
  } else {
    println("chr\tstart\tend\tn_cpgs\tmean_t\tmean_beta\tp_value")
    let mut j = 0
    while j < dmrs.length() {
      let d = dmrs[j]
      println(
        "\{d.chr()}\t\{d.start()}\t\{d.end_()}\t\{d.n_cpgs()}\t\{d.mean_t()}\t\{d.mean_beta()}\t\{d.p_value()}",
      )
      j = j + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Use example data from the module
  // ------------------------------------------------------------------
  println("=== 5. Run on built-in example data ===")
  let example_input = @src.dmrcate_example_data()
  println(
    "Example data: \{example_input.n_cpgs()} CpGs, \{example_input.n_samples()} samples",
  )
  let example_result = @src.dmrcate_test(example_input)
  println(@src.dmrcate_summary(example_result))
  println(@src.dmrcate_results_to_string(example_result))

  // ------------------------------------------------------------------
  // 6. Filter DMRs and extract CpG names
  // ------------------------------------------------------------------
  println("=== 6. Filter DMRs (min 2 CpGs, p < 0.05) ===")
  let example_dmrs = example_result.dmrs()
  let filtered = @src.dmrcate_filter_dmrs(example_dmrs, 2, 0.05)
  println("Filtered DMRs: \{filtered.length()}")
  if filtered.length() > 0 {
    let cpg_names = @src.dmrcate_dmr_cpg_names(filtered)
    println("CpGs in filtered DMRs: \{cpg_names.length()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Spatial kernel demonstration
  // ------------------------------------------------------------------
  println("=== 7. Spatial kernel weights ===")
  let lambda = 1000.0
  let mut d = 0
  while d <= 3000 {
    let w = @src.dmrcate_spatial_kernel(d.to_double(), lambda)
    println("  d=\{d} bp -> w=\{w}")
    d = d + 500
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Spatial autocorrelation
  // ------------------------------------------------------------------
  println("=== 8. Spatial autocorrelation ===")
  let example_cpgs = example_result.cpg_results()
  let autocorr = @src.dmrcate_spatial_autocorrelation(example_cpgs, 1000.0)
  println("Spatial autocorrelation of t-statistics: \{autocorr}")

  println("=== Demo complete ===")
}
