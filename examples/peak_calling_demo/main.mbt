///|
/// MACS2-like ChIP-seq peak calling demo.
///
/// Generates a synthetic ChIP-seq dataset (treatment + matched control),
/// runs the peak calling pipeline with default parameters, displays the
/// called peaks, applies FDR correction and fold-enrichment filtering,
/// and prints a summary.
fn main {
  println("=== MACS2-like Peak Calling Demo ===")
  println("Inspired by the MACS2 algorithm (Zhang et al. 2008, Genome Biol.)")
  println("")

  // 1. Generate sample data.
  println("1. Generating sample ChIP-seq data...")
  let (treatment, control) = @bio.peak_calling_sample_data()
  println("   Treatment reads: \{treatment.length()}")
  println("   Control reads:   \{control.length()}")
  println("")

  // 2. Display read distribution in known enriched regions.
  println("2. Read distribution in candidate enriched regions:")
  let regions = [("chr1", 1000, 1500), ("chr1", 5000, 5500), ("chr1", 9000, 9500)]
  for region in regions {
    let t = @bio.count_reads_in_window(treatment, region.0, region.1, region.2)
    let c = @bio.count_reads_in_window(control, region.0, region.1, region.2)
    println("   \{region.0}:\{region.1}-\{region.2}  treatment=\{t}  control=\{c}")
  }
  println("")

  // 3. Call peaks with default parameters.
  println("3. Calling peaks with default parameters...")
  let default_params = @bio.PeakCallingParams::default()
  println("   window_size=\{default_params.window_size}")
  println("   step_size=\{default_params.step_size}")
  println("   local_lambda_size=\{default_params.local_lambda_size}")
  println("   pvalue_threshold=\{default_params.pvalue_threshold}")
  println("   fdr_threshold=\{default_params.fdr_threshold}")
  println("   min_fold_enrichment=\{default_params.min_fold_enrichment}")
  let peaks_default = @bio.call_peaks(treatment, control, default_params)
  println("   Peaks called (default): \{peaks_default.length()}")
  println("")

  // 4. Call peaks with relaxed parameters to detect more peaks.
  println("4. Calling peaks with relaxed parameters...")
  let relaxed_params = @bio.PeakCallingParams::new(
    window_size=200,
    step_size=50,
    local_lambda_size=2000,
    pvalue_threshold=0.01,
    fdr_threshold=0.5,
    min_fold_enrichment=1.5,
  )
  let peaks_relaxed = @bio.call_peaks(treatment, control, relaxed_params)
  println("   Peaks called (relaxed): \{peaks_relaxed.length()}")
  println("")

  // 5. Display detailed peak information.
  println("5. Detailed peak information (relaxed):")
  let mut i = 0
  while i < peaks_relaxed.length() {
    let p = peaks_relaxed[i]
    println(
      "   Peak \{i + 1}: \{p.chr}:\{p.start}-\{p.end} " +
      "summit=\{p.summit} reads=\{p.read_count} ctrl=\{p.control_count} " +
      "fold=\{p.fold_enrichment} p=\{p.p_value} fdr=\{p.fdr} " +
      "sig=\{p.significant}",
    )
    i = i + 1
  }
  println("")

  // 6. Show significant peaks only.
  println("6. Significant peaks (fdr <= 0.05):")
  let mut n_sig = 0
  for p in peaks_relaxed {
    if p.significant {
      n_sig = n_sig + 1
      println(
        "   \{p.chr}:\{p.start}-\{p.end}  reads=\{p.read_count}  fdr=\{p.fdr}",
      )
    }
  }
  if n_sig == 0 {
    println("   (no significant peaks at FDR <= 0.05)")
  } else {
    println("   Total significant: \{n_sig}")
  }
  println("")

  // 7. Summary statistics.
  println("7. Summary statistics:")
  let summary = @bio.peak_calling_summary(peaks_relaxed)
  println(summary)

  // 8. Demonstrate Poisson p-value computation.
  println("8. Poisson p-value examples:")
  let examples = [(1, 1.0), (3, 2.0), (5, 1.0), (10, 5.0), (20, 1.0)]
  for ex in examples {
    let p = @bio.poisson_pvalue(ex.0, ex.1)
    println("   P(X >= \{ex.0} | lambda=\{ex.1}) = \{p}")
  }
  println("")

  // 9. Demonstrate fold enrichment calculation.
  println("9. Fold enrichment examples:")
  let fold_examples = [(10, 5.0), (20, 2.0), (50, 10.0), (100, 1.0)]
  for fe in fold_examples {
    let f = @bio.compute_fold_enrichment(fe.0, fe.1)
    println("   observed=\{fe.0}  expected=\{fe.1}  fold=\{f}")
  }
  println("")

  println("=== MACS2-like Peak Calling Demo Complete ===")
}
