///|
/// hmmcopy module demo - HMM-based copy number variation calling.
///
/// Demonstrates GC bias correction, Viterbi decoding, CNV calling,
/// and segmentation of genomic bins into copy number segments.
///
/// Reference: Lai D et al. (2018) HMMcopy. Bioconductor package.
/// Run: moon run examples/hmmcopy_demo
fn main {
  println("=== hmmcopy Module Demo ===")
  println("HMM-based copy number variation calling")
  println("Reference: Bioconductor HMMcopy")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let bins = @bio.hmmcopy_sample_bins()
  println("Total bins: " + bins.length().to_string())
  println("Chromosome: " + bins[0].chr)
  println("Bin width: " + bins[0].width().to_string() + " bp")
  println("")

  // Show coverage profile
  println("Coverage profile (first 20 bins):")
  let mut i = 0
  while i < 20 && i < bins.length() {
    let bar_len = bins[i].reads / 5
    let mut bar = ""
    let mut j = 0
    while j < bar_len {
      bar = bar + "#"
      j = j + 1
    }
    println("  bin " + i.to_string() + ": " + bins[i].reads.to_string() + " reads " + bar)
    i = i + 1
  }
  println("  ... (bins 10-19 are CN0 region with low reads)")
  println("")

  // ------------------------------------------------------------------
  // 2. GC bias correction
  // ------------------------------------------------------------------
  println("=== 2. GC Bias Correction ===")
  let corrected = @bio.hmmcopy_correct_gc(bins)
  println("Corrected " + corrected.length().to_string() + " bins")
  // Show before/after for a few bins
  println("Before/after correction (first 5 bins):")
  i = 0
  while i < 5 {
    println("  bin " + i.to_string() + ": " + bins[i].reads.to_string() + " -> " + corrected[i].corrected_reads.to_string())
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. HMM model configuration
  // ------------------------------------------------------------------
  println("=== 3. HMM Model Configuration ===")
  let model = @bio.hmmcopy_sample_model()
  println("Number of states: " + model.num_states.to_string())
  println("Ploidy: " + model.ploidy.to_string())
  println("Emission means (CN2 units):")
  i = 0
  while i < model.num_states {
    println("  CN" + i.to_string() + ": " + model.emissions[i].to_string())
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. CNV calling
  // ------------------------------------------------------------------
  println("=== 4. CNV Calling (Viterbi) ===")
  let called = @bio.hmmcopy_call_cnv(corrected, model)
  println("Called copy numbers for " + called.length().to_string() + " bins")
  println("")

  // Show state distribution
  let state_counts = Array::make(7, 0)
  for b in called {
    state_counts[b.copy] = state_counts[b.copy] + 1
  }
  println("Copy number distribution:")
  i = 0
  while i < 7 {
    if state_counts[i] > 0 {
      println("  CN" + i.to_string() + ": " + state_counts[i].to_string() + " bins")
    }
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Segmentation
  // ------------------------------------------------------------------
  println("=== 5. Segmentation ===")
  let segments = @bio.hmmcopy_segment(called)
  println("Number of segments: " + segments.length().to_string())
  println("")
  println("Segments (showing first 10):")
  let mut shown = 0
  for seg in segments {
    if shown < 10 {
      println(
        "  " +
          seg.chr +
          ":" +
          seg.start.to_string() +
          "-" +
          seg.end.to_string() +
          " CN" +
          seg.copy.to_string() +
          " (" +
          seg.num_bins.to_string() +
          " bins, mean=" +
          seg.mean_reads.to_string() +
          ")",
      )
      shown = shown + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Log-likelihood
  // ------------------------------------------------------------------
  println("=== 6. Model Log-Likelihood ===")
  let ll = @bio.hmmcopy_log_likelihood(corrected, model)
  println("Log-likelihood: " + ll.to_string())
  println("")

  // ------------------------------------------------------------------
  // 7. Identify CNV regions
  // ------------------------------------------------------------------
  println("=== 7. Identified CNV Regions ===")
  for seg in segments {
    if seg.copy != 2 {
      let cnv_type = if seg.copy < 2 { "LOSS" } else { "GAIN" }
      println(
        "  " +
          cnv_type +
          " CN" +
          seg.copy.to_string() +
          " at " +
          seg.chr +
          ":" +
          seg.start.to_string() +
          "-" +
          seg.end.to_string() +
          " (" +
          seg.num_bins.to_string() +
          " bins)",
      )
    }
  }
  println("")

  println("=== Demo complete ===")
}
