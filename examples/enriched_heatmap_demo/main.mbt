///|
/// EnrichedHeatmap module demo - genomic signal normalization.
///
/// Demonstrates normalization of genomic signals (e.g. CpG methylation,
/// histone marks) to target regions (e.g. TSS, promoters) for visualization
/// as enriched heatmaps.
///
/// Reference: Gu, Eils & Schlesner (2018) BMC Genomics 19:234.
/// Run: moon run examples/enriched_heatmap_demo
fn main {
  println("=== EnrichedHeatmap Module Demo ===")
  println("Genomic signal normalization for enriched heatmaps")
  println("Reference: Bioconductor EnrichedHeatmap")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Sample Data ===")
  let (signals, targets) = @bio.enriched_heatmap_sample()
  println("Signals: " + signals.length().to_string())
  println("Targets: " + targets.length().to_string())
  println("Target 0: " + targets[0].chr + ":" + targets[0].start.to_string() + "-" + targets[0].end_.to_string() + " (" + targets[0].strand + ")")
  println("")

  // ------------------------------------------------------------------
  // 2. Basic normalization
  // ------------------------------------------------------------------
  println("=== 2. Basic Normalization ===")
  let config = @bio.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@bio.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @bio.normalize_to_matrix(signals, targets, config)
  println("Matrix dimensions: " + mat.nRows.to_string() + " x " + mat.nCols.to_string())
  println("Upstream windows: " + mat.upstreamIndex.length().to_string())
  println("Downstream windows: " + mat.downstreamIndex.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Normalization with target body
  // ------------------------------------------------------------------
  println("=== 3. With Target Body (k=5) ===")
  let config_t = @bio.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=5,
    meanMode=@bio.mean_mode_absolute(),
    includeTarget=true,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat_t = @bio.normalize_to_matrix(signals, targets, config_t)
  println("Matrix dimensions: " + mat_t.nRows.to_string() + " x " + mat_t.nCols.to_string())
  println("Upstream windows: " + mat_t.upstreamIndex.length().to_string())
  println("Target windows: " + mat_t.targetIndex.length().to_string())
  println("Downstream windows: " + mat_t.downstreamIndex.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Mean mode comparison
  // ------------------------------------------------------------------
  println("=== 4. Mean Mode Comparison ===")
  let modes = [
    ("absolute", @bio.mean_mode_absolute()),
    ("weighted", @bio.mean_mode_weighted()),
    ("w0", @bio.mean_mode_w0()),
    ("coverage", @bio.mean_mode_coverage()),
  ]
  let mut m = 0
  while m < modes.length() {
    let (name, mode) = modes[m]
    let cfg = @bio.EnrichedHeatmapConfig::new(
      extendUp=1000,
      extendDown=1000,
      w=100,
      k=0,
      meanMode=mode,
      includeTarget=false,
      smooth=false,
      smoothWindow=5,
      background=0.0,
      flipUpstream=false,
      keepLow=0.0,
      keepHigh=1.0,
    )
    let res = @bio.normalize_to_matrix(signals, targets, cfg)
    let profile = @bio.enrichment_profile(res)
    let profile_mean = array_mean(profile)
    println("  " + name + ": mean profile value = " + profile_mean.to_string())
    m = m + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Smoothing
  // ------------------------------------------------------------------
  println("=== 5. Row Smoothing ===")
  let config_s = @bio.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@bio.mean_mode_absolute(),
    includeTarget=false,
    smooth=true,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat_s = @bio.normalize_to_matrix(signals, targets, config_s)
  println("Smoothed: " + mat_s.smooth.to_string())
  println("Failed rows: " + mat_s.failedRows.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. Row means
  // ------------------------------------------------------------------
  println("=== 6. Row Means (per-target enrichment) ===")
  let means = @bio.eh_row_means(mat)
  println("Row means (first 5 targets):")
  let mut r = 0
  while r < 5 && r < means.length() {
    println("  Target " + r.to_string() + ": " + means[r].to_string())
    r = r + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Enrichment profile
  // ------------------------------------------------------------------
  println("=== 7. Enrichment Profile (column-wise mean) ===")
  let profile = @bio.enrichment_profile(mat)
  println("Profile length: " + profile.length().to_string())
  println("Profile values (first 5):")
  let mut p = 0
  while p < 5 && p < profile.length() {
    println("  Window " + p.to_string() + ": " + profile[p].to_string())
    p = p + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Minus strand handling
  // ------------------------------------------------------------------
  println("=== 8. Minus Strand Handling ===")
  let m_signals = [
    @bio.GenomicSignal::new("chr1", 9800, 9900, 0.5),
    @bio.GenomicSignal::new("chr1", 10100, 10200, 1.0),
  ]
  let m_targets = [@bio.TargetRegion::new("chr1", 10000, 10000, "-")]
  let m_config = @bio.EnrichedHeatmapConfig::new(
    extendUp=500,
    extendDown=500,
    w=500,
    k=0,
    meanMode=@bio.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let m_mat = @bio.normalize_to_matrix(m_signals, m_targets, m_config)
  println("Minus strand target: " + m_targets[0].chr + ":" + m_targets[0].start.to_string() + " (" + m_targets[0].strand + ")")
  println("Matrix: " + m_mat.nRows.to_string() + " x " + m_mat.nCols.to_string())
  println("Window 0 (upstream, reversed): " + m_mat.matrix[0][0].to_string())
  println("Window 1 (downstream, reversed): " + m_mat.matrix[0][1].to_string())
  println("")

  // ------------------------------------------------------------------
  // 9. Summary
  // ------------------------------------------------------------------
  println("=== 9. Summary ===")
  let summary = @bio.enriched_heatmap_summary(mat)
  println(summary)
  println("=== Demo complete ===")
}

fn array_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}
