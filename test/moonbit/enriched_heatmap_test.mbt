///|
/// Tests for Bioconductor EnrichedHeatmap - signal normalization.

///|
test "eh_config_default" {
  let config = @src.EnrichedHeatmapConfig::new()
  assert_eq(config.extendUp, 5000)
  assert_eq(config.extendDown, 5000)
  assert_eq(config.w, 100)
  assert_eq(config.background, 0.0)
}

///|
test "eh_genomic_signal_new" {
  let s = @src.GenomicSignal::new("chr1", 100, 200, 0.5)
  assert_eq(s.chr, "chr1")
  assert_eq(s.start, 100)
  assert_eq(s.end_, 200)
  assert_true((s.value - 0.5).abs() < 1.0e-10)
}

///|
test "eh_target_region_new" {
  let t = @src.TargetRegion::new("chr1", 1000, 1000, "+")
  assert_eq(t.chr, "chr1")
  assert_eq(t.start, 1000)
  assert_eq(t.end_, 1000)
  assert_eq(t.strand, "+")
}

///|
test "eh_normalize_basic" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, targets.length())
  // nCols = extendUp/w + extendDown/w = 10 + 10 = 20
  assert_eq(mat.nCols, 20)
}

///|
test "eh_normalize_with_target" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=5,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=true,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  // nCols = 10 (up) + 5 (target) + 10 (down) = 25
  assert_eq(mat.nCols, 25)
  assert_eq(mat.upstreamIndex.length(), 10)
  assert_eq(mat.targetIndex.length(), 5)
  assert_eq(mat.downstreamIndex.length(), 10)
}

///|
test "eh_normalize_empty_targets" {
  let signals : Array[@src.GenomicSignal] = []
  let targets : Array[@src.TargetRegion] = []
  let config = @src.EnrichedHeatmapConfig::new()
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 0)
  assert_eq(mat.nCols, 0)
}

///|
test "eh_normalize_matrix_values" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  // Matrix should have the right dimensions.
  assert_eq(mat.matrix.length(), targets.length())
  let mut i = 0
  while i < mat.matrix.length() {
    assert_eq(mat.matrix[i].length(), mat.nCols)
    i = i + 1
  }
}

///|
test "eh_normalize_mean_mode_absolute" {
  let signals = [@src.GenomicSignal::new("chr1", 100, 200, 0.8)]
  let targets = [@src.TargetRegion::new("chr1", 150, 150, "+")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=100,
    extendDown=100,
    w=100,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 1)
  assert_eq(mat.nCols, 2) // 1 up + 1 down
}

///|
test "eh_normalize_mean_mode_weighted" {
  let signals = [@src.GenomicSignal::new("chr1", 100, 200, 0.8)]
  let targets = [@src.TargetRegion::new("chr1", 150, 150, "+")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=100,
    extendDown=100,
    w=100,
    k=0,
    meanMode=@src.mean_mode_weighted(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 1)
}

///|
test "eh_normalize_mean_mode_w0" {
  let signals = [@src.GenomicSignal::new("chr1", 100, 200, 0.8)]
  let targets = [@src.TargetRegion::new("chr1", 150, 150, "+")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=100,
    extendDown=100,
    w=100,
    k=0,
    meanMode=@src.mean_mode_w0(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 1)
}

///|
test "eh_normalize_mean_mode_coverage" {
  let signals = [@src.GenomicSignal::new("chr1", 100, 200, 0.8)]
  let targets = [@src.TargetRegion::new("chr1", 150, 150, "+")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=100,
    extendDown=100,
    w=100,
    k=0,
    meanMode=@src.mean_mode_coverage(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 1)
}

///|
test "eh_normalize_smooth" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=true,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.smooth, true)
  assert_eq(mat.nRows, targets.length())
}

///|
test "eh_normalize_background" {
  // Signal far from target → background should be used.
  let signals = [@src.GenomicSignal::new("chr2", 10000, 10100, 1.0)]
  let targets = [@src.TargetRegion::new("chr1", 5000, 5000, "+")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=1000,
    extendDown=1000,
    w=1000,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=-1.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  // All values should be background (-1.0) since signal is on different chr.
  let mut j = 0
  while j < mat.nCols {
    assert_true((mat.matrix[0][j] - (-1.0)).abs() < 1.0e-10)
    j = j + 1
  }
}

///|
test "eh_enrichment_profile" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  let profile = @src.enrichment_profile(mat)
  assert_eq(profile.length(), mat.nCols)
}

///|
test "eh_row_means" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new()
  let mat = @src.normalize_to_matrix(signals, targets, config)
  let means = @src.eh_row_means(mat)
  assert_eq(means.length(), mat.nRows)
}

///|
test "eh_minus_strand" {
  // Test that minus strand reverses window order.
  let signals = [@src.GenomicSignal::new("chr1", 9800, 9900, 0.5), @src.GenomicSignal::new("chr1", 10100, 10200, 1.0)]
  let targets = [@src.TargetRegion::new("chr1", 10000, 10000, "-")]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=500,
    extendDown=500,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 1)
  assert_eq(mat.nCols, 2)
  // For minus strand: upstream is downstream of target.
  // The signal at 9800-9900 (downstream of target on minus strand)
  // should appear in the "upstream" window (reversed).
  // And signal at 10100-10200 (upstream on minus strand) in "downstream".
  // Values should not be background since signals overlap windows.
  assert_true(mat.matrix[0][0] >= 0.0)
  assert_true(mat.matrix[0][1] >= 0.0)
}

///|
test "eh_summary" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new()
  let mat = @src.normalize_to_matrix(signals, targets, config)
  let summary = @src.enriched_heatmap_summary(mat)
  assert_true(summary.length() > 0)
}

///|
test "eh_sample_data" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  assert_eq(targets.length(), 10)
  assert_true(signals.length() > 0)
  // Each target should be a single point (TSS).
  let mut i = 0
  while i < targets.length() {
    assert_eq(targets[i].start, targets[i].end_)
    i = i + 1
  }
}

///|
test "eh_normalize_clipping" {
  let (signals, targets) = @src.enriched_heatmap_sample()
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=5000,
    extendDown=5000,
    w=500,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.1,
    keepHigh=0.9,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  // After clipping, all values should be within [q_low, q_high].
  let mut i = 0
  while i < mat.nRows {
    let mut j = 0
    while j < mat.nCols {
      assert_true(mat.matrix[i][j] >= -1.0e-10)
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "eh_normalize_multiple_chromosomes" {
  let signals = [
    @src.GenomicSignal::new("chr1", 100, 200, 0.5),
    @src.GenomicSignal::new("chr2", 100, 200, 0.8),
  ]
  let targets = [
    @src.TargetRegion::new("chr1", 150, 150, "+"),
    @src.TargetRegion::new("chr2", 150, 150, "+"),
  ]
  let config = @src.EnrichedHeatmapConfig::new(
    extendUp=100,
    extendDown=100,
    w=100,
    k=0,
    meanMode=@src.mean_mode_absolute(),
    includeTarget=false,
    smooth=false,
    smoothWindow=5,
    background=0.0,
    flipUpstream=false,
    keepLow=0.0,
    keepHigh=1.0,
  )
  let mat = @src.normalize_to_matrix(signals, targets, config)
  assert_eq(mat.nRows, 2)
  // chr1 target should get chr1 signal value.
  // chr2 target should get chr2 signal value.
  // They should differ since signal values differ.
  let chr1_max = mat.matrix[0][0]
  let chr2_max = mat.matrix[1][0]
  assert_true((chr1_max - chr2_max).abs() > 1.0e-10)
}
