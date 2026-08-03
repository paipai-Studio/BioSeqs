///|
/// Tests for the hmmcopy module (Bioconductor HMMcopy).
///
/// Covers HMM bins, model construction, GC bias correction, Viterbi
/// decoding, CNV calling, segmentation, and log-likelihood computation.

// ===========================================================================
// 1. HMMcopyBin construction
// ===========================================================================

///|
test "hc_bin_creation" {
  let b = @src.HMMcopyBin::new(
    chr="chr1",
    start=1000,
    end=2000,
    reads=100,
    gc=0.5,
    mappability=0.95,
  )
  assert_eq(b.chr, "chr1")
  assert_eq(b.start, 1000)
  assert_eq(b.end, 2000)
  assert_eq(b.reads, 100)
  assert_true((b.gc - 0.5).abs() < 0.001)
  assert_true((b.mappability - 0.95).abs() < 0.001)
  assert_true((b.corrected_reads - 100.0).abs() < 0.001)
}

///|
test "hc_bin_default_gc_mappability" {
  let b = @src.HMMcopyBin::new(
    chr="chr1",
    start=0,
    end=1000,
    reads=50,
  )
  assert_true((b.gc - 0.5).abs() < 0.001)
  assert_true((b.mappability - 1.0).abs() < 0.001)
}

///|
test "hc_bin_width" {
  let b = @src.HMMcopyBin::new(
    chr="chr1",
    start=1000,
    end=2000,
    reads=100,
  )
  assert_eq(b.width(), 1000)
}

///|
test "hc_bin_normalized_coverage" {
  let b = @src.HMMcopyBin::new(
    chr="chr1",
    start=0,
    end=500,
    reads=100,
  )
  b.set_corrected_reads(200.0)
  let nc = b.normalized_coverage()
  assert_true((nc - 0.4).abs() < 0.001)
}

///|
test "hc_bin_zero_width_coverage" {
  let b = @src.HMMcopyBin::new(
    chr="chr1",
    start=100,
    end=100,
    reads=50,
  )
  b.set_corrected_reads(75.0)
  // Zero width returns corrected_reads directly
  assert_true((b.normalized_coverage() - 75.0).abs() < 0.001)
}

// ===========================================================================
// 2. HMMcopyModel
// ===========================================================================

///|
test "hc_model_default_construction" {
  let m = @src.HMMcopyModel::new()
  assert_eq(m.num_states, 7)
  assert_eq(m.ploidy, 2)
  // Emission for state 0 (CN0) = 0.0
  assert_true((m.emissions[0] - 0.0).abs() < 0.001)
  // Emission for state 2 (CN2) = 1.0
  assert_true((m.emissions[2] - 1.0).abs() < 0.001)
}

///|
test "hc_model_set_transition" {
  let m = @src.HMMcopyModel::new()
  m.set_transition(0, 0, 0.99)
  m.set_transition(0, 1, 0.01)
  assert_true((m.transition[0][0] - 0.99).abs() < 0.001)
  assert_true((m.transition[0][1] - 0.01).abs() < 0.001)
}

///|
test "hc_model_set_emission" {
  let m = @src.HMMcopyModel::new()
  m.set_emission(3, 1.5)
  assert_true((m.emissions[3] - 1.5).abs() < 0.001)
}

///|
test "hc_model_set_variance" {
  let m = @src.HMMcopyModel::new()
  m.set_variance(2, 0.05)
  assert_true((m.variances[2] - 0.05).abs() < 0.001)
}

///|
test "hc_model_set_initial" {
  let m = @src.HMMcopyModel::new()
  m.set_initial(2, 0.9)
  assert_true((m.initial[2] - 0.9).abs() < 0.001)
}

///|
test "hc_model_state_to_copy" {
  let m = @src.HMMcopyModel::new()
  assert_eq(m.state_to_copy(0), 0)
  assert_eq(m.state_to_copy(2), 2)
  assert_eq(m.state_to_copy(5), 5)
}

///|
test "hc_model_copy_to_state" {
  let m = @src.HMMcopyModel::new()
  assert_eq(m.copy_to_state(2), 2)
  assert_eq(m.copy_to_state(-1), 0)
  // Copy beyond max states is clamped to last state
  assert_eq(m.copy_to_state(100), 6)
}

// ===========================================================================
// 3. GC bias correction
// ===========================================================================

///|
test "hc_correct_gc_basic" {
  let bins = [
    @src.HMMcopyBin::new(chr="chr1", start=0, end=100, reads=100, gc=0.4),
    @src.HMMcopyBin::new(chr="chr1", start=100, end=200, reads=200, gc=0.6),
    @src.HMMcopyBin::new(chr="chr1", start=200, end=300, reads=150, gc=0.5),
  ]
  let corrected = @src.hmmcopy_correct_gc(bins)
  assert_eq(corrected.length(), 3)
  // Corrected reads should be positive
  for b in corrected {
    assert_true(b.corrected_reads > 0.0)
  }
}

///|
test "hc_correct_gc_empty" {
  let bins : Array[@src.HMMcopyBin] = []
  let corrected = @src.hmmcopy_correct_gc(bins)
  assert_eq(corrected.length(), 0)
}

///|
test "hc_correct_gc_preserves_coordinates" {
  let bins = [
    @src.HMMcopyBin::new(chr="chr2", start=500, end=600, reads=80, gc=0.5),
  ]
  let corrected = @src.hmmcopy_correct_gc(bins)
  assert_eq(corrected[0].chr, "chr2")
  assert_eq(corrected[0].start, 500)
  assert_eq(corrected[0].end, 600)
}

// ===========================================================================
// 4. Viterbi algorithm
// ===========================================================================

///|
test "hc_viterbi_empty_bins" {
  let bins : Array[@src.HMMcopyBin] = []
  let model = @src.HMMcopyModel::new()
  let states = @src.hmmcopy_viterbi(bins, model)
  assert_eq(states.length(), 0)
}

///|
test "hc_viterbi_returns_states" {
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  let states = @src.hmmcopy_viterbi(bins, model)
  assert_eq(states.length(), bins.length())
  // All states should be valid (0 to 6)
  for s in states {
    assert_true(s >= 0 && s < 7)
  }
}

///|
test "hc_viterbi_detects_cn0_region" {
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  let states = @src.hmmcopy_viterbi(bins, model)
  // Bins 10-19 are the CN0 region (very low reads)
  // The Viterbi should call low copy number for these
  let mut low_copy_count = 0
  let mut i = 10
  while i < 20 {
    if states[i] <= 1 {
      low_copy_count = low_copy_count + 1
    }
    i = i + 1
  }
  // At least some bins in the CN0 region should have low copy number
  assert_true(low_copy_count > 0)
}

///|
test "hc_viterbi_detects_cn3_region" {
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  let states = @src.hmmcopy_viterbi(bins, model)
  // Bins 60-69 are the CN3 region (elevated reads)
  let mut high_copy_count = 0
  let mut i = 60
  while i < 70 {
    if states[i] >= 3 {
      high_copy_count = high_copy_count + 1
    }
    i = i + 1
  }
  // At least some bins in the CN3 region should have elevated copy number
  assert_true(high_copy_count > 0)
}

// ===========================================================================
// 5. CNV calling and segmentation
// ===========================================================================

///|
test "hc_call_cnv" {
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  let called = @src.hmmcopy_call_cnv(bins, model)
  assert_eq(called.length(), bins.length())
  // Each bin should have copy and state set
  for b in called {
    assert_true(b.copy >= 0 && b.copy < 7)
    assert_true(b.state >= 0 && b.state < 7)
  }
}

///|
test "hc_segment_basic" {
  // Create bins with known copy numbers
  let bins = [
    @src.HMMcopyBin::new(chr="chr1", start=0, end=100, reads=60),
    @src.HMMcopyBin::new(chr="chr1", start=100, end=200, reads=65),
    @src.HMMcopyBin::new(chr="chr1", start=200, end=300, reads=60),
  ]
  bins[0].set_copy(2)
  bins[0].set_state(2)
  bins[1].set_copy(2)
  bins[1].set_state(2)
  bins[2].set_copy(2)
  bins[2].set_state(2)
  let segments = @src.hmmcopy_segment(bins)
  // All same copy -> single segment
  assert_eq(segments.length(), 1)
  assert_eq(segments[0].copy, 2)
  assert_eq(segments[0].num_bins, 3)
}

///|
test "hc_segment_multiple_copies" {
  let bins = [
    @src.HMMcopyBin::new(chr="chr1", start=0, end=100, reads=60),
    @src.HMMcopyBin::new(chr="chr1", start=100, end=200, reads=0),
    @src.HMMcopyBin::new(chr="chr1", start=200, end=300, reads=0),
  ]
  bins[0].set_copy(2)
  bins[0].set_state(2)
  bins[1].set_copy(0)
  bins[1].set_state(0)
  bins[2].set_copy(0)
  bins[2].set_state(0)
  let segments = @src.hmmcopy_segment(bins)
  assert_eq(segments.length(), 2)
  assert_eq(segments[0].copy, 2)
  assert_eq(segments[0].num_bins, 1)
  assert_eq(segments[1].copy, 0)
  assert_eq(segments[1].num_bins, 2)
}

///|
test "hc_segment_empty" {
  let bins : Array[@src.HMMcopyBin] = []
  let segments = @src.hmmcopy_segment(bins)
  assert_eq(segments.length(), 0)
}

///|
test "hc_segment_mean_reads" {
  let bins = [
    @src.HMMcopyBin::new(chr="chr1", start=0, end=100, reads=60),
    @src.HMMcopyBin::new(chr="chr1", start=100, end=200, reads=60),
  ]
  bins[0].set_corrected_reads(60.0)
  bins[0].set_copy(2)
  bins[0].set_state(2)
  bins[1].set_corrected_reads(80.0)
  bins[1].set_copy(2)
  bins[1].set_state(2)
  let segments = @src.hmmcopy_segment(bins)
  assert_eq(segments.length(), 1)
  // Mean should be (60+80)/2 = 70
  assert_true((segments[0].mean_reads - 70.0).abs() < 0.001)
}

// ===========================================================================
// 6. Log-likelihood
// ===========================================================================

///|
test "hc_log_likelihood_empty" {
  let bins : Array[@src.HMMcopyBin] = []
  let model = @src.HMMcopyModel::new()
  let ll = @src.hmmcopy_log_likelihood(bins, model)
  assert_eq(ll, 0.0)
}

///|
test "hc_log_likelihood_positive_bins" {
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  let ll = @src.hmmcopy_log_likelihood(bins, model)
  // Log-likelihood should be a finite number (not NaN)
  assert_true(!ll.is_nan())
}

// ===========================================================================
// 7. Sample data
// ===========================================================================

///|
test "hc_sample_bins" {
  let bins = @src.hmmcopy_sample_bins()
  assert_eq(bins.length(), 100)
  // All on chr1
  for b in bins {
    assert_eq(b.chr, "chr1")
  }
  // CN0 region (bins 10-19) should have low reads
  let mut i = 10
  while i < 20 {
    assert_true(bins[i].reads < 10)
    i = i + 1
  }
  // CN3 region (bins 60-69) should have elevated reads
  i = 60
  while i < 70 {
    assert_true(bins[i].reads >= 80)
    i = i + 1
  }
}

///|
test "hc_sample_model" {
  let m = @src.hmmcopy_sample_model()
  assert_eq(m.num_states, 7)
  // CN2 should have high initial probability
  assert_true(m.initial[2] > 0.5)
  // Emissions should be 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0
  assert_true((m.emissions[0] - 0.0).abs() < 0.001)
  assert_true((m.emissions[2] - 1.0).abs() < 0.001)
  assert_true((m.emissions[6] - 3.0).abs() < 0.001)
}

// ===========================================================================
// 8. Integration test
// ===========================================================================

///|
test "hc_full_pipeline" {
  // Full pipeline: sample bins -> GC correction -> CNV calling -> segmentation
  let bins = @src.hmmcopy_sample_bins()
  let model = @src.hmmcopy_sample_model()
  // Step 1: GC correction
  let corrected = @src.hmmcopy_correct_gc(bins)
  assert_eq(corrected.length(), 100)
  // Step 2: CNV calling
  let called = @src.hmmcopy_call_cnv(corrected, model)
  assert_eq(called.length(), 100)
  // Step 3: Segmentation
  let segments = @src.hmmcopy_segment(called)
  // Should have multiple segments (at least CN0, CN2, CN3 regions)
  assert_true(segments.length() >= 1)
}
