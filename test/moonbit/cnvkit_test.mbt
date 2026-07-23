///|
test "cnvkit_create_probe" {
  let probe = @src.CNVProbe::new(
    "probe_001",
    "chr1",
    1000000,
    1001000,
    -0.5,
    1.0
  )
  
  assert_eq(probe.probe_id, "probe_001")
  assert_eq(probe.chromosome, "chr1")
  assert_eq(probe.start, 1000000)
  assert_eq(probe.end, 1001000)
  assert_true((probe.log2_ratio - (-0.5)).abs() < 0.001)
}

///|
test "cnvkit_create_segment" {
  let segment = @src.CNVSegment::new(
    "chr1",
    1000000,
    2000000,
    50,
    -0.8,
    "deletion"
  )
  
  assert_eq(segment.chromosome, "chr1")
  assert_eq(segment.start, 1000000)
  assert_eq(segment.end, 2000000)
  assert_eq(segment.probes, 50)
  assert_eq(segment.state, "deletion")
}

///|
test "cnvkit_dataset_operations" {
  let mut dataset = @src.CNVDataset::new()
  
  assert_eq(dataset.count_probes(), 0)
  
  let probe = @src.CNVProbe::new(
    "probe_001",
    "chr1",
    1000000,
    1001000,
    -0.5,
    1.0
  )
  
  dataset = dataset.add_probe(probe)
  
  assert_eq(dataset.count_probes(), 1)
  assert_eq(dataset.chromosomes.length(), 1)
  assert_eq(dataset.chromosomes[0], "chr1")
}

///|
test "cnvkit_filter_chromosome" {
  let mut dataset = @src.CNVDataset::new()
  
  let probe1 = @src.CNVProbe::new("probe_001", "chr1", 1000000, 1001000, -0.5, 1.0)
  let probe2 = @src.CNVProbe::new("probe_002", "chr2", 2000000, 2001000, 0.5, 1.0)
  
  dataset = dataset.add_probe(probe1)
  dataset = dataset.add_probe(probe2)
  
  let filtered = dataset.filter_chromosome("chr1")
  
  assert_eq(filtered.count_probes(), 1)
  assert_eq(filtered.chromosomes.length(), 1)
}

///|
test "cnvkit_cbs_segmentation" {
  let probes : Array[@src.CNVProbe] = Array::new()
  
  // Create probes with a clear change point
  let mut i = 0
  while i < 20 {
    let ratio = if i < 10 { -0.8 } else { 0.1 }
    probes.push(@src.CNVProbe::new(
      "probe_" + i.to_string(),
      "chr1",
      1000000 + i * 1000,
      1000000 + i * 1000 + 100,
      ratio,
      1.0
    ))
    i = i + 1
  }
  
  let result = @src.cbs_segment(probes, 0.05)
  
  assert_true(result.n_segments > 0)
  assert_true(result.breakpoints.length() > 0)
}

///|
test "cnvkit_smooth_log2_ratios" {
  let probes : Array[@src.CNVProbe] = Array::new()
  
  let mut i = 0
  while i < 10 {
    probes.push(@src.CNVProbe::new(
      "probe_" + i.to_string(),
      "chr1",
      1000000 + i * 1000,
      1000000 + i * 1000 + 100,
      (i.to_double() - 5.0) * 0.1,
      1.0
    ))
    i = i + 1
  }
  
  let smoothed = @src.smooth_log2_ratios(probes, 3)
  
  assert_eq(smoothed.length(), 10)
  
  // Check that smoothed values are reasonable
  let val = smoothed[0].log2_ratio
  assert_true(val.abs() < 10.0)
}

///|
test "cnvkit_call_state" {
  assert_eq(@src.cbs_call_state(-0.5), "deletion")
  assert_eq(@src.cbs_call_state(0.5), "amplification")
  assert_eq(@src.cbs_call_state(0.0), "neutral")
}

///|
test "cnvkit_detect_breakpoints" {
  let segments : Array[@src.CNVSegment] = [
    @src.CNVSegment::new("chr1", 1000000, 2000000, 50, -0.8, "deletion"),
    @src.CNVSegment::new("chr1", 2000000, 3000000, 50, 0.1, "neutral"),
    @src.CNVSegment::new("chr1", 3000000, 4000000, 50, 0.9, "amplification")
  ]
  
  let breakpoints = @src.detect_breakpoints(segments, 0.3)
  
  assert_true(breakpoints.length() > 0)
}

///|
test "cnvkit_call_copy_numbers" {
  let segments : Array[@src.CNVSegment] = [
    @src.CNVSegment::new("chr1", 1000000, 2000000, 50, -0.8, "deletion"),
    @src.CNVSegment::new("chr1", 2000000, 3000000, 50, 0.0, "neutral"),
    @src.CNVSegment::new("chr1", 3000000, 4000000, 50, 0.8, "amplification")
  ]
  
  let calls = @src.call_copy_numbers(segments, 2)
  
  assert_eq(calls.length(), 3)
  assert_eq(calls[0].state, "deletion")
  assert_eq(calls[1].state, "neutral")
  assert_eq(calls[2].state, "amplification")
}

///|
test "cnvkit_summarize_dataset" {
  let dataset = @src.create_example_cnv_dataset()
  let summary = @src.summarize_cnv(dataset)
  
  assert_true(summary.total_probes > 0)
  assert_true(summary.total_segments >= 0)
}

///|
test "cnvkit_create_example" {
  let dataset = @src.create_example_cnv_dataset()
  
  assert_true(dataset.count_probes() > 0)
  assert_true(dataset.chromosomes.length() >= 2)
}

///|
test "cnvkit_summarize_string" {
  let dataset = @src.create_example_cnv_dataset()
  let summary_str = @src.cbs_summarize_dataset(dataset)
  
  assert_true(summary_str.contains("CNV Dataset Summary"))
  assert_true(summary_str.contains("Total probes:"))
}

///|
test "cnvkit_empty_segmentation" {
  let empty_probes : Array[@src.CNVProbe] = Array::new()
  let result = @src.cbs_segment(empty_probes, 0.05)
  
  assert_eq(result.n_segments, 0)
  assert_eq(result.breakpoints.length(), 0)
}

///|
test "cnvkit_single_probe" {
  let probe = @src.CNVProbe::new("p1", "chr1", 100, 200, 0.5, 1.0)
  let probes : Array[@src.CNVProbe] = [probe]
  
  let result = @src.cbs_segment(probes, 0.05)
  
  assert_eq(result.n_segments, 1)
}

///|
test "cnvkit_smooth_edge_cases" {
  let empty_probes : Array[@src.CNVProbe] = Array::new()
  let smoothed = @src.smooth_log2_ratios(empty_probes, 3)
  
  assert_eq(smoothed.length(), 0)
}
