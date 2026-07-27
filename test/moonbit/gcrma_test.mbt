///|
/// Tests for the GCRMA (GC-content adjusted RMA) module.

///|
test "gcrma_gcrma_config_new" {
  let config = @src.GCRMAConfig::new()
  assert_eq(config.background_method, "IdealMM")
  assert_eq(config.normalize, true)
}

///|
test "gcrma_probe_info_new" {
  let probe = @src.ProbeInfo::new("probe_1_at", "ACGTACGTACGTACGT")
  assert_eq(probe.probe_id, "probe_1_at")
  assert_eq(probe.sequence.length(), 16)
  assert_eq(probe.gc_count > 0, true)
}

///|
test "gcrma_probe_info_with_values" {
  let probe = @src.ProbeInfo::with_values(
    "probe_test", 8, "ACGTACGT", 2.5,
  )
  assert_eq(probe.probe_id, "probe_test")
  assert_eq(probe.gc_count, 8)
  assert_eq(probe.affinity, 2.5)
}

///|
test "gcrma_oligo_gc_count_empty" {
  let result = @src.oligo_gc_count("")
  assert_eq(result, 0)
}

///|
test "gcrma_oligo_gc_count_all_gc" {
  let result = @src.oligo_gc_count("GCGCGCGCGCGC")
  assert_eq(result, 12)
}

///|
test "gcrma_oligo_gc_count_all_at" {
  let result = @src.oligo_gc_count("ATATATATATAT")
  assert_eq(result, 0)
}

///|
test "gcrma_oligo_gc_count_mixed" {
  let result = @src.oligo_gc_count("ACGTACGT")
  assert_eq(result, 4)
}

///|
test "gcrma_oligo_gc_count_lowercase" {
  let result = @src.oligo_gc_count("acgtacgt")
  assert_eq(result, 4)
}

///|
test "gcrma_oligo_gc_fraction" {
  let result = @src.oligo_gc_fraction("ACGTACGT")
  assert_eq(result, 0.5)
}

///|
test "gcrma_oligo_gc_fraction_empty" {
  let result = @src.oligo_gc_fraction("")
  assert_eq(result, 0.0)
}

///|
test "gcrma_estimate_affinity_empty" {
  let result = @src.estimate_affinity("")
  assert_eq(result, 0.0)
}

///|
test "gcrma_estimate_affinity_gc_rich" {
  let result = @src.estimate_affinity("GCGCGCGCGCGC")
  assert_eq(result > 2.0, true)
}

///|
test "gcrma_estimate_affinity_at_rich" {
  let result = @src.estimate_affinity("ATATATATATAT")
  assert_eq(result <= 2.0, true)
}

///|
test "gcrma_estimate_affinity_mixed" {
  let result = @src.estimate_affinity("ACGTACGTACGT")
  assert_eq(result > 0.0, true)
}

///|
test "gcrma_compute_gc_lookup_table_empty" {
  let table = @src.compute_gc_lookup_table([])
  assert_eq(table.keys().length(), 0)
}

///|
test "gcrma_compute_gc_lookup_table_basic" {
  let probes = [
    @src.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
    @src.ProbeInfo::new("p3", "ATATATATATAT"),
  ]
  let table = @src.compute_gc_lookup_table(probes)
  assert_eq(table.keys().length() > 0, true)
}

///|
test "gcrma_background_correction_express" {
  let pm = [100.0, 200.0, 150.0, 300.0]
  let mm = [50.0, 80.0, 70.0, 120.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.5, "Express",
  )
  assert_eq(result.length(), 4)
  assert_eq(result[0], 50.0)
  assert_eq(result[1], 120.0)
}

///|
test "gcrma_background_correction_express_negative" {
  let pm = [30.0, 50.0]
  let mm = [80.0, 100.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.5, "Express",
  )
  assert_eq(result.length(), 2)
  assert_eq(result[0], 1.0)
}

///|
test "gcrma_background_correction_idealm" {
  let pm = [100.0, 200.0, 150.0, 300.0]
  let mm = [50.0, 80.0, 70.0, 120.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.5, "IdealMM",
  )
  assert_eq(result.length(), 4)
  assert_eq(result[0] > 1.0, true)
}

///|
test "gcrma_background_correction_empty" {
  let result = @src.gcrma_background_correction(
    [], [], 0.5, "IdealMM",
  )
  assert_eq(result.length(), 0)
}

///|
test "gcrma_gc_correction_basic" {
  let pm_values = [100.0, 200.0, 150.0]
  let gc_counts = [4, 6, 8]
  let probes = [
    @src.ProbeInfo::new("p1", "ACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCATAT"),
    @src.ProbeInfo::new("p3", "GCGCGCGC"),
  ]
  let result = @src.gcrma_gc_correction(pm_values, gc_counts, probes)
  assert_eq(result.length(), 3)
}

///|
test "gcrma_gc_correction_empty" {
  let result = @src.gcrma_gc_correction([], [], [])
  assert_eq(result.length(), 0)
}

///|
test "gcrma_normalize_single_row" {
  let data = [[100.0, 200.0, 150.0]]
  let result = @src.gcrma_normalize(data)
  assert_eq(result.length(), 1)
  assert_eq(result[0].length(), 3)
}

///|
test "gcrma_normalize_multiple_rows" {
  let data = [
    [100.0, 200.0],
    [150.0, 250.0],
    [120.0, 180.0],
  ]
  let result = @src.gcrma_normalize(data)
  assert_eq(result.length(), 3)
  assert_eq(result[0].length(), 2)
}

///|
test "gcrma_normalize_empty" {
  let result = @src.gcrma_normalize([])
  assert_eq(result.length(), 0)
}

///|
test "gcrma_summarize_basic" {
  let expression = [
    [100.0, 200.0],
    [150.0, 250.0],
    [120.0, 180.0],
    [180.0, 220.0],
  ]
  let probe_sets : Map[String, Array[Int]] = Map([], capacity=0)
  probe_sets.set("gene1_at", [0, 1, 2, 3])
  let result = @src.gcrma_summarize(expression, probe_sets)
  assert_eq(result.keys().length(), 1)
  let vals = result.get_or_default("gene1_at", [])
  assert_eq(vals.length(), 2)
}

///|
test "gcrma_summarize_empty_sets" {
  let expression = [[100.0, 200.0]]
  let probe_sets : Map[String, Array[Int]] = Map([], capacity=0)
  probe_sets.set("gene1_at", [])
  let result = @src.gcrma_summarize(expression, probe_sets)
  let vals = result.get_or_default("gene1_at", [])
  assert_eq(vals.length(), 0)
}

///|
test "gcrma_summarize_multiple_probesets" {
  let expression = [
    [100.0, 200.0],
    [150.0, 250.0],
    [120.0, 180.0],
    [180.0, 220.0],
    [90.0, 160.0],
  ]
  let probe_sets : Map[String, Array[Int]] = Map([], capacity=0)
  probe_sets.set("gene1_at", [0, 1, 2])
  probe_sets.set("gene2_at", [3, 4])
  let result = @src.gcrma_summarize(expression, probe_sets)
  assert_eq(result.keys().length(), 2)
  let g1 = result.get_or_default("gene1_at", [])
  let g2 = result.get_or_default("gene2_at", [])
  assert_eq(g1.length(), 2)
  assert_eq(g2.length(), 2)
}

///|
test "gcrma_process_basic" {
  let cel_data = [
    [100.0, 200.0, 150.0],
    [150.0, 250.0, 180.0],
    [120.0, 180.0, 160.0],
    [180.0, 220.0, 190.0],
  ]
  let probe_info = [
    @src.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
    @src.ProbeInfo::new("p3", "ATATATATATAT"),
    @src.ProbeInfo::new("p4", "ACACACACACAC"),
  ]
  let result = @src.gcrma_process(cel_data, probe_info)
  assert_eq(result.expression_matrix.length(), 4)
  assert_eq(result.expression_matrix[0].length(), 3)
}

///|
test "gcrma_process_empty" {
  let result = @src.gcrma_process([], [])
  assert_eq(result.expression_matrix.length(), 0)
}

///|
test "gcrma_process_with_config" {
  let cel_data = [
    [100.0, 200.0],
    [150.0, 250.0],
  ]
  let probe_info = [
    @src.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
  ]
  let config = @src.GCRMAConfig::new()
  let result = @src.gcrma_process_with_config(cel_data, probe_info, config)
  assert_eq(result.expression_matrix.length(), 2)
  assert_eq(result.expression_matrix[0].length(), 2)
}

///|
test "gcrma_process_no_normalize" {
  let cel_data = [
    [100.0, 200.0],
    [150.0, 250.0],
  ]
  let probe_info = [
    @src.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
  ]
  let config = @src.GCRMAConfig::new(normalize=false)
  let result = @src.gcrma_process_with_config(cel_data, probe_info, config)
  assert_eq(result.expression_matrix.length(), 2)
}

///|
test "gcrma_background_correction_mm_longer_than_pm" {
  let pm = [100.0, 200.0]
  let mm = [50.0, 80.0, 90.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.5, "Express",
  )
  assert_eq(result.length(), 2)
  assert_eq(result[0], 50.0)
  assert_eq(result[1], 120.0)
}

///|
test "gcrma_background_correction_mm_shorter_than_pm" {
  let pm = [100.0, 200.0, 150.0]
  let mm = [50.0, 80.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.5, "Express",
  )
  assert_eq(result.length(), 3)
  assert_eq(result[0], 50.0)
  assert_eq(result[2], 150.0)
}

///|
test "gcrma_background_correction_idealm_clamps_to_minimum" {
  let pm = [10.0, 20.0]
  let mm = [8.0, 15.0]
  let result = @src.gcrma_background_correction(
    pm, mm, 0.9, "IdealMM",
  )
  assert_eq(result.length(), 2)
  assert_eq(result[0] >= 1.0, true)
  assert_eq(result[1] >= 1.0, true)
}

///|
test "gcrma_gc_correction_all_same_gc" {
  let pm_values = [100.0, 200.0, 150.0]
  let gc_counts = [4, 4, 4]
  let probes = [
    @src.ProbeInfo::new("p1", "ACGTACGT"),
    @src.ProbeInfo::new("p2", "TGCAACGT"),
    @src.ProbeInfo::new("p3", "GCATACGT"),
  ]
  let result = @src.gcrma_gc_correction(pm_values, gc_counts, probes)
  assert_eq(result.length(), 3)
}

///|
test "gcrma_multiple_probeset_summary" {
  let expression = [
    [100.0, 200.0, 300.0],
    [150.0, 250.0, 350.0],
    [120.0, 220.0, 320.0],
  ]
  let probe_sets : Map[String, Array[Int]] = Map([], capacity=0)
  probe_sets.set("geneA", [0, 1, 2])
  let result = @src.gcrma_summarize(expression, probe_sets)
  let vals = result.get_or_default("geneA", [])
  assert_eq(vals.length(), 3)
  assert_eq(vals[0] > 0.0, true)
}

///|
test "gcrma_config_background_method" {
  let config = @src.GCRMAConfig::new(background_method="Express")
  assert_eq(config.background_method, "Express")
}

///|
test "gcrma_full_pipeline_express" {
  let cel_data = [
    [100.0, 200.0],
    [150.0, 250.0],
    [120.0, 180.0],
    [180.0, 220.0],
    [90.0, 160.0],
    [110.0, 190.0],
  ]
  let probe_info = [
    @src.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @src.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
    @src.ProbeInfo::new("p3", "ATATATATATAT"),
    @src.ProbeInfo::new("p4", "ACACACACACAC"),
    @src.ProbeInfo::new("p5", "TGTGTGTGTGTG"),
    @src.ProbeInfo::new("p6", "CGCGCGCGCGCG"),
  ]
  let config = @src.GCRMAConfig::new(background_method="Express", normalize=true)
  let result = @src.gcrma_process_with_config(cel_data, probe_info, config)
  assert_eq(result.expression_matrix.length(), 6)
  assert_eq(result.expression_matrix[0].length(), 2)
  let params = result.gc_correction_params
  assert_eq(params.get_or_default("n_probes", 0.0), 6.0)
}

///|
test "gcrma_oligo_gc_count_different_lengths" {
  let short = @src.oligo_gc_count("ACGT")
  let long = @src.oligo_gc_count("ACGTACGTACGTACGT")
  assert_eq(short, 2)
  assert_eq(long, 8)
}

///|
test "gcrma_estimate_affinity_dinucleotide" {
  let result = @src.estimate_affinity("GC")
  assert_eq(result > 0.0, true)
}

///|
test "gcrma_background_correction_high_gc" {
  let pm = [100.0, 200.0, 150.0]
  let mm = [60.0, 100.0, 80.0]
  let result_low = @src.gcrma_background_correction(
    pm, mm, 0.2, "IdealMM",
  )
  let result_high = @src.gcrma_background_correction(
    pm, mm, 0.8, "IdealMM",
  )
  assert_eq(result_low.length(), 3)
  assert_eq(result_high.length(), 3)
}

///|
test "gcrma_gc_correction_preserves_length" {
  let pm_values = [50.0, 60.0, 70.0, 80.0, 90.0]
  let gc_counts = [2, 3, 4, 5, 6]
  let probes = [
    @src.ProbeInfo::new("p1", "ACGTACGT"),
    @src.ProbeInfo::new("p2", "ACGTACGC"),
    @src.ProbeInfo::new("p3", "GCGCGCGC"),
    @src.ProbeInfo::new("p4", "GCGCGCAC"),
    @src.ProbeInfo::new("p5", "GCGCGCGG"),
  ]
  let result = @src.gcrma_gc_correction(pm_values, gc_counts, probes)
  assert_eq(result.length(), 5)
}