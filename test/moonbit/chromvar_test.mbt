///|
test "chromvar_create_example_peak_set" {
  let peak_set = @src.chromvar_create_example_peak_set()
  assert_eq(peak_set.peaks.length(), 100)
  assert_eq(peak_set.names.length(), 100)
}

///|
test "chromvar_create_example_motifs" {
  let motifs = @src.chromvar_create_example_motifs()
  assert_eq(motifs.length(), 3)
}

///|
test "chromvar_create_example_data" {
  let data = @src.chromvar_create_example_data()
  assert_eq(data.counts.length(), 50)
  assert_eq(data.cell_names.length(), 50)
  assert_eq(data.gc_content.length(), 100)
}

///|
test "chromvar_compute_deviations" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let bias = @src.BiasCorrection::new()
  let deviations = @src.chromvar_compute_deviations(data, motifs, bias)
  assert_eq(deviations.deviations.length(), 3)
  assert_eq(deviations.deviations[0].length(), 50)
}

///|
test "chromvar_z_score" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let bias = @src.BiasCorrection::new()
  let deviations = @src.chromvar_compute_deviations(data, motifs, bias)
  let z_scores = @src.chromvar_z_score(deviations.deviations)
  assert_eq(z_scores.length(), 3)
}

///|
test "chromvar_variability" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let bias = @src.BiasCorrection::new()
  let deviations = @src.chromvar_compute_deviations(data, motifs, bias)
  let variability = @src.chromvar_variability(deviations.deviations)
  assert_eq(variability.length(), 3)
}

///|
test "chromvar_differential_deviation" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let bias = @src.BiasCorrection::new()
  let deviations = @src.chromvar_compute_deviations(data, motifs, bias)
  let groups = Array::make(50, 0)
  let mut i = 25
  while i < 50 {
    groups[i] = 1
    i = i + 1
  }
  let diff = @src.chromvar_differential_deviation(deviations.deviations, groups)
  assert_eq(diff.length(), 3)
}

///|
test "chromvar_enrichment_analysis" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let enrichment = @src.chromvar_enrichment_analysis(data, motifs, [
    "Group A", "Group B",
  ])
  assert_eq(enrichment.length(), 3)
}

///|
test "chromvar_filter_peaks" {
  let data = @src.chromvar_create_example_data()
  let filtered = @src.chromvar_filter_peaks(data, 2, 1)
  assert_eq(filtered.counts.length(), 50)
}

///|
test "chromvar_normalize" {
  let data = @src.chromvar_create_example_data()
  let normalized = @src.chromvar_normalize(data.counts, "cpm")
  assert_eq(normalized.length(), 50)
}

///|
test "chromvar_cluster_cells" {
  let data = @src.chromvar_create_example_data()
  let motifs = @src.chromvar_create_example_motifs()
  let bias = @src.BiasCorrection::new()
  let deviations = @src.chromvar_compute_deviations(data, motifs, bias)
  let clusters = @src.chromvar_cluster_cells(deviations.deviations, "kmeans")
  assert_eq(clusters.length(), 50)
}
