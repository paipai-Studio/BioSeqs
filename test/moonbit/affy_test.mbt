///|
test "affy_probe_set_create" {
  let pm = [[100.0, 200.0], [150.0, 250.0], [120.0, 180.0]]
  let mm = [[50.0, 80.0], [70.0, 90.0], [60.0, 75.0]]
  let probe_set = @src.ProbeSet::new("gene1_at", pm, mm)

  assert_eq(probe_set.probe_set_id, "gene1_at")
  assert_eq(probe_set.pm_intensities.length(), 3)
}

///|
test "affy_affy_batch_create" {
  let probe_set_ids = ["gene1_at", "gene2_at"]
  let pm_matrix = [
    [[100.0, 200.0], [150.0, 250.0]],
    [[180.0, 220.0], [200.0, 240.0]],
  ]
  let sample_names = ["sample1", "sample2"]
  let batch = @src.affy_batch_from_matrix(
    probe_set_ids, pm_matrix, sample_names,
  )

  assert_eq(batch.n_probe_sets, 2)
  assert_eq(batch.n_samples, 2)
}

///|
test "affy_background_correct_rma" {
  let pm = [[60.0, 150.0], [40.0, 200.0]]
  let corrected = @src.background_correct_rma(pm)

  assert_eq(corrected.length(), 2)
  assert_eq(corrected[0][0] > 0.0, true)
}

///|
test "affy_compute_pm_mm_difference" {
  let pm = [[100.0, 200.0], [150.0, 250.0]]
  let mm = [[50.0, 80.0], [70.0, 90.0]]
  let probe_set = @src.ProbeSet::new("gene1_at", pm, mm)
  let diff = @src.compute_pm_mm_difference(probe_set)

  assert_eq(diff.length(), 2)
  assert_eq(diff[0][0], 50.0)
}

///|
test "affy_summarize_probeset_median" {
  let pm = [[100.0, 200.0], [150.0, 250.0], [120.0, 180.0]]
  let mm : Array[Array[Double]] = []
  let probe_set = @src.ProbeSet::new("gene1_at", pm, mm)
  let summary = @src.summarize_probeset(probe_set, "median")

  assert_eq(summary.length(), 2)
}

///|
test "affy_summarize_probeset_mean" {
  let pm = [[100.0, 200.0], [150.0, 250.0], [120.0, 180.0]]
  let mm : Array[Array[Double]] = []
  let probe_set = @src.ProbeSet::new("gene1_at", pm, mm)
  let summary = @src.summarize_probeset(probe_set, "mean")

  assert_eq(summary.length(), 2)
}

///|
test "affy_rma_normalize" {
  let probe_set_ids = ["gene1_at", "gene2_at"]
  let pm_matrix = [
    [[100.0, 200.0], [150.0, 250.0], [120.0, 180.0]],
    [[180.0, 220.0], [200.0, 240.0], [190.0, 210.0]],
  ]
  let sample_names = ["sample1", "sample2"]
  let batch = @src.affy_batch_from_matrix(
    probe_set_ids, pm_matrix, sample_names,
  )
  let normalized = @src.rma_normalize(batch)

  assert_eq(normalized.length(), 2)
  assert_eq(normalized[0].length(), 2)
}
