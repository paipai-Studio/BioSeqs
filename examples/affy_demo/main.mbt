///|
/// Affymetrix microarray data analysis example.
/// Inspired by Biopython Bio.Affy.
///
/// Run: moon run examples/affy_demo/main.mbt
fn main {
  println("=== Affy Demo ===")
  println("Inspired by Biopython Bio.Affy package")
  println("")

  println("1. Create AffyBatch from probe intensities")
  let probe_set_ids = ["gene1_at", "gene2_at", "gene3_at"]
  let pm_matrix = [
    [[100.0, 200.0, 150.0], [150.0, 250.0, 180.0], [120.0, 180.0, 220.0]],
    [[180.0, 220.0, 200.0], [200.0, 240.0, 210.0], [190.0, 210.0, 230.0]],
    [[250.0, 300.0, 280.0], [270.0, 320.0, 290.0], [260.0, 310.0, 275.0]],
  ]
  let sample_names = ["control1", "treatment1", "treatment2"]
  let batch = @bio.affy_batch_from_matrix(
    probe_set_ids, pm_matrix, sample_names,
  )
  println(
    "  Created AffyBatch with \{batch.n_probe_sets} probe sets and \{batch.n_samples} samples",
  )
  println("")

  println("2. Background correction (RMA)")
  let bg_corrected = @bio.background_correct_rma(pm_matrix[0])
  println("  Background corrected probe set 0:")
  let mut i = 0
  while i < bg_corrected.length() {
    println("    Probe \{i}: \{bg_corrected[i]}")
    i = i + 1
  }
  println("")

  println("3. PM-MM difference calculation")
  let pm = [[100.0, 200.0], [150.0, 250.0]]
  let mm = [[50.0, 80.0], [70.0, 90.0]]
  let probe_set = @bio.ProbeSet::new("gene1_at", pm, mm)
  let diff = @bio.compute_pm_mm_difference(probe_set)
  println("  PM-MM differences:")
  i = 0
  while i < diff.length() {
    println("    Probe \{i}: \{diff[i]}")
    i = i + 1
  }
  println("")

  println("4. Probe set summarization (median)")
  let median_summary = @bio.summarize_probeset(probe_set, "median")
  println("  Median summary: \{median_summary}")
  let mean_summary = @bio.summarize_probeset(probe_set, "mean")
  println("  Mean summary: \{mean_summary}")
  println("")

  println("5. Full RMA normalization")
  let normalized = @bio.rma_normalize(batch)
  println("  RMA normalized expression values:")
  i = 0
  while i < normalized.length() {
    println("    \{probe_set_ids[i]}: \{normalized[i]}")
    i = i + 1
  }
}
