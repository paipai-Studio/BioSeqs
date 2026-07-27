///|
/// GCRMA (GC-content adjusted Robust Multi-array Average) Demo.
///
/// This example demonstrates the complete GCRMA pipeline for Affymetrix
/// microarray analysis, including:
///   1. Creating probe info with GC content
///   2. Computing a GC lookup table
///   3. Background correction with GC adjustment
///   4. GC correction on probe intensities
///   5. Quantile normalization
///   6. Full GCRMA pipeline
///   7. Probe set summarization
///
/// Inspired by Bioconductor's gcrma package.

///|
fn main {
  println("=== GCRMA Demo ===")
  println("GC-content adjusted Robust Multi-array Average")
  println("Inspired by Bioconductor's gcrma package\n")

  println("=== 1. Create Probe Info with GC Content ===")
  let probes = [
    @bio.ProbeInfo::new("p1", "ACGTACGTACGT"),
    @bio.ProbeInfo::new("p2", "GCGCGCGCGCGC"),
    @bio.ProbeInfo::new("p3", "ATATATATATAT"),
    @bio.ProbeInfo::new("p4", "ACACACACACAC"),
    @bio.ProbeInfo::new("p5", "TGTGTGTGTGTG"),
    @bio.ProbeInfo::new("p6", "CGCGCGCGCGCG"),
  ]
  let mut i = 0
  while i < probes.length() {
    let p = probes[i]
    let gc_frac = @bio.oligo_gc_fraction(p.sequence)
    println(
      "  \{p.probe_id}: seq=\{p.sequence}, gc_count=\{p.gc_count}, " +
      "gc_frac=\{gc_frac}, affinity=\{p.affinity}",
    )
    i = i + 1
  }

  println("\n=== 2. Compute GC Lookup Table ===")
  let gc_lookup = @bio.compute_gc_lookup_table(probes)
  println("  GC lookup table size: \{gc_lookup.keys().length()}")
  let keys = gc_lookup.keys().collect()
  keys.sort()
  let mut ki = 0
  while ki < keys.length() {
    let k = keys[ki]
    let v = gc_lookup.get_or_default(k, 0.0)
    println("    GC count \{k}: median affinity = \{v}")
    ki = ki + 1
  }

  println("\n=== 3. Background Correction with GC Adjustment ===")
  let pm = [100.0, 200.0, 150.0, 300.0, 250.0, 180.0]
  let mm = [50.0, 80.0, 70.0, 120.0, 110.0, 90.0]
  println("  PM intensities:")
  i = 0
  while i < pm.length() {
    println("    PM[\{i}] = \{pm[i]}")
    i = i + 1
  }
  println("  MM intensities:")
  i = 0
  while i < mm.length() {
    println("    MM[\{i}] = \{mm[i]}")
    i = i + 1
  }

  let bg_idealmm = @bio.gcrma_background_correction(
    pm, mm, 0.5, "IdealMM",
  )
  println("\n  Background correction (IdealMM, gc=0.5):")
  i = 0
  while i < bg_idealmm.length() {
    println("    Corrected[\{i}] = \{bg_idealmm[i]}")
    i = i + 1
  }

  let bg_express = @bio.gcrma_background_correction(
    pm, mm, 0.5, "Express",
  )
  println("\n  Background correction (Express, gc=0.5):")
  i = 0
  while i < bg_express.length() {
    println("    Corrected[\{i}] = \{bg_express[i]}")
    i = i + 1
  }

  println("\n=== 4. GC Correction on Probe Intensities ===")
  let gc_counts = [6, 12, 0, 6, 6, 12]
  let gc_corrected = @bio.gcrma_gc_correction(pm, gc_counts, probes)
  println("  Original PM values -> GC-corrected values:")
  i = 0
  while i < pm.length() {
    println("    PM[\{i}] = \{pm[i]}  ->  corrected = \{gc_corrected[i]}")
    i = i + 1
  }

  println("\n=== 5. Quantile Normalization ===")
  let expression_matrix = [
    [100.0, 200.0, 150.0],
    [150.0, 250.0, 180.0],
    [120.0, 180.0, 160.0],
    [180.0, 220.0, 190.0],
    [90.0, 160.0, 140.0],
    [110.0, 190.0, 170.0],
  ]
  println("  Before normalization:")
  i = 0
  while i < expression_matrix.length() {
    println("    Probe \{i}: \{expression_matrix[i]}")
    i = i + 1
  }

  let normalized = @bio.gcrma_normalize(expression_matrix)
  println("\n  After quantile normalization:")
  i = 0
  while i < normalized.length() {
    println("    Probe \{i}: \{normalized[i]}")
    i = i + 1
  }

  println("\n=== 6. Full GCRMA Pipeline ===")
  let cel_data = [
    [100.0, 200.0, 150.0],
    [150.0, 250.0, 180.0],
    [120.0, 180.0, 160.0],
    [180.0, 220.0, 190.0],
    [90.0, 160.0, 140.0],
    [110.0, 190.0, 170.0],
  ]
  println("  Input CEL data (6 probes x 3 samples):")
  i = 0
  while i < cel_data.length() {
    println("    Probe \{i}: \{cel_data[i]}")
    i = i + 1
  }

  let config = @bio.GCRMAConfig::new(
    background_method="IdealMM",
    normalize=true,
  )
  println("\n  Configuration: background_method=\{config.background_method}, " +
    "normalize=\{config.normalize}")

  let result = @bio.gcrma_process_with_config(cel_data, probes, config)
  println("\n  Processed expression matrix:")
  i = 0
  while i < result.expression_matrix.length() {
    println("    Probe \{i}: \{result.expression_matrix[i]}")
    i = i + 1
  }

  let params = result.gc_correction_params
  println("\n  GC correction parameters:")
  println("    n_probes: \{params.get_or_default("n_probes", 0.0)}")
  println("    n_samples: \{params.get_or_default("n_samples", 0.0)}")
  println(
    "    gc_correction_applied: \{params.get_or_default("gc_correction_applied", 0.0)}",
  )

  println("\n=== 7. Probe Set Summarization ===")
  let expr_for_summary = [
    [100.0, 200.0, 150.0],
    [150.0, 250.0, 180.0],
    [120.0, 180.0, 160.0],
    [180.0, 220.0, 190.0],
    [90.0, 160.0, 140.0],
    [110.0, 190.0, 170.0],
  ]
  let probe_sets : Map[String, Array[Int]] = Map([], capacity=0)
  probe_sets.set("gene1_at", [0, 1, 2])
  probe_sets.set("gene2_at", [3, 4, 5])
  println("  Probe set definition:")
  println("    gene1_at: probes [0, 1, 2]")
  println("    gene2_at: probes [3, 4, 5]")

  let summarized = @bio.gcrma_summarize(expr_for_summary, probe_sets)
  println("\n  Summarized expression (median across probes):")
  let skeys = summarized.keys().collect()
  skeys.sort()
  let mut si = 0
  while si < skeys.length() {
    let gene_id = skeys[si]
    let vals = summarized.get_or_default(gene_id, [])
    println("    \{gene_id}: \{vals}")
    si = si + 1
  }

  println("\n=== GCRMA Demo Complete ===")
}
