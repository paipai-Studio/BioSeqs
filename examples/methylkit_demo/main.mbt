///|
/// methylKit module demo - bisulfite sequencing methylation analysis.
///
/// Run: moon run examples/methylkit_demo
fn main {
  println("=== methylKit Module Demo ===")
  println("Bisulfite sequencing methylation analysis")
  println("Reference: Akalin A et al. (2012) Genome Biol 13:R87")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data (control, treatment)
  // ------------------------------------------------------------------
  println("=== 1. Load sample data ===")
  let (control, treatment) = @bio.methylkit_sample_data()
  println("control:   \{control.sample_id()}, \{control.n_cpgs()} CpGs")
  println("treatment: \{treatment.sample_id()}, \{treatment.n_cpgs()} CpGs")
  println("")

  // ------------------------------------------------------------------
  // 2. Sample summaries (n_cpgs, mean_methylation)
  // ------------------------------------------------------------------
  println("=== 2. Sample summaries ===")
  println(@bio.methylkit_sample_summary(control))
  println(@bio.methylkit_sample_summary(treatment))
  println("")

  // ------------------------------------------------------------------
  // 3. Create config and filter samples by coverage
  // ------------------------------------------------------------------
  println("=== 3. Filter samples by coverage/context ===")
  let config = @bio.MethylKitConfig::new()
  println(
    "config: min_coverage=\{config.min_coverage}, max_coverage=\{config.max_coverage}, context=\{config.context_filter}, q_threshold=\{config.q_threshold}",
  )
  let samples = [control, treatment]
  let filtered = @bio.methylkit_filter_samples(samples, config)
  let mut i = 0
  while i < filtered.length() {
    println(
      "  filtered \{filtered[i].sample_id()}: \{filtered[i].n_cpgs()} CpGs retained",
    )
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Compute correlation between samples
  // ------------------------------------------------------------------
  println("=== 4. Sample correlation ===")
  let corr = @bio.methylkit_sample_correlation(control, treatment)
  println("Pearson correlation (control vs treatment): \{corr}")
  let matrix = @bio.methylkit_correlation_matrix(filtered)
  println("Correlation matrix (\{matrix.length()} x \{matrix[0].length()}):")
  let mut r = 0
  while r < matrix.length() {
    let row = matrix[r]
    let mut s = "  ["
    let mut c = 0
    while c < row.length() {
      s = s + "\{row[c]}"
      if c < row.length() - 1 {
        s = s + ", "
      }
      c = c + 1
    }
    s = s + "]"
    println(s)
    r = r + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Differential methylation analysis
  // ------------------------------------------------------------------
  println("=== 5. Differential methylation analysis ===")
  let diffs = @bio.methylkit_diff_meth(filtered[0], filtered[1], config)
  println(@bio.methylkit_diff_summary(diffs))
  let mut n_sig = 0
  let mut j = 0
  while j < diffs.length() {
    if diffs[j].significant() {
      n_sig = n_sig + 1
    }
    j = j + 1
  }
  println("Number of significant DMPs: \{n_sig}")
  println("")

  // ------------------------------------------------------------------
  // 6. Top 5 DMPs (position, meth_diff, q_value)
  // ------------------------------------------------------------------
  println("=== 6. Top 5 DMPs ===")
  println("chr\tposition\tmeth_diff\tq_value\tsignificant")
  let limit = if diffs.length() < 5 { diffs.length() } else { 5 }
  let mut t = 0
  while t < limit {
    let d = diffs[t]
    println(
      "\{d.chr()}\t\{d.position()}\t\{d.meth_diff()}\t\{d.q_value()}\t\{d.significant()}",
    )
    t = t + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Identify DMRs
  // ------------------------------------------------------------------
  println("=== 7. Differentially methylated regions ===")
  let dmrs = @bio.methylkit_get_dmr(diffs)
  println("Found \{dmrs.length()} DMR(s):")
  let mut d = 0
  while d < dmrs.length() {
    let (chr, start, end) = dmrs[d]
    println("  \{chr}:\{start}-\{end}")
    d = d + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Export first DMP to BED format
  // ------------------------------------------------------------------
  println("=== 8. Export first DMP to BED format ===")
  if diffs.length() > 0 {
    let first_chr = diffs[0].chr()
    let first_pos = diffs[0].position()
    match filtered[0].get_cytosine(first_chr, first_pos) {
      Some(mc) =>
        println("BED: " + @bio.methylkit_cpg_to_bed(mc))
      None => println("Cytosine not found in control sample")
    }
  } else {
    println("No DMPs available to export")
  }

  println("=== Demo complete ===")
}
