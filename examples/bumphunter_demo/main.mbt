///|
/// Bumphunter module demo - methylation/signal region (bump) detection.
///
/// Run: moon run examples/bumphunter_demo
fn main {
  println("=== Bumphunter Module Demo ===")
  println("Differentially methylated region (DMR) detection")
  println("Reference: Jaffe AE et al. (2012) Int J Epidemiol 41(1):200-9")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample methylation data with a known bump
  // ------------------------------------------------------------------
  println("=== 1. Load sample methylation data ===")
  let positions = @bio.bump_sample_data()
  println("Positions: \{positions.length()} CpG sites on chr1")
  println("Samples: \{positions[0].values().length()} (4 cases + 4 controls)")
  println("Bump region: positions 10-20 (cases elevated)")
  println("")

  // ------------------------------------------------------------------
  // 2. Compute per-position t-statistics
  // ------------------------------------------------------------------
  println("=== 2. Per-position t-statistics ===")
  let stats = @bio.bump_compute_t_statistics_test(positions, 4)
  let mut i = 0
  while i < stats.length() {
    let marker = if i >= 10 && i <= 20 { " [bump]" } else { "" }
    if i % 5 == 0 {
      println(
        "  pos \{i}: t=\{stats[i]}\{marker}",
      )
    }
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Smooth the signal
  // ------------------------------------------------------------------
  println("=== 3. Smoothed signal (running mean, k=2) ===")
  let smoothed = @bio.bump_smooth_test(stats, 2)
  let mut si = 0
  while si < smoothed.length() {
    if si % 5 == 0 {
      println(
        "  pos \{si}: smoothed=\{smoothed[si]}",
      )
    }
    si = si + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Find candidate bumps
  // ------------------------------------------------------------------
  println("=== 4. Candidate bumps (cutoff=2.0) ===")
  let candidates = @bio.bump_find_candidates_test(positions, smoothed, 2.0)
  println("Found \{candidates.length()} candidate bump(s)")
  let mut ci = 0
  while ci < candidates.length() {
    let b = candidates[ci]
    println(
      "  chr=\{b.chrom()} start=\{b.start()} end=\{b.end_()} value=\{b.value()} area=\{b.area()}",
    )
    ci = ci + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Run full bumphunter with permutation testing
  // ------------------------------------------------------------------
  println("=== 5. Full bumphunter analysis (20 permutations) ===")
  let results = @bio.bump_hunt(positions, 4, 2.0, 2, 20)
  println("Significant bumps: \{results.length()}")
  let mut ri = 0
  while ri < results.length() {
    let r = results[ri]
    println(
      "  chr=\{r.chrom()} start=\{r.start()} end=\{r.end_()} value=\{r.value()} p=\{r.p_value()} fdr=\{r.fdr()}",
    )
    ri = ri + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. High cutoff finds nothing
  // ------------------------------------------------------------------
  println("=== 6. High cutoff (100.0) finds no bumps ===")
  let none = @bio.bump_hunt(positions, 4, 100.0, 2, 10)
  println("Bumps found: \{none.length()}")
  println("")

  // ------------------------------------------------------------------
  // 7. Full results table
  // ------------------------------------------------------------------
  println("=== 7. Results table ===")
  println(@bio.bump_results_to_string(results))

  println("=== Demo complete ===")
}
