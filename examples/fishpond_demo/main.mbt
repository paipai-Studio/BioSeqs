///|
/// Fishpond (Swish) module demo - non-parametric differential expression.
///
/// Run: moon run examples/fishpond_demo
fn main {
  println("=== Fishpond (Swish) Module Demo ===")
  println("Non-parametric differential expression analysis")
  println("Reference: Zhu A et al. (2019) Nucleic Acids Res 47(18):e105")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample count data
  // ------------------------------------------------------------------
  println("=== 1. Load sample count data ===")
  let counts = @bio.fish_sample_counts()
  println("Transcripts: \{counts.transcript_names().length()}")
  println("Samples: \{counts.sample_names().length()} (8 case + 8 control)")
  println("Transcripts: \{counts.transcript_names()}")
  println("")

  // ------------------------------------------------------------------
  // 2. Mann-Whitney-Wilcoxon statistic example
  // ------------------------------------------------------------------
  println("=== 2. Mann-Whitney-Wilcoxon statistic example ===")
  let w_stat = @bio.fish_mann_whitney_test(
    [10.0, 20.0, 30.0],
    [1.0, 2.0, 3.0],
  )
  println("W statistic for clear separation: \{w_stat}")
  let w_tie = @bio.fish_mann_whitney_test(
    [5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0],
  )
  println("W statistic for all-tied data: \{w_tie}")
  println("")

  // ------------------------------------------------------------------
  // 3. Log2 fold change example
  // ------------------------------------------------------------------
  println("=== 3. Log2 fold change example ===")
  let lfc = @bio.fish_log2fc_test([8.0, 8.0], [2.0, 2.0])
  println("log2FC for case=8, control=2: \{lfc}")
  println("")

  // ------------------------------------------------------------------
  // 4. Run full Swish analysis
  // ------------------------------------------------------------------
  println("=== 4. Run Swish analysis (500 permutations) ===")
  let results = @bio.fish_swish(counts, "case", 500)
  println("Analyzed \{results.length()} transcripts")
  println("")

  // ------------------------------------------------------------------
  // 5. Display results per transcript
  // ------------------------------------------------------------------
  println("=== 5. Per-transcript results ===")
  for r in results {
    let marker = if r.direction() == "up" {
      " ++"
    } else if r.direction() == "down" {
      " --"
    } else {
      ""
    }
    println(
      "  \{r.transcript()}: log2FC=\{r.log2_fold_change()}, p=\{r.p_value()}, fdr=\{r.fdr()}, dir=\{r.direction()}\{marker}",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Filter significant results
  // ------------------------------------------------------------------
  println("=== 6. Significant transcripts (FDR < 0.05) ===")
  let sig = @bio.fish_significant(results, 0.05)
  println("Significant: \{sig.length()} transcripts")
  for r in sig {
    println(
      "  \{r.transcript()} (\{r.direction()}): log2FC=\{r.log2_fold_change()}, fdr=\{r.fdr()}",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Count up/down/ns
  // ------------------------------------------------------------------
  println("=== 7. Direction summary ===")
  let mut n_up = 0
  let mut n_down = 0
  let mut n_ns = 0
  for r in results {
    if r.direction() == "up" {
      n_up = n_up + 1
    } else if r.direction() == "down" {
      n_down = n_down + 1
    } else {
      n_ns = n_ns + 1
    }
  }
  println("  Up-regulated: \{n_up}")
  println("  Down-regulated: \{n_down}")
  println("  Not significant: \{n_ns}")
  println("")

  // ------------------------------------------------------------------
  // 8. Full results table
  // ------------------------------------------------------------------
  println("=== 8. Full results table ===")
  println(@bio.fish_results_to_string(results))

  println("=== Demo complete ===")
}
