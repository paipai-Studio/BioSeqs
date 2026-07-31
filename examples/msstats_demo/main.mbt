///|
/// MSstats module example - Protein significance analysis for mass spectrometry.
///
/// Run: moon run examples/msstats_demo
fn main {
  println("=== MSstats Module Demo ===")
  println("Protein Significance Analysis for Mass Spectrometry Data")
  println("")

  // Generate sample data
  println("=== Sample MS Data ===")
  let features = @bio.msstats_sample_data()
  println("Total features: \{features.length()}")
  let mut n_ctrl = 0
  let mut n_trt = 0
  let mut i = 0
  while i < features.length() {
    if features[i].group == "Ctrl" {
      n_ctrl = n_ctrl + 1
    } else {
      n_trt = n_trt + 1
    }
    i = i + 1
  }
  println("  Control features: \{n_ctrl}")
  println("  Treatment features: \{n_trt}")
  println("")

  // Data processing
  println("=== Data Processing (Median Normalization) ===")
  let processed = @bio.ms_data_process(features, normalization=@bio.ms_norm_median())
  println("Processed features: \{processed.length()}")
  println("  First feature intensity (log2): \{processed[0].intensity}")
  println("")

  // Summarization
  println("=== Protein Summarization (Tukey) ===")
  let summarized = @bio.ms_summarize(processed, method=@bio.ms_summary_tukey())
  println("Protein-level results: \{summarized.length()}")
  let mut j = 0
  while j < summarized.length() {
    println("  \{summarized[j].protein} [\{summarized[j].group}]: log2_abundance=\{summarized[j].log2_abundance}")
    j = j + 1
  }
  println("")

  // Group comparison
  println("=== Group Comparison ===")
  let results = @bio.ms_group_comparison(summarized, fdr_threshold=0.1, log2fc_threshold=0.5)
  println(results.summary())
  println("")

  // Top proteins
  println("=== Top 5 Proteins ===")
  let top = results.get_top_proteins(5)
  let mut k = 0
  while k < top.length() {
    let r = top[k]
    println("  \{r.protein}: log2FC=\{r.log2fc}, p_val=\{r.p_val}, p_adj=\{r.p_adj}, sig=\{r.significant}")
    k = k + 1
  }
  println("")

  // Sample size calculation
  println("=== Sample Size Design ===")
  let ss = @bio.ms_design_sample_size(results.results, desired_power=0.8)
  println("  Desired power: \{ss.desired_power}")
  println("  Required samples: \{ss.n_samples}")
  println("  Detectable proteins: \{ss.n_proteins_detectable}")
  println("  Effect size (delta): \{ss.delta}")
  println("")

  // Normalization methods comparison
  println("=== Normalization Methods ===")
  let norm_none = @bio.ms_data_process(features, normalization=@bio.ms_norm_none())
  let norm_med = @bio.ms_data_process(features, normalization=@bio.ms_norm_median())
  let norm_qu = @bio.ms_data_process(features, normalization=@bio.ms_norm_quantile())
  println("None: first intensity=\{norm_none[0].intensity}")
  println("Median: first intensity=\{norm_med[0].intensity}")
  println("Quantile: first intensity=\{norm_qu[0].intensity}")
}
