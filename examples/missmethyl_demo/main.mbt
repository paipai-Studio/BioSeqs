///|
/// missMethyl demo - Differential methylation analysis with GO enrichment.
///
/// Demonstrates beta/M-value conversion, DMP analysis with BH-FDR correction,
/// gene-level GO enrichment, probe-bias-corrected GO enrichment, probe
/// filtering, and textual summaries.
///
/// Reference: Phipson B, Maksimovic J, Oshlack A (2016) missMethyl: an R
/// package for analyzing data from Illumina's HumanMethylation450 platform.
/// Bioinformatics 32(2):286-288.
/// Run: moon run examples/missmethyl_demo
fn main {
  println("=== missMethyl Demo ===")
  println("Differential methylation analysis with GO enrichment")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample annotations, GO database, and beta-value matrix
  // ------------------------------------------------------------------
  println("=== 1. Load sample data ===")
  let annotations = @bio.mm_sample_annotations()
  let go_db = @bio.mm_sample_go_database()
  let beta = @bio.mm_sample_beta()
  let group : Array[Bool] = [
    false, false, false, false, true, true, true, true,
  ]
  println("Loaded \{annotations.length()} probe annotations")
  println("Loaded \{beta.length()} probes x \{beta[0].length()} samples beta matrix")
  println("Group: 4 control (false) + 4 treatment (true)")
  println("")

  // ------------------------------------------------------------------
  // 2. Beta/M value conversion examples
  // ------------------------------------------------------------------
  println("=== 2. Beta/M value conversions ===")
  let betas : Array[Double] = [0.1, 0.25, 0.5, 0.75, 0.9]
  for b in betas {
    let m = @bio.mm_beta_to_m(b)
    let b_back = @bio.mm_m_to_beta(m)
    println("  beta=" + b.to_string() + " -> M=" + m.to_string() + " -> beta(back)=" + b_back.to_string())
  }
  // Edge cases: clamping
  let m0 = @bio.mm_beta_to_m(0.0)
  let m1 = @bio.mm_beta_to_m(1.0)
  println("  beta=0.0 (clamped) -> M=" + m0.to_string())
  println("  beta=1.0 (clamped) -> M=" + m1.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. DMP analysis (limma-style t-test on M-values + BH-FDR)
  // ------------------------------------------------------------------
  println("=== 3. Differential Methylation Position (DMP) analysis ===")
  let dmp_results = @bio.mm_dmp_analysis(beta, group, annotations)
  let mut n_sig_05 = 0
  let mut n_sig_01 = 0
  for r in dmp_results {
    if r.is_significant(0.05) {
      n_sig_05 = n_sig_05 + 1
    }
    if r.is_significant(0.01) {
      n_sig_01 = n_sig_01 + 1
    }
  }
  println("Total probes tested: \{dmp_results.length()}")
  println("Significant (adj p < 0.05): \{n_sig_05}")
  println("Significant (adj p < 0.01): \{n_sig_01}")
  println("")

  // Manual bubble sort by adjusted p-value (ascending) for display.
  let sorted = dmp_results.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      if sorted[j].adj_p_value() > sorted[j + 1].adj_p_value() {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  println("Top 5 DMPs by significance:")
  let top_n = if n < 5 { n } else { 5 }
  for k in 0..<top_n {
    let r = sorted[k]
    println(
      "  " +
        r.probe_id() +
        " " +
        r.chr() +
        ":" +
        r.pos().to_string() +
        " gene=" +
        r.gene() +
        " logFC=" +
        r.log_fc().to_string() +
        " t=" +
        r.t_stat().to_string() +
        " adjP=" +
        r.adj_p_value().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Probe filtering by chromosome and CpG island relation
  // ------------------------------------------------------------------
  println("=== 4. Probe filtering ===")
  let chr1 = @bio.mm_filter_probes(annotations, "chr1", "all")
  let islands = @bio.mm_filter_probes(annotations, "all", "Island")
  let chr1_islands = @bio.mm_filter_probes(annotations, "chr1", "Island")
  println("chr1 probes: \{chr1.length()}")
  println("Island probes (all chr): \{islands.length()}")
  println("chr1 + Island probes: \{chr1_islands.length()}")
  for ann in chr1_islands {
    println(
      "  " +
        ann.probe_id() +
        " gene=" +
        ann.gene() +
        " pos=" +
        ann.pos().to_string() +
        " type=" +
        ann.probe_type(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 5. GO enrichment at the gene level (Fisher's exact + BH-FDR)
  // ------------------------------------------------------------------
  println("=== 5. GO enrichment (gene-level) ===")
  let sig_probes : Array[String] = Array::new()
  for r in sorted {
    if r.is_significant(0.05) {
      sig_probes.push(r.probe_id())
    }
  }
  let all_probes : Array[String] = Array::new()
  for ann in annotations {
    all_probes.push(ann.probe_id())
  }
  println("Significant probes: \{sig_probes.length()} / \{all_probes.length()} total")
  let go_results = @bio.mm_go_enrichment(
    sig_probes, all_probes, annotations, go_db,
  )
  println("GO terms tested: \{go_results.length()}")
  for r in go_results {
    println(
      "  " +
        r.go_id() +
        " [" +
        r.ontology() +
        "] " +
        r.go_term() +
        " sig=" +
        r.n_sig_probes().to_string() +
        "/" +
        r.n_probes().to_string() +
        " expected=" +
        r.expected().to_string() +
        " p=" +
        r.p_value().to_string() +
        " adjP=" +
        r.adj_p_value().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Probe-bias-corrected GO enrichment (the core missMethyl method)
  // ------------------------------------------------------------------
  println("=== 6. Probe-bias-corrected GO enrichment ===")
  let corrected_results = @bio.mm_probe_bias_correction(
    sig_probes, all_probes, annotations, go_db,
  )
  println("GO terms tested: \{corrected_results.length()}")
  for r in corrected_results {
    println(
      "  " +
        r.go_id() +
        " [" +
        r.ontology() +
        "] " +
        r.go_term() +
        " sig=" +
        r.n_sig_probes().to_string() +
        "/" +
        r.n_probes().to_string() +
        " expected=" +
        r.expected().to_string() +
        " p=" +
        r.p_value().to_string() +
        " adjP=" +
        r.adj_p_value().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Textual summaries
  // ------------------------------------------------------------------
  println("=== 7. Summaries ===")
  println(@bio.mm_summary(dmp_results))
  println(@bio.mm_go_summary(go_results))
  println(@bio.mm_go_summary(corrected_results))

  println("=== missMethyl Demo complete ===")
}
