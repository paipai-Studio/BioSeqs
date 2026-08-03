///|
/// Primer3 module demo - PCR primer design and analysis.
///
/// Run: moon run examples/primer3_demo
fn main {
  println("=== Primer3 Module Demo ===")
  println("PCR primer design with Tm calculation, GC content, and complementarity scoring")
  println("Reference: Rozen S, Skaletsky H (2000) Methods Mol Biol 132:365-86")
  println("")

  // ------------------------------------------------------------------
  // 1. GC content calculation
  // ------------------------------------------------------------------
  println("=== 1. GC Content Calculation ===")
  let seq1 = "ATGCGTACGTAGATGCGTAC"
  let gc1 = @bio.primer3_gc_content(seq1)
  println("Sequence: \{seq1}")
  println("GC content: \{gc1 * 100.0}%")
  println("")

  let seq2 = "GGCCGGCCGGCCGGCCGGCC"
  let gc2 = @bio.primer3_gc_content(seq2)
  println("Sequence: \{seq2}")
  println("GC content: \{gc2 * 100.0}%")
  println("")

  // ------------------------------------------------------------------
  // 2. Melting temperature (Tm) calculation
  // ------------------------------------------------------------------
  println("=== 2. Melting Temperature (Tm) ===")
  let short_seq = "AATTGC"
  let tm_wallace = @bio.primer3_tm_wallace(short_seq)
  println("Short oligo (\{short_seq}):")
  println("  Wallace rule Tm: \{tm_wallace} C")
  println("")

  let long_seq = "ATGCGTACGTAGATGCGTAC"
  let tm_salt = @bio.primer3_tm_salt(long_seq, 0.05)
  let tm_auto = @bio.primer3_tm(long_seq, 0.05)
  println("Longer primer (\{long_seq}):")
  println("  Salt-adjusted Tm (50mM Na+): \{tm_salt} C")
  println("  Auto-selected Tm: \{tm_auto} C")
  println("")

  // ------------------------------------------------------------------
  // 3. Self-complementarity and cross-dimer scoring
  // ------------------------------------------------------------------
  println("=== 3. Complementarity Scoring ===")
  let primer1 = "GAATTC"
  let sc = @bio.primer3_self_complement(primer1)
  println("Primer: \{primer1}")
  println("  Self-complementarity score: \{sc} (lower is better)")
  let hp = @bio.primer3_hairpin_3prime(primer1, 4)
  println("  3' hairpin potential (last 4 nt): \{hp}")
  println("")

  let primer2 = "ATGCATGC"
  let cd = @bio.primer3_cross_dimer(primer1, primer2)
  println("Cross-dimer between \{primer1} and \{primer2}: \{cd}")
  println("")

  // ------------------------------------------------------------------
  // 4. Primer design pipeline
  // ------------------------------------------------------------------
  println("=== 4. Primer Design Pipeline ===")
  let cfg = @bio.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @bio.primer3_sample_template()
  println("Template length: \{template.length()} bp")
  println("Config: len=\{cfg.primer_length_min()}-\{cfg.primer_length_max()}, Tm=\{cfg.tm_min()}-\{cfg.tm_max()}, GC=\{cfg.gc_min() * 100.0}-\{cfg.gc_max() * 100.0}%")
  println("")

  let design = @bio.primer3_design_primers(template, cfg)
  println("Forward candidates: \{design.n_forward()}")
  println("Reverse candidates: \{design.n_reverse()}")
  println("Has pair: \{design.has_pair()}")
  println("")

  // ------------------------------------------------------------------
  // 5. Display best primers
  // ------------------------------------------------------------------
  println("=== 5. Best Primer Candidates ===")
  match design.best_forward() {
    Some(f) => println("Best forward: \{f.to_string()}")
    None => println("Best forward: (none)")
  }
  match design.best_reverse() {
    Some(r) => println("Best reverse: \{r.to_string()}")
    None => println("Best reverse: (none)")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Amplicon length
  // ------------------------------------------------------------------
  println("=== 6. Amplicon Information ===")
  match design.best_pair() {
    Some((f, r)) => {
      let amp_len = @bio.primer3_pair_amplicon_length(f, r)
      println("Best pair amplicon length: \{amp_len} bp")
    }
    None => println("No primer pair available")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Display top 5 forward and reverse candidates
  // ------------------------------------------------------------------
  println("=== 7. Top Primer Candidates ===")
  let fwds = design.forward()
  let n_fwd_show = if fwds.length() < 5 { fwds.length() } else { 5 }
  println("Top \{n_fwd_show} forward primers:")
  for i in 0..<n_fwd_show {
    println("  \{i + 1}. \{fwds[i].to_string()}")
  }
  println("")

  let revs = design.reverse()
  let n_rev_show = if revs.length() < 5 { revs.length() } else { 5 }
  println("Top \{n_rev_show} reverse primers:")
  for i in 0..<n_rev_show {
    println("  \{i + 1}. \{revs[i].to_string()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Design summary
  // ------------------------------------------------------------------
  println("=== 8. Design Summary ===")
  let summary = @bio.primer3_design_summary(design)
  println(summary)

  println("=== Demo complete ===")
}
