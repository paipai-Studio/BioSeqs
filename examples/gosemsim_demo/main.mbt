/// Demo: GOSemSim — GO semantic similarity (Bioconductor).
///
/// Ported from GOSemSim R package (Yu et al., 2010).
/// Five measures: Resnik, Lin, Rel (Schlicker), Jiang, Wang.
/// Gene-level combine: max, avg, rcmax, BMA.

fn main {
  let g = @src.gosemsim_example_graph()

  println("=== GOSemSim Demo (Bioconductor) ===")
  println("GO semantic similarity — 5 measures × gene-level combine")
  println("")

  println("Example GO DAG (root=GO:0008150):")
  println("  GO:0008150(ic=0.0) biological_process")
  println("    ├── GO:0050896(0.3) response_to_stimulus")
  println("    │    ├── GO:0009987(0.5) cellular_process → GO:0006810(0.7) transport")
  println("    │    └── GO:0023046(0.4) signaling → GO:0007165(0.6) signal_transduction")
  println("    └── GO:0009058(0.2) biosynthetic_process → GO:0010467(0.35) gene_expression")
  println("")

  // ── Term-level comparison ──────────────────────────────────────────
  println("=== Term-level similarity (5 measures) ===")

  fn fmt(v : Double) -> String {
    v.to_string()
  }

  let pair : Array[(String, String, String)] = Array(capacity=0)
  pair.push(("GO:0006810", "GO:0006810", "self"))
  pair.push(("GO:0006810", "GO:0009987", "parent-child"))
  pair.push(("GO:0006810", "GO:0007165", "sibling_subtree"))
  pair.push(("GO:0009987", "GO:0023046", "siblings"))
  pair.push(("GO:0006810", "GO:0010467", "no_common"))
  pair.push(("GO:0007165", "GO:0010467", "distant_root"))

  println("Pair\tResnik\tLin\tRel\tJiang\tWang")
  let mut i = 0
  while i < pair.length() {
    let (t1, t2, desc) = pair[i]
    let resnik = @src.gosemsim_term_sim(g, t1, t2, @src.resnik_measure())
    let lin = @src.gosemsim_term_sim(g, t1, t2, @src.lin_measure())
    let rel = @src.gosemsim_term_sim(g, t1, t2, @src.rel_measure())
    let jiang = @src.gosemsim_term_sim(g, t1, t2, @src.jiang_measure())
    let wang = @src.gosemsim_term_sim(g, t1, t2, @src.wang_measure())
    println(desc + "\t" + fmt(resnik) + "\t" + fmt(lin) + "\t" +
            fmt(rel) + "\t" + fmt(jiang) + "\t" + fmt(wang))
    i = i + 1
  }
  println("")

  // ── Gene-level combine methods ──────────────────────────────────────
  println("=== Gene-level similarity (Resnik × combine) ===")
  let go1 : Array[String] = ["GO:0006810", "GO:0007165"]
  let go2 : Array[String] = ["GO:0006810", "GO:0009987"]

  // Resnik 2x2 matrix:
  //   GO:0006810 vs GO:0006810 = 0.7
  //   GO:0006810 vs GO:0009987 = 0.5
  //   GO:0007165 vs GO:0006810 = 0.3
  //   GO:0007165 vs GO:0009987 = 0.3
  let max_v = @src.gosemsim_gen_sim(g, go1, go2, @src.resnik_measure(), combine="max")
  let avg_v = @src.gosemsim_gen_sim(g, go1, go2, @src.resnik_measure(), combine="avg")
  let rcmax_v = @src.gosemsim_gen_sim(g, go1, go2, @src.resnik_measure(), combine="rcmax")
  let bma_v = @src.gosemsim_gen_sim(g, go1, go2, @src.resnik_measure(), combine="BMA")

  println("  Gene1: transport(0.7) + signal_transduction(0.6)")
  println("  Gene2: transport(0.7) + cellular_process(0.5)")
  println("  Resnik matrix: [[0.7, 0.5], [0.3, 0.3]]")
  println("  combine=max    → " + fmt(max_v))
  println("  combine=avg    → " + fmt(avg_v))
  println("  combine=rcmax  → " + fmt(rcmax_v))
  println("  combine=BMA    → " + fmt(bma_v))
  println("")

  println("=== Done ===")
}
