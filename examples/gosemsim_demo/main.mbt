///|
/// GOSemSim GO semantic similarity demo.
/// Demonstrates computing semantic similarity between GO terms.

fn main {
  println("=== GOSemSim GO Semantic Similarity Demo ===")
  println("")

  // Build example GO graph
  let g = @bio.gosemsim_example_graph()

  println("GO graph built with root: biological_process")
  println("Terms in graph: 9")

  // Define test term pairs
  println("")
  println("Term Similarity Measures:")
  println("")

  // Compute similarity for each pair
  let resnik1 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @bio.resnik_measure())
  let lin1 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @bio.lin_measure())
  let rel1 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @bio.rel_measure())
  let jiang1 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @bio.jiang_measure())
  let wang1 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @bio.wang_measure())
  println("  transport vs cellular_process:")
  println("    Resnik=\{resnik1.to_string()}, Lin=\{lin1.to_string()}, Rel=\{rel1.to_string()}")
  println("    Jiang=\{jiang1.to_string()}, Wang=\{wang1.to_string()}")

  let resnik2 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @bio.resnik_measure())
  let lin2 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @bio.lin_measure())
  let rel2 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @bio.rel_measure())
  println("  transport vs signal_transduction:")
  println("    Resnik=\{resnik2.to_string()}, Lin=\{lin2.to_string()}, Rel=\{rel2.to_string()}")

  let resnik3 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0006810", @bio.resnik_measure())
  let lin3 = @bio.gosemsim_term_sim(g, "GO:0006810", "GO:0006810", @bio.lin_measure())
  println("  transport vs transport (self):")
  println("    Resnik=\{resnik3.to_string()}, Lin=\{lin3.to_string()}")

  let resnik4 = @bio.gosemsim_term_sim(g, "GO:0009987", "GO:0009058", @bio.resnik_measure())
  let lin4 = @bio.gosemsim_term_sim(g, "GO:0009987", "GO:0009058", @bio.lin_measure())
  println("  cellular_process vs biosynthetic_process:")
  println("    Resnik=\{resnik4.to_string()}, Lin=\{lin4.to_string()}")

  // Gene-level similarity
  println("")
  println("=== Gene-Level Similarity ===")
  println("")

  let go1 = ["GO:0006810", "GO:0007165", "GO:0009987"]
  let go2 = ["GO:0009987", "GO:0023046", "GO:0009058"]

  println("Gene 1 GO terms: \{go1}")
  println("Gene 2 GO terms: \{go2}")

  let sim_max = @bio.gosemsim_gen_sim(g, go1, go2, @bio.resnik_measure(), combine = "max")
  let sim_avg = @bio.gosemsim_gen_sim(g, go1, go2, @bio.lin_measure(), combine = "avg")
  let sim_rcmax = @bio.gosemsim_gen_sim(g, go1, go2, @bio.rel_measure(), combine = "rcmax")

  println("")
  println("Combined similarity:")
  println("  Max (Resnik): \{sim_max.to_string()}")
  println("  Avg (Lin):    \{sim_avg.to_string()}")
  println("  RCMax (Rel):  \{sim_rcmax.to_string()}")

  println("")
  println("Done!")
}
