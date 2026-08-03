///|
/// MCPcounter module demo - marker-based cell population counting.
///
/// Run: moon run examples/mcp_counter_demo
fn main {
  println("=== MCPcounter Module Demo ===")
  println("Microenvironment Cell Population counter")
  println("Reference: Becht E et al. (2016) eLife 5:e18516")
  println("")

  // ------------------------------------------------------------------
  // 1. Display built-in cell populations
  // ------------------------------------------------------------------
  println("=== 1. Built-in cell populations ===")
  let pops = @bio.mcp_default_populations()
  println("Total populations: \{pops.length()}")
  for pop in pops {
    println("  \{pop.name}: \{pop.markers.length()} markers")
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Load sample data
  // ------------------------------------------------------------------
  println("=== 2. Load sample expression data ===")
  let (genes, samples, matrix) = @bio.mcp_sample_data()
  println("Genes: \{genes.length()}")
  println("Samples: \{samples}")
  println("")

  // ------------------------------------------------------------------
  // 3. Run MCPcounter
  // ------------------------------------------------------------------
  println("=== 3. Run MCPcounter ===")
  let result = @bio.mcp_run(genes, samples, matrix)
  println("Populations scored: \{result.populations().length()}")
  println("")

  // ------------------------------------------------------------------
  // 4. Display scores per population
  // ------------------------------------------------------------------
  println("=== 4. Per-population abundance scores ===")
  let mut pi = 0
  while pi < result.populations().length() {
    let pop_name = result.populations()[pi]
    println("--- \{pop_name} ---")
    let mut si = 0
    while si < result.sample_names().length() {
      let score = result.scores()[pi][si]
      let bar_len = if score > 0.0 { (score * 2.0).to_int() } else { 0 }
      let mut bar = ""
      let mut bi = 0
      while bi < bar_len {
        bar = bar + "#"
        bi = bi + 1
      }
      println("  \{result.sample_names()[si]}: \{score} \{bar}")
      si = si + 1
    }
    pi = pi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Compare tumor compositions
  // ------------------------------------------------------------------
  println("=== 5. Tumor composition comparison ===")
  let t_scores = result.get_population_scores("T cells")
  let b_scores = result.get_population_scores("B lineage")
  let f_scores = result.get_population_scores("Fibroblasts")
  let mut i = 0
  while i < samples.length() {
    println(
      "\{samples[i]}: T=\{t_scores[i]}, B=\{b_scores[i]}, Fibro=\{f_scores[i]}",
    )
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Query specific score
  // ------------------------------------------------------------------
  println("=== 6. Query specific scores ===")
  let t_a = result.get_score("T cells", "Tumor_A")
  let b_a = result.get_score("B lineage", "Tumor_A")
  println("Tumor_A - T cells: \{t_a}")
  println("Tumor_A - B lineage: \{b_a}")
  if t_a > b_a {
    println("-> Tumor_A is T-cell rich (as designed)")
  } else {
    println("-> Tumor_A is B-cell rich (unexpected)")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Custom population example
  // ------------------------------------------------------------------
  println("=== 7. Custom cell population example ===")
  let custom_pop = @bio.McpCellPopulation::new(
    "Custom lineage",
    ["CD3D", "CD3E", "CD2", "TRAC"],
  )
  let custom_result = @bio.mcp_run(
    genes, samples, matrix,
    populations=[custom_pop],
  )
  println("Custom population scores:")
  let mut ci = 0
  while ci < samples.length() {
    let s = custom_result.scores()[0][ci]
    println("  \{samples[ci]}: \{s}")
    ci = ci + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Format full result table
  // ------------------------------------------------------------------
  println("=== 8. Full result table ===")
  println(result.to_string())

  println("=== Demo Complete ===")
}
