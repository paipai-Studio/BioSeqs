///|
fn main {
  println("=== FGSEA (Fast Gene Set Enrichment Analysis) Demo ===")

  let gene_sets = @bio.sample_gene_sets()
  let stats = @bio.sample_gene_stats()

  println("")
  println("Input:")
  println("  Gene sets: " + gene_sets.length().to_string())
  println("  Ranked genes: " + stats.length().to_string())

  println("")
  println("Top 10 ranked genes:")
  let mut i = 0
  while i < 10 && i < stats.length() {
    let (gene, stat) = stats[i]
    println("  " + (i + 1).to_string() + ". " + gene + " (stat=" + stat.to_string() + ")")
    i = i + 1
  }

  println("")
  println("Running FGSEA with 100 permutations...")
  let results = @bio.fgsea(gene_sets, stats, nperm = 100)

  println("")
  println("FGSEA Results:")
  println("==============")

  i = 0
  while i < results.length() {
    let r = results[i]
    println("")
    println("Pathway: " + r.pathway)
    println("  ES: " + r.es.to_string())
    println("  NES: " + r.nes.to_string())
    println("  p-value: " + r.pval.to_string())
    println("  adj.pval: " + r.padj.to_string())
    println("  Size: " + r.size.to_string())
    println("  Leading edge genes: " + r.leadingEdge.length().to_string())
    i = i + 1
  }

  let top_up = @bio.fgsea_top_up(results, n = 3)
  println("")
  println("Top 3 Up-regulated Pathways:")
  i = 0
  while i < top_up.length() {
    let r = top_up[i]
    println("  " + (i + 1).to_string() + ". " + r.pathway + " (NES=" + r.nes.to_string() + ")")
    i = i + 1
  }

  let top_down = @bio.fgsea_top_down(results, n = 3)
  println("")
  println("Top 3 Down-regulated Pathways:")
  i = 0
  while i < top_down.length() {
    let r = top_down[i]
    println("  " + (i + 1).to_string() + ". " + r.pathway + " (NES=" + r.nes.to_string() + ")")
    i = i + 1
  }

  println("")
  println("=== Demo Complete ===")
}
