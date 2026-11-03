///|
/// edgeR Differential Expression Analysis Demo.
fn main {
  println("=== edgeR Differential Expression Analysis Demo ===\n")

  let counts = [
    [10, 20, 30, 40],
    [50, 60, 70, 80],
    [90, 100, 110, 120],
    [130, 140, 150, 160],
  ]
  let group = ["A", "A", "B", "B"]
  let genes = ["GeneA", "GeneB", "GeneC", "GeneD"]

  println("1. Creating DGEList...")
  let dge = @src.dge_list(counts, group, genes)
  println("   Genes: \{genes.join(", ")}")
  println("   Groups: \{group.join(", ")}")

  println("\n2. Calculating normalization factors...")
  let dge_norm = @src.calc_norm_factors(dge)
  let norm_str = dge_norm.norm_factors.map(fn(x) { x.to_string() }).join(", ")
  println("   Norm factors: \{norm_str}")

  println("\n3. Estimating dispersion...")
  let disp = @src.estimate_dispersion(dge_norm)
  println("   Dispersion: \{disp}")

  println("\n4. Performing exact test (A vs B)...")
  let results = @src.exact_test(dge_norm, "A", "B")
  println("   Results for \{results.genes.length()} genes")

  println("\n5. Getting top 2 genes...")
  let top = @src.top_tags(results, 2)
  for i = 0; i < top.genes.length(); i = i + 1 {
    println("   \{top.genes[i]}: logFC=\{top.log_fc[i]}, p=\{top.p_value[i]}")
  }

  println("\n=== edgeR Demo Complete ===")
}
