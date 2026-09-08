///|
/// goseq Demo - GO enrichment analysis with length bias correction.
/// Demonstrates Wallenius non-central hypergeometric distribution
/// for gene length bias correction in GO analysis.

fn main {
  println("=== goseq Demo ===")
  println("")

  // 1. Create synthetic gene data with length bias
  // Longer genes are more likely to be differentially expressed
  println("--- Step 1: Create Gene Data ---")

  let gene_names : Array[String] = Array::new()
  let gene_lengths : Array[Double] = Array::new()
  let de_status : Array[Int] = Array::new()

  // Short genes (500-1500 bp): 10% DE
  let mut i = 0
  while i < 50 {
    gene_names.push("ShortGene_" + i.to_string())
    gene_lengths.push(500.0 + (i * 20).to_double())
    de_status.push(if i < 5 { 1 } else { 0 })
    i = i + 1
  }

  // Medium genes (2000-3500 bp): 30% DE
  i = 0
  while i < 50 {
    gene_names.push("MedGene_" + i.to_string())
    gene_lengths.push(2000.0 + (i * 30).to_double())
    de_status.push(if i < 15 { 1 } else { 0 })
    i = i + 1
  }

  // Long genes (4000-6000 bp): 60% DE
  i = 0
  while i < 50 {
    gene_names.push("LongGene_" + i.to_string())
    gene_lengths.push(4000.0 + (i * 40).to_double())
    de_status.push(if i < 30 { 1 } else { 0 })
    i = i + 1
  }

  let n_genes = gene_names.length()
  let n_de = {
    let mut count = 0
    let mut j = 0
    while j < n_genes {
      if de_status[j] == 1 { count = count + 1 }
      j = j + 1
    }
    count
  }

  println("Total genes: \{n_genes.to_string()}")
  println("DE genes: \{n_de.to_string()}")

  // 2. Calculate PWF (Probability Weighting Function)
  println("")
  println("--- Step 2: Calculate PWF ---")

  let pwf = @src.goseq_nullp(gene_lengths, de_status)
  println("PWF calculated for \{pwf.length()} genes")

  // Show PWF values for first few genes
  println("")
  println("Sample PWF values (first 5 genes):")
  let mut j = 0
  while j < 5 {
    println(
      "  \{gene_names[j]}: length=\{gene_lengths[j].to_string()} bp, DE=\{de_status[j].to_string()}, PWF=\{pwf[j].to_string()}",
    )
    j = j + 1
  }

  // 3. Build GO term sets
  println("")
  println("--- Step 3: Build GO Term Sets ---")

  let go_sets = [
    @src.GoseqGeneSet::new(
      "GO:0006915",
      "apoptotic process",
      "BP",
      // Contains mostly short genes (few DE)
      ["ShortGene_0", "ShortGene_1", "ShortGene_2", "ShortGene_3", "ShortGene_4",
       "ShortGene_5", "ShortGene_10", "ShortGene_20", "ShortGene_30", "ShortGene_40"],
    ),
    @src.GoseqGeneSet::new(
      "GO:0007049",
      "cell cycle",
      "BP",
      // Contains mostly medium genes (moderate DE)
      ["MedGene_0", "MedGene_1", "MedGene_2", "MedGene_3", "MedGene_4",
       "MedGene_10", "MedGene_20", "MedGene_30", "MedGene_40"],
    ),
    @src.GoseqGeneSet::new(
      "GO:0008150",
      "biological process",
      "BP",
      // Contains mostly long genes (high DE)
      ["LongGene_0", "LongGene_1", "LongGene_2", "LongGene_3", "LongGene_4",
       "LongGene_10", "LongGene_20", "LongGene_30", "LongGene_40",
       "LongGene_5", "LongGene_15", "LongGene_25"],
    ),
    @src.GoseqGeneSet::new(
      "GO:0003674",
      "molecular function",
      "MF",
      // Mixed composition
      ["ShortGene_0", "MedGene_0", "LongGene_0", "ShortGene_10", "MedGene_10",
       "LongGene_10", "ShortGene_20", "MedGene_20", "LongGene_20"],
    ),
    @src.GoseqGeneSet::new(
      "GO:0005575",
      "cellular component",
      "CC",
      // Only non-DE genes
      ["ShortGene_40", "ShortGene_41", "ShortGene_42", "ShortGene_43",
       "MedGene_40", "MedGene_41", "MedGene_42"],
    ),
  ]

  println("Built \{go_sets.length().to_string()} GO term sets")

  // 4. Run goseq enrichment analysis
  println("")
  println("--- Step 4: Run goseq Analysis ---")

  let results = @src.goseq(gene_names, gene_lengths, de_status, go_sets)
  println("Results: \{results.length().to_string()} GO terms")

  // 5. Display results
  println("")
  println("--- Step 5: Results ---")

  let mut k = 0
  while k < results.length() {
    let r = results[k]
    println(
      "  \{r.go_id} (\{r.ontology}) \{r.term_name}:",
    )
    println(
      "    p-value=\{r.p_value.to_string()}, adjusted=\{r.adjusted_p_value.to_string()}",
    )
    println(
      "    overlap=\{r.overlap.to_string()}, DE in term=\{r.num_de_in_term.to_string()}, term size=\{r.num_in_term.to_string()}",
    )
    k = k + 1
  }

  // 6. Filter significant results
  println("")
  println("--- Step 6: Filter Significant Results (p.adj <= 0.05, overlap >= 2) ---")

  let sig = @src.goseq_filter_results(results, 0.05, 2)
  if sig.length() == 0 {
    println("No significant terms at this threshold")
    // Relaxed threshold
    let relaxed = @src.goseq_filter_results(results, 0.5, 1)
    println("Relaxed threshold (p.adj <= 0.5, overlap >= 1): \{relaxed.length().to_string()} terms")
    let mut m = 0
    while m < relaxed.length() {
      let r = relaxed[m]
      println(
        "  \{r.go_id} \{r.term_name}: p=\{r.p_value.to_string()}, adj=\{r.adjusted_p_value.to_string()}, overlap=\{r.overlap.to_string()}",
      )
      m = m + 1
    }
  } else {
    let mut m = 0
    while m < sig.length() {
      let r = sig[m]
      println(
        "  \{r.go_id} \{r.term_name}: p=\{r.p_value.to_string()}, adj=\{r.adjusted_p_value.to_string()}, overlap=\{r.overlap.to_string()}",
      )
      m = m + 1
    }
  }

  // 7. Demonstrate length bias with Wallenius
  println("")
  println("--- Step 7: Wallenius vs Hypergeometric Comparison ---")

  // Compare p-values for a term with bias
  // GO:0008150 (biological_process) has mostly long genes
  let term_with_bias = go_sets[2]
  let mut term_de = 0
  let mut term_total = 0
  let mut gi = 0
  while gi < term_with_bias.gene_ids.length() {
    let mut found = false
    let mut g = 0
    while g < n_genes {
      if gene_names[g] == term_with_bias.gene_ids[gi] {
        if de_status[g] == 1 { term_de = term_de + 1 }
        term_total = term_total + 1
        found = true
        g = n_genes // break
      }
      g = g + 1
    }
    gi = gi + 1
  }

  let wallenius_p = @src.wallenius_pvalue(term_de, term_total, n_de, n_genes, 1.5)
  let hypergeo_p = @src.hypergeo_pvalue(term_de, term_total, n_de, n_genes)
  println("Term: \{term_with_bias.term_name} (GO term with length bias)")
  println("  DE in term = \{term_de.to_string()}, term size = \{term_total.to_string()}")
  println("  Wallenius p-value (odds=1.5): \{wallenius_p.to_string()}")
  println("  Hypergeometric p-value:       \{hypergeo_p.to_string()}")
  println("  Difference: \{(hypergeo_p - wallenius_p).to_string()}")

  println("")
  println("Done!")
}