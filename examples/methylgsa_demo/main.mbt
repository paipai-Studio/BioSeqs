///|
fn main {
  println("=== methylGSA Demo ===")

  let input = @src.methylgsa_sample_input()
  println("\n1. Input: " + input.cpg_ids.length().to_string() + " CpGs")
  for i in 0..<input.cpg_ids.length() {
    if i < input.cpg_pvalues.length() {
      println("   " + input.cpg_ids[i] + ": p=" + input.cpg_pvalues[i].to_string())
    }
  }

  println("\n2. Per-gene minimal CpG p-values")
  let min_p = @src.methylgsa_gene_min_p(input)
  for g in min_p.keys() {
    match min_p.get(g) {
      Some(p) => println("   " + g + ": " + p.to_string())
      None => ()
    }
  }

  println("\n3. Gene sets")
  let sets = @src.methylgsa_sample_gene_sets()
  for name in sets.keys() {
    match sets.get(name) {
      Some(genes) => {
        let mut line = name + ": "
        for g in genes {
          line = line + g + " "
        }
        println("   " + line)
      }
      None => ()
    }
  }

  println("\n4. methylRRA-OLS results")
  let ols_results = @src.methylrra_ols(input, gene_sets=sets)
  for r in ols_results {
    println(
      "   " + r.set_name +
      " | size=" + r.set_size.to_string() +
      " | stat=" + r.statistic.to_string() +
      " | p=" + r.pvalue.to_string() +
      " | q=" + r.qvalue.to_string(),
    )
  }

  println("\n5. methylRRA-Robust results")
  let robust_results = @src.methylrra_robust(input, gene_sets=sets)
  for r in robust_results {
    println(
      "   " + r.set_name +
      " | size=" + r.set_size.to_string() +
      " | stat=" + r.statistic.to_string() +
      " | p=" + r.pvalue.to_string() +
      " | q=" + r.qvalue.to_string(),
    )
  }

  println("\n6. BH-FDR correction")
  let results : Array[@src.MethylgsaResult] = Array::new()
  results.push(
    @src.methylgsa_result_new(
      set_name="SetA", set_size=10, statistic=5.0,
      pvalue=0.001, qvalue=1.0,
    ),
  )
  results.push(
    @src.methylgsa_result_new(
      set_name="SetB", set_size=10, statistic=3.0,
      pvalue=0.01, qvalue=1.0,
    ),
  )
  results.push(
    @src.methylgsa_result_new(
      set_name="SetC", set_size=10, statistic=1.0,
      pvalue=0.05, qvalue=1.0,
    ),
  )
  let adjusted = @src.methylgsa_bh(results)
  for r in adjusted {
    println("   " + r.set_name + ": p=" + r.pvalue.to_string() + " -> q=" + r.qvalue.to_string())
  }

  println("\n7. Empty set handling")
  let empty_sets : Map[String, Array[String]] = Map([], capacity=1)
  empty_sets.set("EMPTY", ["NO_SUCH_GENE"])
  let empty_res = @src.methylrra_robust(input, gene_sets=empty_sets)
  for r in empty_res {
    println("   " + r.set_name + ": size=" + r.set_size.to_string() + ", p=" + r.pvalue.to_string())
  }
}