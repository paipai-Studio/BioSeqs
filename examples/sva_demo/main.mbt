///|
fn main {
  println("=== SVA (Surrogate Variable Analysis) Demo ===")

  let (dat, genes, samples, batch) = @bio.sample_expression_matrix()

  println("")
  println("Input Data Summary:")
  println("  Genes: " + genes.length().to_string())
  println("  Samples: " + samples.length().to_string())
  println("  Batches: 2 (0 and 1)")

  println("")
  println("Sample expression values (first 5 genes, first 6 samples):")
  let mut i = 0
  while i < 5 {
    let mut j = 0
    let mut line = genes[i] + ": "
    while j < 6 {
      line = line + dat[i][j].to_string() + " "
      j = j + 1
    }
    println("  " + line)
    i = i + 1
  }

  println("")
  println("=== Running ComBat Batch Correction ===")
  let combat_result = @bio.combat(dat, batch)

  println("")
  println("After ComBat correction (first 5 genes, first 6 samples):")
  i = 0
  while i < 5 {
    let mut j = 0
    let mut line = genes[i] + ": "
    while j < 6 {
      line = line + combat_result.adjusted[i][j].to_string() + " "
      j = j + 1
    }
    println("  " + line)
    i = i + 1
  }

  println("")
  println("=== Running Simple Batch Correction ===")
  let _simple_adj = @bio.combat_simple(dat, batch)
  println("Simple batch correction completed.")

  println("")
  println("=== Running SVA (Surrogate Variable Analysis) ===")
  let mod = [
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
    [1.0],
  ]
  let sva_result = @bio.sva(dat, mod, n_sv=2)

  println("")
  println("SVA Results:")
  println("  Number of surrogate variables: " + sva_result.n_sv.to_string())
  println(
    "  Number of genes tested: " + sva_result.pprob_gam.length().to_string(),
  )

  let mut n_high = 0
  i = 0
  while i < sva_result.pprob_gam.length() {
    if sva_result.pprob_gam[i] > 0.5 {
      n_high = n_high + 1
    }
    i = i + 1
  }
  println(
    "  Genes with high probability of being affected: " + n_high.to_string(),
  )

  println("")
  println("=== Demo Complete ===")
}
