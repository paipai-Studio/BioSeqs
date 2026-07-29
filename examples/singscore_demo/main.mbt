///|
/// Singscore gene set scoring demo.
/// Demonstrates rank-based single-sample gene set scoring.

fn main {
  println("=== Singscore Gene Set Scoring Demo ===")
  println("")

  // Create example data
  let (sample, spec) = @bio.singscore_create_example()

  println("Sample: \{sample.sample_id}")
  println("Genes: \{sample.gene_names.length().to_string()}")
  println("Gene set: \{spec.name}")
  println("  Up-regulated: \{spec.up_genes.length().to_string()} genes")
  println("  Down-regulated: \{spec.down_genes.length().to_string()} genes")

  // Score the sample
  let result = @bio.singscore_score(sample, spec, sample.sample_id)

  println("")
  println("Results:")
  println("  Sample: \{result.sample_id}")
  println("  Score: \{result.score.to_string()}")
  println("  Up score: \{result.up_score.to_string()}")
  println("  Down score: \{result.down_score.to_string()}")
  println("  Up dispersion: \{result.up_dispersion.to_string()}")
  println("  Down dispersion: \{result.down_dispersion.to_string()}")
  println("  Genes found (up): \{result.n_up_found.to_string()}")
  println("  Genes found (down): \{result.n_down_found.to_string()}")

  // Create multiple samples with different expression patterns
  println("")
  println("=== Multiple Sample Scoring ===")
  println("")

  let genes = sample.gene_names
  let n_genes = genes.length()
  let samples : Array[@bio.SampleExpression] = Array::new()

  // Activated immune cell (high expression)
  {
    let expr : Array[Double] = Array::make(n_genes, 5.0)
    samples.push(@bio.SampleExpression::new("activated", genes, expr))
  }

  // Resting cell (low expression)
  {
    let expr : Array[Double] = Array::make(n_genes, 0.5)
    samples.push(@bio.SampleExpression::new("resting", genes, expr))
  }

  // Mixed expression
  {
    let expr : Array[Double] = Array::make(n_genes, 2.0)
    samples.push(@bio.SampleExpression::new("mixed", genes, expr))
  }

  // Score each sample
  for s in samples {
    let r = @bio.singscore_score(s, spec, s.sample_id)
    println("  \{r.sample_id}: score=\{r.score.to_string()}, up=\{r.up_score.to_string()}, down=\{r.down_score.to_string()}")
  }

  println("")
  println("Done!")
}
