///|
/// MOFA2 multi-omics factor analysis demo.
/// Demonstrates sparse Bayesian factor analysis on multi-omics data.

fn main {
  println("=== MOFA2 Multi-Omics Factor Analysis Demo ===")
  println("")

  // Create example data (RNA-seq + ATAC-seq)
  let n_samples = 50
  let views = @bio.mofa_create_example(n_samples, 5)

  println("Input: \{n_samples.to_string()} samples, \{views.length().to_string()} omics views")
  {
    let mut i = 0
    while i < views.length() {
      let view = views[i]
      println("  View \{i.to_string()}: \{view.view_name} (\{view.n_features().to_string()} features)")
      i = i + 1
    }
  }

  // Configure MOFA2
  let params = @bio.MofaParams::create(
    n_factors = 5,
    max_iterations = 100,
    tolerance = 0.01,
    seed = 42,
  )

  // Run MOFA2
  println("")
  println("Running MOFA2 (EM algorithm)...")
  let result = @bio.mofa_run(views, params)

  println("")
  println("Results:")
  println("  Active factors: \{result.active_factors.to_string()}/5")
  println("  Iterations: \{result.n_iterations.to_string()}")
  println("  Converged: \{result.converged.to_string()}")

  // Print factor variance explained
  println("")
  println("Variance explained per factor:")
  {
    let mut f = 0
    while f < result.factor_vars.length() {
      let pct = result.factor_vars[f] * 100.0
      println("  Factor \{(f + 1).to_string()}: \{pct.to_string()}%")
      f = f + 1
    }
  }

  // Print variance explained per view
  println("")
  println("Total variance explained per view:")
  {
    let mut vi = 0
    while vi < result.view_vars.length() {
      let view = views[vi]
      let mut sum = 0.0
      for v in result.view_vars[vi] {
        sum = sum + v
      }
      let pct = sum * 100.0
      println("  \{view.view_name}: \{pct.to_string()}%")
      vi = vi + 1
    }
  }

  // Print sample factors (first 5 samples)
  println("")
  println("Factor activations (first 5 samples):")
  {
    let mut s = 0
    while s < 5 && s < result.factors.length() {
      let mut line = "  Sample " + (s + 1).to_string() + ": ["
      let mut f = 0
      while f < result.factors[s].length() {
        if f > 0 {
          line = line + ", "
        }
        line = line + result.factors[s][f].to_string()
        f = f + 1
      }
      line = line + "]"
      println(line)
      s = s + 1
    }
  }

  println("")
  println("Done!")
}
