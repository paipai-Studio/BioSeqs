///|
/// t-SNE Dimensionality Reduction Demo
///
/// This demo shows how to use the Rtsne module for
/// t-Distributed Stochastic Neighbor Embedding (t-SNE) analysis.
fn main {
  println("=== t-SNE Dimensionality Reduction Demo ===")
  println("")

  println("1. Default t-SNE configuration:")
  let config = @src.TsneConfig::new()
  println("   Perplexity: \{config.perplexity}")
  println("   Max iterations: \{config.max_iter}")
  println("   Dimensions: \{config.dims}")
  println("   Learning rate (eta): \{config.eta}")
  println("   Theta (Barnes-Hut): \{config.theta}")
  println("")

  println("2. Distance matrix calculation:")
  let sample_data = [
    [0.0, 0.0, 0.0],
    [3.0, 4.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let dist_matrix = @src.calc_distance_matrix(sample_data)
  println("   Samples: \{sample_data.length()}")
  println("   Features per sample: \{sample_data[0].length()}")
  println("   Distance matrix:")
  let mut i = 0
  while i < dist_matrix.length() {
    println("     \{dist_matrix[i]}")
    i = i + 1
  }
  println("")

  println("3. t-SNE on clustered test data:")
  let test_data = @src.create_tsne_test_data(20, 5)
  println("   Test data samples: \{test_data.length()}")
  println("   Test data features: \{test_data[0].length()}")
  println("   Data has 3 clusters")
  println("")

  println("4. Running t-SNE (reduced iterations for demo):")
  let test_config = @src.TsneConfig::new_custom(5.0, 200, 2, 42)
  let tsne_result = @src.tsne(test_data, test_config)
  println("   Embedding shape: \{tsne_result.embedding.length()} x \{tsne_result.embedding[0].length()}")
  println("   Iterations run: \{tsne_result.n_iter}")
  println("   Final cost: \{tsne_result.costs[tsne_result.costs.length() - 1]}")
  println("")

  println("5. First few embedding coordinates:")
  let mut j = 0
  while j < 5 && j < tsne_result.embedding.length() {
    println("   Sample \{j}: [\{tsne_result.embedding[j][0]}, \{tsne_result.embedding[j][1]}]")
    j = j + 1
  }
  println("")

  println("=== Demo Complete ===")
}
