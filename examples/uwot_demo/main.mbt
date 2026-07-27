///|
/// UMAP Dimensionality Reduction Demo
///
/// This demo shows how to use the uwot module for
/// Uniform Manifold Approximation and Projection (UMAP) analysis.
fn main {
  println("=== UMAP Dimensionality Reduction Demo ===")
  println("")

  println("1. Default UMAP configuration:")
  let config = @src.UmapConfig::new()
  println("   n_neighbors: \{config.n_neighbors}")
  println("   n_components: \{config.n_components}")
  println("   n_epochs: \{config.n_epochs}")
  println("   min_dist: \{config.min_dist}")
  println("   spread: \{config.spread}")
  println("   learning_rate: \{config.learning_rate}")
  println("   metric: \{config.metric}")
  println("")

  println("2. Distance matrix calculation:")
  let sample_data = [
    [0.0, 0.0, 0.0],
    [3.0, 4.0, 0.0],
    [1.0, 1.0, 1.0],
  ]
  let dist_matrix = @src.uwot_distance_matrix(sample_data)
  println("   Samples: \{sample_data.length()}")
  println("   Features per sample: \{sample_data[0].length()}")
  println("   Distance matrix:")
  let mut i = 0
  while i < dist_matrix.length() {
    println("     \{dist_matrix[i]}")
    i = i + 1
  }
  println("")

  println("3. UMAP on clustered test data:")
  let test_data = @src.create_umap_test_data(20, 5)
  println("   Test data samples: \{test_data.length()}")
  println("   Test data features: \{test_data[0].length()}")
  println("   Data has 4 clusters")
  println("")

  println("4. Running UMAP (reduced epochs for demo):")
  let test_config = @src.UmapConfig::new_custom(5, 2, 100, 42)
  let umap_result = @src.uwot_umap(test_data, test_config)
  println("   Embedding shape: \{umap_result.embedding.length()} x \{umap_result.embedding[0].length()}")
  println("   Epochs run: \{umap_result.n_epochs}")
  println("")

  println("5. First few embedding coordinates:")
  let mut j = 0
  while j < 5 && j < umap_result.embedding.length() {
    println("   Sample \{j}: [\{umap_result.embedding[j][0]}, \{umap_result.embedding[j][1]}]")
    j = j + 1
  }
  println("")

  println("6. UMAP with 3 components:")
  let config_3d = @src.UmapConfig::new_custom(5, 3, 100, 42)
  let result_3d = @src.uwot_umap(test_data, config_3d)
  println("   Embedding shape: \{result_3d.embedding.length()} x \{result_3d.embedding[0].length()}")
  println("   Sample 0: [\{result_3d.embedding[0][0]}, \{result_3d.embedding[0][1]}, \{result_3d.embedding[0][2]}]")
  println("")

  println("=== Demo Complete ===")
}
