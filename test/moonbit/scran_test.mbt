///|
@test "scran_compute_sum_factors" {
  let (counts, gene_names) = @src.create_example_scran_data()
  let sizes = [10]
  let sf = @src.compute_sum_factors(counts, sizes)
  assert_eq(sf.factors.length(), 100)
}

@test "scran_log_normalize" {
  let (counts, gene_names) = @src.create_example_scran_data()
  let size_factors = Array::make(100, 1.0)
  let normalized = @src.log_normalize(counts, size_factors)
  assert_eq(normalized.length(), 100)
  assert_eq(normalized[0].length(), 500)
}

@test "scran_build_snn_graph" {
  let (counts, gene_names) = @src.create_example_scran_data()
  let size_factors = Array::make(100, 1.0)
  let normalized = @src.log_normalize(counts, size_factors)
  let (adjacency, weights) = @src.build_snn_graph(normalized, 10)
  assert_eq(adjacency.length(), 100)
  assert_eq(weights.length(), 100)
}

@test "scran_leiden_clustering" {
  let (counts, gene_names) = @src.create_example_scran_data()
  let size_factors = Array::make(100, 1.0)
  let normalized = @src.log_normalize(counts, size_factors)
  let (adjacency, weights) = @src.build_snn_graph(normalized, 10)
  let clusters = @src.leiden_clustering(adjacency, weights, 1.0)
  assert_eq(clusters.clusters.length(), 100)
  assert_true(clusters.n_clusters >= 1)
}

@test "scran_score_markers" {
  let (counts, gene_names) = @src.create_example_scran_data()
  let clusters = Array::make(100, 0)
  let mut i = 0
  while i < 100 {
    clusters[i] = i / 20
    i = i + 1
  }
  let markers = @src.score_markers(counts, clusters, gene_names)
  assert_true(markers.gene.length() > 0)
  assert_eq(markers.p_value.length(), markers.gene.length())
}

@test "scran_create_example_data" {
  let (counts, gene_names) = @src.create_example_scran_data()
  assert_eq(counts.length(), 100)
  assert_eq(gene_names.length(), 500)
  assert_eq(counts[0].length(), 500)
}
