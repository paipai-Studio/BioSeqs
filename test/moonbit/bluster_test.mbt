///|
fn bluster_test_close(
  left : Double,
  right : Double,
  tolerance : Double,
) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn bluster_test_some(value : Double?) -> Double {
  match value {
    Some(actual) => actual
    None => abort("expected a numeric value")
  }
}

///|
fn bluster_test_unique_count(values : Array[Int]) -> Int {
  let unique : Array[Int] = []
  for value in values {
    let mut found = false
    for existing in unique {
      if existing == value {
        found = true
        break
      }
    }
    if !found {
      unique.push(value)
    }
  }
  unique.length()
}

///|
fn bluster_test_two_groups() -> Array[Array[Double]] {
  [[0.0, 0.0], [0.1, 0.0], [0.2, 0.1], [10.0, 10.0], [10.1, 10.0], [10.2, 10.1]]
}

///|
fn bluster_test_sce() -> @src.SingleCellExperiment {
  let experiment = @src.SingleCellExperiment::new(
    [[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]],
    ["G1", "G2"],
    ["C1", "C2", "C3", "C4"],
  )
  experiment.row_data["symbol"] = ["A", "B"]
  experiment.col_data["batch"] = ["A", "A", "B", "B"]
  experiment.reduced_dims["PCA"] = [
    [0.0, 0.0],
    [0.1, 0.0],
    [10.0, 10.0],
    [10.1, 10.0],
  ]
  experiment.metadata["project"] = "bluster-test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[5.0, 6.0, 7.0, 8.0]],
    ["Spike1"],
    ["C1", "C2", "C3", "C4"],
  )
  experiment
}

///|
test "bluster: distance and weighting accessors expose enum values" {
  assert_true(
    @src.bluster_euclidean_distance() != @src.bluster_manhattan_distance(),
  )
  assert_true(
    @src.bluster_manhattan_distance() != @src.bluster_cosine_distance(),
  )
  assert_true(@src.bluster_rank_weight() != @src.bluster_number_weight())
  assert_true(@src.bluster_number_weight() != @src.bluster_jaccard_weight())
}

///|
test "bluster: link denominator accessors expose enum values" {
  assert_true(@src.bluster_link_minimum() != @src.bluster_link_maximum())
  assert_true(@src.bluster_link_maximum() != @src.bluster_link_union())
}

///|
test "bluster: k-means configuration preserves explicit settings" {
  let config = @src.BlusterKmeansConfig::create(
    3,
    max_iterations=25,
    starts=4,
    tolerance=1.0e-6,
    seed=17,
  ) catch {
    _ => abort("valid k-means configuration should build")
  }
  assert_eq(config.centers, 3)
  assert_eq(config.max_iterations, 25)
  assert_eq(config.starts, 4)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.seed, 17)
}

///|
test "bluster: k-means configuration rejects invalid controls" {
  let mut failures = 0
  ignore(@src.BlusterKmeansConfig::create(0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterKmeansConfig::create(2, starts=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterKmeansConfig::create(2, max_iterations=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterKmeansConfig::create(2, tolerance=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterKmeansConfig::create(2, seed=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 5)
}

///|
test "bluster: graph configuration preserves explicit settings" {
  let config = @src.BlusterGraphConfig::create(
    k=7,
    shared=false,
    snn_weight=@src.bluster_jaccard_weight(),
    distance=@src.bluster_manhattan_distance(),
    resolution=0.5,
    max_iterations=30,
    seed=19,
  ) catch {
    _ => abort("valid graph configuration should build")
  }
  assert_eq(config.k, 7)
  assert_false(config.shared)
  assert_true(config.snn_weight == @src.bluster_jaccard_weight())
  assert_true(config.distance == @src.bluster_manhattan_distance())
  assert_eq(config.resolution, 0.5)
  assert_eq(config.max_iterations, 30)
  assert_eq(config.seed, 19)
}

///|
test "bluster: graph configuration rejects invalid controls" {
  let mut failures = 0
  ignore(@src.BlusterGraphConfig::create(k=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterGraphConfig::create(resolution=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterGraphConfig::create(max_iterations=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterGraphConfig::create(seed=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "bluster: two-step configuration validates both stages" {
  let valid = @src.BlusterTwoStepConfig::create(
    first_centers=4,
    first_starts=2,
    second_k=3,
    seed=23,
  ) catch {
    _ => abort("valid two-step configuration should build")
  }
  assert_eq(valid.first_centers, 4)
  assert_eq(valid.first_starts, 2)
  assert_eq(valid.second_k, 3)
  let mut failures = 0
  ignore(@src.BlusterTwoStepConfig::create(first_centers=-1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterTwoStepConfig::create(first_starts=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.BlusterTwoStepConfig::create(second_k=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "bluster: directed KNN uses exact nearest neighbors" {
  let graph = @src.bluster_make_knn_graph(
    [[0.0], [1.0], [3.0], [10.0]],
    k=1,
    directed=true,
  ) catch {
    _ => abort("valid directed KNN graph should build")
  }
  assert_true(graph.directed)
  assert_eq(graph.k, 1)
  assert_eq(graph.graph_kind, "KNN")
  assert_eq(graph.node_count(), 4)
  assert_eq(graph.edge_count(), 4)
  assert_eq(graph.edge_weight(0, 1), 1.0)
  assert_eq(graph.edge_weight(1, 0), 1.0)
  assert_eq(graph.edge_weight(2, 1), 1.0)
  assert_eq(graph.edge_weight(3, 2), 1.0)
  assert_eq(graph.edge_weight(1, 2), 0.0)
}

///|
test "bluster: KNN tie ordering is deterministic by observation index" {
  let graph = @src.bluster_make_knn_graph(
    [[0.0], [2.0], [-2.0]],
    k=1,
    directed=true,
  ) catch {
    _ => abort("tied KNN graph should build")
  }
  assert_eq(graph.edge_weight(0, 1), 1.0)
  assert_eq(graph.edge_weight(0, 2), 0.0)
}

///|
test "bluster: undirected KNN symmetrizes one-way relationships" {
  let graph = @src.bluster_make_knn_graph([[0.0], [1.0], [3.0], [10.0]], k=1) catch {
    _ => abort("valid undirected KNN graph should build")
  }
  assert_false(graph.directed)
  assert_eq(graph.edge_count(), 3)
  assert_eq(graph.edge_weight(1, 2), 1.0)
  assert_eq(graph.edge_weight(2, 1), 1.0)
}

///|
test "bluster: KNN caps k at the available observation count" {
  let graph = @src.bluster_make_knn_graph([[0.0], [1.0], [2.0]], k=20) catch {
    _ => abort("oversized k should be capped")
  }
  assert_eq(graph.k, 2)
  assert_eq(graph.edge_count(), 3)
}

///|
test "bluster: rank SNN weighting includes self at rank zero" {
  let graph = @src.bluster_make_snn_graph(
    [[0.0], [1.0], [10.0]],
    k=1,
    weighting=@src.bluster_rank_weight(),
  ) catch {
    _ => abort("rank SNN graph should build")
  }
  assert_eq(graph.graph_kind, "SNN")
  assert_eq(graph.edge_count(), 3)
  assert_true(bluster_test_close(graph.edge_weight(0, 1), 0.5, 1.0e-12))
  assert_true(bluster_test_close(graph.edge_weight(0, 2), 1.0e-6, 1.0e-15))
  assert_true(bluster_test_close(graph.edge_weight(1, 2), 0.5, 1.0e-12))
}

///|
test "bluster: number SNN weighting counts shared neighborhoods" {
  let graph = @src.bluster_make_snn_graph(
    [[0.0], [1.0], [10.0]],
    k=1,
    weighting=@src.bluster_number_weight(),
  ) catch {
    _ => abort("number SNN graph should build")
  }
  assert_eq(graph.edge_weight(0, 1), 2.0)
  assert_eq(graph.edge_weight(0, 2), 1.0)
  assert_eq(graph.edge_weight(1, 2), 1.0)
}

///|
test "bluster: Jaccard SNN weighting uses neighborhood unions" {
  let graph = @src.bluster_make_snn_graph(
    [[0.0], [1.0], [10.0]],
    k=1,
    weighting=@src.bluster_jaccard_weight(),
  ) catch {
    _ => abort("Jaccard SNN graph should build")
  }
  assert_eq(graph.edge_weight(0, 1), 1.0)
  assert_true(bluster_test_close(graph.edge_weight(0, 2), 1.0 / 3.0, 1.0e-12))
  assert_true(bluster_test_close(graph.edge_weight(1, 2), 1.0 / 3.0, 1.0e-12))
}

///|
test "bluster: SNN caps k and forms complete full-neighborhood graph" {
  let graph = @src.bluster_make_snn_graph(
    [[0.0], [1.0], [2.0]],
    k=20,
    weighting=@src.bluster_number_weight(),
  ) catch {
    _ => abort("oversized SNN k should be capped")
  }
  assert_eq(graph.k, 2)
  assert_eq(graph.edge_count(), 3)
  assert_eq(graph.edge_weight(0, 2), 3.0)
}

///|
test "bluster: graph constructors reject invalid matrices and k" {
  let mut failures = 0
  ignore(@src.bluster_make_snn_graph([], k=1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.bluster_make_knn_graph([[1.0], [1.0, 2.0]], k=1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.bluster_make_snn_graph([[1.0e301]], k=1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.bluster_make_knn_graph([[1.0]], k=0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "bluster: k-means separates clear groups and reports WSS" {
  let config = @src.BlusterKmeansConfig::create(2, starts=5, seed=31) catch {
    _ => abort("valid k-means configuration should build")
  }
  let result = @src.bluster_cluster_rows_kmeans(
    [[0.0, 0.0], [0.0, 1.0], [10.0, 10.0], [10.0, 11.0]],
    config,
  ) catch {
    _ => abort("k-means should cluster valid data")
  }
  assert_eq(result.clusters.length(), 4)
  assert_eq(result.centers.length(), 2)
  assert_true(result.clusters[0] == result.clusters[1])
  assert_true(result.clusters[2] == result.clusters[3])
  assert_true(result.clusters[0] != result.clusters[2])
  assert_true(bluster_test_close(result.total_within_sum_squares, 1.0, 1.0e-12))
  assert_true(
    bluster_test_close(
      result.within_sum_squares[0] + result.within_sum_squares[1],
      result.total_within_sum_squares,
      1.0e-12,
    ),
  )
}

///|
test "bluster: k-means is reproducible for a fixed seed" {
  let config = @src.BlusterKmeansConfig::create(2, starts=4, seed=37) catch {
    _ => abort("valid k-means configuration should build")
  }
  let data = bluster_test_two_groups()
  let first = @src.bluster_cluster_rows_kmeans(data, config) catch {
    _ => abort("first k-means run should work")
  }
  let second = @src.bluster_cluster_rows_kmeans(data, config) catch {
    _ => abort("second k-means run should work")
  }
  assert_eq(first.clusters, second.clusters)
  assert_eq(first.centers, second.centers)
  assert_eq(first.total_within_sum_squares, second.total_within_sum_squares)
}

///|
test "bluster: multiple starts cannot worsen the selected WSS" {
  let one = @src.BlusterKmeansConfig::create(3, starts=1, seed=41) catch {
    _ => abort("single-start configuration should build")
  }
  let many = @src.BlusterKmeansConfig::create(3, starts=8, seed=41) catch {
    _ => abort("multi-start configuration should build")
  }
  let data = [
    [0.0, 0.0],
    [0.2, 0.0],
    [5.0, 5.0],
    [5.2, 5.0],
    [10.0, 0.0],
    [10.2, 0.0],
  ]
  let first = @src.bluster_cluster_rows_kmeans(data, one) catch {
    _ => abort("single-start k-means should work")
  }
  let second = @src.bluster_cluster_rows_kmeans(data, many) catch {
    _ => abort("multi-start k-means should work")
  }
  assert_true(
    second.total_within_sum_squares <= first.total_within_sum_squares + 1.0e-12,
  )
}

///|
test "bluster: k-means rejects more centers than observations" {
  let config = @src.BlusterKmeansConfig::create(3) catch {
    _ => abort("configuration validation is independent of data")
  }
  let mut raised = false
  ignore(@src.bluster_cluster_rows_kmeans([[0.0], [1.0]], config)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "bluster: graph clustering separates disconnected populations" {
  let config = @src.BlusterGraphConfig::create(k=1, resolution=0.5, seed=43) catch {
    _ => abort("valid graph configuration should build")
  }
  let result = @src.bluster_cluster_rows_graph(
    bluster_test_two_groups(),
    config,
  ) catch {
    _ => abort("graph clustering should work")
  }
  assert_eq(result.clusters.length(), 6)
  assert_true(result.clusters[0] == result.clusters[1])
  assert_true(result.clusters[1] == result.clusters[2])
  assert_true(result.clusters[3] == result.clusters[4])
  assert_true(result.clusters[4] == result.clusters[5])
  assert_true(result.clusters[0] != result.clusters[3])
  assert_eq(bluster_test_unique_count(result.clusters), 2)
}

///|
test "bluster: KNN graph clustering uses the non-shared path" {
  let config = @src.BlusterGraphConfig::create(
    k=1,
    shared=false,
    resolution=0.5,
    seed=47,
  ) catch {
    _ => abort("valid KNN graph configuration should build")
  }
  let result = @src.bluster_cluster_rows_graph(
    bluster_test_two_groups(),
    config,
  ) catch {
    _ => abort("KNN graph clustering should work")
  }
  assert_eq(result.graph.graph_kind, "KNN")
  assert_true(result.clusters[0] != result.clusters[3])
}

///|
test "bluster: graph clustering is reproducible for a fixed seed" {
  let config = @src.BlusterGraphConfig::create(k=1, seed=53) catch {
    _ => abort("valid graph configuration should build")
  }
  let first = @src.bluster_cluster_rows_graph(bluster_test_two_groups(), config) catch {
    _ => abort("first graph clustering should work")
  }
  let second = @src.bluster_cluster_rows_graph(
    bluster_test_two_groups(),
    config,
  ) catch {
    _ => abort("second graph clustering should work")
  }
  assert_eq(first.clusters, second.clusters)
  assert_eq(first.modularity, second.modularity)
}

///|
test "bluster: higher resolution does not coarsen this graph" {
  let low = @src.BlusterGraphConfig::create(k=1, resolution=0.2, seed=59) catch {
    _ => abort("low-resolution configuration should build")
  }
  let high = @src.BlusterGraphConfig::create(k=1, resolution=3.0, seed=59) catch {
    _ => abort("high-resolution configuration should build")
  }
  let data = bluster_test_two_groups()
  let low_result = @src.bluster_cluster_rows_graph(data, low) catch {
    _ => abort("low-resolution clustering should work")
  }
  let high_result = @src.bluster_cluster_rows_graph(data, high) catch {
    _ => abort("high-resolution clustering should work")
  }
  assert_true(
    bluster_test_unique_count(high_result.clusters) >=
    bluster_test_unique_count(low_result.clusters),
  )
}

///|
test "bluster: a one-node graph remains one cluster" {
  let config = @src.BlusterGraphConfig::create(k=10) catch {
    _ => abort("valid graph configuration should build")
  }
  let result = @src.bluster_cluster_rows_graph([[1.0, 2.0]], config) catch {
    _ => abort("one-node graph clustering should work")
  }
  assert_eq(result.clusters, [0])
  assert_eq(result.graph.k, 0)
  assert_eq(result.modularity, 0.0)
}

///|
test "bluster: two-step clustering maps centroid labels to observations" {
  let config = @src.BlusterTwoStepConfig::create(
    first_centers=4,
    first_starts=4,
    second_k=1,
    resolution=0.5,
    seed=61,
  ) catch {
    _ => abort("valid two-step configuration should build")
  }
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.1],
    [0.3, 0.1],
    [10.0, 10.0],
    [10.1, 10.0],
    [10.2, 10.1],
    [10.3, 10.1],
  ]
  let result = @src.bluster_cluster_rows_two_step(data, config) catch {
    _ => abort("two-step clustering should work")
  }
  assert_eq(result.clusters.length(), data.length())
  assert_eq(result.centroids.length(), 4)
  assert_eq(result.centroids, result.first.centers)
  assert_eq(result.centroid_clusters.length(), 4)
  for observation in 0..<data.length() {
    assert_eq(
      result.clusters[observation],
      result.centroid_clusters[result.first.clusters[observation]],
    )
  }
}

///|
test "bluster: automatic two-step center count is rounded square root" {
  let config = @src.BlusterTwoStepConfig::create(second_k=1, seed=67) catch {
    _ => abort("automatic two-step configuration should build")
  }
  let result = @src.bluster_cluster_rows_two_step(
    [[0.0], [1.0], [2.0], [10.0], [11.0], [12.0], [20.0], [21.0], [22.0]],
    config,
  ) catch {
    _ => abort("automatic two-step clustering should work")
  }
  assert_eq(result.centroids.length(), 3)
}

///|
test "bluster: unadjusted pairwise Rand decomposes correct pairs" {
  let result = @src.bluster_pairwise_rand(
    [0, 0, 1, 1],
    [0, 1, 0, 1],
    adjusted=false,
  ) catch {
    _ => abort("pairwise Rand should work")
  }
  assert_eq(result.cluster_names, [0, 1])
  assert_eq(bluster_test_some(result.correct[0][0]), 0.0)
  assert_eq(bluster_test_some(result.totals[0][0]), 1.0)
  assert_eq(bluster_test_some(result.correct[0][1]), 2.0)
  assert_eq(bluster_test_some(result.totals[0][1]), 4.0)
  assert_true(
    bluster_test_close(bluster_test_some(result.index), 1.0 / 3.0, 1.0e-12),
  )
  assert_true(result.ratios[1][0] is None)
}

///|
test "bluster: adjusted pairwise Rand matches the standard ARI" {
  let result = @src.bluster_pairwise_rand([0, 0, 1, 1], [0, 1, 0, 1]) catch {
    _ => abort("adjusted pairwise Rand should work")
  }
  assert_true(
    bluster_test_close(bluster_test_some(result.index), -0.5, 1.0e-12),
  )
  assert_true(
    bluster_test_close(bluster_test_some(result.ratios[0][0]), -0.5, 1.0e-12),
  )
  assert_true(
    bluster_test_close(bluster_test_some(result.ratios[0][1]), -0.5, 1.0e-12),
  )
}

///|
test "bluster: identical clusterings have Rand index one" {
  let result = @src.bluster_pairwise_rand([2, 2, 5, 5, 9, 9], [2, 2, 5, 5, 9, 9]) catch {
    _ => abort("identical pairwise Rand should work")
  }
  assert_true(bluster_test_close(bluster_test_some(result.index), 1.0, 1.0e-12))
  for right in 0..<result.cluster_names.length() {
    for left in 0..<(right + 1) {
      assert_true(
        bluster_test_close(
          bluster_test_some(result.ratios[left][right]),
          1.0,
          1.0e-12,
        ),
      )
    }
  }
}

///|
test "bluster: clustering comparison is symmetric with unit diagonal" {
  let compared = @src.bluster_compare_clusterings([
    [0, 0, 1, 1],
    [0, 1, 0, 1],
    [0, 0, 1, 1],
  ]) catch {
    _ => abort("clustering comparison should work")
  }
  assert_eq(bluster_test_some(compared[0][0]), 1.0)
  assert_eq(bluster_test_some(compared[2][2]), 1.0)
  assert_eq(compared[0][1], compared[1][0])
  assert_eq(compared[0][2], compared[2][0])
  assert_eq(bluster_test_some(compared[0][2]), 1.0)
}

///|
test "bluster: Rand and comparison validate equal lengths" {
  let mut failures = 0
  ignore(@src.bluster_pairwise_rand([0], [0, 1])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.bluster_compare_clusterings([[0], [0, 1]])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "bluster: approximate silhouette is one for pure duplicates" {
  let result = @src.bluster_approx_silhouette([[0.0], [0.0], [10.0], [10.0]], [
    0, 0, 1, 1,
  ]) catch {
    _ => abort("silhouette approximation should work")
  }
  assert_eq(result.other_clusters, [1, 1, 0, 0])
  for width in result.widths {
    assert_eq(width, 1.0)
  }
}

///|
test "bluster: approximate silhouette is zero for identical distributions" {
  let result = @src.bluster_approx_silhouette([[0.0], [10.0], [0.0], [10.0]], [
    0, 0, 1, 1,
  ]) catch {
    _ => abort("silhouette approximation should work")
  }
  for width in result.widths {
    assert_true(bluster_test_close(width, 0.0, 1.0e-12))
  }
}

///|
test "bluster: approximate silhouette handles one cluster" {
  let result = @src.bluster_approx_silhouette([[0.0], [1.0], [2.0]], [7, 7, 7]) catch {
    _ => abort("one-cluster silhouette should work")
  }
  assert_eq(result.other_clusters, [7, 7, 7])
  assert_eq(result.widths, [0.0, 0.0, 0.0])
}

///|
test "bluster: cluster RMSD uses sample variances" {
  let result = @src.bluster_cluster_rmsd([[0.0], [2.0], [10.0]], [0, 0, 1]) catch {
    _ => abort("cluster RMSD should work")
  }
  assert_eq(result.cluster_names, [0, 1])
  assert_true(
    bluster_test_close(bluster_test_some(result.values[0]), 2.0.sqrt(), 1.0e-12),
  )
  assert_true(result.values[1] is None)
}

///|
test "bluster: cluster RMSD sum mode returns within-cluster sums" {
  let result = @src.bluster_cluster_rmsd(
    [[0.0], [2.0], [10.0], [14.0]],
    [0, 0, 1, 1],
    sum_squares=true,
  ) catch {
    _ => abort("cluster RMSD sum mode should work")
  }
  assert_eq(bluster_test_some(result.values[0]), 2.0)
  assert_eq(bluster_test_some(result.values[1]), 8.0)
}

///|
test "bluster: neighbor purity is one for separated duplicates" {
  let result = @src.bluster_neighbor_purity(
    [[0.0], [0.0], [10.0], [10.0]],
    [0, 0, 1, 1],
    k=1,
  ) catch {
    _ => abort("neighbor purity should work")
  }
  assert_eq(result.radius, 0.0)
  assert_eq(result.purity, [1.0, 1.0, 1.0, 1.0])
  assert_eq(result.maximum_clusters, [0, 0, 1, 1])
}

///|
test "bluster: balanced purity removes cluster frequency effects" {
  let result = @src.bluster_neighbor_purity(
    [[0.0], [0.0], [0.0]],
    [0, 0, 1],
    k=1,
  ) catch {
    _ => abort("balanced neighbor purity should work")
  }
  assert_eq(result.purity, [0.5, 0.5, 0.5])
  assert_eq(result.maximum_clusters, [0, 0, 0])
}

///|
test "bluster: unbalanced purity retains cluster frequency effects" {
  let result = @src.bluster_neighbor_purity(
    [[0.0], [0.0], [0.0]],
    [0, 0, 1],
    k=1,
    balanced=false,
  ) catch {
    _ => abort("unbalanced neighbor purity should work")
  }
  assert_true(bluster_test_close(result.purity[0], 2.0 / 3.0, 1.0e-12))
  assert_true(bluster_test_close(result.purity[2], 1.0 / 3.0, 1.0e-12))
}

///|
test "bluster: custom purity weights override frequency balancing" {
  let result = @src.bluster_neighbor_purity(
    [[0.0], [0.0], [0.0]],
    [0, 0, 1],
    k=1,
    custom_weights=Some([1.0, 1.0, 2.0]),
  ) catch {
    _ => abort("custom-weight neighbor purity should work")
  }
  assert_eq(result.purity, [0.5, 0.5, 0.5])
}

///|
test "bluster: neighbor purity validates labels and custom weights" {
  let mut failures = 0
  ignore(@src.bluster_neighbor_purity([[0.0]], [], k=1)) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.bluster_neighbor_purity(
      [[0.0], [1.0]],
      [0, 1],
      custom_weights=Some([1.0]),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.bluster_neighbor_purity(
      [[0.0], [1.0]],
      [0, 1],
      custom_weights=Some([1.0, -1.0]),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "bluster: pairwise modularity decomposes observed and expected weights" {
  let graph = @src.bluster_make_snn_graph([[0.0], [1.0], [10.0], [11.0]], k=1) catch {
    _ => abort("SNN graph should build")
  }
  let result = @src.bluster_pairwise_modularity(graph, [0, 0, 1, 1]) catch {
    _ => abort("pairwise modularity should work")
  }
  assert_true(bluster_test_close(result.total_weight, 1.0, 1.0e-12))
  assert_true(
    bluster_test_close(bluster_test_some(result.observed[0][0]), 0.5, 1.0e-12),
  )
  assert_eq(bluster_test_some(result.observed[0][1]), 0.0)
  assert_true(
    bluster_test_close(bluster_test_some(result.expected[0][0]), 0.25, 1.0e-12),
  )
  assert_true(
    bluster_test_close(bluster_test_some(result.expected[0][1]), 0.5, 1.0e-12),
  )
  assert_true(
    bluster_test_close(
      bluster_test_some(result.modularity[0][0]),
      0.25,
      1.0e-12,
    ),
  )
  assert_eq(bluster_test_some(result.ratios[0][0]), 2.0)
  assert_eq(bluster_test_some(result.ratios[0][1]), 0.0)
}

///|
test "bluster: pairwise modularity upper triangle sums to zero" {
  let graph = @src.bluster_make_knn_graph([[0.0], [1.0], [3.0], [10.0]], k=1) catch {
    _ => abort("KNN graph should build")
  }
  let result = @src.bluster_pairwise_modularity(graph, [0, 0, 1, 1]) catch {
    _ => abort("pairwise modularity should work")
  }
  let mut total = 0.0
  for right in 0..<result.cluster_names.length() {
    for left in 0..<(right + 1) {
      total = total + bluster_test_some(result.modularity[left][right])
    }
  }
  assert_true(bluster_test_close(total, 0.0, 1.0e-12))
}

///|
test "bluster: community merging joins connected singleton pairs" {
  let graph = @src.bluster_make_snn_graph([[0.0], [1.0], [10.0], [11.0]], k=1) catch {
    _ => abort("SNN graph should build")
  }
  let merged = @src.bluster_merge_communities(graph, [10, 11, 20, 21], 2) catch {
    _ => abort("community merging should work")
  }
  assert_eq(bluster_test_unique_count(merged), 2)
  assert_true(merged[0] == merged[1])
  assert_true(merged[2] == merged[3])
  assert_true(merged[0] != merged[2])
}

///|
test "bluster: community merging preserves labels when no merge is needed" {
  let graph = @src.bluster_make_knn_graph([[0.0], [1.0], [10.0], [11.0]], k=1) catch {
    _ => abort("KNN graph should build")
  }
  let labels = [10, 10, 20, 20]
  let merged = @src.bluster_merge_communities(graph, labels, 5) catch {
    _ => abort("no-op community merging should work")
  }
  assert_eq(merged, labels)
}

///|
test "bluster: nested clusters identify perfect one-to-many mappings" {
  let result = @src.bluster_nested_clusters([0, 0, 1, 1], [10, 11, 20, 20]) catch {
    _ => abort("nested cluster mapping should work")
  }
  assert_eq(result.reference_names, [0, 1])
  assert_eq(result.alternative_names, [10, 11, 20])
  assert_eq(result.proportions[0], [1.0, 0.0])
  assert_eq(result.proportions[1], [1.0, 0.0])
  assert_eq(result.proportions[2], [0.0, 1.0])
  assert_eq(result.alternative_mapping, [0, 0, 1])
  assert_eq(result.alternative_maximum, [1.0, 1.0, 1.0])
  assert_eq(result.reference_scores, [1.0, 1.0])
}

///|
test "bluster: nested cluster mapping reports imperfect mixtures" {
  let result = @src.bluster_nested_clusters([0, 0, 1, 1], [10, 10, 10, 10]) catch {
    _ => abort("mixed nested cluster mapping should work")
  }
  assert_eq(result.proportions[0], [0.5, 0.5])
  assert_eq(result.alternative_mapping, [0])
  assert_eq(result.reference_scores, [0.5, 0.5])
}

///|
test "bluster: minimum link denominator favors nested clusters" {
  let result = @src.bluster_link_clusters_matrix(
    [0, 0, 1, 1],
    [0, 1, 1, 1],
    denominator=@src.bluster_link_minimum(),
  ) catch {
    _ => abort("minimum-denominator linking should work")
  }
  assert_eq(result, [[1.0, 0.5], [0.0, 1.0]])
}

///|
test "bluster: maximum link denominator requires similar sizes" {
  let result = @src.bluster_link_clusters_matrix(
    [0, 0, 1, 1],
    [0, 1, 1, 1],
    denominator=@src.bluster_link_maximum(),
  ) catch {
    _ => abort("maximum-denominator linking should work")
  }
  assert_true(bluster_test_close(result[0][0], 0.5, 1.0e-12))
  assert_true(bluster_test_close(result[0][1], 1.0 / 3.0, 1.0e-12))
  assert_true(bluster_test_close(result[1][1], 2.0 / 3.0, 1.0e-12))
}

///|
test "bluster: union link denominator computes Jaccard correspondence" {
  let result = @src.bluster_link_clusters_matrix([0, 0, 1, 1], [0, 1, 1, 1]) catch {
    _ => abort("union-denominator linking should work")
  }
  assert_true(bluster_test_close(result[0][0], 0.5, 1.0e-12))
  assert_true(bluster_test_close(result[0][1], 0.25, 1.0e-12))
  assert_true(bluster_test_close(result[1][1], 2.0 / 3.0, 1.0e-12))
}

///|
test "bluster: nested and linked mappings validate equal lengths" {
  let mut failures = 0
  ignore(@src.bluster_nested_clusters([0], [0, 1])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.bluster_link_clusters_matrix([0], [0, 1])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "bluster: bootstrap k-means stability is reproducible" {
  let config = @src.BlusterKmeansConfig::create(2, starts=4, seed=71) catch {
    _ => abort("bootstrap k-means configuration should build")
  }
  let data = [[0.0], [0.1], [0.2], [0.3], [10.0], [10.1], [10.2], [10.3]]
  let first = @src.bluster_bootstrap_kmeans_stability(
    data,
    config,
    iterations=8,
    seed=73,
  ) catch {
    _ => abort("first bootstrap stability run should work")
  }
  let second = @src.bluster_bootstrap_kmeans_stability(
    data,
    config,
    iterations=8,
    seed=73,
  ) catch {
    _ => abort("second bootstrap stability run should work")
  }
  assert_eq(first.cluster_names, [0, 1])
  assert_eq(first.ratios, second.ratios)
  assert_eq(first.adjusted_indices, second.adjusted_indices)
  assert_eq(first.adjusted_indices.length(), 8)
  assert_true(bluster_test_some(first.ratios[0][0]) > 0.5)
  assert_true(bluster_test_some(first.ratios[0][1]) > 0.5)
  assert_true(bluster_test_some(first.ratios[1][1]) > 0.5)
}

///|
test "bluster: bootstrap mean aggregation returns the full matrix shape" {
  let config = @src.BlusterKmeansConfig::create(2, starts=3, seed=79) catch {
    _ => abort("bootstrap k-means configuration should build")
  }
  let result = @src.bluster_bootstrap_kmeans_stability(
    [[0.0], [0.1], [0.2], [10.0], [10.1], [10.2]],
    config,
    iterations=5,
    use_mean=true,
    seed=83,
  ) catch {
    _ => abort("mean bootstrap stability should work")
  }
  assert_eq(result.ratios.length(), 2)
  assert_eq(result.ratios[0].length(), 2)
  assert_true(result.ratios[1][0] is None)
}

///|
test "bluster: bootstrap stability rejects non-positive iterations" {
  let config = @src.BlusterKmeansConfig::create(2) catch {
    _ => abort("valid k-means configuration should build")
  }
  let mut raised = false
  ignore(
    @src.bluster_bootstrap_kmeans_stability(
      [[0.0], [1.0]],
      config,
      iterations=0,
    ),
  ) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "bluster: SCE clustering writes labels without mutating the input" {
  let experiment = bluster_test_sce()
  let config = @src.BlusterGraphConfig::create(k=1, resolution=0.5, seed=89) catch {
    _ => abort("valid SCE graph configuration should build")
  }
  let result = @src.bluster_cluster_sce(
    experiment,
    config,
    output_column="community",
  ) catch {
    _ => abort("SCE graph clustering should work")
  }
  assert_false(experiment.col_data.contains("community"))
  assert_true(result.experiment.col_data.contains("community"))
  assert_eq(result.experiment.col_data["community"].length(), 4)
  assert_eq(result.experiment.metadata["bluster.reduced_dim"], "PCA")
  assert_eq(result.experiment.metadata["bluster.graph"], "SNN")
  assert_eq(result.clustering.clusters.length(), 4)
}

///|
test "bluster: SCE clustering deep-copies nested containers" {
  let experiment = bluster_test_sce()
  let config = @src.BlusterGraphConfig::create(k=1, seed=97) catch {
    _ => abort("valid SCE graph configuration should build")
  }
  let result = @src.bluster_cluster_sce(experiment, config) catch {
    _ => abort("SCE graph clustering should work")
  }
  result.experiment.assays["counts"][0][0] = 999.0
  result.experiment.row_data["symbol"][0] = "changed"
  result.experiment.reduced_dims["PCA"][0][0] = 999.0
  result.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 999.0
  assert_eq(experiment.assays["counts"][0][0], 1.0)
  assert_eq(experiment.row_data["symbol"][0], "A")
  assert_eq(experiment.reduced_dims["PCA"][0][0], 0.0)
  assert_eq(
    experiment.alternative_experiments["spike"].assays["counts"][0][0],
    5.0,
  )
}

///|
test "bluster: SCE clustering validates reduced dimensions" {
  let config = @src.BlusterGraphConfig::create(k=1) catch {
    _ => abort("valid SCE graph configuration should build")
  }
  let missing = @src.SingleCellExperiment::new([[1.0, 2.0]], ["G1"], [
    "C1", "C2",
  ])
  let mut failures = 0
  ignore(@src.bluster_cluster_sce(missing, config)) catch {
    _ => failures = failures + 1
  }
  missing.reduced_dims["PCA"] = [[0.0]]
  ignore(@src.bluster_cluster_sce(missing, config)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}
