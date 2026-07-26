///|
/// BiocNeighbors Tests - Nearest neighbor search for single-cell analysis

///|
test "knn_options_new" {
  let opts = @src.knn_options_new()
  assert_eq(opts.k, 10)
  assert_eq(opts.distance, "euclidean")
  assert_eq(opts.n_trees, 50)
  assert_eq(opts.search_k, -1)
}

///|
test "knn_options_with_method" {
  let opts = @src.knn_options_with_method(@src.kmknn_method(), 5)
  assert_eq(opts.k, 5)
  assert_true(opts.method == @src.kmknn_method())
  assert_eq(opts.distance, "euclidean")
}

///|
test "knn_options_full" {
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    8,
    "manhattan",
    30,
    100,
  )
  assert_true(opts.method == @src.annoy_method())
  assert_eq(opts.k, 8)
  assert_eq(opts.distance, "manhattan")
  assert_eq(opts.n_trees, 30)
  assert_eq(opts.search_k, 100)
}

///|
test "knn_compute_distance_euclidean" {
  let p1 = [1.0, 2.0, 3.0]
  let p2 = [4.0, 6.0, 3.0]
  let dist = @src.knn_compute_distance(p1, p2, "euclidean")
  assert_eq(dist, 5.0)
}

///|
test "knn_compute_distance_manhattan" {
  let p1 = [1.0, 2.0, 3.0]
  let p2 = [4.0, 6.0, 3.0]
  let dist = @src.knn_compute_distance(p1, p2, "manhattan")
  assert_eq(dist, 7.0)
}

///|
test "knn_compute_distance_cosine" {
  let p1 = [1.0, 0.0]
  let p2 = [0.0, 1.0]
  let dist = @src.knn_compute_distance(p1, p2, "cosine")
  assert_eq(dist, 1.0)
}

///|
test "knn_brute_force_basic" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
  ]
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 4, 2, query, 1, 2, "euclidean")
  assert_eq(result.n_query, 1)
  assert_eq(result.k, 2)
  assert_eq(result.indices.length(), 1)
  assert_eq(result.indices[0][0], 0)
  assert_eq(result.indices[0].length(), 2)
}

///|
test "knn_brute_force_k_equals_n" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 3, 2, query, 1, 5, "euclidean")
  assert_eq(result.k, 3)
  assert_eq(result.indices[0].length(), 3)
}

///|
test "knn_brute_force_empty_data" {
  let data : Array[Array[Double]] = []
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 0, 2, query, 1, 2, "euclidean")
  assert_eq(result.k, 0)
  assert_eq(result.indices.length(), 0)
}

///|
test "knn_brute_force_zero_k" {
  let data = [[0.0, 0.0], [1.0, 0.0]]
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 2, 2, query, 1, 0, "euclidean")
  assert_eq(result.k, 0)
  assert_eq(result.indices.length(), 0)
}

///|
test "knn_brute_force_multiple_queries" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let query = [[0.0, 0.0], [1.0, 1.0]]
  let result = @src.knn_brute_force(data, 4, 2, query, 2, 2, "euclidean")
  assert_eq(result.n_query, 2)
  assert_eq(result.indices.length(), 2)
}

///|
test "knn_brute_force_manhattan" {
  let data = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]]
  let query = [[1.0, 0.0]]
  let result = @src.knn_brute_force(data, 3, 2, query, 1, 3, "manhattan")
  assert_eq(result.k, 3)
  assert_eq(result.distances[0][0], 1.0)
}

///|
test "knn_find_k_smallest" {
  let dists = [5.0, 2.0, 8.0, 1.0, 3.0]
  let indices = @src.knn_find_k_smallest(dists, 3)
  assert_eq(indices.length(), 3)
  assert_eq(indices[0], 3)
  assert_eq(indices[1], 1)
  assert_eq(indices[2], 4)
}

///|
test "knn_find_k_smallest_k_equals_n" {
  let dists = [3.0, 1.0, 2.0]
  let indices = @src.knn_find_k_smallest(dists, 5)
  assert_eq(indices.length(), 3)
}

///|
test "build_knn_index_brute_force" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.brute_force_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_knn_index(data, 3, 2, opts)
  assert_eq(index.n_points, 3)
  assert_eq(index.n_dim, 2)
}

///|
test "build_knn_index_kmknn" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_knn_index(data, 8, 2, opts)
  assert_true(index.method == @src.kmknn_method())
  assert_true(index.centroids.length() > 0)
}

///|
test "build_knn_index_annoy" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    5,
    -1,
  )
  let index = @src.build_knn_index(data, 8, 2, opts)
  assert_true(index.method == @src.annoy_method())
  assert_true(index.tree_nodes.length() > 0)
}

///|
test "build_kmknn_index_empty" {
  let data : Array[Array[Double]] = []
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 0, 2, opts)
  assert_eq(index.n_points, 0)
  assert_eq(index.centroids.length(), 0)
}

///|
test "build_annoy_index_empty" {
  let data : Array[Array[Double]] = []
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_annoy_index(data, 0, 2, opts)
  assert_eq(index.n_points, 0)
  assert_eq(index.tree_nodes.length(), 0)
}

///|
test "run_knn_brute_force" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.brute_force_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_knn_index(data, 4, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.run_knn(index, query, 1, opts)
  assert_eq(result.n_query, 1)
  assert_eq(result.indices[0][0], 0)
}

///|
test "run_knn_kmknn" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    3,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_knn_index(data, 8, 2, opts)
  let query = [[0.0, 0.0], [10.0, 10.0]]
  let result = @src.run_knn(index, query, 2, opts)
  assert_eq(result.n_query, 2)
  assert_eq(result.indices.length(), 2)
}

///|
test "run_knn_annoy" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    3,
    "euclidean",
    5,
    -1,
  )
  let index = @src.build_knn_index(data, 8, 2, opts)
  let query = [[0.0, 0.0], [10.0, 10.0]]
  let result = @src.run_knn(index, query, 2, opts)
  assert_eq(result.n_query, 2)
  assert_eq(result.indices.length(), 2)
}

///|
test "knn_annoy_cosine_distance" {
  let data = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [-1.0, 0.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "cosine",
    3,
    -1,
  )
  let index = @src.build_annoy_index(data, 4, 2, opts)
  let query = [[1.0, 0.0]]
  let result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_kmknn_cosine_distance" {
  let data = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [-1.0, 0.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "cosine",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 4, 2, opts)
  let query = [[1.0, 0.0]]
  let result = @src.knn_kmknn(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_brute_force_cosine" {
  let data = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let query = [[1.0, 0.0]]
  let result = @src.knn_brute_force(data, 3, 2, query, 1, 3, "cosine")
  assert_eq(result.k, 3)
  assert_eq(result.distances[0][0], 0.0)
}

///|
test "neighbor_search_method_eq" {
  let m1 = @src.kmknn_method()
  let m2 = @src.annoy_method()
  let m3 = @src.brute_force_method()
  assert_true(m1 == @src.kmknn_method())
  assert_true(m2 == @src.annoy_method())
  assert_true(m3 == @src.brute_force_method())
}

///|
test "knn_options_eq" {
  let opts1 = @src.knn_options_new()
  let opts2 = @src.knn_options_new()
  assert_eq(opts1, opts2)
}

///|
test "knn_brute_force_single_dim" {
  let data = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let query = [[1.0], [5.0]]
  let result = @src.knn_brute_force(data, 5, 1, query, 2, 2, "euclidean")
  assert_eq(result.n_query, 2)
  assert_eq(result.indices[0][0], 0)
  assert_eq(result.indices[1][0], 4)
}

///|
test "knn_kmknn_single_point" {
  let data = [[1.0, 2.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    1,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 1, 2, opts)
  let query = [[1.0, 2.0]]
  let result = @src.knn_kmknn(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_annoy_single_point" {
  let data = [[1.0, 2.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    1,
    "euclidean",
    3,
    -1,
  )
  let index = @src.build_annoy_index(data, 1, 2, opts)
  let query = [[1.0, 2.0]]
  let result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_brute_force_vs_kmknn" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 2.0],
    [2.0, 3.0],
    [3.0, 3.0],
  ]
  let query = [[0.0, 0.0]]
  let brute_result = @src.knn_brute_force(data, 8, 2, query, 1, 3, "euclidean")
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    3,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 8, 2, opts)
  let kmknn_result = @src.knn_kmknn(index, query, 1, opts)
  assert_eq(brute_result.k, kmknn_result.k)
}

///|
test "knn_brute_force_vs_annoy" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 2.0],
    [2.0, 3.0],
    [3.0, 3.0],
  ]
  let query = [[0.0, 0.0]]
  let brute_result = @src.knn_brute_force(data, 8, 2, query, 1, 3, "euclidean")
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    3,
    "euclidean",
    5,
    -1,
  )
  let index = @src.build_annoy_index(data, 8, 2, opts)
  let annoy_result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(brute_result.k, 3)
  assert_true(annoy_result.k >= 1)
}

///|
test "knn_brute_force_no_query" {
  let data = [[0.0, 0.0], [1.0, 0.0]]
  let query : Array[Array[Double]] = []
  let result = @src.knn_brute_force(data, 2, 2, query, 0, 2, "euclidean")
  assert_eq(result.n_query, 0)
  assert_eq(result.indices.length(), 0)
}

///|
test "build_kmknn_index_small" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 3, 2, opts)
  assert_eq(index.n_points, 3)
  assert_true(index.centroids.length() >= 1)
}

///|
test "build_annoy_index_small" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    3,
    -1,
  )
  let index = @src.build_annoy_index(data, 3, 2, opts)
  assert_eq(index.n_points, 3)
  assert_true(index.tree_nodes.length() > 0)
}

///|
test "knn_kmknn_empty_query" {
  let data = [[0.0, 0.0], [1.0, 0.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 2, 2, opts)
  let query : Array[Array[Double]] = []
  let result = @src.knn_kmknn(index, query, 0, opts)
  assert_eq(result.n_query, 0)
}

///|
test "knn_annoy_empty_query" {
  let data = [[0.0, 0.0], [1.0, 0.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    3,
    -1,
  )
  let index = @src.build_annoy_index(data, 2, 2, opts)
  let query : Array[Array[Double]] = []
  let result = @src.knn_annoy(index, query, 0, opts)
  assert_eq(result.n_query, 0)
}

///|
test "knn_annoy_search_k_auto" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    3,
    -1,
  )
  let index = @src.build_annoy_index(data, 4, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_annoy_search_k_custom" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    3,
    1,
  )
  let index = @src.build_annoy_index(data, 4, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_brute_force_3d" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
    [2.0, 3.0, 4.0],
  ]
  let query = [[1.5, 2.5, 3.5]]
  let result = @src.knn_brute_force(data, 4, 3, query, 1, 2, "euclidean")
  assert_eq(result.indices[0][0], 0)
  assert_eq(result.indices[0][1], 3)
}

///|
test "knn_options_full_annoy_defaults" {
  let opts = @src.knn_options_with_method(@src.annoy_method(), 5)
  assert_eq(opts.n_trees, 50)
  assert_eq(opts.search_k, -1)
}

///|
test "knn_quick_sort_sorted" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let indices = [0, 1, 2, 3, 4]
  let (sorted_idx, _) = @src.knn_quick_sort(data, indices, 0, 4)
  assert_eq(sorted_idx[0], 0)
  assert_eq(sorted_idx[4], 4)
}

///|
test "knn_quick_sort_reverse" {
  let data = [5.0, 4.0, 3.0, 2.0, 1.0]
  let indices = [0, 1, 2, 3, 4]
  let (sorted_idx, _) = @src.knn_quick_sort(data, indices, 0, 4)
  assert_eq(sorted_idx[0], 4)
  assert_eq(sorted_idx[4], 0)
}

///|
test "knn_quick_sort_single" {
  let data = [3.0]
  let indices = [0]
  let (sorted_idx, _) = @src.knn_quick_sort(data, indices, 0, 0)
  assert_eq(sorted_idx.length(), 1)
}

///|
test "knn_brute_force_large_k" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 4, 2, query, 1, 10, "euclidean")
  assert_eq(result.k, 4)
  assert_eq(result.indices[0].length(), 4)
}

///|
test "knn_brute_force_different_k_values" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 2.0],
  ]
  let query = [[0.0, 0.0]]
  let result_k1 = @src.knn_brute_force(data, 5, 2, query, 1, 1, "euclidean")
  assert_eq(result_k1.k, 1)
  assert_eq(result_k1.indices[0].length(), 1)
  assert_eq(result_k1.indices[0][0], 0)
  let result_k3 = @src.knn_brute_force(data, 5, 2, query, 1, 3, "euclidean")
  assert_eq(result_k3.k, 3)
  assert_eq(result_k3.indices[0].length(), 3)
}

///|
test "knn_brute_force_manhattan_distance" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let query = [[0.5, 0.5]]
  let result = @src.knn_brute_force(data, 4, 2, query, 1, 2, "manhattan")
  assert_eq(result.k, 2)
  assert_true(result.distances[0][0] <= result.distances[0][1])
}

///|
test "knn_brute_force_cosine_normalized" {
  let data = [[1.0, 0.0], [0.0, 1.0], [0.7071, 0.7071], [-1.0, 0.0]]
  let query = [[1.0, 0.0]]
  let result = @src.knn_brute_force(data, 4, 2, query, 1, 4, "cosine")
  assert_eq(result.k, 4)
  assert_true((result.distances[0][0] - 0.0).abs() < 1.0e-9)
  assert_eq(result.indices[0][0], 0)
}

///|
test "knn_brute_force_vs_kmknn_same_results" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 2.0],
    [2.0, 3.0],
    [3.0, 3.0],
  ]
  let query = [[0.0, 0.0], [3.0, 3.0]]
  let brute_result = @src.knn_brute_force(data, 8, 2, query, 2, 3, "euclidean")
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    3,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 8, 2, opts)
  let kmknn_result = @src.knn_kmknn(index, query, 2, opts)
  assert_eq(brute_result.n_query, kmknn_result.n_query)
  assert_eq(brute_result.k, kmknn_result.k)
}

///|
test "knn_annoy_approximately_correct" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [10.0, 11.0],
    [11.0, 11.0],
  ]
  let query = [[0.0, 0.0]]
  let brute_result = @src.knn_brute_force(data, 8, 2, query, 1, 3, "euclidean")
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    3,
    "euclidean",
    10,
    -1,
  )
  let index = @src.build_annoy_index(data, 8, 2, opts)
  let annoy_result = @src.knn_annoy(index, query, 1, opts)
  assert_true(annoy_result.k >= 1)
  assert_true(annoy_result.k <= 3)
}

///|
test "knn_brute_force_cosine_3d" {
  let data = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 1.0, 0.0],
  ]
  let query = [[1.0, 0.0, 0.0]]
  let result = @src.knn_brute_force(data, 4, 3, query, 1, 2, "cosine")
  assert_eq(result.k, 2)
  assert_eq(result.indices[0][0], 0)
  assert_true((result.distances[0][0] - 0.0).abs() < 1.0e-9)
}

///|
test "knn_kmknn_manhattan_distance" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "manhattan",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 4, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.knn_kmknn(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_brute_force_single_point_data" {
  let data = [[1.0, 2.0]]
  let query = [[1.0, 2.0]]
  let result = @src.knn_brute_force(data, 1, 2, query, 1, 1, "euclidean")
  assert_eq(result.k, 1)
  assert_eq(result.indices[0][0], 0)
  assert_true((result.distances[0][0] - 0.0).abs() < 1.0e-9)
}

///|
test "knn_brute_force_k_exceeds_points" {
  let data = [[0.0, 0.0], [1.0, 0.0]]
  let query = [[0.0, 0.0]]
  let result = @src.knn_brute_force(data, 2, 2, query, 1, 100, "euclidean")
  assert_eq(result.k, 2)
  assert_eq(result.indices[0].length(), 2)
}

///|
test "knn_brute_force_identical_points" {
  let data = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
  let query = [[1.0, 1.0]]
  let result = @src.knn_brute_force(data, 3, 2, query, 1, 3, "euclidean")
  assert_eq(result.k, 3)
  for i = 0; i < 3; i = i + 1 {
    assert_true((result.distances[0][i] - 0.0).abs() < 1.0e-9)
  }
}

///|
test "knn_quick_sort_duplicates" {
  let data = [3.0, 1.0, 3.0, 2.0]
  let indices = [0, 1, 2, 3]
  let (sorted_idx, sorted_dist) = @src.knn_quick_sort(data, indices, 0, 3)
  assert_eq(sorted_idx.length(), 4)
  assert_true(sorted_dist[0] <= sorted_dist[1])
  assert_true(sorted_dist[1] <= sorted_dist[2])
  assert_true(sorted_dist[2] <= sorted_dist[3])
}

///|
test "knn_kmknn_empty_data" {
  let data : Array[Array[Double]] = []
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    2,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_kmknn_index(data, 0, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.knn_kmknn(index, query, 1, opts)
  assert_eq(result.n_query, 1)
  assert_eq(result.k, 0)
}

///|
test "knn_annoy_manhattan_distance" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "manhattan",
    5,
    -1,
  )
  let index = @src.build_annoy_index(data, 4, 2, opts)
  let query = [[0.0, 0.0]]
  let result = @src.knn_annoy(index, query, 1, opts)
  assert_eq(result.n_query, 1)
}

///|
test "knn_compute_distance_zero" {
  let p1 = [1.0, 2.0, 3.0]
  let dist = @src.knn_compute_distance(p1, p1, "euclidean")
  assert_true(dist.abs() < 1.0e-15)
  let dist_m = @src.knn_compute_distance(p1, p1, "manhattan")
  assert_true(dist_m.abs() < 1.0e-15)
  let dist_c = @src.knn_compute_distance(p1, p1, "cosine")
  assert_true(dist_c.abs() < 1.0e-15)
}

///|
test "knn_options_roundtrip" {
  let opts = @src.knn_options_full(
    @src.brute_force_method(),
    7,
    "cosine",
    25,
    50,
  )
  assert_true(opts.method == @src.brute_force_method())
  assert_eq(opts.k, 7)
  assert_eq(opts.distance, "cosine")
  assert_eq(opts.n_trees, 25)
  assert_eq(opts.search_k, 50)
}

///|
test "knn_find_k_smallest_duplicates" {
  let dists = [2.0, 1.0, 1.0, 3.0]
  let indices = @src.knn_find_k_smallest(dists, 2)
  assert_eq(indices.length(), 2)
}

///|
test "knn_brute_force_5d" {
  let data = [
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 1.0],
  ]
  let query = [[1.0, 0.0, 0.0, 0.0, 0.0]]
  let result = @src.knn_brute_force(data, 5, 5, query, 1, 2, "euclidean")
  assert_eq(result.indices[0][0], 0)
}

///|
test "build_knn_index_kmknn_small_data" {
  let data = [[1.0, 2.0]]
  let opts = @src.knn_options_full(
    @src.kmknn_method(),
    1,
    "euclidean",
    50,
    -1,
  )
  let index = @src.build_knn_index(data, 1, 2, opts)
  assert_true(index.method == @src.kmknn_method())
  assert_eq(index.n_points, 1)
}

///|
test "run_knn_annoy_multiple_queries" {
  let data = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
  let opts = @src.knn_options_full(
    @src.annoy_method(),
    2,
    "euclidean",
    3,
    -1,
  )
  let index = @src.build_knn_index(data, 4, 2, opts)
  let query = [[0.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
  let result = @src.run_knn(index, query, 3, opts)
  assert_eq(result.n_query, 3)
}