///|
/// K-means Clustering Tests
test "kmeans_create" {
  let kmeans = @src.KMeans::new(3, 100, 0.000001)
  assert_eq(kmeans.k, 3)
  assert_eq(kmeans.max_iter, 100)
}

///|
test "euclidean_distance" {
  let a = [1.0, 2.0]
  let b = [4.0, 6.0]
  let dist = @src.euclidean_distance(a, b)
  assert_eq(dist, 5.0)
}

///|
test "manhattan_distance" {
  let a = [1.0, 2.0]
  let b = [4.0, 6.0]
  let dist = @src.manhattan_distance(a, b)
  assert_eq(dist, 7.0)
}

///|
test "cosine_distance" {
  let a = [1.0, 0.0]
  let b = [0.0, 1.0]
  let dist = @src.cosine_distance(a, b)
  assert_eq(dist, 1.0)
}

///|
test "kmeans_fit" {
  let data = [
    [1.0, 2.0],
    [2.0, 1.0],
    [2.0, 3.0],
    [8.0, 7.0],
    [9.0, 8.0],
    [7.0, 9.0],
    [4.0, 5.0],
    [5.0, 4.0],
    [5.0, 6.0],
  ]
  let kmeans = @src.KMeans::new(3, 100, 0.000001).fit(data)
  assert_eq(kmeans.labels.length(), 9)
  assert_eq(kmeans.centroids.length(), 3)
}

///|
test "kmeans_predict" {
  let data = [[1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]]
  let kmeans = @src.KMeans::new(2, 100, 0.000001).fit(data)
  let new_data = [[1.5, 1.5], [8.5, 8.5]]
  let predictions = kmeans.predict(new_data)
  assert_eq(predictions.length(), 2)
}

///|
test "kmeans_predict_single" {
  let data = [[1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]]
  let kmeans = @src.KMeans::new(2, 100, 0.000001).fit(data)
  let prediction = kmeans.predict_single([1.5, 1.5])
  assert_true(prediction >= 0 && prediction < 2)
}

///|
test "kmeans_get_centroids" {
  let data = [[1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]]
  let kmeans = @src.KMeans::new(2, 100, 0.000001).fit(data)
  let centroids = kmeans.get_centroids()
  assert_eq(centroids.length(), 2)
}

///|
test "kmeans_get_labels" {
  let data = [[1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]]
  let kmeans = @src.KMeans::new(2, 100, 0.000001).fit(data)
  let labels = kmeans.get_labels()
  assert_eq(labels.length(), 4)
}

///|
test "kmeans_inertia" {
  let data = [[1.0, 2.0], [2.0, 1.0], [8.0, 7.0], [9.0, 8.0]]
  let kmeans = @src.KMeans::new(2, 100, 0.000001).fit(data)
  let inertia = kmeans.inertia(data)
  assert_true(inertia >= 0.0)
}

///|
test "kmeans_silhouette_score" {
  let data = [
    [1.0, 2.0],
    [2.0, 1.0],
    [2.0, 3.0],
    [8.0, 7.0],
    [9.0, 8.0],
    [7.0, 9.0],
  ]
  let labels = [0, 0, 0, 1, 1, 1]
  let score = @src.silhouette_score(data, labels)
  assert_true(score > 0.0)
}

///|
test "kmeans_cluster_gene_expression" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [8.0, 7.0, 6.0],
    [9.0, 8.0, 7.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let kmeans = @src.cluster_gene_expression(data, 3)
  assert_eq(kmeans.labels.length(), 6)
}

///|
test "kmeans_find_optimal_k" {
  let data = [
    [1.0, 2.0],
    [2.0, 1.0],
    [8.0, 7.0],
    [9.0, 8.0],
    [4.0, 5.0],
    [5.0, 4.0],
  ]
  let k = @src.find_optimal_k(data, 2, 4)
  assert_true(k >= 2 && k <= 4)
}

///|
test "kmeans_empty_data" {
  let data = @src.KMeans::new(3, 100, 0.000001).centroids
  let kmeans = @src.KMeans::new(3, 100, 0.000001).fit(data)
  assert_eq(kmeans.labels.length(), 0)
}

///|
test "kmeans_single_point" {
  let data = [[1.0, 2.0]]
  let kmeans = @src.KMeans::new(1, 100, 0.000001).fit(data)
  assert_eq(kmeans.labels.length(), 1)
}
