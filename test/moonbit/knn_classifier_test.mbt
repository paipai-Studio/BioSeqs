///|
/// Tests for the Biopython Bio.kNN-compatible classifier.

///|
test "knn manhattan distance matches reference" {
  assert_eq(@src.knn_distance([1.0, 2.0, 3.0], [4.0, 6.0, 3.0]), 7.0)
  assert_eq(@src.knn_distance([0.0, 0.0], [0.0, 0.0]), 0.0)
}

///|
test "knn euclidean distance" {
  let d = @src.knn_euclidean_distance([3.0, 0.0], [0.0, 4.0])
  assert_true((d - 5.0).abs() < 1.0e-9)
}

///|
test "knn classify exact match wins" {
  let train = [
    @src.KnnSample::new(features=[1.0, 1.0], label="A"),
    @src.KnnSample::new(features=[5.0, 5.0], label="B"),
  ]
  assert_eq(@src.knn_classify(train, 1, [1.0, 1.0]), "A")
  assert_eq(@src.knn_classify(train, 1, [5.0, 5.0]), "B")
}

///|
test "knn classify weighted vote" {
  // Query near A (dist 1) and B (dist 2): 1/1 > 1/2 so A wins.
  let train = [
    @src.KnnSample::new(features=[0.0], label="A"),
    @src.KnnSample::new(features=[3.0], label="B"),
  ]
  assert_eq(@src.knn_classify(train, 2, [1.0]), "A")
}

///|
test "knn classify full returns vote table" {
  let train = [
    @src.KnnSample::new(features=[0.0], label="A"),
    @src.KnnSample::new(features=[3.0], label="B"),
  ]
  let result = @src.knn_classify_full(train, 2, [1.0])
  assert_eq(result.label, "A")
  assert_true(result.votes.length() == 2)
  assert_true((result.vote_of("A") - 1.0).abs() < 1.0e-9)
  assert_true((result.vote_of("B") - 0.5).abs() < 1.0e-9)
}

///|
test "knn classify probs sum to one" {
  let train = [
    @src.KnnSample::new(features=[0.0], label="A"),
    @src.KnnSample::new(features=[3.0], label="B"),
    @src.KnnSample::new(features=[6.0], label="A"),
  ]
  let probs = @src.knn_classify_probs(train, 3, [2.0])
  let mut total = 0.0
  for vote in probs {
    total += vote.weight
  }
  assert_true((total - 1.0).abs() < 1.0e-9)
}

///|
test "knn leave-one-out accuracy" {
  let train = [
    @src.KnnSample::new(features=[0.0], label="A"),
    @src.KnnSample::new(features=[1.0], label="A"),
    @src.KnnSample::new(features=[10.0], label="B"),
    @src.KnnSample::new(features=[11.0], label="B"),
  ]
  assert_true(@src.knn_loo_accuracy(train, 1) == 1.0)
}

///|
test "knn confusion matrix" {
  let train = [
    @src.KnnSample::new(features=[0.0], label="A"),
    @src.KnnSample::new(features=[1.0], label="A"),
    @src.KnnSample::new(features=[10.0], label="B"),
    @src.KnnSample::new(features=[11.0], label="B"),
  ]
  let (labels, matrix) = @src.knn_confusion_matrix(train, 1)
  assert_eq(labels, ["A", "B"])
  assert_eq(matrix[0][0], 2)
  assert_eq(matrix[1][1], 2)
  assert_eq(matrix[0][1], 0)
  assert_eq(matrix[1][0], 0)
}

///|
test "knn minmax normalize column" {
  let normalized = @src.knn_minmax_normalize_column([0.0, 5.0, 10.0])
  assert_true((normalized[0] - 0.0).abs() < 1.0e-9)
  assert_true((normalized[1] - 0.5).abs() < 1.0e-9)
  assert_true((normalized[2] - 1.0).abs() < 1.0e-9)
}

///|
test "knn zscore normalize column" {
  let normalized = @src.knn_zscore_normalize_column([1.0, 2.0, 3.0])
  let mut sum = 0.0
  for v in normalized {
    sum += v
  }
  assert_true(sum.abs() < 1.0e-9)
}

///|
test "knn normalize samples preserves labels" {
  let train = [
    @src.KnnSample::new(features=[0.0, 10.0], label="A"),
    @src.KnnSample::new(features=[10.0, 0.0], label="B"),
  ]
  let normalized = @src.knn_normalize_samples(train, mode="minmax")
  assert_eq(normalized.length(), 2)
  assert_eq(normalized[0].label, "A")
  assert_eq(normalized[1].label, "B")
  assert_true((normalized[0].features[0] - 0.0).abs() < 1.0e-9)
  assert_true((normalized[0].features[1] - 1.0).abs() < 1.0e-9)
}

///|
test "knn sample operon prediction" {
  let train = @src.knn_sample_train()
  let query = @src.knn_sample_query()
  let label = @src.knn_classify(train, 3, query)
  assert_true(label == "SameOperon" || label == "DifferentOperon")
  // Query [3, -1, 0.12] is closest to the DifferentOperon samples at
  // features [3, -1, 0.1] and [3, -1, 0.15]: predict DifferentOperon.
  assert_eq(label, "DifferentOperon")
}