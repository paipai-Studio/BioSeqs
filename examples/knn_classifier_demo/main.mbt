///|
fn main {
  println("=== Bio.kNN Demo ===")

  let train = @src.knn_sample_train()
  let query = @src.knn_sample_query()

  println("\n1. Training set (operon-prediction tutorial data)")
  for sample in train {
    println("  " + sample.label + " <- " + sample.features.to_string())
  }

  println("\n2. Classify query " + query.to_string())
  let result = @src.knn_classify_full(train, 3, query)
  println("  Predicted label: " + result.label)
  println("  Vote table:")
  for vote in result.votes {
    println("    " + vote.label + ": " + vote.weight.to_string())
  }

  println("\n3. Leave-one-out cross validation")
  let accuracy = @src.knn_loo_accuracy(train, 3)
  println("  LOO accuracy (k=3): " + accuracy.to_string())

  println("\n4. Confusion matrix (k=1)")
  let (labels, matrix) = @src.knn_confusion_matrix(train, 1)
  println("  Labels: " + labels.to_string())
  for i in 0..<matrix.length() {
    println("  " + labels[i] + " -> " + matrix[i].to_string())
  }

  println("\n5. Feature normalization")
  let normalized = @src.knn_normalize_samples(train, mode="minmax")
  println("  Min-max normalized first sample: " + normalized[0].features.to_string())
}