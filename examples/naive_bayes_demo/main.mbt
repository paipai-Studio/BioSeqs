///|
/// NaiveBayes module example - k-mer based biological sequence classifier.
///
/// Run: moon run examples/naive_bayes_demo
fn main {
  println("=== Bio.NaiveBayes Module Demo ===")
  println("k-mer based Biological Sequence Classification")
  println("")

  // ---- Section 1: Train on sample data ----
  println("=== 1. Load and Inspect Training Data ===")
  let (train_seqs, train_labels) = @bio.naive_bayes_sample_data()
  println("Total sequences: " + train_seqs.length().to_string())
  let mut at_count = 0
  let mut gc_count = 0
  let mut i = 0
  while i < train_labels.length() {
    if train_labels[i] == "AT_rich" { at_count = at_count + 1 }
    else if train_labels[i] == "GC_rich" { gc_count = gc_count + 1 }
    i = i + 1
  }
  println("Class distribution: AT_rich=" + at_count.to_string() +
    ", GC_rich=" + gc_count.to_string())
  println("")

  // ---- Section 2: Train classifier ----
  println("=== 2. Train Naive Bayes Classifier (k=3, Laplace smoothing alpha=1.0) ===")
  let clf_base = @bio.NaiveBayesClassifier::new()
  let clf = @bio.naive_bayes_train(clf_base, train_seqs, train_labels)
  println("Vocabulary size (unique 3-mers): " + clf.vocabulary.length().to_string())
  println("Class labels: " + clf.class_labels.length().to_string())
  println("")

  // ---- Section 3: Predict new sequences ----
  println("=== 3. Predict Classes ===")
  let at_seq = "ATATATATATATATATATAT"
  let gc_seq = "GCGCGCGCGCGCGCGCGCGC"
  let mixed = "AGCTAGCTAGCTAGCTAGCT"

  let pred_at = @bio.naive_bayes_predict(clf, at_seq)
  let pred_gc = @bio.naive_bayes_predict(clf, gc_seq)
  let pred_mixed = @bio.naive_bayes_predict(clf, mixed)
  println("Sequence: " + at_seq)
  println("  Predicted class: " + pred_at.0 + "  (log-prob: " + pred_at.1.to_string() + ")")
  println("Sequence: " + gc_seq)
  println("  Predicted class: " + pred_gc.0 + "  (log-prob: " + pred_gc.1.to_string() + ")")
  println("Sequence: " + mixed)
  println("  Predicted class: " + pred_mixed.0 + "  (log-prob: " + pred_mixed.1.to_string() + ")")
  println("")

  // ---- Section 4: Probabilistic predictions (softmax-normalized) ----
  println("=== 4. Probabilistic Predictions (posterior probabilities) ===")
  let probs = @bio.naive_bayes_predict_proba(clf, mixed)
  println("Sequence: " + mixed)
  let mut pi = 0
  while pi < probs.length() {
    let (label, p) = probs[pi]
    println("  P(" + label + ") = " + p.to_string())
    pi = pi + 1
  }
  println("")

  // ---- Section 5: Top-K predictions ----
  println("=== 5. Top-K Predictions ===")
  let top2 = @bio.naive_bayes_top_k(clf, mixed, 2)
  let mut ti = 0
  while ti < top2.length() {
    let (label, score) = top2[ti]
    println("  Rank #" + (ti + 1).to_string() + ": " + label + " (log-prob: " + score.to_string() + ")")
    ti = ti + 1
  }
  println("")

  // ---- Section 6: Effect of k-mer size ----
  println("=== 6. Effect of k-mer Size on Prediction ===")
  let test_seq = "AAATTTAAATTTAAATTT"
  println("Test sequence: " + test_seq)
  for k in [2, 3, 4, 5] {
    let clf_k = @bio.NaiveBayesClassifier::new().set_kmer_size(k)
    let trained_k = @bio.naive_bayes_train(clf_k, train_seqs, train_labels)
    let pred_k = @bio.naive_bayes_predict(trained_k, test_seq)
    println("  k=" + k.to_string() + " -> class=" + pred_k.0 +
      ", log-prob=" + pred_k.1.to_string() +
      ", vocab=" + trained_k.vocabulary.length().to_string())
  }
  println("")

  // ---- Section 7: Accuracy on training set ----
  println("=== 7. Training Set Accuracy ===")
  let acc = @bio.naive_bayes_accuracy(clf, train_seqs, train_labels)
  println("Accuracy on training data: " + (acc * 100.0).to_string() + "%")
  println("")

  println("=== Demo complete ===")
}
