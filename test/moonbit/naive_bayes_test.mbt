///|
/// Tests for Bio.NaiveBayes sequence classifier module.

test "classifier_new_defaults" {
  let clf = @src.NaiveBayesClassifier::new()
  assert_eq(clf.kmer_size, 3)
  assert_eq(clf.alpha, 1.0)
  assert_eq(clf.vocabulary.length(), 0)
  assert_eq(clf.class_labels.length(), 0)
}

test "classifier_set_kmer_size" {
  let clf = @src.NaiveBayesClassifier::new().set_kmer_size(5)
  assert_eq(clf.kmer_size, 5)
  assert_eq(clf.alpha, 1.0)
}

test "classifier_set_alpha" {
  let clf = @src.NaiveBayesClassifier::new().set_alpha(0.5)
  assert_eq(clf.kmer_size, 3)
  assert_eq(clf.alpha, 0.5)
}

test "classifier_setters_chain" {
  let clf = @src.NaiveBayesClassifier::new().set_kmer_size(4).set_alpha(0.1)
  assert_eq(clf.kmer_size, 4)
  assert_eq(clf.alpha, 0.1)
}

test "extract_kmers_basic" {
  let kmers = @src.naive_bayes_extract_kmers("ABCDE", 2)
  assert_eq(kmers.length(), 4)
  assert_eq(kmers.get("AB").unwrap(), 1)
  assert_eq(kmers.get("BC").unwrap(), 1)
  assert_eq(kmers.get("CD").unwrap(), 1)
  assert_eq(kmers.get("DE").unwrap(), 1)
}

test "extract_kmers_repeated" {
  let kmers = @src.naive_bayes_extract_kmers("AAAA", 2)
  assert_eq(kmers.length(), 1)
  assert_eq(kmers.get("AA").unwrap(), 3)
}

test "extract_kmers_k_larger_than_length" {
  let kmers = @src.naive_bayes_extract_kmers("ABC", 5)
  assert_eq(kmers.length(), 0)
}

test "extract_kmers_k_equals_length" {
  let kmers = @src.naive_bayes_extract_kmers("ABC", 3)
  assert_eq(kmers.length(), 1)
  assert_eq(kmers.get("ABC").unwrap(), 1)
}

test "extract_kmers_empty_sequence" {
  let kmers = @src.naive_bayes_extract_kmers("", 3)
  assert_eq(kmers.length(), 0)
}

test "train_basic" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  assert_eq(trained.class_labels.length(), 2)
  assert_true(trained.vocabulary.length() > 0)
  assert_eq(trained.models.length(), 2)
}

test "predict_at_rich_correct" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let pred = @src.naive_bayes_predict(trained, "ATATATATATATATATATAT")
  assert_eq(pred.0, "AT_rich")
}

test "predict_gc_rich_correct" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let pred = @src.naive_bayes_predict(trained, "GCGCGCGCGCGCGCGCGCGC")
  assert_eq(pred.0, "GC_rich")
}

test "predict_proba_sums_to_one" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let probs = @src.naive_bayes_predict_proba(trained, "ATATATATATATATATATAT")
  let mut sum = 0.0
  for p in probs {
    sum = sum + p.1
  }
  assert_true(sum > 0.99 && sum < 1.01)
}

test "predict_log_probs_returns_all_classes" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let log_probs = @src.naive_bayes_predict_log_probs(trained, "ATATATATATATATATATAT")
  assert_eq(log_probs.length(), 2)
}

test "top_k_returns_k_items" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let top = @src.naive_bayes_top_k(trained, "ATATATATATATATATATAT", 1)
  assert_eq(top.length(), 1)
}

test "top_k_descending_order" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let top = @src.naive_bayes_top_k(trained, "ATATATATATATATATATAT", 2)
  assert_eq(top.length(), 2)
  assert_true(top[0].1 >= top[1].1)
}

test "accuracy_on_sample_data" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let acc = @src.naive_bayes_accuracy(trained, data.0, data.1)
  assert_true(acc > 0.5)
}

test "unknown_sequence_reasonable_defaults" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let probs = @src.naive_bayes_predict_proba(trained, "NNNNNNNNNNNNNNNNNNNN")
  assert_eq(probs.length(), 2)
  let mut sum = 0.0
  for p in probs {
    sum = sum + p.1
    assert_true(p.1 >= 0.0 && p.1 <= 1.0)
  }
  assert_true(sum > 0.99 && sum < 1.01)
}

test "kmer_size_affects_results" {
  let data = @src.naive_bayes_sample_data()
  let clf3 = @src.NaiveBayesClassifier::new()
  let trained3 = @src.naive_bayes_train(clf3, data.0, data.1)
  let clf2 = @src.NaiveBayesClassifier::new().set_kmer_size(2)
  let trained2 = @src.naive_bayes_train(clf2, data.0, data.1)
  let log_probs3 = @src.naive_bayes_predict_log_probs(trained3, "ATATATATATATATATATAT")
  let log_probs2 = @src.naive_bayes_predict_log_probs(trained2, "ATATATATATATATATATAT")
  let diff = (log_probs3[0].1 - log_probs2[0].1).abs()
  assert_true(diff > 0.0 || trained3.vocabulary.length() != trained2.vocabulary.length())
}

test "class_labels_array_after_training" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  assert_eq(trained.class_labels.length(), 2)
  let has_at = trained.class_labels.contains("AT_rich")
  let has_gc = trained.class_labels.contains("GC_rich")
  assert_true(has_at && has_gc)
}

test "sample_data_has_8_sequences" {
  let data = @src.naive_bayes_sample_data()
  assert_eq(data.0.length(), 8)
  assert_eq(data.1.length(), 8)
}

test "sample_data_labels_count" {
  let data = @src.naive_bayes_sample_data()
  let mut at_count = 0
  let mut gc_count = 0
  for label in data.1 {
    if label == "AT_rich" {
      at_count = at_count + 1
    } else if label == "GC_rich" {
      gc_count = gc_count + 1
    }
  }
  assert_eq(at_count, 4)
  assert_eq(gc_count, 4)
}

test "train_empty_sequences" {
  let clf = @src.NaiveBayesClassifier::new()
  let seqs : Array[String] = []
  let labels : Array[String] = []
  let trained = @src.naive_bayes_train(clf, seqs, labels)
  assert_eq(trained.class_labels.length(), 0)
  assert_eq(trained.vocabulary.length(), 0)
}

test "predict_empty_sequence" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let pred = @src.naive_bayes_predict(trained, "")
  assert_true(trained.class_labels.contains(pred.0))
}

test "accuracy_empty_data" {
  let data = @src.naive_bayes_sample_data()
  let clf = @src.NaiveBayesClassifier::new()
  let trained = @src.naive_bayes_train(clf, data.0, data.1)
  let empty_seqs : Array[String] = []
  let empty_labels : Array[String] = []
  let acc = @src.naive_bayes_accuracy(trained, empty_seqs, empty_labels)
  assert_eq(acc, 0.0)
}
