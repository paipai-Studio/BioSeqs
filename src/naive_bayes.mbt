///|
/// Bio.NaiveBayes - Naive Bayes sequence classifier using k-mer frequency features,
/// inspired by Biopython Bio.NaiveBayes.

///|
/// Class-specific model for Naive Bayes classifier.
pub struct NBClassModel {
  class_label : String
  n_sequences : Int
  total_kmers : Int
  kmer_counts : Map[String, Double]
  class_prior : Double
}

///|
/// Naive Bayes classifier for biological sequences based on k-mer frequencies.
pub struct NaiveBayesClassifier {
  kmer_size : Int
  alpha : Double
  vocabulary : Array[String]
  class_labels : Array[String]
  models : Map[String, NBClassModel]
}

///|
/// Create a new NaiveBayesClassifier with default kmer_size=3 and alpha=1.0 (Laplace smoothing).
pub fn NaiveBayesClassifier::new() -> NaiveBayesClassifier {
  NaiveBayesClassifier::{
    kmer_size: 3,
    alpha: 1.0,
    vocabulary: [],
    class_labels: [],
    models: Map([], capacity=10)
  }
}

///|
/// Set the k-mer size and return a new classifier instance.
pub fn NaiveBayesClassifier::set_kmer_size(self : NaiveBayesClassifier, val : Int) -> NaiveBayesClassifier {
  NaiveBayesClassifier::{
    kmer_size: val,
    alpha: self.alpha,
    vocabulary: self.vocabulary,
    class_labels: self.class_labels,
    models: self.models
  }
}

///|
/// Set the Laplace smoothing alpha parameter and return a new classifier instance.
pub fn NaiveBayesClassifier::set_alpha(self : NaiveBayesClassifier, val : Double) -> NaiveBayesClassifier {
  NaiveBayesClassifier::{
    kmer_size: self.kmer_size,
    alpha: val,
    vocabulary: self.vocabulary,
    class_labels: self.class_labels,
    models: self.models
  }
}

///|
/// Extract k-mers from a sequence using a sliding window of size k.
/// Returns a map of k-mer to count. If k > sequence length, returns empty map.
/// For protein sequences, any characters are kept.
pub fn naive_bayes_extract_kmers(sequence : String, k : Int) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=100)
  let n = sequence.length()
  if k <= 0 || k > n {
    return counts
  }
  let num_windows = n - k + 1
  let mut i = 0
  while i < num_windows {
    let kmer = substring(sequence, i, k)
    if counts.contains(kmer) {
      counts.set(kmer, counts.get(kmer).unwrap() + 1)
    } else {
      counts.set(kmer, 1)
    }
    i = i + 1
  }
  counts
}

///|
/// Train the Naive Bayes classifier on sequences with their class labels.
/// Builds per-class k-mer counts, total k-mers per class, class priors, and vocabulary.
pub fn naive_bayes_train(
  classifier : NaiveBayesClassifier,
  sequences : Array[String],
  labels : Array[String]
) -> NaiveBayesClassifier {
  let n = sequences.length()
  if n == 0 || n != labels.length() {
    return NaiveBayesClassifier::{
      kmer_size: classifier.kmer_size,
      alpha: classifier.alpha,
      vocabulary: [],
      class_labels: [],
      models: Map([], capacity=10)
    }
  }

  let k = classifier.kmer_size

  let class_seq_counts : Map[String, Int] = Map([], capacity=10)
  let class_total_kmers : Map[String, Int] = Map([], capacity=10)
  let class_kmer_counts : Map[String, Map[String, Double]] = Map([], capacity=10)
  let vocab_set : Map[String, Bool] = Map([], capacity=1000)

  for idx = 0; idx < n; idx = idx + 1 {
    let seq = sequences[idx]
    let label = labels[idx]
    let kmers = naive_bayes_extract_kmers(seq, k)

    if class_seq_counts.contains(label) {
      class_seq_counts.set(label, class_seq_counts.get(label).unwrap() + 1)
    } else {
      class_seq_counts.set(label, 1)
      class_total_kmers.set(label, 0)
      class_kmer_counts.set(label, Map([], capacity=1000))
    }

    let mut total_for_class = class_total_kmers.get(label).unwrap()
    let counts_for_class = class_kmer_counts.get(label).unwrap()

    for kmer in kmers.keys() {
      let kmer_count = kmers.get(kmer).unwrap()
      total_for_class = total_for_class + kmer_count
      if counts_for_class.contains(kmer) {
        counts_for_class.set(kmer, counts_for_class.get(kmer).unwrap() + kmer_count.to_double())
      } else {
        counts_for_class.set(kmer, kmer_count.to_double())
      }
      if !vocab_set.contains(kmer) {
        vocab_set.set(kmer, true)
      }
    }
    class_total_kmers.set(label, total_for_class)
  }

  let class_labels_arr : Array[String] = class_seq_counts.keys().collect()
  let vocabulary_arr : Array[String] = vocab_set.keys().collect()

  let models_map : Map[String, NBClassModel] = Map([], capacity=class_labels_arr.length())
  for label in class_labels_arr {
    let n_seq = class_seq_counts.get(label).unwrap()
    let prior = n_seq.to_double() / n.to_double()
    let total_k = class_total_kmers.get(label).unwrap()
    let kmer_counts_copy = class_kmer_counts.get(label).unwrap()
    let model = NBClassModel::{
      class_label: label,
      n_sequences: n_seq,
      total_kmers: total_k,
      kmer_counts: kmer_counts_copy,
      class_prior: prior
    }
    models_map.set(label, model)
  }

  NaiveBayesClassifier::{
    kmer_size: classifier.kmer_size,
    alpha: classifier.alpha,
    vocabulary: vocabulary_arr,
    class_labels: class_labels_arr,
    models: models_map
  }
}

///|
/// Predict log probabilities for each class given a sequence.
/// log P(class|seq) proportional to log(P(class)) + sum_i log(P(kmer_i|class))
/// Uses Laplace smoothing for k-mer probabilities.
pub fn naive_bayes_predict_log_probs(
  classifier : NaiveBayesClassifier,
  sequence : String
) -> Array[(String, Double)] {
  let result : Array[(String, Double)] = Array::new()
  let vocab_size = classifier.vocabulary.length().to_double()
  let alpha = classifier.alpha
  let kmers = naive_bayes_extract_kmers(sequence, classifier.kmer_size)

  for class_label in classifier.class_labels {
    let model_option = classifier.models.get(class_label)
    if model_option is None {
      continue
    }
    let model = model_option.unwrap()
    let mut log_prob = @math.ln(model.class_prior)
    let denom = model.total_kmers.to_double() + alpha * vocab_size

    for kmer in kmers.keys() {
      let kmer_count_in_seq = kmers.get(kmer).unwrap()
      let count_in_model = if model.kmer_counts.contains(kmer) {
        model.kmer_counts.get(kmer).unwrap()
      } else {
        0.0
      }
      let smoothed_prob = (count_in_model + alpha) / denom
      log_prob = log_prob + kmer_count_in_seq.to_double() * @math.ln(smoothed_prob)
    }

    result.push((class_label, log_prob))
  }

  result
}

///|
/// Predict the most likely class for a sequence.
/// Returns (best_label, log_probability) = argmax of predict_log_probs.
pub fn naive_bayes_predict(
  classifier : NaiveBayesClassifier,
  sequence : String
) -> (String, Double) {
  let log_probs = naive_bayes_predict_log_probs(classifier, sequence)
  if log_probs.length() == 0 {
    return ("", 0.0)
  }
  let mut best_label = log_probs[0].0
  let mut best_prob = log_probs[0].1
  for i = 1; i < log_probs.length(); i = i + 1 {
    if log_probs[i].1 > best_prob {
      best_label = log_probs[i].0
      best_prob = log_probs[i].1
    }
  }
  (best_label, best_prob)
}

///|
/// Predict class probabilities using softmax normalization of log probabilities.
/// Returns array of (class_label, probability) where probabilities sum to ~1.0.
pub fn naive_bayes_predict_proba(
  classifier : NaiveBayesClassifier,
  sequence : String
) -> Array[(String, Double)] {
  let log_probs = naive_bayes_predict_log_probs(classifier, sequence)
  let result : Array[(String, Double)] = Array::new()
  if log_probs.length() == 0 {
    return result
  }

  let mut max_log = log_probs[0].1
  for i = 1; i < log_probs.length(); i = i + 1 {
    if log_probs[i].1 > max_log {
      max_log = log_probs[i].1
    }
  }

  let exps : Array[Double] = Array::new()
  let mut sum = 0.0
  for lp in log_probs {
    let e = @math.exp(lp.1 - max_log)
    exps.push(e)
    sum = sum + e
  }

  for i = 0; i < log_probs.length(); i = i + 1 {
    let prob = if sum == 0.0 { 1.0 / log_probs.length().to_double() } else { exps[i] / sum }
    result.push((log_probs[i].0, prob))
  }

  result
}

///|
/// Return the top-k most probable classes for a sequence, sorted by probability descending.
pub fn naive_bayes_top_k(
  classifier : NaiveBayesClassifier,
  sequence : String,
  k : Int
) -> Array[(String, Double)] {
  let probs = naive_bayes_predict_proba(classifier, sequence)

  let mut i = 0
  while i < probs.length() {
    let mut j = i + 1
    while j < probs.length() {
      if probs[j].1 > probs[i].1 {
        let temp = probs[i]
        probs[i] = probs[j]
        probs[j] = temp
      }
      j = j + 1
    }
    i = i + 1
  }

  let result : Array[(String, Double)] = Array::new()
  let limit = if k < probs.length() { k } else { probs.length() }
  let mut idx = 0
  while idx < limit {
    result.push(probs[idx])
    idx = idx + 1
  }
  result
}

///|
/// Compute classification accuracy on a test set of sequences with true labels.
/// Returns the fraction of correctly classified sequences (0.0 to 1.0).
pub fn naive_bayes_accuracy(
  classifier : NaiveBayesClassifier,
  sequences : Array[String],
  true_labels : Array[String]
) -> Double {
  let n = sequences.length()
  if n == 0 || n != true_labels.length() {
    return 0.0
  }
  let mut correct = 0
  for i = 0; i < n; i = i + 1 {
    let pred = naive_bayes_predict(classifier, sequences[i])
    if pred.0 == true_labels[i] {
      correct = correct + 1
    }
  }
  correct.to_double() / n.to_double()
}

///|
/// Create sample training data: 8 DNA sequences.
/// 4 sequences with label "AT_rich" (high A/T content) and
/// 4 sequences with label "GC_rich" (high G/C content).
pub fn naive_bayes_sample_data() -> (Array[String], Array[String]) {
  let sequences : Array[String] = [
    "ATATATATATATATATATAT",
    "AATTAAAATTTTAAATTTTA",
    "TATATAAATTATATTAATTT",
    "ATTTAAATAAATTTTAAAAT",
    "GCGCGCGCGCGCGCGCGCGC",
    "GGCCGGCCCGGGCCCGGCCC",
    "CGGCCGGGCCCCGGGCCTGG",
    "GCGCCGCGGGGCCCGGCCGC"
  ]
  let labels : Array[String] = [
    "AT_rich",
    "AT_rich",
    "AT_rich",
    "AT_rich",
    "GC_rich",
    "GC_rich",
    "GC_rich",
    "GC_rich"
  ]
  (sequences, labels)
}
