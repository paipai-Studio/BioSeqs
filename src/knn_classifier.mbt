///|
/// k-Nearest Neighbors supervised classifier, ported from Biopython's
/// `Bio.kNN` module.
///
/// The Biopython semantics are preserved exactly:
/// - distance is the Manhattan (city block) distance between feature vectors;
/// - classification uses inverse-distance weighted voting over the k nearest
///   training samples;
/// - an exact match (distance 0) wins immediately.
///
/// Beyond the original API this module adds min-max / z-score normalization,
/// leave-one-out cross validation and a confusion matrix helper.

///|

/// A labeled training sample with an n-dimensional feature vector.
pub struct KnnSample {
  features : Array[Double]
  label : String
} derive(Eq)
pub extend KnnSample with Eq::{equal, not_equal}

///|
pub fn KnnSample::new(
  features~ : Array[Double],
  label~ : String,
) -> KnnSample {
  KnnSample::{ features, label, }
}

///|
/// One (label, weight) entry produced by the voting stage.
pub struct KnnVote {
  label : String
  weight : Double
} derive(Eq)
pub extend KnnVote with Eq::{equal, not_equal}

///|
/// Result of a weighted classification: winner plus the full vote table.
pub struct KnnClassification {
  label : String
  votes : Array[KnnVote]
} derive(Eq)
pub extend KnnClassification with Eq::{equal, not_equal}

///|
/// Create a new classification result.
pub fn KnnClassification::new(
  label~ : String,
  votes~ : Array[KnnVote],
) -> KnnClassification {
  KnnClassification::{ label, votes, }
}

///|
/// Look up the vote weight of a label, or 0.0 when absent.
pub fn KnnClassification::vote_of(
  self : KnnClassification,
  label : String,
) -> Double {
  for vote in self.votes {
    if vote.label == label {
      return vote.weight
    }
  }
  0.0
}

///|
/// Manhattan (city block) distance, identical to `Bio.kNN.distance`.
pub fn knn_distance(u : Array[Double], v : Array[Double]) -> Double {
  let n = if u.length() < v.length() { u.length() } else { v.length() }
  let mut total = 0.0
  for i in 0..<n {
    let diff = u[i] - v[i]
    total += if diff < 0.0 { 0.0 - diff } else { diff }
  }
  total
}

///|
/// Euclidean distance between two feature vectors.
pub fn knn_euclidean_distance(u : Array[Double], v : Array[Double]) -> Double {
  let n = if u.length() < v.length() { u.length() } else { v.length() }
  let mut total = 0.0
  for i in 0..<n {
    let diff = u[i] - v[i]
    total += diff * diff
  }
  total.sqrt()
}

///|
/// Majority label among exact matches (distance 0), if any.
fn knn_exact_match_label(
  train : Array[KnnSample],
  query : Array[Double],
) -> String? {
  let winners : Array[String] = Array(capacity=0)
  for sample in train {
    if knn_distance(sample.features, query) == 0.0 {
      if !(winners.contains(sample.label)) {
        winners.push(sample.label)
      }
    }
  }
  // Deterministic choice: the first distinct exact-match label, mirroring
  // the "first hit wins" behaviour of the reference implementation.
  if winners.is_empty() {
    None
  } else {
    Some(winners[0])
  }
}

///|
/// Inverse-distance weighted vote counts of the k nearest neighbors,
/// sorted by descending weight (ties broken by first occurrence),
/// mirroring `Bio.kNN.classify` internals.
pub fn knn_classify_votes(
  train : Array[KnnSample],
  k : Int,
  query : Array[Double],
) -> Array[KnnVote] {
  let distances : Array[(Double, String)] = Array(capacity=0)
  for sample in train {
    distances.push((knn_distance(sample.features, query), sample.label))
  }
  distances.sort_by((a, b) => compare_double(a.0, b.0))
  let take = if k < distances.length() { k } else { distances.length() }
  let votes : Map[String, Double] = Map([], capacity=take)
  let order : Array[String] = Array(capacity=0)
  for i in 0..<take {
    let (d, label) = distances[i]
    if d == 0.0 {
      // Exact match: keep weight 0 so the caller can short-circuit,
      // matching the reference "return immediately" semantics.
      continue
    }
    if !(votes.contains(label)) {
      votes.set(label, 0.0)
      order.push(label)
    }
    votes.set(label, votes.get(label).unwrap_or(0.0) + 1.0 / d)
  }
  let result : Array[KnnVote] = Array(capacity=0)
  for label in order {
    result.push(KnnVote::{ label, weight: votes.get(label).unwrap_or(0.0), })
  }
  result.sort_by((a, b) => compare_double(b.weight, a.weight))
  result
}

///|
/// Classify a query vector against the training set, returning the
/// winning label (the equivalent of `Bio.kNN.classify`).
pub fn knn_classify(
  train : Array[KnnSample],
  k : Int,
  query : Array[Double],
) -> String {
  match knn_exact_match_label(train, query) {
    Some(label) => return label
    None => ()
  }
  let votes = knn_classify_votes(train, k, query)
  if votes.is_empty() {
    ""
  } else {
    votes[0].label
  }
}

///|
/// Full classification result including the vote table.
pub fn knn_classify_full(
  train : Array[KnnSample],
  k : Int,
  query : Array[Double],
) -> KnnClassification {
  match knn_exact_match_label(train, query) {
    Some(label) =>
      return KnnClassification::new(label~, votes=[KnnVote::{ label, weight: 1.0, }])
    None => ()
  }
  let votes = knn_classify_votes(train, k, query)
  if votes.is_empty() {
    KnnClassification::new(label="", votes=[])
  } else {
    KnnClassification::new(label=votes[0].label, votes~, )
  }
}

///|
/// Normalized vote probabilities (weights divided by their sum).
pub fn knn_classify_probs(
  train : Array[KnnSample],
  k : Int,
  query : Array[Double],
) -> Array[KnnVote] {
  match knn_exact_match_label(train, query) {
    Some(label) => return [KnnVote::{ label, weight: 1.0, }]
    None => ()
  }
  let votes = knn_classify_votes(train, k, query)
  let mut total = 0.0
  for vote in votes {
    total += vote.weight
  }
  if total == 0.0 {
    return votes
  }
  votes.map(vote => KnnVote::{ label: vote.label, weight: vote.weight / total, })
}

///|
/// Leave-one-out cross validation accuracy on the training set.
pub fn knn_loo_accuracy(train : Array[KnnSample], k : Int) -> Double {
  if train.is_empty() {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<train.length() {
    let rest : Array[KnnSample] = Array(capacity=0)
    for j in 0..<train.length() {
      if j != i {
        rest.push(train[j])
      }
    }
    if knn_classify(rest, k, train[i].features) == train[i].label {
      correct += 1
    }
  }
  correct.to_double() / train.length().to_double()
}

///|
/// Distinct labels in order of first appearance.
pub fn knn_labels(train : Array[KnnSample]) -> Array[String] {
  let labels : Array[String] = Array(capacity=0)
  for sample in train {
    if !(labels.contains(sample.label)) {
      labels.push(sample.label)
    }
  }
  labels
}

///|
/// Confusion matrix over the distinct labels of the training set.
/// Returns the ordered label list and an n x n matrix where the entry
/// `[i][j]` counts samples with true label `labels[i]` predicted as
/// `labels[j]`.
pub fn knn_confusion_matrix(
  train : Array[KnnSample],
  k : Int,
) -> (Array[String], Array[Array[Int]]) {
  let labels = knn_labels(train)
  let index : Map[String, Int] = Map([], capacity=labels.length())
  for i in 0..<labels.length() {
    index.set(labels[i], i)
  }
  let matrix : Array[Array[Int]] = Array(capacity=0)
  for _i in 0..<labels.length() {
    let row : Array[Int] = Array(capacity=0)
    for _j in 0..<labels.length() {
      row.push(0)
    }
    matrix.push(row)
  }
  for sample in train {
    let predicted = knn_classify(train, k, sample.features)
    match (index.get(sample.label), index.get(predicted)) {
      (Some(i), Some(j)) => matrix[i][j] += 1
      _ => ()
    }
  }
  (labels, matrix)
}

///|
/// Min-max normalize a feature column to [0, 1]. Constant columns map to 0.
pub fn knn_minmax_normalize_column(values : Array[Double]) -> Array[Double] {
  if values.is_empty() {
    return []
  }
  let mut lo = values[0]
  let mut hi = values[0]
  for v in values {
    if v < lo {
      lo = v
    }
    if v > hi {
      hi = v
    }
  }
  let span = hi - lo
  values.map(v => if span == 0.0 { 0.0 } else { (v - lo) / span })
}

///|
/// Z-score normalize a feature column (mean 0, standard deviation 1).
pub fn knn_zscore_normalize_column(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  if n == 0 {
    return []
  }
  let mut mean = 0.0
  for v in values {
    mean += v
  }
  mean = mean / n.to_double()
  let mut variance = 0.0
  for v in values {
    let d = v - mean
    variance += d * d
  }
  variance = variance / n.to_double()
  let sd = variance.sqrt()
  values.map(v => if sd == 0.0 { 0.0 } else { (v - mean) / sd })
}

///|
/// Normalize every column (feature dimension) of the samples.
/// `mode` is `"minmax"` or `"zscore"`.
pub fn knn_normalize_samples(
  train : Array[KnnSample],
  mode~ : String = "minmax",
) -> Array[KnnSample] {
  if train.is_empty() {
    return []
  }
  let n_dims = train[0].features.length()
  let normalized : Array[Array[Double]] = Array(capacity=0)
  for dim in 0..<n_dims {
    let column : Array[Double] = Array(capacity=0)
    for sample in train {
      column.push(sample.features[dim])
    }
    normalized.push(
      if mode == "zscore" {
        knn_zscore_normalize_column(column)
      } else {
        knn_minmax_normalize_column(column)
      },
    )
  }
  let result : Array[KnnSample] = Array(capacity=0)
  for i in 0..<train.length() {
    let feats : Array[Double] = Array(capacity=0)
    for dim in 0..<n_dims {
      feats.push(normalized[dim][i])
    }
    result.push(KnnSample::new(features=feats, label=train[i].label))
  }
  result
}

///|
fn compare_double(a : Double, b : Double) -> Int {
  if a < b {
    -1
  } else if a > b {
    1
  } else {
    0
  }
}

///|
/// Sample data in the spirit of the Biopython operon-prediction tutorial:
/// gene-pair features followed by a filled `-1.0` column, labeled by
/// whether the pair lies in the same operon.
pub fn knn_sample_train() -> Array[KnnSample] {
  [
    KnnSample::new(features=[0, -1, 0.9], label="SameOperon"),
    KnnSample::new(features=[1, -1, 0.2], label="DifferentOperon"),
    KnnSample::new(features=[2, 0, 0.8], label="SameOperon"),
    KnnSample::new(features=[3, -1, 0.1], label="DifferentOperon"),
    KnnSample::new(features=[2, 0, 0.75], label="SameOperon"),
    KnnSample::new(features=[3, -1, 0.15], label="DifferentOperon"),
  ]
}

///|
/// A held-out query pair used in the reference tutorial.
pub fn knn_sample_query() -> Array[Double] {
  [3, -1, 0.12]
}