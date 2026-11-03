/// Unified sequence distance interface.
///
/// scikit-bio equivalent: `sequence_distance(seq1, seq2, method=...)`
///
/// Provides a single API for computing distances between sequences
/// using various methods: hamming, levenshtein, kmer, normalized.

///|
/// Compute Hamming distance between two sequences.
/// Returns -1 if sequences have different lengths.
pub fn skbio_hamming_distance(seq1 : String, seq2 : String) -> Int32 {
  if seq1.length() != seq2.length() {
    return -1
  }
  let mut count : Int32 = 0
  for i = 0; i < seq1.length(); i = i + 1 {
    if seq1.unsafe_get(i).unsafe_to_char() !=
      seq2.unsafe_get(i).unsafe_to_char() {
      count = count + 1
    }
  }
  count
}

///|
/// Compute normalized Hamming distance (0.0 to 1.0).
/// Returns -1.0 if sequences have different lengths.
pub fn skbio_hamming_normalized(seq1 : String, seq2 : String) -> Double {
  if seq1.length() != seq2.length() {
    return -1.0
  }
  if seq1.length() == 0 {
    return 0.0
  }
  skbio_hamming_distance(seq1, seq2).to_double() / seq1.length().to_double()
}

///|
/// Compute Levenshtein (edit) distance.
/// Allows insertions, deletions, and substitutions.
pub fn skbio_levenshtein_distance(seq1 : String, seq2 : String) -> Int32 {
  let m = seq1.length()
  let n = seq2.length()
  if m == 0 {
    return n
  }
  if n == 0 {
    return m
  }
  // Use two rows for space efficiency
  let mut prev : Array[Int32] = Array::make(n + 1, 0)
  let mut curr : Array[Int32] = Array::make(n + 1, 0)
  for j = 0; j <= n; j = j + 1 {
    prev[j] = j
  }
  for i = 1; i <= m; i = i + 1 {
    curr[0] = i
    for j = 1; j <= n; j = j + 1 {
      let cost = if seq1.unsafe_get(i - 1).unsafe_to_char() ==
        seq2.unsafe_get(j - 1).unsafe_to_char() {
        0
      } else {
        1
      }
      let del = prev[j] + 1
      let ins = curr[j - 1] + 1
      let sub = prev[j - 1] + cost
      curr[j] = if del <= ins && del <= sub {
        del
      } else if ins <= sub {
        ins
      } else {
        sub
      }
    }
    let tmp = prev
    prev = curr
    curr = tmp
  }
  prev[n]
}

///|
/// Compute normalized Levenshtein distance (0.0 to 1.0).
/// Normalized by the length of the longer sequence.
pub fn skbio_levenshtein_normalized(seq1 : String, seq2 : String) -> Double {
  let max_len = if seq1.length() > seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  if max_len == 0 {
    return 0.0
  }
  skbio_levenshtein_distance(seq1, seq2).to_double() / max_len.to_double()
}

///|
/// Compute k-mer Jaccard distance between two sequences.
/// Distance = 1 - |intersection| / |union|
pub fn skbio_kmer_jaccard_distance(
  seq1 : String,
  seq2 : String,
  k : Int32,
) -> Double {
  let kmers1 = skbio_extract_kmers(seq1, k)
  let kmers2 = skbio_extract_kmers(seq2, k)
  let set1 : Set[String] = Set::from_array(kmers1)
  let set2 : Set[String] = Set::from_array(kmers2)
  let mut intersection = 0
  let mut union = 0
  let all_keys : Set[String] = Set::from_array(kmers1)
  for s in kmers2 {
    all_keys.add(s)
  }
  union = all_keys.size()
  for s in set1.to_array() {
    if set2.contains(s) {
      intersection = intersection + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  1.0 - intersection.to_double() / union.to_double()
}

///|
/// Extract all k-mers from a sequence.
fn skbio_extract_kmers(seq : String, k : Int32) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let n = seq.length()
  if k <= 0 || k > n {
    return result
  }
  for i = 0; i <= n - k; i = i + 1 {
    result.push(seq.substring(start=i, end=i + k))
  }
  result
}

///|
/// Compute sequence distance using a specified method.
///
/// Methods:
/// - "hamming": Hamming distance (requires equal length)
/// - "hamming_norm": Normalized Hamming (0-1)
/// - "levenshtein": Edit distance
/// - "levenshtein_norm": Normalized edit distance (0-1)
/// - "kmer_jaccard": 1 - Jaccard similarity of k-mer sets
pub fn skbio_sequence_distance(
  seq1 : String,
  seq2 : String,
  method~ : String = "levenshtein",
  k~ : Int32 = 3,
) -> Double {
  match method {
    "hamming" => skbio_hamming_distance(seq1, seq2).to_double()
    "hamming_norm" => skbio_hamming_normalized(seq1, seq2)
    "levenshtein" => skbio_levenshtein_distance(seq1, seq2).to_double()
    "levenshtein_norm" => skbio_levenshtein_normalized(seq1, seq2)
    "kmer_jaccard" => skbio_kmer_jaccard_distance(seq1, seq2, k)
    _ => skbio_levenshtein_distance(seq1, seq2).to_double()
  }
}

///|
/// Compute pairwise distance matrix for a set of sequences.
pub fn skbio_sequence_distance_matrix(
  sequences : Array[String],
  method~ : String = "levenshtein",
  k~ : Int32 = 3,
) -> Array[Array[Double]] {
  let n = sequences.length()
  let result : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    result[i] = Array::make(n, 0.0)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let d = skbio_sequence_distance(sequences[i], sequences[j], method~, k~)
      result[i][j] = d
      result[j][i] = d
    }
  }
  result
}
