///|
/// K-mer similarity table — fast sequence comparison via k-mer counting.
///
/// Instead of doing a full alignment, we can compare sequences by counting
/// their k-mer compositions and computing a similarity score. This is much
/// faster (O(n) per sequence) and useful for pre-filtering, clustering,
/// or estimating evolutionary distance.

///|
/// A k-mer count table for a sequence.
pub struct KmerTable {
  k : Int
  counts : Map[String, Int]  // k-mer -> count
  total : Int                 // total k-mers (sequence length - k + 1)
  alphabet : String           // "DNA", "RNA", or "protein"
}

///|
/// Build a k-mer count table from a sequence.
pub fn kmer_table_build(
  seq : String,
  k : Int,
  alphabet : String,
) -> KmerTable {
  let counts : Map[String, Int] = {}
  let n = seq.length()
  let mut total = 0
  if k <= 0 || n < k {
    return KmerTable::{ k, counts, total: 0, alphabet }
  }
  let valid_chars = kmer_valid_set(alphabet)
  for i = 0; i <= n - k; i = i + 1 {
    let mut valid = true
    for j = 0; j < k; j = j + 1 {
      let c = seq.unsafe_get(i + j).unsafe_to_char()
      let key = c.to_string()
      if not(valid_chars.contains(key)) {
        valid = false
        break
      }
    }
    if valid {
      let kmer = seq.substring(start=i, end=i + k)
      let cur = match counts.get(kmer) {
        Some(v) => v + 1
        None => 1
      }
      counts[kmer] = cur
      total = total + 1
    }
  }
  KmerTable::{ k, counts, total, alphabet }
}

///|
/// Compute the Jaccard similarity between two k-mer tables.
/// Jaccard = |A ∩ B| / |A ∪ B|
pub fn kmer_jaccard(a : KmerTable, b : KmerTable) -> Double {
  if a.k != b.k {
    return 0.0
  }
  let mut intersection = 0
  let mut union = 0
  let keys_a = a.counts.keys().collect()
  for i = 0; i < keys_a.length(); i = i + 1 {
    let key = keys_a[i]
    let ca = a.counts.get(key).unwrap()
    let cb = match b.counts.get(key) {
      Some(v) => v
      None => 0
    }
    intersection = intersection + min_int(ca, cb)
    union = union + kmer_max_int(ca, cb)
  }
  let keys_b = b.counts.keys().collect()
  for i = 0; i < keys_b.length(); i = i + 1 {
    let key = keys_b[i]
    if not(a.counts.contains(key)) {
      union = union + b.counts.get(key).unwrap()
    }
  }
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}

///|
/// Compute the cosine similarity between two k-mer tables.
pub fn kmer_cosine(a : KmerTable, b : KmerTable) -> Double {
  if a.k != b.k {
    return 0.0
  }
  let mut dot = 0.0
  let mut norm_a = 0.0
  let mut norm_b = 0.0
  let keys_a = a.counts.keys().collect()
  for i = 0; i < keys_a.length(); i = i + 1 {
    let key = keys_a[i]
    let ca = a.counts.get(key).unwrap().to_double()
    norm_a = norm_a + ca * ca
    let cb = match b.counts.get(key) {
      Some(v) => v.to_double()
      None => 0.0
    }
    dot = dot + ca * cb
  }
  let keys_b = b.counts.keys().collect()
  for i = 0; i < keys_b.length(); i = i + 1 {
    let key = keys_b[i]
    let cb = b.counts.get(key).unwrap().to_double()
    norm_b = norm_b + cb * cb
  }
  if norm_a == 0.0 || norm_b == 0.0 {
    return 0.0
  }
  dot / (kmer_sqrt(norm_a) * kmer_sqrt(norm_b))
}

///|
/// Compute the Bray-Curtis dissimilarity between two k-mer tables.
pub fn kmer_bray_curtis(a : KmerTable, b : KmerTable) -> Double {
  if a.k != b.k {
    return 1.0
  }
  let mut diff_sum = 0.0
  let mut total_sum = 0.0
  let keys_a = a.counts.keys().collect()
  for i = 0; i < keys_a.length(); i = i + 1 {
    let key = keys_a[i]
    let ca = a.counts.get(key).unwrap().to_double()
    let cb = match b.counts.get(key) {
      Some(v) => v.to_double()
      None => 0.0
    }
    diff_sum = diff_sum + kmer_abs(ca - cb)
    total_sum = total_sum + ca + cb
  }
  let keys_b = b.counts.keys().collect()
  for i = 0; i < keys_b.length(); i = i + 1 {
    let key = keys_b[i]
    if not(a.counts.contains(key)) {
      let cb = b.counts.get(key).unwrap().to_double()
      diff_sum = diff_sum + cb
      total_sum = total_sum + cb
    }
  }
  if total_sum == 0.0 {
    return 1.0
  }
  diff_sum / total_sum
}

///|
/// Estimate the Mash distance between two sequences.
/// d = -1/k * ln(jaccard)
pub fn kmer_mash_distance(a : KmerTable, b : KmerTable) -> Double {
  if a.k != b.k {
    return 1.0
  }
  let j = kmer_jaccard(a, b)
  if j <= 0.0 {
    return 1.0
  }
  if j >= 1.0 {
    return 0.0
  }
  let d = -1.0 / a.k.to_double() * @math.ln(j)
  if d > 1.0 {
    1.0
  } else {
    d
  }
}

///|
/// Build a k-mer similarity matrix for a set of sequences.
pub fn kmer_similarity_matrix(
  sequences : Array[String],
  k : Int,
  alphabet : String,
  metric : String,
) -> Array[Array[Double]] {
  let n = sequences.length()
  let tables : Array[KmerTable] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    tables.push(kmer_table_build(sequences[i], k, alphabet))
  }
  let matrix : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      let sim = if metric == "jaccard" {
        kmer_jaccard(tables[i], tables[j])
      } else if metric == "cosine" {
        kmer_cosine(tables[i], tables[j])
      } else if metric == "bray_curtis" {
        1.0 - kmer_bray_curtis(tables[i], tables[j])
      } else if metric == "mash" {
        1.0 - kmer_mash_distance(tables[i], tables[j])
      } else {
        kmer_jaccard(tables[i], tables[j])
      }
      row.push(sim)
    }
    matrix.push(row)
  }
  matrix
}

///|
/// Get the number of unique k-mers in the table.
pub fn kmer_unique_count(table : KmerTable) -> Int {
  table.counts.length()
}

///|
/// Get the total number of k-mers (including duplicates).
pub fn kmer_total_count(table : KmerTable) -> Int {
  table.total
}

///|
/// Get the k-mer size.
pub fn kmer_k(table : KmerTable) -> Int {
  table.k
}

// === Internal helpers ===

///|
fn kmer_valid_set(alphabet : String) -> Map[String, Bool] {
  let chars = if alphabet == "DNA" {
    ["A", "C", "G", "T"]
  } else if alphabet == "RNA" {
    ["A", "C", "G", "U"]
  } else if alphabet == "protein" {
    ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
     "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
  } else {
    ["A", "C", "G", "T"]
  }
  let m : Map[String, Bool] = {}
  for i = 0; i < chars.length(); i = i + 1 {
    m[chars[i]] = true
  }
  m
}

///|
fn kmer_max_int(a : Int, b : Int) -> Int {
  if a > b { a } else { b }
}

///|
fn kmer_sqrt(x : Double) -> Double {
  @math.exp(@math.ln(x) / 2.0)
}

///|
fn kmer_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}
