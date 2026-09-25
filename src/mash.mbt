///|
/// Mash / MinHash: genomic sketching for fast ANI / distance / clustering.
///
/// References:
/// - Ondov BD et al. (2016) Mash: fast genome and metagenome distance
///   estimation using MinHash. Genome Biology 17:132.
/// - Broder AZ (1997) On the resemblance and containment of documents.
///   Compression and Complexity of Sequences.
///
/// Core idea: for each genome, compute hashes of all canonical k-mers,
/// keep the smallest `s` hashes (the "sketch"). Two sketches can be
/// compared in O(s) time via Jaccard similarity, then converted to
/// Mash distance and ANI.

// ============================================================
// MinHash hash function
// ============================================================

///|
/// 64-bit FNV-1a hash for a k-mer string.
/// Returns a UInt64 — used as the hash value for MinHash sketching.
/// We use unsigned arithmetic so that "smallest hash" comparison
/// matches the original Mash implementation.
fn mash_hash64(s : String) -> UInt64 {
  let mut hash : UInt64 = 14695981039346656037UL // FNV offset basis (64-bit)
  let prime : UInt64 = 1099511628211UL // FNV prime (64-bit)
  for i = 0; i < s.length(); i = i + 1 {
    hash = hash ^ (s.unsafe_get(i).to_uint64())
    hash = hash * prime
  }
  hash
}

// ============================================================
// Canonical k-mer
// ============================================================

///|
/// Compute the reverse complement of a DNA string (ACGT only).
/// For non-ACGT characters, the complement is the character itself
/// (matching Mash behavior — non-ACGT k-mers are skipped at a higher level).
fn mash_reverse_complement(s : String) -> String {
  let n = s.length()
  let buf = StringBuilder()
  // Build in reverse order
  let mut i = n
  while i > 0 {
    i = i - 1
    let c = s.unsafe_get(i)
    buf.write_char(
      match c {
        65 => 'T' // A
        84 => 'A' // T
        71 => 'C' // G
        67 => 'G' // C
        97 => 't' // a
        116 => 'a' // t
        103 => 'c' // g
        99 => 'g' // c
        _ => c.to_char().unwrap_or('?')
      },
    )
  }
  buf.to_string()
}

///|
/// Return the canonical form of a k-mer: the lexicographically smaller
/// of the k-mer and its reverse complement.
/// Mash uses canonical k-mers so that a sequence and its reverse strand
/// produce the same sketch.
fn canonical_kmer(kmer : String) -> String {
  let rc = mash_reverse_complement(kmer)
  if kmer <= rc {
    kmer
  } else {
    rc
  }
}

///|
/// Check if a string contains only ACGT characters (upper or lower case).
/// Non-ACGT k-mers are skipped during sketching.
fn is_acgt(s : String) -> Bool {
  for i = 0; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i)
    match c {
      65 | 67 | 71 | 84 | 97 | 99 | 103 | 116 => continue // ACGTacgt
      _ => return false
    }
  }
  true
}

// ============================================================
// MashSketch
// ============================================================

///|
/// A MinHash sketch of a genome: the `s` smallest canonical k-mer hashes.
pub struct MashSketch {
  name : String
  kmer_size : Int
  sketch_size : Int
  // Sorted ascending — the s smallest hash values
  hashes : Array[UInt64]
  // Total number of valid canonical k-mers hashed (before taking min-s)
  total_kmers : Int
}

///|
/// Get the name of the sketched genome.
pub fn MashSketch::name(self : MashSketch) -> String {
  self.name
}

///|
/// Get the k-mer size used for sketching.
pub fn MashSketch::kmer_size(self : MashSketch) -> Int {
  self.kmer_size
}

///|
/// Get the sketch size (number of hashes retained).
pub fn MashSketch::sketch_size(self : MashSketch) -> Int {
  self.sketch_size
}

///|
/// Get the total number of valid k-mers seen (before min-s selection).
pub fn MashSketch::total_kmers(self : MashSketch) -> Int {
  self.total_kmers
}

// ============================================================
// Sketch construction
// ============================================================

///|
/// Build a MinHash sketch from a DNA sequence.
///
/// Parameters:
/// - `name`: genome/sequence identifier
/// - `seq`: DNA sequence (FASTA content, single line, uppercase recommended)
/// - `k`: k-mer size (Mash default: 21)
/// - `s`: sketch size — number of smallest hashes to retain (Mash default: 1000)
///
/// Non-ACGT k-mers are skipped. Canonical k-mers (min of k-mer and its
/// reverse complement) are used, so the sketch is strand-independent.
pub fn mash_sketch(
  name : String,
  seq : String,
  k : Int,
  s : Int
) -> MashSketch {
  let n = seq.length()
  let mut total = 0
  // Use a max-heap-like approach: keep the s smallest hashes
  // For simplicity, collect all hashes then sort and take the smallest s
  // (Mash uses a priority queue; for MoonBit we use sort which is O(n log n))
  let all_hashes : Array[UInt64] = []

  if n >= k && k > 0 && s > 0 {
    let mut i = 0
    while i <= n - k {
      let kmer = seq[i:i + k].to_owned()
      if is_acgt(kmer) {
        let canon = canonical_kmer(kmer)
        let h = mash_hash64(canon)
        all_hashes.push(h)
        total = total + 1
      }
      i = i + 1
    }
  }

  // Sort ascending and take the smallest s
  all_hashes.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })

  let sketch_hashes : Array[UInt64] = []
  let limit = if all_hashes.length() < s { all_hashes.length() } else { s }
  for i = 0; i < limit; i = i + 1 {
    sketch_hashes.push(all_hashes[i])
  }

  MashSketch::{
    name,
    kmer_size: k,
    sketch_size: limit,
    hashes: sketch_hashes,
    total_kmers: total,
  }
}

// ============================================================
// Distance estimation
// ============================================================

///|
/// Count the number of shared hashes between two sketches
/// (intersection size). Both sketch arrays are sorted ascending.
fn count_shared(a : Array[UInt64], b : Array[UInt64]) -> Int {
  let mut i = 0
  let mut j = 0
  let mut shared = 0
  while i < a.length() && j < b.length() {
    if a[i] == b[j] {
      shared = shared + 1
      i = i + 1
      j = j + 1
    } else if a[i] < b[j] {
      i = i + 1
    } else {
      j = j + 1
    }
  }
  shared
}

///|
/// Jaccard similarity between two sketches: J = |intersection| / |union|.
/// If both sketches are empty, returns 0.0.
pub fn mash_jaccard(a : MashSketch, b : MashSketch) -> Double {
  if a.hashes.length() == 0 && b.hashes.length() == 0 {
    return 0.0
  }
  let shared = count_shared(a.hashes, b.hashes)
  let union = a.hashes.length() + b.hashes.length() - shared
  if union == 0 {
    return 0.0
  }
  shared.to_double() / union.to_double()
}

///|
/// Mash distance between two sketches.
///
/// Formula (Ondov et al. 2016):
///   d = -1/k * ln(2J / (1 + J))
/// where J = Jaccard similarity, k = k-mer size.
///
/// The distance is clipped to [0, 1]. For J = 0, distance = 1 (maximally
/// distant). For J = 1, distance = 0 (identical).
pub fn mash_distance(a : MashSketch, b : MashSketch) -> Double {
  let j = mash_jaccard(a, b)
  if j == 0.0 {
    return 1.0
  }
  if j == 1.0 {
    return 0.0
  }
  // d = -1/k * ln(2J / (1+J))
  let k = a.kmer_size.to_double()
  let ratio = 2.0 * j / (1.0 + j)
  // ln(x) for x in (0, 1) — ratio is always < 1 when j < 1
  let d = -1.0 / k * @math.ln(ratio)
  if d < 0.0 {
    0.0
  } else if d > 1.0 {
    1.0
  } else {
    d
  }
}

///|
/// Estimate Average Nucleotide Identity (ANI) from Mash distance.
///
/// ANI ≈ 1 - d (valid for small distances; Mash paper shows this
/// approximation holds well for d < 0.1, i.e., ANI > 90%).
pub fn mash_ani(a : MashSketch, b : MashSketch) -> Double {
  let d = mash_distance(a, b)
  1.0 - d
}

///|
/// P-value for the Mash distance, testing whether two sketches come from
/// the same genome. Uses a simple binomial model:
/// n = shared hashes, N = sketch size, p0 = 1 / (total possible k-mers)
/// For practical purposes, we use the Mash approximation:
///   p ≈ (1 - J)^s  (probability of s or fewer matches by chance)
/// This is a conservative upper bound.
pub fn mash_pvalue(a : MashSketch, b : MashSketch) -> Double {
  let j = mash_jaccard(a, b)
  let s = a.hashes.length().to_double()
  // Conservative: (1-J)^s
  @math.pow(1.0 - j, s)
}

// ============================================================
// Batch operations: distance matrix, clustering
// ============================================================

///|
/// Compute a pairwise Mash distance matrix for an array of sketches.
/// Returns a flat array of size n*n, where element [i*n + j] is the
/// distance between sketch i and sketch j.
/// The diagonal is 0.0.
pub fn mash_distance_matrix(sketches : Array[MashSketch]) -> Array[Double] {
  let n = sketches.length()
  let matrix : Array[Double] = Array(capacity=n * n)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i == j {
        matrix.push(0.0)
      } else {
        matrix.push(mash_distance(sketches[i], sketches[j]))
      }
    }
  }
  matrix
}

///|
/// A neighbor in a simple clustering result.
pub struct MashClusterMember {
  name : String
  distance : Double
}

///|
/// Simple single-linkage clustering from a Mash distance matrix.
/// Given a threshold, two genomes are in the same cluster if their
/// Mash distance is <= threshold. This is equivalent to connected
/// components in the graph where edges connect genomes within threshold.
///
/// Returns an array of clusters, each cluster is an array of member names.
/// The first member of each cluster is the "seed" (first encountered).
pub fn mash_cluster(
  sketches : Array[MashSketch],
  threshold : Double
) -> Array[Array[String]] {
  let n = sketches.length()
  let matrix = mash_distance_matrix(sketches)

  // Union-Find for connected components
  let parent : Array[Int] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    parent.push(i)
  }

  fn find(parent : Array[Int], x : Int) -> Int {
    let mut root = x
    while parent[root] != root {
      root = parent[root]
    }
    // Path compression
    let mut curr = x
    while parent[curr] != root {
      let next = parent[curr]
      parent[curr] = root
      curr = next
    }
    root
  }

  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if matrix[i * n + j] <= threshold {
        let ri = find(parent, i)
        let rj = find(parent, j)
        if ri != rj {
          parent[ri] = rj
        }
      }
    }
  }

  // Group by root
  let groups : Map[Int, Array[String]] = {}
  for i = 0; i < n; i = i + 1 {
    let root = find(parent, i)
    match groups.get(root) {
      Some(arr) => arr.push(sketches[i].name)
      None => groups[root] = [sketches[i].name]
    }
  }

  // Convert to array of arrays
  let clusters : Array[Array[String]] = []
  for _, members in groups {
    clusters.push(members)
  }
  clusters
}

// ============================================================
// Convenience: sketch from multiple sequences
// ============================================================

///|
/// Build sketches for multiple sequences at once.
/// Returns an array of MashSketch, one per input sequence.
pub fn mash_sketch_all(
  names : Array[String],
  seqs : Array[String],
  k : Int,
  s : Int
) -> Array[MashSketch] {
  let sketches : Array[MashSketch] = []
  let n = if names.length() < seqs.length() { names.length() } else { seqs.length() }
  for i = 0; i < n; i = i + 1 {
    sketches.push(mash_sketch(names[i], seqs[i], k, s))
  }
  sketches
}
