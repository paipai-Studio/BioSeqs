///|
/// Bloom filter and hash-based k-mer counting.
///
/// References:
/// - Bloom BH (1970). "Space/time trade-offs in hash coding with allowable
///   errors." Communications of the ACM. 13(7):422-426.
/// - Marcais G, Kingsford C (2011). "A fast, lock-free approach for efficient
///   parallel counting of occurrences of k-mers." Bioinformatics. 27(6):764-70.

// ---------------------------------------------------------------------------
// Hash functions for Bloom filter
// ---------------------------------------------------------------------------

///|
/// Compute a hash using FNV-1a variant.
/// Simple, deterministic, and fast for short strings (k-mers).
fn fnv1a_hash(s : String) -> Int {
  let mut hash = 16777619 // FNV prime as offset (fits Int32)
  for i = 0; i < s.length(); i = i + 1 {
    hash = hash ^ s.unsafe_get(i).unsafe_to_char().to_int()
    hash = hash * 16777619
  }
  if hash < 0 {
    hash = -hash
  }
  hash
}

///|
/// Compute a second independent hash using a polynomial rolling hash.
fn poly_hash(s : String) -> Int {
  let mut hash = 0
  for i = 0; i < s.length(); i = i + 1 {
    hash = hash * 131 + s.unsafe_get(i).unsafe_to_char().to_int()
  }
  if hash < 0 {
    hash = -hash
  }
  if hash == 0 {
    hash = 1 // ensure non-zero for double-hashing
  }
  hash
}

///|
/// Generate `num_hashes` independent bit positions from two base hashes
/// using the double-hashing scheme: h_i = (h1 + i * h2) % bits.
fn bloom_positions(
  h1 : Int,
  h2 : Int,
  num_hashes : Int,
  bits : Int,
) -> Array[Int] {
  let positions : Array[Int] = Array::new()
  for i = 0; i < num_hashes; i = i + 1 {
    let combined = h1 + i * h2
    let pos = combined % bits
    positions.push(if pos < 0 { -pos } else { pos })
  }
  positions
}

// ---------------------------------------------------------------------------
// BitArray — compact bit storage using Int array
// ---------------------------------------------------------------------------

///|
/// Compact bit array backed by an Int array.
struct BitArray {
  bits : FixedArray[Int]
  num_bits : Int
}

///|
fn BitArray::new(num_bits : Int) -> BitArray {
  BitArray::{ bits: FixedArray::make(num_bits, 0), num_bits }
}

///|
fn BitArray::set(self : BitArray, idx : Int) -> Unit {
  self.bits[idx] = 1
}

///|
fn BitArray::get(self : BitArray, idx : Int) -> Bool {
  self.bits[idx] != 0
}

///|
fn BitArray::count_set(self : BitArray) -> Int {
  let mut count = 0
  for i = 0; i < self.num_bits; i = i + 1 {
    if self.bits[i] != 0 {
      count = count + 1
    }
  }
  count
}

// ---------------------------------------------------------------------------
// BloomFilter
// ---------------------------------------------------------------------------

///|
/// A Bloom filter: a space-efficient probabilistic data structure for
/// approximate set membership testing.
///
/// - No false negatives: if `contains` returns `false`, the element was
///   definitely never inserted.
/// - Possible false positives: if `contains` returns `true`, the element
///   *might* have been inserted (with probability `false_positive_rate`).
pub struct BloomFilter {
  bit_array : BitArray
  num_hashes : Int
  num_bits : Int
  mut count : Int // number of items inserted
}

///|
/// Create a Bloom filter sized for `expected_items` elements with a target
/// false-positive rate of `false_positive_rate` (0.0 – 1.0).
///
/// The optimal number of bits and hash functions are computed from the
/// standard formulas:
///   m = -(n * ln(p)) / (ln(2)^2)
///   k = (m / n) * ln(2)
///
/// Since MoonBit's `Double` lacks `log`/`exp`, we use a precomputed
/// table for common FPR values and a fallback linear approximation.
pub fn BloomFilter::new(
  expected_items : Int,
  false_positive_rate : Double,
) -> BloomFilter {
  let n = expected_items.to_double()
  let p = false_positive_rate
  // Approximate m/n ratio: m/n ≈ -ln(p) / (ln2)^2 = -1.44 * ln(p)
  // For common FPR values:
  //   p=0.01 → m/n ≈ 9.6, k ≈ 7
  //   p=0.05 → m/n ≈ 6.2, k ≈ 4
  //   p=0.10 → m/n ≈ 4.8, k ≈ 3
  let bits_per_item = if p <= 0.001 {
    19.2
  } else if p <= 0.01 {
    9.6
  } else if p <= 0.05 {
    6.2
  } else if p <= 0.10 {
    4.8
  } else {
    3.0 // fallback for high FPR
  }
  let k = if p <= 0.001 {
    10
  } else if p <= 0.01 {
    7
  } else if p <= 0.05 {
    4
  } else if p <= 0.10 {
    3
  } else {
    2
  }
  let num_bits_raw = (n * bits_per_item).to_int()
  let num_bits = if num_bits_raw < 1 { 1 } else { num_bits_raw }
  BloomFilter::{
    bit_array: BitArray::new(num_bits),
    num_hashes: k,
    num_bits,
    count: 0,
  }
}

///|
/// Create a Bloom filter with a fixed number of bits and hash functions.
pub fn BloomFilter::with_params(
  num_bits : Int,
  num_hashes : Int,
) -> BloomFilter {
  BloomFilter::{
    bit_array: BitArray::new(num_bits),
    num_hashes,
    num_bits,
    count: 0,
  }
}

///|
/// Insert an element into the Bloom filter.
pub fn BloomFilter::insert(self : BloomFilter, key : String) -> Unit {
  let h1 = fnv1a_hash(key)
  let h2 = poly_hash(key)
  let positions = bloom_positions(h1, h2, self.num_hashes, self.num_bits)
  for i = 0; i < positions.length(); i = i + 1 {
    self.bit_array.set(positions[i])
  }
  self.count = self.count + 1
}

///|
/// Test whether an element *might* be in the set.
///
/// Returns `true` if the element is possibly present (could be a false
/// positive), or `false` if it is definitely not present.
pub fn BloomFilter::contains(self : BloomFilter, key : String) -> Bool {
  let h1 = fnv1a_hash(key)
  let h2 = poly_hash(key)
  let positions = bloom_positions(h1, h2, self.num_hashes, self.num_bits)
  for i = 0; i < positions.length(); i = i + 1 {
    if not(self.bit_array.get(positions[i])) {
      return false
    }
  }
  true
}

///|
/// Get the number of items inserted.
pub fn BloomFilter::len(self : BloomFilter) -> Int {
  self.count
}

///|
/// Get the number of bits in the filter.
pub fn BloomFilter::capacity(self : BloomFilter) -> Int {
  self.num_bits
}

///|
/// Get the number of hash functions used.
pub fn BloomFilter::hash_count(self : BloomFilter) -> Int {
  self.num_hashes
}

///|
/// Get the number of bits currently set to 1.
pub fn BloomFilter::bits_set(self : BloomFilter) -> Int {
  self.bit_array.count_set()
}

///|
/// Estimate the current false-positive rate based on the fill ratio.
///
/// Formula: p ≈ (1 - e^(-kn/m))^k
/// Approximated using Taylor series: e^x ≈ 1 + x + x^2/2
pub fn BloomFilter::estimated_fpr(self : BloomFilter) -> Double {
  let n = self.count.to_double()
  let k = self.num_hashes.to_double()
  let m = self.num_bits.to_double()
  if m == 0.0 {
    return 0.0
  }
  // ratio of bits set ≈ 1 - e^(-kn/m)
  let x = -k * n / m
  // Taylor approximation of e^x: 1 + x + x*x/2 + x*x*x/6
  let exp_x = 1.0 + x + x * x / 2.0 + x * x * x / 6.0
  let ratio = 1.0 - exp_x
  if ratio <= 0.0 {
    return 0.0
  }
  // p ≈ ratio^k — use @math.pow for integer exponent
  let mut result = 1.0
  let mut i = 0
  while i < self.num_hashes {
    result = result * ratio
    i = i + 1
  }
  result
}

///|
/// Clear all bits, resetting the filter to empty.
pub fn BloomFilter::clear(self : BloomFilter) -> Unit {
  for i = 0; i < self.bit_array.num_bits; i = i + 1 {
    self.bit_array.bits[i] = 0
  }
  self.count = 0
}

// ---------------------------------------------------------------------------
// KmerCounter — exact counting via Map
// ---------------------------------------------------------------------------

///|
/// Result entry for a single k-mer count.
pub struct KmerCount {
  kmer : String
  count : Int
}

///|
/// Get the k-mer string.
pub fn KmerCount::kmer(self : KmerCount) -> String {
  self.kmer
}

///|
/// Get the count value.
pub fn KmerCount::count(self : KmerCount) -> Int {
  self.count
}

///|
/// Exact k-mer counter using a hash map.
///
/// Suitable for small-to-medium datasets where exact counts are required.
/// For large datasets where approximate membership suffices, use
/// `BloomFilter` instead.
pub struct KmerCounter {
  k : Int
  counts : Map[String, Int]
  mut total_kmers : Int
}

///|
/// Create a k-mer counter for k-mer size `k`.
pub fn KmerCounter::new(k : Int) -> KmerCounter {
  KmerCounter::{ k, counts: Map([], capacity=0), total_kmers: 0 }
}

///|
/// Count all k-mers in a sequence.
pub fn KmerCounter::count(self : KmerCounter, seq : String) -> Unit {
  let n = seq.length()
  if self.k <= 0 || self.k > n {
    return
  }
  let total = n - self.k + 1
  for i = 0; i < total; i = i + 1 {
    let kmer = seq[i:i + self.k].to_owned()
    match self.counts.get(kmer) {
      Some(v) => self.counts[kmer] = v + 1
      None => self.counts[kmer] = 1
    }
    self.total_kmers = self.total_kmers + 1
  }
}

///|
/// Get the count for a specific k-mer. Returns 0 if not found.
pub fn KmerCounter::get(self : KmerCounter, kmer : String) -> Int {
  match self.counts.get(kmer) {
    Some(v) => v
    None => 0
  }
}

///|
/// Get the total number of k-mers counted (including duplicates).
pub fn KmerCounter::total(self : KmerCounter) -> Int {
  self.total_kmers
}

///|
/// Get the number of unique k-mers.
pub fn KmerCounter::unique(self : KmerCounter) -> Int {
  self.counts.size()
}

///|
/// Get the k value.
pub fn KmerCounter::k(self : KmerCounter) -> Int {
  self.k
}

///|
/// Get all k-mer counts as an array of `KmerCount` entries,
/// sorted by count (descending).
pub fn KmerCounter::top_k(self : KmerCounter, n : Int) -> Array[KmerCount] {
  let entries : Array[KmerCount] = Array::new()
  let keys = self.counts.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    let kmer = keys[i]
    let cnt = match self.counts.get(kmer) {
      Some(v) => v
      None => 0
    }
    entries.push(KmerCount::{ kmer, count: cnt })
  }
  // Simple selection sort for top-n (efficient for small n).
  let result : Array[KmerCount] = Array::new()
  let len = entries.length()
  let mut remaining = entries
  let limit = if n < len { n } else { len }
  for iter = 0; iter < limit; iter = iter + 1 {
    let mut max_idx = 0
    let mut max_val = -1
    for j = 0; j < remaining.length(); j = j + 1 {
      if remaining[j].count > max_val {
        max_val = remaining[j].count
        max_idx = j
      }
    }
    result.push(remaining[max_idx])
    remaining = remove_at(remaining, max_idx)
  }
  result
}

///|
/// Remove element at index from array, returning a new array.
fn remove_at(arr : Array[KmerCount], idx : Int) -> Array[KmerCount] {
  let result : Array[KmerCount] = Array::new()
  for i = 0; i < arr.length(); i = i + 1 {
    if i != idx {
      result.push(arr[i])
    }
  }
  result
}

// ---------------------------------------------------------------------------
// BloomKmerFilter — approximate k-mer presence using Bloom filter
// ---------------------------------------------------------------------------

///|
/// Approximate k-mer presence filter using a Bloom filter.
///
/// Useful for quickly testing whether a k-mer has been seen before,
/// without storing exact counts. Memory-efficient for large genomes.
pub struct BloomKmerFilter {
  bloom : BloomFilter
  k : Int
}

///|
/// Create a Bloom-based k-mer filter for k-mer size `k`,
/// sized for `expected_kmers` unique k-mers.
pub fn BloomKmerFilter::new(
  k : Int,
  expected_kmers : Int,
  false_positive_rate : Double,
) -> BloomKmerFilter {
  BloomKmerFilter::{
    bloom: BloomFilter::new(expected_kmers, false_positive_rate),
    k,
  }
}

///|
/// Add all k-mers from a sequence into the Bloom filter.
pub fn BloomKmerFilter::add_sequence(
  self : BloomKmerFilter,
  seq : String,
) -> Unit {
  let n = seq.length()
  if self.k <= 0 || self.k > n {
    return
  }
  let total = n - self.k + 1
  for i = 0; i < total; i = i + 1 {
    let kmer = seq[i:i + self.k].to_owned()
    self.bloom.insert(kmer)
  }
}

///|
/// Test whether a k-mer might be present.
/// Returns `true` if possibly present (may be false positive),
/// `false` if definitely absent.
pub fn BloomKmerFilter::contains(self : BloomKmerFilter, kmer : String) -> Bool {
  self.bloom.contains(kmer)
}

///|
/// Get the k value.
pub fn BloomKmerFilter::k(self : BloomKmerFilter) -> Int {
  self.k
}

///|
/// Get the number of items inserted.
pub fn BloomKmerFilter::len(self : BloomKmerFilter) -> Int {
  self.bloom.len()
}

///|
/// Get the underlying Bloom filter's estimated false-positive rate.
pub fn BloomKmerFilter::estimated_fpr(self : BloomKmerFilter) -> Double {
  self.bloom.estimated_fpr()
}

///|
/// Get the number of bits used by the underlying Bloom filter.
pub fn BloomKmerFilter::capacity(self : BloomKmerFilter) -> Int {
  self.bloom.capacity()
}

///|
/// Get the number of hash functions used.
pub fn BloomKmerFilter::hash_count(self : BloomKmerFilter) -> Int {
  self.bloom.hash_count()
}

///|
/// Count approximate distinct k-mers present in a sequence.
/// Uses the Bloom filter to test each k-mer: if not seen, it's new;
/// if seen, it's possibly a duplicate (may have false positives).
///
/// Returns `(approx_distinct, approx_duplicate)`.
pub fn BloomKmerFilter::count_distinct(
  self : BloomKmerFilter,
  seq : String,
) -> (Int, Int) {
  let n = seq.length()
  if self.k <= 0 || self.k > n {
    return (0, 0)
  }
  let total = n - self.k + 1
  let mut distinct = 0
  let mut duplicate = 0
  for i = 0; i < total; i = i + 1 {
    let kmer = seq[i:i + self.k].to_owned()
    if self.bloom.contains(kmer) {
      duplicate = duplicate + 1
    } else {
      distinct = distinct + 1
      self.bloom.insert(kmer)
    }
  }
  (distinct, duplicate)
}
