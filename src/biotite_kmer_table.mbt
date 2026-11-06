/// biotite_kmer_table.mbt — K-mer lookup table with perfect/space-efficient hashing
/// Inspired by biotite.sequence.align.KmerTable, BucketKmerTable, KmerAlphabet

///|
/// K-mer alphabet: treats k-mers as individual symbols in an extended alphabet.
/// For DNA with k=3, alphabet size = 4^3 = 64.
pub struct BiotiteKmerAlphabet {
  base_size : Int32 // size of base alphabet (e.g. 4 for DNA)
  k : Int32 // k-mer size
  alphabet_size : Int32 // base_size^k
}

///|
/// Create a k-mer alphabet.
pub fn BiotiteKmerAlphabet::new(
  base_size : Int32,
  k : Int32,
) -> BiotiteKmerAlphabet {
  let mut size = 1
  for _i = 0; _i < k; _i = _i + 1 {
    size = size * base_size
  }
  { base_size, k, alphabet_size: size, }
}

///|
/// Encode a k-mer (array of base symbol codes) into a single integer code.
pub fn BiotiteKmerAlphabet::encode(
  self : BiotiteKmerAlphabet,
  kmer : Array[Int32],
) -> Int32 {
  let mut code = 0
  for i = 0; i < self.k; i = i + 1 {
    code = code * self.base_size + kmer[i]
  }
  code
}

///|
/// Decode an integer code back to k-mer symbol codes.
pub fn BiotiteKmerAlphabet::decode(
  self : BiotiteKmerAlphabet,
  code : Int32,
) -> Array[Int32] {
  let result : Array[Int32] = Array::make(self.k, 0)
  let mut c = code
  for i = self.k - 1; i >= 0; i = i - 1 {
    result[i] = c % self.base_size
    c = c / self.base_size
  }
  result
}

///|
/// K-mer table using perfect hashing (direct array indexing).
/// Maps each k-mer code to a list of positions in a reference sequence.
pub struct BiotiteKmerTable {
  alphabet : BiotiteKmerAlphabet
  table : Array[Array[Int32]] // kmer_code -> positions
}

///|
/// Build a k-mer table from a sequence of symbol codes.
pub fn BiotiteKmerTable::build(
  alphabet : BiotiteKmerAlphabet,
  sequence : Array[Int32],
) -> BiotiteKmerTable {
  let table : Array[Array[Int32]] = Array::make(alphabet.alphabet_size, [])
  let n = sequence.length()
  if n >= alphabet.k {
    for pos = 0; pos <= n - alphabet.k; pos = pos + 1 {
      let kmer : Array[Int32] = Array(capacity=0)
      for i = 0; i < alphabet.k; i = i + 1 {
        kmer.push(sequence[pos + i])
      }
      let code = alphabet.encode(kmer)
      table[code].push(pos)
    }
  }
  { alphabet, table, }
}

///|
/// Lookup positions of a k-mer in the table.
pub fn BiotiteKmerTable::lookup(
  self : BiotiteKmerTable,
  kmer : Array[Int32],
) -> Array[Int32] {
  let code = self.alphabet.encode(kmer)
  if code >= 0 && code < self.table.length() {
    self.table[code]
  } else {
    []
  }
}

///|
/// Lookup by code directly.
pub fn BiotiteKmerTable::lookup_code(
  self : BiotiteKmerTable,
  code : Int32,
) -> Array[Int32] {
  if code >= 0 && code < self.table.length() {
    self.table[code]
  } else {
    []
  }
}

///|
/// Get total number of k-mers stored.
pub fn BiotiteKmerTable::size(self : BiotiteKmerTable) -> Int32 {
  let mut total = 0
  for i = 0; i < self.table.length(); i = i + 1 {
    total = total + self.table[i].length()
  }
  total
}

///|
/// Bucket k-mer table: space-efficient hashing using buckets.
/// Instead of one entry per possible k-mer, uses hash buckets.
pub struct BiotiteBucketKmerTable {
  alphabet : BiotiteKmerAlphabet
  n_buckets : Int32
  buckets : Array[Array[(Int32, Array[Int32])]] // (kmer_code, positions) per bucket
}

///|
/// Hash function for bucket assignment.
fn biotite_kmer_hash(code : Int32, n_buckets : Int32) -> Int32 {
  // Simple multiplicative hash
  let golden = 265443576
  let h = code * golden % n_buckets
  if h < 0 {
    -h
  } else {
    h
  }
}

///|
/// Build a bucket k-mer table.
pub fn BiotiteBucketKmerTable::build(
  alphabet : BiotiteKmerAlphabet,
  sequence : Array[Int32],
  n_buckets~ : Int32 = 1024,
) -> BiotiteBucketKmerTable {
  let buckets : Array[Array[(Int32, Array[Int32])]] = Array::make(n_buckets, [])
  let n = sequence.length()
  if n >= alphabet.k {
    for pos = 0; pos <= n - alphabet.k; pos = pos + 1 {
      let kmer : Array[Int32] = Array(capacity=0)
      for i = 0; i < alphabet.k; i = i + 1 {
        kmer.push(sequence[pos + i])
      }
      let code = alphabet.encode(kmer)
      let bucket_idx = biotite_kmer_hash(code, n_buckets)
      // Find or create entry
      let mut found = false
      for j = 0; j < buckets[bucket_idx].length(); j = j + 1 {
        if buckets[bucket_idx][j].0 == code {
          buckets[bucket_idx][j].1.push(pos)
          found = true
          break
        }
      }
      if !found {
        buckets[bucket_idx].push((code, [pos]))
      }
    }
  }
  { alphabet, n_buckets, buckets, }
}

///|
/// Lookup in bucket k-mer table.
pub fn BiotiteBucketKmerTable::lookup(
  self : BiotiteBucketKmerTable,
  kmer : Array[Int32],
) -> Array[Int32] {
  let code = self.alphabet.encode(kmer)
  let bucket_idx = biotite_kmer_hash(code, self.n_buckets)
  for j = 0; j < self.buckets[bucket_idx].length(); j = j + 1 {
    if self.buckets[bucket_idx][j].0 == code {
      return self.buckets[bucket_idx][j].1
    }
  }
  []
}

///|
/// Calculate optimal bucket number for a given sequence length and k.
pub fn biotite_bucket_number(
  sequence_length : Int32,
  k : Int32,
  _base_size : Int32,
) -> Int32 {
  // Number of distinct k-mers expected
  let n_kmers = sequence_length - k + 1
  if n_kmers <= 0 {
    return 1
  }
  // Target load factor ~0.75
  let n = n_kmers * 4 / 3
  // Round up to next power of 2
  let mut result = 1
  while result < n {
    result = result * 2
  }
  if result < 16 {
    16
  } else {
    result
  }
}
