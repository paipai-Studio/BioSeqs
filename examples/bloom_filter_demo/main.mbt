///|
/// Bloom filter and k-mer counting demo.
///
/// Demonstrates:
/// 1. Bloom filter for approximate set membership
/// 2. Exact k-mer counting with hash map
/// 3. Bloom-based approximate k-mer presence filtering
///
/// Run: moon run examples/bloom_filter_demo/main.mbt
fn main {
  println("=== Bloom Filter & K-mer Counting Demo ===")
  println("")

  // --- 1. Bloom Filter Basics ---
  println("1. Bloom Filter Basics:")
  let bf = @bio.BloomFilter::new(100, 0.01)
  println("   Capacity (bits): \{@bio.BloomFilter::capacity(bf)}")
  println("   Hash functions: \{@bio.BloomFilter::hash_count(bf)}")
  @bio.BloomFilter::insert(bf, "ACGTACGT")
  @bio.BloomFilter::insert(bf, "TTTTGGGG")
  @bio.BloomFilter::insert(bf, "CCCCAAAA")
  println("   Items inserted: \{@bio.BloomFilter::len(bf)}")
  println("   Bits set: \{@bio.BloomFilter::bits_set(bf)}")
  println(
    "   Contains 'ACGTACGT': \{@bio.BloomFilter::contains(bf, "ACGTACGT")}",
  )
  println(
    "   Contains 'GTGTCACA': \{@bio.BloomFilter::contains(bf, "GTGTCACA")}",
  )
  println("   Estimated FPR: \{@bio.BloomFilter::estimated_fpr(bf)}")
  println("")

  // --- 2. No False Negatives ---
  println("2. No False Negatives (all inserted items found):")
  let bf2 = @bio.BloomFilter::new(500, 0.05)
  let kmers = ["ACGT", "CGTA", "GTAC", "TACG", "AAAA", "TTTT", "GGGG", "CCCC"]
  for i = 0; i < kmers.length(); i = i + 1 {
    @bio.BloomFilter::insert(bf2, kmers[i])
  }
  let mut all_found = true
  for i = 0; i < kmers.length(); i = i + 1 {
    if not(@bio.BloomFilter::contains(bf2, kmers[i])) {
      all_found = false
    }
  }
  println("   All \{kmers.length()} items found: \{all_found}")
  println("   Estimated FPR: \{@bio.BloomFilter::estimated_fpr(bf2)}")
  println("")

  // --- 3. Clear and Reuse ---
  println("3. Clear and Reuse:")
  let bf3 = @bio.BloomFilter::with_params(1000, 3)
  @bio.BloomFilter::insert(bf3, "test1")
  @bio.BloomFilter::insert(bf3, "test2")
  println(
    "   Before clear: \{@bio.BloomFilter::len(bf3)} items, \{@bio.BloomFilter::bits_set(bf3)} bits set",
  )
  @bio.BloomFilter::clear(bf3)
  println(
    "   After clear: \{@bio.BloomFilter::len(bf3)} items, \{@bio.BloomFilter::bits_set(bf3)} bits set",
  )
  @bio.BloomFilter::insert(bf3, "new_item")
  println(
    "   After re-insert: \{@bio.BloomFilter::len(bf3)} items, \{@bio.BloomFilter::bits_set(bf3)} bits set",
  )
  println("")

  // --- 4. Exact K-mer Counting ---
  println("4. Exact K-mer Counting (k=3):")
  let kc = @bio.KmerCounter::new(3)
  @bio.KmerCounter::count(kc, "ACGTACGTACGT")
  println("   Sequence: ACGTACGTACGT")
  println("   k = \{@bio.KmerCounter::k(kc)}")
  println("   Total k-mers: \{@bio.KmerCounter::total(kc)}")
  println("   Unique k-mers: \{@bio.KmerCounter::unique(kc)}")
  println("   Count of 'ACG': \{@bio.KmerCounter::get(kc, "ACG")}")
  println("   Count of 'CGT': \{@bio.KmerCounter::get(kc, "CGT")}")
  println("   Count of 'GTA': \{@bio.KmerCounter::get(kc, "GTA")}")
  println("   Count of 'TTT': \{@bio.KmerCounter::get(kc, "TTT")}")
  println("")

  // --- 5. Top-k K-mers ---
  println("5. Top K-mers:")
  let top = @bio.KmerCounter::top_k(kc, 3)
  println("   Top 3 k-mers by count:")
  for i = 0; i < top.length(); i = i + 1 {
    println("     \{top[i].kmer()}: \{top[i].count()}")
  }
  println("")

  // --- 6. Multiple Sequence Counting ---
  println("6. Multiple Sequence Counting:")
  let kc2 = @bio.KmerCounter::new(4)
  @bio.KmerCounter::count(kc2, "ACGTACGT")
  @bio.KmerCounter::count(kc2, "ACGTAAAA")
  @bio.KmerCounter::count(kc2, "TTTTACGT")
  println("   Sequences: ACGTACGT, ACGTAAAA, TTTTACGT")
  println("   k = 4")
  println("   Total k-mers: \{@bio.KmerCounter::total(kc2)}")
  println("   Unique k-mers: \{@bio.KmerCounter::unique(kc2)}")
  println("   Count of 'ACGT': \{@bio.KmerCounter::get(kc2, "ACGT")}")
  println("")

  // --- 7. Bloom-based K-mer Filter ---
  println("7. Bloom-based K-mer Filter:")
  let bkf = @bio.BloomKmerFilter::new(3, 1000, 0.01)
  @bio.BloomKmerFilter::add_sequence(bkf, "ACGTACGTACGT")
  println("   Sequence indexed: ACGTACGTACGT")
  println("   k = \{@bio.BloomKmerFilter::k(bkf)}")
  println("   Items inserted: \{@bio.BloomKmerFilter::len(bkf)}")
  println("   Filter capacity: \{@bio.BloomKmerFilter::capacity(bkf)} bits")
  println("   Hash functions: \{@bio.BloomKmerFilter::hash_count(bkf)}")
  println("   Contains 'ACG': \{@bio.BloomKmerFilter::contains(bkf, "ACG")}")
  println("   Contains 'TTT': \{@bio.BloomKmerFilter::contains(bkf, "TTT")}")
  println("")

  // --- 8. Count Distinct K-mers ---
  println("8. Count Distinct K-mers (approximate):")
  let bkf2 = @bio.BloomKmerFilter::new(3, 10000, 0.001)
  let (distinct, duplicate) = @bio.BloomKmerFilter::count_distinct(
    bkf2, "ACGTACGT",
  )
  println("   Sequence: ACGTACGT")
  println("   Distinct k-mers: \{distinct}")
  println("   Duplicate k-mers: \{duplicate}")
  println("   Total k-mers: \{distinct + duplicate}")
  println("")

  // --- 9. FPR Comparison ---
  println("9. False Positive Rate Comparison:")
  let fpr_values = [0.001, 0.01, 0.05, 0.10]
  for i = 0; i < fpr_values.length(); i = i + 1 {
    let bf_fpr = @bio.BloomFilter::new(100, fpr_values[i])
    println(
      "   FPR=\{fpr_values[i]}: bits=\{@bio.BloomFilter::capacity(bf_fpr)}, hashes=\{@bio.BloomFilter::hash_count(bf_fpr)}",
    )
  }
  println("")

  // --- 10. Large-scale Simulation ---
  println("10. Large-scale Simulation (1000 items):")
  let bf_large = @bio.BloomFilter::new(1000, 0.01)
  for i = 0; i < 100; i = i + 1 {
    @bio.BloomFilter::insert(bf_large, "kmer_\{i}")
  }
  println("   Inserted: \{@bio.BloomFilter::len(bf_large)} items")
  println("   Bits set: \{@bio.BloomFilter::bits_set(bf_large)}")
  println("   Estimated FPR: \{@bio.BloomFilter::estimated_fpr(bf_large)}")
  // Verify no false negatives
  let mut found = 0
  for i = 0; i < 100; i = i + 1 {
    if @bio.BloomFilter::contains(bf_large, "kmer_\{i}") {
      found = found + 1
    }
  }
  println("   Found \{found}/100 inserted items (should be 100)")
  println("")

  println("=== Bloom Filter & K-mer Counting Demo Complete ===")
}
