///|
/// Tests for Bloom filter and k-mer counting.

// ---------------------------------------------------------------------------
// BloomFilter tests
// ---------------------------------------------------------------------------

///|
test "bloom_insert_and_contains" {
  let bf = @src.BloomFilter::new(100, 0.01)
  @src.BloomFilter::insert(bf, "ACGT")
  @src.BloomFilter::insert(bf, "GGTA")
  @src.BloomFilter::insert(bf, "CCGC")
  assert_true(@src.BloomFilter::contains(bf, "ACGT"))
  assert_true(@src.BloomFilter::contains(bf, "GGTA"))
  assert_true(@src.BloomFilter::contains(bf, "CCGC"))
}

///|
test "bloom_no_false_negative" {
  // Inserted items must always be found (no false negatives).
  let bf = @src.BloomFilter::new(50, 0.05)
  let items = ["ACGT", "TTTT", "GGGG", "CCCC", "ATAT"]
  for i = 0; i < items.length(); i = i + 1 {
    @src.BloomFilter::insert(bf, items[i])
  }
  for i = 0; i < items.length(); i = i + 1 {
    assert_true(@src.BloomFilter::contains(bf, items[i]))
  }
}

///|
test "bloom_absent_items" {
  // Items never inserted should (mostly) return false.
  // With a properly sized filter, false positive rate is low.
  let bf = @src.BloomFilter::new(1000, 0.01)
  @src.BloomFilter::insert(bf, "ACGT")
  @src.BloomFilter::insert(bf, "GGTA")
  // "TTTT" was never inserted — should be absent.
  assert_false(@src.BloomFilter::contains(bf, "TTTT"))
}

///|
test "bloom_len" {
  let bf = @src.BloomFilter::new(100, 0.01)
  assert_eq(@src.BloomFilter::len(bf), 0)
  @src.BloomFilter::insert(bf, "A")
  @src.BloomFilter::insert(bf, "B")
  @src.BloomFilter::insert(bf, "C")
  assert_eq(@src.BloomFilter::len(bf), 3)
}

///|
test "bloom_capacity" {
  let bf = @src.BloomFilter::new(100, 0.01)
  let cap = @src.BloomFilter::capacity(bf)
  assert_true(cap > 0)
}

///|
test "bloom_hash_count" {
  let bf = @src.BloomFilter::new(100, 0.01)
  let hc = @src.BloomFilter::hash_count(bf)
  assert_true(hc >= 1)
}

///|
test "bloom_bits_set" {
  let bf = @src.BloomFilter::with_params(1000, 3)
  assert_eq(@src.BloomFilter::bits_set(bf), 0)
  @src.BloomFilter::insert(bf, "ACGT")
  // After one insert, at most 3 bits should be set.
  assert_true(@src.BloomFilter::bits_set(bf) <= 3)
  assert_true(@src.BloomFilter::bits_set(bf) >= 1)
}

///|
test "bloom_clear" {
  let bf = @src.BloomFilter::with_params(1000, 3)
  @src.BloomFilter::insert(bf, "ACGT")
  @src.BloomFilter::insert(bf, "GGTA")
  assert_true(@src.BloomFilter::bits_set(bf) > 0)
  @src.BloomFilter::clear(bf)
  assert_eq(@src.BloomFilter::bits_set(bf), 0)
  assert_eq(@src.BloomFilter::len(bf), 0)
  assert_false(@src.BloomFilter::contains(bf, "ACGT"))
}

///|
test "bloom_estimated_fpr_empty" {
  let bf = @src.BloomFilter::new(100, 0.01)
  // Empty filter should have FPR ≈ 0
  let fpr = @src.BloomFilter::estimated_fpr(bf)
  assert_true(fpr < 0.001)
}

///|
test "bloom_estimated_fpr_increases" {
  // As we insert more items, FPR should increase.
  let bf = @src.BloomFilter::new(100, 0.01)
  let fpr_empty = @src.BloomFilter::estimated_fpr(bf)
  for i = 0; i < 50; i = i + 1 {
    @src.BloomFilter::insert(bf, "item\{i}")
  }
  let fpr_half = @src.BloomFilter::estimated_fpr(bf)
  assert_true(fpr_half > fpr_empty)
}

///|
test "bloom_with_params" {
  let bf = @src.BloomFilter::with_params(500, 4)
  assert_eq(@src.BloomFilter::capacity(bf), 500)
  assert_eq(@src.BloomFilter::hash_count(bf), 4)
  @src.BloomFilter::insert(bf, "test")
  assert_true(@src.BloomFilter::contains(bf, "test"))
}

///|
test "bloom_many_items" {
  // Insert many items and verify no false negatives.
  let bf = @src.BloomFilter::new(500, 0.05)
  let items = [
    "ACGT", "CGTA", "GTAC", "TACG", "AAAA", "TTTT", "GGGG", "CCCC", "ATGC", "TGGC",
  ]
  for i = 0; i < items.length(); i = i + 1 {
    @src.BloomFilter::insert(bf, items[i])
  }
  for i = 0; i < items.length(); i = i + 1 {
    assert_true(@src.BloomFilter::contains(bf, items[i]))
  }
  assert_eq(@src.BloomFilter::len(bf), 10)
}

// ---------------------------------------------------------------------------
// KmerCounter tests
// ---------------------------------------------------------------------------

///|
test "kmer_counter_basic" {
  let kc = @src.KmerCounter::new(3)
  @src.KmerCounter::count(kc, "ACGTACGT")
  // 3-mers: ACG, CGT, GTA, TAC, ACG, CGT → 6 total, 4 unique
  assert_eq(@src.KmerCounter::total(kc), 6)
  assert_eq(@src.KmerCounter::unique(kc), 4)
  assert_eq(@src.KmerCounter::get(kc, "ACG"), 2)
  assert_eq(@src.KmerCounter::get(kc, "CGT"), 2)
  assert_eq(@src.KmerCounter::get(kc, "GTA"), 1)
  assert_eq(@src.KmerCounter::get(kc, "TAC"), 1)
}

///|
test "kmer_counter_k_value" {
  let kc = @src.KmerCounter::new(5)
  assert_eq(@src.KmerCounter::k(kc), 5)
}

///|
test "kmer_counter_absent" {
  let kc = @src.KmerCounter::new(3)
  @src.KmerCounter::count(kc, "ACGT")
  assert_eq(@src.KmerCounter::get(kc, "TTT"), 0)
}

///|
test "kmer_counter_multiple_sequences" {
  let kc = @src.KmerCounter::new(2)
  @src.KmerCounter::count(kc, "ACGT")
  @src.KmerCounter::count(kc, "ACGT")
  // 2-mers: AC, CG, GT × 2 = 6 total, 3 unique
  assert_eq(@src.KmerCounter::total(kc), 6)
  assert_eq(@src.KmerCounter::unique(kc), 3)
  assert_eq(@src.KmerCounter::get(kc, "AC"), 2)
}

///|
test "kmer_counter_short_sequence" {
  // Sequence shorter than k → no k-mers
  let kc = @src.KmerCounter::new(5)
  @src.KmerCounter::count(kc, "ACG")
  assert_eq(@src.KmerCounter::total(kc), 0)
  assert_eq(@src.KmerCounter::unique(kc), 0)
}

///|
test "kmer_counter_top_k" {
  let kc = @src.KmerCounter::new(3)
  @src.KmerCounter::count(kc, "ACGTACGTACGT")
  // ACG appears 3 times, CGT 3 times, GTA 2 times, TAC 2 times
  let top = @src.KmerCounter::top_k(kc, 2)
  assert_eq(top.length(), 2)
  assert_eq(top[0].count(), 3)
  assert_eq(top[1].count(), 3)
}

///|
test "kmer_counter_exact_k" {
  // Sequence length equals k → exactly one k-mer
  let kc = @src.KmerCounter::new(4)
  @src.KmerCounter::count(kc, "ACGT")
  assert_eq(@src.KmerCounter::total(kc), 1)
  assert_eq(@src.KmerCounter::unique(kc), 1)
  assert_eq(@src.KmerCounter::get(kc, "ACGT"), 1)
}

///|
test "kmer_counter_repeated" {
  let kc = @src.KmerCounter::new(1)
  @src.KmerCounter::count(kc, "AAA")
  assert_eq(@src.KmerCounter::total(kc), 3)
  assert_eq(@src.KmerCounter::unique(kc), 1)
  assert_eq(@src.KmerCounter::get(kc, "A"), 3)
}

// ---------------------------------------------------------------------------
// BloomKmerFilter tests
// ---------------------------------------------------------------------------

///|
test "bloom_kmer_filter_basic" {
  let bkf = @src.BloomKmerFilter::new(3, 100, 0.01)
  @src.BloomKmerFilter::add_sequence(bkf, "ACGTACGT")
  // All 3-mers from the sequence should be present (no false negatives)
  assert_true(@src.BloomKmerFilter::contains(bkf, "ACG"))
  assert_true(@src.BloomKmerFilter::contains(bkf, "CGT"))
  assert_true(@src.BloomKmerFilter::contains(bkf, "GTA"))
  assert_true(@src.BloomKmerFilter::contains(bkf, "TAC"))
}

///|
test "bloom_kmer_filter_absent" {
  let bkf = @src.BloomKmerFilter::new(3, 1000, 0.01)
  @src.BloomKmerFilter::add_sequence(bkf, "ACGTACGT")
  // "TTT" was never in the sequence
  assert_false(@src.BloomKmerFilter::contains(bkf, "TTT"))
}

///|
test "bloom_kmer_filter_k" {
  let bkf = @src.BloomKmerFilter::new(5, 100, 0.01)
  assert_eq(@src.BloomKmerFilter::k(bkf), 5)
}

///|
test "bloom_kmer_filter_len" {
  let bkf = @src.BloomKmerFilter::new(3, 100, 0.01)
  assert_eq(@src.BloomKmerFilter::len(bkf), 0)
  @src.BloomKmerFilter::add_sequence(bkf, "ACGTACGT")
  // 6 k-mers inserted
  assert_eq(@src.BloomKmerFilter::len(bkf), 6)
}

///|
test "bloom_kmer_filter_count_distinct" {
  let bkf = @src.BloomKmerFilter::new(3, 1000, 0.01)
  let (distinct, duplicate) = @src.BloomKmerFilter::count_distinct(
    bkf, "ACGTACGT",
  )
  // 3-mers: ACG, CGT, GTA, TAC, ACG, CGT
  // First pass: 4 distinct, then ACG and CGT are duplicates
  assert_eq(distinct, 4)
  assert_eq(duplicate, 2)
}

///|
test "bloom_kmer_filter_count_distinct_all_unique" {
  // Use a large filter to minimize false positives.
  let bkf = @src.BloomKmerFilter::new(3, 10000, 0.001)
  let (distinct, duplicate) = @src.BloomKmerFilter::count_distinct(bkf, "ACGT")
  // "ACGT" has two 3-mers: ACG, CGT — both distinct.
  // With a large filter, both should be counted as distinct.
  assert_eq(distinct, 2)
  assert_eq(duplicate, 0)
}

///|
test "bloom_kmer_filter_capacity" {
  let bkf = @src.BloomKmerFilter::new(3, 100, 0.01)
  let cap = @src.BloomKmerFilter::capacity(bkf)
  assert_true(cap > 0)
}

///|
test "bloom_kmer_filter_hash_count" {
  let bkf = @src.BloomKmerFilter::new(3, 100, 0.01)
  let hc = @src.BloomKmerFilter::hash_count(bkf)
  assert_true(hc >= 1)
}

///|
test "bloom_kmer_filter_estimated_fpr" {
  let bkf = @src.BloomKmerFilter::new(3, 1000, 0.01)
  let fpr = @src.BloomKmerFilter::estimated_fpr(bkf)
  // Empty filter → low FPR
  assert_true(fpr < 0.001)
}

///|
test "bloom_kmer_filter_short_sequence" {
  let bkf = @src.BloomKmerFilter::new(5, 100, 0.01)
  @src.BloomKmerFilter::add_sequence(bkf, "ACG")
  // Sequence too short for k=5
  assert_eq(@src.BloomKmerFilter::len(bkf), 0)
}
