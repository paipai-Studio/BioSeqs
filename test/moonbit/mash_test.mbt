///|
/// Tests for Mash / MinHash genomic sketching.

// ============================================================
// Helper: deterministic test sequences
// ============================================================

///|
/// A simple DNA sequence for sketching (all ACGT, no N).
fn mash_test_seq_a() -> String {
  "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
}

///|
/// A slightly different DNA sequence (a few SNPs relative to seq_a).
fn mash_test_seq_b() -> String {
  "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA"
}

///|
/// A completely different DNA sequence.
fn mash_test_seq_c() -> String {
  "TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAA"
}

// ============================================================
// Sketch construction tests
// ============================================================

test "mash_sketch_basic" {
  let seq = mash_test_seq_a()
  let sketch = @src.mash_sketch("genome_a", seq, 5, 10)
  assert_eq(sketch.name, "genome_a")
  assert_eq(sketch.kmer_size, 5)
  // 64 - 5 + 1 = 60 k-mers, all ACGT
  assert_eq(sketch.total_kmers, 60)
  // Should have up to 10 hashes (sketch_size = 10)
  assert_true(sketch.sketch_size <= 10)
  assert_true(sketch.sketch_size > 0)
  // Hashes should be sorted ascending
  let mut sorted = true
  for i in 1..<sketch.hashes.length() {
    if sketch.hashes[i] < sketch.hashes[i - 1] {
      sorted = false
    }
  }
  assert_true(sorted)
}

test "mash_sketch_empty_sequence" {
  let sketch = @src.mash_sketch("empty", "", 21, 1000)
  assert_eq(sketch.total_kmers, 0)
  assert_eq(sketch.sketch_size, 0)
}

test "mash_sketch_short_sequence" {
  // Sequence shorter than k → no k-mers
  let sketch = @src.mash_sketch("short", "ACGT", 21, 1000)
  assert_eq(sketch.total_kmers, 0)
  assert_eq(sketch.sketch_size, 0)
}

test "mash_sketch_skips_non_acgt" {
  // Contains N → k-mers with N should be skipped
  let seq = "ACGTACGTACGTNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGT"
  let sketch = @src.mash_sketch("has_n", seq, 5, 100)
  // Some k-mers should be generated, but fewer than if all ACGT
  assert_true(sketch.total_kmers > 0)
  assert_true(sketch.total_kmers < 45) // 49 - 5 + 1 = 45, but N's reduce it
}

test "mash_sketch_canonical_strand_independence" {
  // A sequence and its reverse complement should produce the same sketch
  // (because we use canonical k-mers)
  let seq = "ACGTACGTGCAATGCAGGTACGTACGTACGTACGTACGTACGTACGTAC"
  let rc = @src.biostrings_reverse_complement(seq)
  let sketch1 = @src.mash_sketch("forward", seq, 7, 50)
  let sketch2 = @src.mash_sketch("reverse", rc, 7, 50)
  // Same number of k-mers
  assert_eq(sketch1.total_kmers, sketch2.total_kmers)
  // Jaccard should be ~1.0 (identical sketches)
  let j = @src.mash_jaccard(sketch1, sketch2)
  assert_true(j > 0.95)
}

// ============================================================
// Jaccard / distance / ANI tests
// ============================================================

test "mash_jaccard_identical" {
  let seq = mash_test_seq_a()
  let s1 = @src.mash_sketch("a1", seq, 5, 20)
  let s2 = @src.mash_sketch("a2", seq, 5, 20)
  let j = @src.mash_jaccard(s1, s2)
  // Identical sequences → Jaccard = 1.0
  assert_eq(j, 1.0)
}

test "mash_jaccard_different" {
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 5, 20)
  let s3 = @src.mash_sketch("c", mash_test_seq_c(), 5, 20)
  let j = @src.mash_jaccard(s1, s3)
  // Very different sequences → low Jaccard
  assert_true(j < 0.3)
}

test "mash_distance_identical" {
  let seq = mash_test_seq_a()
  let s1 = @src.mash_sketch("a1", seq, 5, 20)
  let s2 = @src.mash_sketch("a2", seq, 5, 20)
  let d = @src.mash_distance(s1, s2)
  // Identical → distance = 0
  assert_eq(d, 0.0)
}

test "mash_distance_similar" {
  // seq_a and seq_b differ by 1 SNP at the end
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 5, 20)
  let s2 = @src.mash_sketch("b", mash_test_seq_b(), 5, 20)
  let d = @src.mash_distance(s1, s2)
  // Similar sequences → small distance
  assert_true(d < 0.1)
}

test "mash_distance_completely_different" {
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 5, 20)
  let s3 = @src.mash_sketch("c", mash_test_seq_c(), 5, 20)
  let d = @src.mash_distance(s1, s3)
  // Very different → large distance
  assert_true(d > 0.5)
}

test "mash_distance_range_01" {
  // Distance should always be in [0, 1]
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 7, 15)
  let s2 = @src.mash_sketch("b", mash_test_seq_b(), 7, 15)
  let s3 = @src.mash_sketch("c", mash_test_seq_c(), 7, 15)
  let d1 = @src.mash_distance(s1, s2)
  let d2 = @src.mash_distance(s1, s3)
  let d3 = @src.mash_distance(s2, s3)
  assert_true(d1 >= 0.0 && d1 <= 1.0)
  assert_true(d2 >= 0.0 && d2 <= 1.0)
  assert_true(d3 >= 0.0 && d3 <= 1.0)
}

test "mash_ani_identical" {
  let seq = mash_test_seq_a()
  let s1 = @src.mash_sketch("a1", seq, 5, 20)
  let s2 = @src.mash_sketch("a2", seq, 5, 20)
  let ani = @src.mash_ani(s1, s2)
  // Identical → ANI = 1.0
  assert_eq(ani, 1.0)
}

test "mash_ani_different" {
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 5, 20)
  let s3 = @src.mash_sketch("c", mash_test_seq_c(), 5, 20)
  let ani = @src.mash_ani(s1, s3)
  // Different → ANI < 1.0
  assert_true(ani < 0.9)
}

// ============================================================
// P-value test
// ============================================================

test "mash_pvalue_identical" {
  let seq = mash_test_seq_a()
  let s1 = @src.mash_sketch("a1", seq, 5, 20)
  let s2 = @src.mash_sketch("a2", seq, 5, 20)
  let p = @src.mash_pvalue(s1, s2)
  // Identical → p-value should be very small (significant)
  assert_true(p < 0.01)
}

test "mash_pvalue_different" {
  let s1 = @src.mash_sketch("a", mash_test_seq_a(), 5, 20)
  let s3 = @src.mash_sketch("c", mash_test_seq_c(), 5, 20)
  let p = @src.mash_pvalue(s1, s3)
  // Different → p-value should be large (not significant)
  assert_true(p > 0.5)
}

// ============================================================
// Distance matrix test
// ============================================================

test "mash_distance_matrix" {
  let names = ["a", "b", "c"]
  let seqs = [mash_test_seq_a(), mash_test_seq_b(), mash_test_seq_c()]
  let sketches = @src.mash_sketch_all(names, seqs, 5, 20)
  let matrix = @src.mash_distance_matrix(sketches)
  // 3x3 matrix = 9 elements
  assert_eq(matrix.length(), 9)
  // Diagonal should be 0
  assert_eq(matrix[0], 0.0)
  assert_eq(matrix[4], 0.0)
  assert_eq(matrix[8], 0.0)
  // Symmetric
  assert_eq(matrix[1], matrix[3])
  assert_eq(matrix[2], matrix[6])
  assert_eq(matrix[5], matrix[7])
}

// ============================================================
// Clustering test
// ============================================================

test "mash_cluster_groups_similar" {
  let names = ["genome_a", "genome_a2", "genome_c"]
  let seqs = [mash_test_seq_a(), mash_test_seq_a(), mash_test_seq_c()]
  let sketches = @src.mash_sketch_all(names, seqs, 5, 20)
  // Threshold 0.1: a and a2 (identical) should cluster, c separate
  let clusters = @src.mash_cluster(sketches, 0.1)
  // Should have 2 clusters
  assert_eq(clusters.length(), 2)
  // One cluster should have 2 members (a, a2), one should have 1 (c)
  let mut has_size_2 = false
  let mut has_size_1 = false
  for c in clusters {
    if c.length() == 2 {
      has_size_2 = true
    }
    if c.length() == 1 {
      has_size_1 = true
    }
  }
  assert_true(has_size_2)
  assert_true(has_size_1)
}

test "mash_cluster_all_same" {
  let names = ["g1", "g2", "g3"]
  let seqs = [mash_test_seq_a(), mash_test_seq_a(), mash_test_seq_a()]
  let sketches = @src.mash_sketch_all(names, seqs, 5, 20)
  let clusters = @src.mash_cluster(sketches, 0.5)
  // All identical → 1 cluster with 3 members
  assert_eq(clusters.length(), 1)
  assert_eq(clusters[0].length(), 3)
}

test "mash_cluster_all_different" {
  let names = ["g1", "g2", "g3"]
  let seq_a = mash_test_seq_a()
  // Make 3 truly different sequences
  let seq_b = mash_test_seq_b()
  let seq_c = mash_test_seq_c()
  let seqs = [seq_a, seq_b, seq_c]
  let sketches = @src.mash_sketch_all(names, seqs, 5, 20)
  // Very low threshold → all separate
  let clusters = @src.mash_cluster(sketches, 0.0)
  // With threshold 0, only identical (distance=0) cluster together
  // a and b are similar but not identical, c is different
  // So a, b, c are all separate with threshold 0
  assert_eq(clusters.length(), 3)
}

// ============================================================
// Edge cases
// ============================================================

test "mash_distance_empty_sketches" {
  let s1 = @src.mash_sketch("empty1", "", 21, 1000)
  let s2 = @src.mash_sketch("empty2", "", 21, 1000)
  let d = @src.mash_distance(s1, s2)
  // Both empty → Jaccard = 0 → distance = 1.0
  assert_eq(d, 1.0)
}

test "mash_sketch_all_parallel" {
  let names = ["s1", "s2"]
  let seqs = [mash_test_seq_a(), mash_test_seq_c()]
  let sketches = @src.mash_sketch_all(names, seqs, 7, 15)
  assert_eq(sketches.length(), 2)
  assert_eq(sketches[0].name, "s1")
  assert_eq(sketches[1].name, "s2")
}
