///|
/// Tests for Kmer module.

test "BioKmerCounter creation" {
  let counter = @src.create_example_kmer_counter()
  assert_eq(counter.k, 3)
}

test "BioKmerCounter count_kmers" {
  let counter = @src.BioKmerCounter::new(3)
  let counter = counter.count_kmers("ABCDEFGH")
  assert_eq(counter.unique_kmers(), 6)
}

test "BioKmerCounter count" {
  let counter = @src.BioKmerCounter::new(2)
  let counter = counter.count_kmers("AAAA")
  assert_eq(counter.count("AA"), 3)
}

test "BioKmerCounter total_kmers" {
  let counter = @src.BioKmerCounter::new(2)
  let counter = counter.count_kmers("ATCG")
  assert_eq(counter.total_kmers(), 3)
}

test "BioKmerCounter frequency" {
  let counter = @src.BioKmerCounter::new(2)
  let counter = counter.count_kmers("AAAA")
  let freq = counter.frequency("AA")
  assert_eq(freq, 1.0)
}

test "BioKmerCounter get_top_kmers" {
  let counter = @src.BioKmerCounter::new(2)
  let counter = counter.count_kmers("AATAA")
  let top = counter.get_top_kmers(2)
  assert_eq(top.length(), 2)
}

test "BioKmerCounter jaccard" {
  let counter1 = @src.BioKmerCounter::new(2)
  let counter1 = counter1.count_kmers("ATCG")
  let counter2 = @src.BioKmerCounter::new(2)
  let counter2 = counter2.count_kmers("ATGG")
  let jaccard = counter1.jaccard(counter2)
  assert_eq(jaccard, 0.2)
}

test "BioKmerCounter hamming_distance" {
  let counter1 = @src.BioKmerCounter::new(2)
  let counter1 = counter1.count_kmers("AAAA")
  let counter2 = @src.BioKmerCounter::new(2)
  let counter2 = counter2.count_kmers("AAAT")
  let dist = counter1.hamming_distance(counter2)
  assert_eq(dist, 2)
}

test "kmer_spectrum" {
  let spectrum = @src.kmer_spectrum("ATCG", 2, 3)
  assert_eq(spectrum.length(), 2)
}