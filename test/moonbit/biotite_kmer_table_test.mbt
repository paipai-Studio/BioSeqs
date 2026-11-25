///|
/// Tests for biotite_kmer_table.mbt — k-mer alphabets and tables.
test "biotite_kmer_alphabet_encode_decode" {
  let alphabet = @src.BiotiteKmerAlphabet::new(4, 3)
  assert_eq(alphabet.alphabet_size, 64)
  let code = alphabet.encode([1, 2, 3])
  assert_eq(code, 27)
  let decoded = alphabet.decode(27)
  assert_eq(decoded.length(), 3)
  assert_eq(decoded[0], 1)
  assert_eq(decoded[1], 2)
  assert_eq(decoded[2], 3)
}

///|
test "biotite_kmer_table_lookup" {
  let alphabet = @src.BiotiteKmerAlphabet::new(4, 3)
  let table = @src.BiotiteKmerTable::build(alphabet, [0, 1, 2, 3, 0, 1, 2])
  assert_eq(table.size(), 5)
  let hits = table.lookup([0, 1, 2])
  assert_eq(hits.length(), 2)
  assert_eq(hits[0], 0)
  assert_eq(hits[1], 4)
  assert_eq(table.lookup_code(27).length(), 1)
}

///|
test "biotite_bucket_kmer_table_lookup" {
  let alphabet = @src.BiotiteKmerAlphabet::new(4, 3)
  let table = @src.BiotiteBucketKmerTable::build(
    alphabet,
    [0, 1, 2, 3, 0, 1, 2],
    n_buckets=64,
  )
  let hits = table.lookup([0, 1, 2])
  assert_eq(hits.length(), 2)
  assert_eq(hits[0], 0)
  assert_eq(hits[1], 4)
}

///|
test "biotite_kmer_table_short_sequence" {
  let alphabet = @src.BiotiteKmerAlphabet::new(4, 3)
  let table = @src.BiotiteKmerTable::build(alphabet, [0, 1])
  assert_eq(table.size(), 0)
  assert_eq(table.lookup([0, 1, 2]).length(), 0)
}

///|
test "biotite_bucket_number" {
  let n = @src.biotite_bucket_number(1000, 3, 4)
  assert_true(n > 0)
}
