test "kmer_build_dna" {
  let t = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  assert_eq(@src.kmer_k(t), 3)
  // ACGTACGT has 6 3-mers: ACG, CGT, GTA, TAC, ACG, CGT
  assert_eq(@src.kmer_total_count(t), 6)
  // Unique: ACG, CGT, GTA, TAC = 4
  assert_eq(@src.kmer_unique_count(t), 4)
}

test "kmer_build_protein" {
  let t = @src.kmer_table_build("MVLSPADK", 2, "protein")
  assert_eq(@src.kmer_k(t), 2)
  // MVLSPADK has 7 2-mers, all unique
  assert_eq(@src.kmer_total_count(t), 7)
  assert_eq(@src.kmer_unique_count(t), 7)
}

test "kmer_jaccard_identical" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  assert_true(@moonbit.assert_close(@src.kmer_jaccard(a, b), 1.0, 1.0e-9))
}

test "kmer_jaccard_different" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("TTTTTTTT", 3, "DNA")
  // No shared 3-mers
  assert_true(@moonbit.assert_close(@src.kmer_jaccard(a, b), 0.0, 1.0e-9))
}

test "kmer_jaccard_partial" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("ACGTAAAA", 3, "DNA")
  // Some shared 3-mers, Jaccard should be between 0 and 1
  let j = @src.kmer_jaccard(a, b)
  assert_true(j > 0.0)
  assert_true(j < 1.0)
}

test "kmer_cosine_identical" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  assert_true(@moonbit.assert_close(@src.kmer_cosine(a, b), 1.0, 1.0e-9))
}

test "kmer_cosine_orthogonal" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("TTTTTTTT", 3, "DNA")
  assert_true(@moonbit.assert_close(@src.kmer_cosine(a, b), 0.0, 1.0e-9))
}

test "kmer_bray_curtis_identical" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  assert_true(@moonbit.assert_close(@src.kmer_bray_curtis(a, b), 0.0, 1.0e-9))
}

test "kmer_bray_curtis_different" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("TTTTTTTT", 3, "DNA")
  assert_true(@moonbit.assert_close(@src.kmer_bray_curtis(a, b), 1.0, 1.0e-9))
}

test "kmer_mash_distance" {
  let a = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  let b = @src.kmer_table_build("ACGTACGT", 3, "DNA")
  // Identical sequences → distance 0
  assert_true(@moonbit.assert_close(@src.kmer_mash_distance(a, b), 0.0, 1.0e-9))
  let c = @src.kmer_table_build("TTTTTTTT", 3, "DNA")
  // Completely different → distance 1
  assert_true(@moonbit.assert_close(@src.kmer_mash_distance(a, c), 1.0, 1.0e-9))
}

test "kmer_similarity_matrix" {
  let seqs = ["ACGTACGT", "ACGTACGA", "TTTTTTTT"]
  let m = @src.kmer_similarity_matrix(seqs, 3, "DNA", "jaccard")
  assert_eq(m.length(), 3)
  assert_eq(m[0].length(), 3)
  // Diagonal should be 1.0
  assert_true(@moonbit.assert_close(m[0][0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_close(m[1][1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_close(m[2][2], 1.0, 1.0e-9))
  // Symmetric
  assert_true(@moonbit.assert_close(m[0][1], m[1][0], 1.0e-9))
  assert_true(@moonbit.assert_close(m[0][2], m[2][0], 1.0e-9))
}

test "kmer_skip_invalid" {
  // N should be skipped in DNA mode
  let t = @src.kmer_table_build("ACGTNACG", 3, "DNA")
  // Valid 3-mers: ACG, CGT, ACG (skipping ones with N)
  assert_eq(@src.kmer_total_count(t), 3)
}
