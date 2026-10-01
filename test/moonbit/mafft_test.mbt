///|
test "mafft_kmer_distance identical sequences" {
  let d = @src.mafft_kmer_distance("ACGTACGTACGT", "ACGTACGTACGT", 4)
  assert_true(d < 0.01)
}

///|
test "mafft_kmer_distance different sequences" {
  let d = @src.mafft_kmer_distance("ACGTACGTACGT", "TTTTGGGGCCCC", 4)
  assert_true(d > 0.5)
}

///|
test "mafft_kmer_distance similar sequences" {
  let d = @src.mafft_kmer_distance("ACGTACGTACGT", "ACGTACGTACCT", 4)
  assert_true(d >= 0.0 && d <= 1.0)
}

///|
test "mafft_distance_matrix basic" {
  let seqs = ["ACGTACGT", "ACGTACGT", "TTTTGGGG"]
  let m = @src.mafft_distance_matrix(seqs, 4)
  assert_eq(m.length(), 3)
  assert_true(m[0][0] < 0.01)
  assert_true(m[0][1] < 0.01)
  assert_true(m[0][2] > 0.0)
  // Symmetric
  assert_true((m[0][2] - m[2][0]).abs() < 0.001)
}

///|
test "mafft_upgma_tree basic" {
  let dist = [
    [0.0, 0.1, 0.8],
    [0.1, 0.0, 0.7],
    [0.8, 0.7, 0.0],
  ]
  let tree = @src.mafft_upgma_tree(dist)
  assert_true(!tree.is_leaf)
  assert_eq(tree.children.length(), 2)
}

///|
test "mafft_upgma_tree single sequence" {
  let dist = [[0.0]]
  let tree = @src.mafft_upgma_tree(dist)
  assert_true(tree.is_leaf)
  assert_eq(tree.seq_index, 0)
}

///|
test "mafft_align_profiles basic" {
  let matrix = @src.blosum62()
  let p1 = ["ACDE"]
  let p2 = ["ACGE"]
  let result = @src.mafft_align_profiles(p1, p2, matrix, -10.0)
  assert_eq(result.length(), 2)
  assert_eq(result[0].length(), result[1].length())
}

///|
test "mafft_align identical sequences" {
  let seqs = ["ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.mafft_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  // All should be the same length
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
  // No gaps for identical sequences
  for s in result.aligned_sequences {
    assert_false(s.contains("-"))
  }
}

///|
test "mafft_align similar sequences" {
  let seqs = [
    "ACDEFGHIK",
    "ACDEAGHIK",
    "ACDEFGHLK",
  ]
  let names = ["s1", "s2", "s3"]
  let result = @src.mafft_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "mafft_align with gaps needed" {
  let seqs = [
    "ACDEFGHIK",
    "ACFGHIK",
    "ACDEFGHIKL",
  ]
  let names = ["s1", "s2", "s3"]
  let result = @src.mafft_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
  // Should have some gaps
  let mut has_gap = false
  for s in result.aligned_sequences {
    if s.contains("-") {
      has_gap = true
    }
  }
  assert_true(has_gap)
}

///|
test "mafft_align single sequence" {
  let result = @src.mafft_align(["ACDEFG"], ["s1"])
  assert_eq(result.aligned_sequences.length(), 1)
  assert_eq(result.aligned_sequences[0], "ACDEFG")
}

///|
test "mafft_align empty" {
  let result = @src.mafft_align([], [])
  assert_eq(result.aligned_sequences.length(), 0)
}

///|
test "mafft_align with params" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.mafft_align_with_params(seqs, names, 3, -8.0, 2)
  assert_eq(result.aligned_sequences.length(), 3)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "mafft_to_fasta output" {
  let seqs = ["ACDE", "ACGE"]
  let names = ["seq1", "seq2"]
  let result = @src.mafft_align(seqs, names)
  let fasta = @src.msa_to_fasta(result)
  assert_true(fasta.contains(">seq1"))
  assert_true(fasta.contains(">seq2"))
}

///|
test "mafft_to_clustal output" {
  let seqs = ["ACDE", "ACGE"]
  let names = ["seq1", "seq2"]
  let result = @src.mafft_align(seqs, names)
  let clustal = @src.msa_to_clustal(result)
  assert_true(clustal.length() > 0)
  assert_true(clustal.contains("seq1"))
  assert_true(clustal.contains("seq2"))
}

///|
test "mafft_align preserves residues" {
  // After removing gaps, each aligned sequence should match the original
  let seqs = ["ACDEFGHIK", "ACFGHIK", "ACDEFGHIKL"]
  let names = ["s1", "s2", "s3"]
  let result = @src.mafft_align(seqs, names)
  for i = 0; i < seqs.length(); i = i + 1 {
    let aligned = result.aligned_sequences[i]
    let mut ungapped = ""
    for j = 0; j < aligned.length(); j = j + 1 {
      let c = aligned[j].unsafe_to_char()
      if c != '-' {
        ungapped += c.to_string()
      }
    }
    assert_eq(ungapped, seqs[i])
  }
}

///|
test "mafft MsaResult fields" {
  let seqs = ["ACDEFG", "ACDEAG"]
  let result = @src.mafft_align(seqs, ["a", "b"])
  assert_eq(result.method, "mafft")
  assert_eq(result.names.length(), 2)
  assert_true(result.score != 0.0)
}
