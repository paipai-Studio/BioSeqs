///|
/// Tests for Bio.Align.AlignClusterer - Progressive Multiple Sequence Alignment

///|
test "align_cluster_new" {
  let c = @src.align_cluster_new()
  assert_eq(c.match_score, 2.0)
  assert_eq(c.mismatch_score, -1.0)
  assert_eq(c.open_gap_penalty, -2.0)
  assert_eq(c.extend_gap_penalty, -0.5)
}

///|
test "align_cluster_with_scoring" {
  let c = @src.align_cluster_with_scoring(3.0, -2.0, -3.0, -1.0)
  assert_eq(c.match_score, 3.0)
  assert_eq(c.mismatch_score, -2.0)
  assert_eq(c.open_gap_penalty, -3.0)
  assert_eq(c.extend_gap_penalty, -1.0)
}

///|
test "compute_distance_identical" {
  let c = @src.align_cluster_new()
  let d = @src.compute_distance(c, "ACGT", "ACGT")
  assert_eq(d, 0.0)
}

///|
test "compute_distance_different" {
  let c = @src.align_cluster_new()
  let d = @src.compute_distance(c, "AAAA", "CCCC")
  assert_true(d > 0.0)
}

///|
test "compute_distance_similar" {
  let c = @src.align_cluster_new()
  let d1 = @src.compute_distance(c, "ACGTACGT", "ACGTACGT")
  let d2 = @src.compute_distance(c, "ACGTACGT", "TTTTTTTT")
  assert_true(d1 < d2)
}

///|
test "build_distance_matrix" {
  let c = @src.align_cluster_new()
  let seqs = ["AAAA", "CCCC", "GGGG"]
  let ids = ["seq1", "seq2", "seq3"]
  let dm = @src.build_distance_matrix(c, seqs, ids)
  assert_eq(dm.sequence_ids.length(), 3)
  assert_eq(dm.matrix.length(), 3)
  assert_eq(dm.matrix[0].length(), 3)
  assert_eq(dm.matrix[0][0], 0.0)
  assert_eq(dm.matrix[1][1], 0.0)
  assert_eq(dm.matrix[2][2], 0.0)
}

///|
test "upgma_tree_two_sequences" {
  let c = @src.align_cluster_new()
  let seqs = ["AAAA", "CCCC"]
  let ids = ["seq1", "seq2"]
  let dm = @src.build_distance_matrix(c, seqs, ids)
  let tree = @src.upgma_tree(dm)
  assert_true(tree.nodes.length() >= 3)
}

///|
test "upgma_tree_identical" {
  let c = @src.align_cluster_new()
  let seqs = ["AAAA", "AAAA", "AAAA"]
  let ids = ["seq1", "seq2", "seq3"]
  let dm = @src.build_distance_matrix(c, seqs, ids)
  let tree = @src.upgma_tree(dm)
  assert_true(tree.nodes.length() >= 5)
}

///|
test "align_pairwise_identical" {
  let (a1, a2, score) = @src.align_pairwise("ACGT", "ACGT")
  assert_eq(a1, "ACGT")
  assert_eq(a2, "ACGT")
  assert_true(score > 0.0)
}

///|
test "align_pairwise_with_gaps" {
  let (a1, a2, _score) = @src.align_pairwise("ACGT", "ACG")
  assert_eq(a1.length(), a2.length())
}

///|
test "progressive_alignment_two_identical" {
  let c = @src.align_cluster_new()
  let seqs = ["ACGT", "ACGT"]
  let ids = ["seq1", "seq2"]
  let msa = @src.progressive_alignment(c, seqs, ids)
  assert_eq(msa.num_records(), 2)
  assert_true(msa.get_alignment_length() >= 4)
}

///|
test "progressive_alignment_single_sequence" {
  let c = @src.align_cluster_new()
  let seqs = ["ACGTACGT"]
  let ids = ["seq1"]
  let msa = @src.progressive_alignment(c, seqs, ids)
  assert_eq(msa.num_records(), 1)
  assert_eq(msa.get_alignment_length(), 8)
}

///|
test "progressive_alignment_empty" {
  let c = @src.align_cluster_new()
  let seqs : Array[String] = []
  let ids : Array[String] = []
  let msa = @src.progressive_alignment(c, seqs, ids)
  assert_eq(msa.num_records(), 0)
}

///|
test "get_aligned_sequences" {
  let c = @src.align_cluster_new()
  let seqs = ["ACGT", "ACGT"]
  let aligned = @src.get_aligned_sequences(c, seqs)
  assert_eq(aligned.length(), 2)
  assert_eq(aligned[0].length(), aligned[1].length())
}

///|
test "align_profiles" {
  let p1 = ["ACGT", "ACGT"]
  let p2 = ["ACGT", "ACGT"]
  let (a1, a2, score) = @src.align_profiles(p1, p2)
  assert_eq(a1.length(), a2.length())
  assert_true(score > 0.0)
}

///|
test "profile_from_sequences" {
  let seqs = ["ACGT", "ACGT"]
  let profile = @src.profile_from_sequences(seqs)
  assert_eq(profile.length(), 4)
  assert_eq(profile[0].length(), 27)
}

///|
test "guide_tree_leaf_count" {
  let c = @src.align_cluster_new()
  let seqs = ["AAAA", "CCCC", "GGGG", "TTTT"]
  let ids = ["s1", "s2", "s3", "s4"]
  let dm = @src.build_distance_matrix(c, seqs, ids)
  let tree = @src.upgma_tree(dm)
  let mut leaf_count = 0
  for node in tree.nodes {
    if node.is_leaf {
      leaf_count = leaf_count + 1
    }
  }
  assert_eq(leaf_count, 4)
}