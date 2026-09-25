///|
/// Tests for DECIPHER module.

///|
test "decipher nt_idx" {
  assert_eq(@src.decipher_nt_idx('A'), 0)
  assert_eq(@src.decipher_nt_idx('C'), 1)
  assert_eq(@src.decipher_nt_idx('G'), 2)
  assert_eq(@src.decipher_nt_idx('T'), 3)
  assert_eq(@src.decipher_nt_idx('N'), -1)
}

///|
test "decipher idx_nt" {
  assert_eq(@src.decipher_idx_nt(0), 'A')
  assert_eq(@src.decipher_idx_nt(1), 'C')
  assert_eq(@src.decipher_idx_nt(2), 'G')
  assert_eq(@src.decipher_idx_nt(3), 'T')
}

///|
test "decipher default dna matrix" {
  let m = @src.default_dna_matrix()
  assert_eq(m.scores[0][0], 5) // A-A match
  assert_eq(m.scores[0][1], -4) // A-C mismatch
  assert_eq(m.gap_open, -12)
  assert_eq(m.gap_extend, -4)
}

///|
test "decipher align identical" {
  let aln = @src.align_seqs("ACGTACGT", "ACGTACGT")
  assert_eq(aln.aligned_a, "ACGTACGT")
  assert_eq(aln.aligned_b, "ACGTACGT")
  assert_true((aln.identity - 1.0).abs() < 0.001)
  assert_eq(aln.mismatches, 0)
  assert_eq(aln.gaps, 0)
}

///|
test "decipher align with mismatch" {
  let aln = @src.align_seqs("ACGTACGT", "ACGTACGA")
  assert_true((aln.identity - 0.875).abs() < 0.001)
  assert_eq(aln.mismatches, 1)
}

///|
test "decipher align with gap" {
  let aln = @src.align_seqs("ACGT", "ACT")
  assert_true(aln.gaps >= 1)
}

///|
test "decipher align batch" {
  let seqs = ["ACGTACGT", "ACGTACGA", "ACGTACGT"]
  let results = @src.align_seqs_batch(seqs, "ACGTACGT")
  assert_eq(results.length(), 3)
  assert_true((results[0].identity - 1.0).abs() < 0.001)
}

///|
test "decipher build profile" {
  let aligned = ["ACGT", "ACGT", "ACGA"]
  let profile = @src.build_profile(aligned)
  assert_eq(profile.length, 4)
  // Position 0: all A, so freq[0] = 1.0
  assert_true((profile.freqs[0][0] - 1.0).abs() < 0.001)
  // Position 3: 2 T + 1 A, so freq T = 2/3, freq A = 1/3
  assert_true((profile.freqs[3][3] - 0.6667).abs() < 0.01)
  assert_true((profile.freqs[3][0] - 0.3333).abs() < 0.01)
}

///|
test "decipher profile to consensus" {
  let aligned = ["ACGT", "ACGT", "ACGA"]
  let profile = @src.build_profile(aligned)
  let cons = @src.profile_to_consensus(profile)
  assert_eq(cons, "ACGT")
}

///|
test "decipher align profiles" {
  let p1 = @src.build_profile(["ACGTACGT", "ACGTACGT"])
  let p2 = @src.build_profile(["ACGTACGT", "ACGTACGA"])
  let matrix = @src.default_dna_matrix()
  let aln = @src.decipher_align_profiles(p1, p2, matrix)
  assert_true(aln.identity > 0.8)
}

///|
test "decipher distance matrix" {
  let seqs = ["ACGTACGT", "ACGTACGT", "TTTTTTTT"]
  let dist = @src.decipher_distance_matrix(seqs)
  assert_eq(dist.length(), 3)
  assert_true((dist[0][0] - 0.0).abs() < 0.001)
  assert_true((dist[0][1] - 0.0).abs() < 0.001) // identical sequences
  assert_true(dist[0][2] > 0.5) // very different sequences
}

///|
test "decipher id clusters" {
  let seqs = [
    "ACGTACGTACGTACGTACGT",
    "ACGTACGTACGTACGTACGA",
    "ACGTACGTACGTACGTACGT",
    "TTGCGCGCATATATATGCGC",
    "TTGCGCGCATATATATGCGT",
  ]
  let result = @src.id_clusters(seqs, threshold=0.9)
  assert_true(result.n_clusters >= 1)
  assert_eq(result.cluster_ids.length(), 5)
}

///|
test "decipher id clusters high threshold" {
  let seqs = ["ACGTACGT", "ACGTACGT", "TTTTTTTT"]
  let result = @src.id_clusters(seqs, threshold=0.99)
  // ACGTACGT and ACGTACGT should be in same cluster
  assert_eq(result.cluster_ids[0], result.cluster_ids[1])
  // TTTTTTTT should be in different cluster
  assert_true(result.cluster_ids[0] != result.cluster_ids[2])
}

///|
test "decipher neighbor joining" {
  let dist = [
    [0.0, 5.0, 9.0, 9.0],
    [5.0, 0.0, 10.0, 10.0],
    [9.0, 10.0, 0.0, 8.0],
    [9.0, 10.0, 8.0, 0.0],
  ]
  let names = ["A", "B", "C", "D"]
  let tree = @src.neighbor_joining(dist, names)
  assert_eq(tree.n_taxa, 4)
  assert_true(tree.root >= 0)
}

///|
test "decipher tree to newick" {
  let dist = [
    [0.0, 5.0, 9.0, 9.0],
    [5.0, 0.0, 10.0, 10.0],
    [9.0, 10.0, 0.0, 8.0],
    [9.0, 10.0, 8.0, 0.0],
  ]
  let names = ["A", "B", "C", "D"]
  let tree = @src.neighbor_joining(dist, names)
  let newick = @src.tree_to_newick(tree)
  assert_true(newick.length() > 5)
  assert_true(newick.contains("("))
  assert_true(newick.contains(";"))
}

///|
test "decipher find chimeras" {
  let parent_a = "AAAAAAATTTTTTT"
  let parent_b = "GGGGGGGCCCCCCC"
  let chimera = "AAAAAAACCCCCCC"
  let seqs = [parent_a, parent_b, chimera]
  let abundances = [100, 100, 10]
  let flags = @src.find_chimeras(seqs, abundances, threshold=0.8)
  assert_true(flags[2]) // chimera detected
  assert_true(!flags[0]) // parent_a not chimera
  assert_true(!flags[1]) // parent_b not chimera
}

///|
test "decipher consensus" {
  let aligned = [
    "ACGTACGT",
    "ACGTACGT",
    "ACGTACGA",
  ]
  let cons = @src.decipher_consensus(aligned)
  assert_eq(cons, "ACGTACGT")
}

///|
test "decipher browse seqs" {
  let seqs = ["ACGT", "TTGC"]
  let result = @src.browse_seqs(seqs, names=["seq1", "seq2"])
  assert_true(result.contains("seq1"))
  assert_true(result.contains("ACGT"))
}

///|
test "decipher database operations" {
  let db = @src.seq_database_new()
  assert_eq(@src.db_size(db), 0)
  let db2 = @src.db_add_seq(db, "ACGTACGT", description="test1")
  assert_eq(@src.db_size(db2), 1)
  let db3 = @src.db_add_seq(db2, "TTTTGGGG", description="test2")
  assert_eq(@src.db_size(db3), 2)
  let entry = @src.db_get(db3, 0)
  assert_true(entry is Some(_))
  assert_eq(entry.unwrap().description, "test1")
}

///|
test "decipher database search" {
  let db = @src.create_example_database()
  assert_eq(@src.db_size(db), 4)
  let (idx, aln) = @src.db_search(db, "ACGTACGTACGTACGTACGTACGTACGTACGT")
  assert_true(idx >= 0)
  assert_true(aln.identity > 0.9)
}

///|
test "decipher database remove" {
  let db = @src.create_example_database()
  let db2 = @src.db_remove(db, 0)
  assert_eq(@src.db_size(db2), 3)
}

///|
test "decipher gc content" {
  assert_true((@src.decipher_gc_content("ACGT") - 0.5).abs() < 0.001)
  assert_true((@src.decipher_gc_content("GGCC") - 1.0).abs() < 0.001)
  assert_true((@src.decipher_gc_content("AATT") - 0.0).abs() < 0.001)
}

///|
test "decipher sequence complexity" {
  // ACGT has max complexity (all 4 nucleotides)
  let c1 = @src.decipher_complexity("ACGTACGT")
  assert_true(c1 > 0.9)
  // AAAA has zero complexity
  let c2 = @src.decipher_complexity("AAAAAAAA")
  assert_true(c2 < 0.1)
}

///|
test "decipher example sequences" {
  let seqs = @src.create_example_sequences()
  assert_eq(seqs.length(), 6)
}

///|
test "decipher example database" {
  let db = @src.create_example_database()
  assert_eq(@src.db_size(db), 4)
  let entry = @src.db_get(db, 0)
  assert_true(entry is Some(_))
  assert_eq(entry.unwrap().description, "E_coli_16S")
}