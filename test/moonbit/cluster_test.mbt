///|
test "cluster_seq_distance" {
  let dist = @src.seq_distance("ACGT", "ACAT")

  assert_true(dist >= 0.0 && dist <= 1.0)
}

///|
test "cluster_distance_matrix" {
  let seqs = ["ACGT", "ACAT", "AGGT"]
  let matrix = @src.distance_matrix(seqs)

  assert_eq(matrix.length(), 3)
}

///|
test "cluster_single_linkage" {
  let seqs = ["ACGT", "ACAT", "AGGT", "TTTT"]
  let distances = @src.distance_matrix(seqs)
  let result = @src.single_linkage_clustering(distances, 0.1)

  assert_true(@src.total_clusters(result) >= 1)
}

///|
test "cluster_complete_linkage" {
  let seqs = ["ACGT", "ACAT", "AGGT", "TTTT"]
  let distances = @src.distance_matrix(seqs)
  let result = @src.complete_linkage_clustering(distances, 0.1)

  assert_true(@src.total_clusters(result) >= 1)
}

///|
test "cluster_average_linkage" {
  let seqs = ["ACGT", "ACAT", "AGGT", "TTTT"]
  let distances = @src.distance_matrix(seqs)
  let result = @src.average_linkage_clustering(distances, 0.1)

  assert_true(@src.total_clusters(result) >= 1)
}

///|
test "cluster_hclust_to_newick" {
  let seqs = ["ACGT", "ACAT", "AGGT"]
  let distances = @src.distance_matrix(seqs)
  let result = @src.single_linkage_clustering(distances, 0.1)
  let newick = @src.hclust_to_newick(result, seqs)

  assert_true(newick.length() > 0)
}

///|
test "cluster_silhouette_score" {
  let seqs = ["ACGT", "ACAT", "AGGT", "TTTT"]
  let labels = [0, 0, 0, 1]
  let score = @src.seq_silhouette_score(seqs, labels)

  assert_true(score >= -1.0 && score <= 1.0)
}

///|
test "cluster_total_clusters" {
  let seqs = ["ACGT", "ACAT", "AGGT", "TTTT"]
  let distances = @src.distance_matrix(seqs)
  let result = @src.single_linkage_clustering(distances, 0.1)

  assert_eq(@src.total_clusters(result), 4)
}
