///|
/// Tests for motif_clusterbuster module (Bio.motifs clusterbuster format).

test "clusterbuster_parse basic" {
  let content = ">cluster1 2.5\nmotifA\n[0.7 0.1 0.1 0.1]\n[0.1 0.1 0.7 0.1]\n"
  let clusters = @src.clusterbuster_parse(content)
  assert_eq(clusters.length(), 1)
  assert_eq(clusters[0].name, "cluster1")
  assert_true((clusters[0].score - 2.5).abs() < 1.0e-9)
  assert_eq(clusters[0].motifs.length(), 1)
  assert_eq(clusters[0].motifs[0].name, "motifA")
  assert_eq(clusters[0].motifs[0].pfm.length(), 2)
  assert_true((clusters[0].motifs[0].pfm[0][0] - 0.7).abs() < 1.0e-9)
}

test "clusterbuster_parse multiple clusters and motifs" {
  let content = ">cluster1 1.0\nmotifA\n[0.7 0.1 0.1 0.1]\nmotifB\n[0.3 0.2 0.2 0.3]\n>cluster2 2.0\nmotifC\n[0.1 0.1 0.1 0.7]\n"
  let clusters = @src.clusterbuster_parse(content)
  assert_eq(clusters.length(), 2)
  assert_eq(clusters[0].motifs.length(), 2)
  assert_eq(clusters[1].motifs.length(), 1)
  assert_eq(clusters[1].name, "cluster2")
}

test "clusterbuster_parse comments and blank lines" {
  let content = "# a comment\n\n>cluster1 1.0\nmotifA\n[0.7 0.1 0.1 0.1]\n"
  let clusters = @src.clusterbuster_parse(content)
  assert_eq(clusters.length(), 1)
}

test "clusterbuster_to_string roundtrip" {
  let content = ">cluster1 2.5\nmotifA\n[0.7 0.1 0.1 0.1]\n"
  let clusters = @src.clusterbuster_parse(content)
  let output = @src.clusterbuster_to_string(clusters)
  assert_true(output.contains(">cluster1"))
  assert_true(output.contains("motifA"))
  assert_true(output.contains("["))
  // Re-parse the output
  let reparsed = @src.clusterbuster_parse(output)
  assert_eq(reparsed.length(), 1)
  assert_true((reparsed[0].score - 2.5).abs() < 1.0e-9)
}

test "clusterbuster_to_motif_pfm conversion" {
  let content = ">cluster1 1.0\nmotifA\n[0.7 0.1 0.1 0.1]\n[0.1 0.1 0.7 0.1]\n"
  let clusters = @src.clusterbuster_parse(content)
  let pfm = @src.clusterbuster_to_motif_pfm(clusters[0].motifs[0])
  assert_eq(pfm.counts.length(), 2)
  // Probability 0.7 * 100 = 70 counts
  assert_true((pfm.counts[0][0] - 70.0).abs() < 1.0e-9)
  assert_true((pfm.counts[0][1] - 10.0).abs() < 1.0e-9)
}

test "motif_pfm_to_clusterbuster conversion" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["AC", "AC", "AG"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let cb = @src.motif_pfm_to_clusterbuster("converted", pfm)
  assert_eq(cb.name, "converted")
  assert_eq(cb.pfm.length(), 2)
  // Position 0: A=3/3=1.0
  assert_true((cb.pfm[0][0] - 1.0).abs() < 1.0e-9)
}

test "clusterbuster_cluster_score" {
  let content = ">cluster1 1.0\nmotifA\n[0.7 0.1 0.1 0.1]\n[0.1 0.1 0.7 0.1]\n"
  let clusters = @src.clusterbuster_parse(content)
  let score = @src.clusterbuster_cluster_score(clusters[0])
  // max probs: 0.7 + 0.7 = 1.4
  assert_true((score - 1.4).abs() < 1.0e-9)
}

test "clusterbuster_example" {
  let clusters = @src.clusterbuster_example()
  assert_eq(clusters.length(), 1)
  assert_eq(clusters[0].name, "cluster1")
}