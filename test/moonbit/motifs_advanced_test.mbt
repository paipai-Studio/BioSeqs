///|
/// Tests for motifs_advanced module.

// ============================================================================
// JASPAR format tests
// ============================================================================

///|
test "jaspar_motif_creation" {
  let pwm : Array[Array[Double]] = [
    [0.1, 0.2, 0.3, 0.4],
    [0.4, 0.3, 0.2, 0.1],
  ]
  let motif = @src.JasparMotif::new(
    "MA0001.1",
    "AGL3",
    "CORE",
    "MADS-box",
    "AGL3",
    "9606",
    "vertebrates",
    ["P12345"],
    pwm,
    ["A", "C", "G", "T"],
  )
  assert_eq(motif.matrix_id, "MA0001.1")
  assert_eq(motif.name, "AGL3")
  assert_eq(motif.collection, "CORE")
  assert_eq(motif.pwm.length(), 2)
}

///|
test "jaspar_parse_simple_pfm" {
  let pfm_content =
    ">MA0001.1 AGL3\n" +
    "A  [ 0  3 79 40 66 48 65 11 65  0 ]\n" +
    "C  [ 94 75  4  3  1  2  5  2  7  0 ]\n" +
    "G  [  1  0  3  4  1  0  5  3  4  0 ]\n" +
    "T  [  2 20 11 50 30 48 23 82 22 99 ]\n"
  let motifs = @src.parse_jaspar_pfm(pfm_content)
  assert_eq(motifs.length(), 1)
  assert_eq(motifs[0].matrix_id, "MA0001.1")
  assert_eq(motifs[0].name, "AGL3")
  assert_eq(motifs[0].pwm.length(), 10)
}

///|
test "jaspar_parse_multiple_motifs" {
  let pfm_content =
    ">MA0001.1 AGL3\n" +
    "A [ 10 20 30 40 ]\n" +
    "C [ 40 30 20 10 ]\n" +
    "G [ 20 20 20 20 ]\n" +
    "T [ 30 30 30 30 ]\n" +
    ">MA0002.1 TF2\n" +
    "A [ 5 5 5 5 ]\n" +
    "C [ 5 5 5 5 ]\n" +
    "G [ 5 5 5 5 ]\n" +
    "T [ 5 5 5 5 ]\n"
  let motifs = @src.parse_jaspar_pfm(pfm_content)
  assert_eq(motifs.length(), 2)
  assert_eq(motifs[0].name, "AGL3")
  assert_eq(motifs[1].name, "TF2")
}

///|
test "jaspar_pfm_to_pwm" {
  let pfm : Array[Array[Double]] = [
    [10.0, 20.0, 30.0, 40.0],
    [40.0, 30.0, 20.0, 10.0],
  ]
  let pwm = @src.jaspar_pfm_to_pwm(pfm, 0.5)
  assert_eq(pwm.length(), 2)
  assert_eq(pwm[0].length(), 4)
  // Check that each row sums to approximately 1
  let mut sum0 = 0.0
  for i = 0; i < 4; i = i + 1 {
    sum0 = sum0 + pwm[0][i]
  }
  assert_true(sum0 > 0.99 && sum0 < 1.01)
}

///|
test "jaspar_pwm_to_jaspar_format" {
  let pwm : Array[Array[Double]] = [
    [0.1, 0.2, 0.3, 0.4],
    [0.4, 0.3, 0.2, 0.1],
  ]
  let jaspar_str = @src.pwm_to_jaspar(pwm, "MA0001.1", "TestMotif")
  assert_true(jaspar_str.has_prefix(">MA0001.1 TestMotif"))
  assert_true(jaspar_str.contains("A ["))
  assert_true(jaspar_str.contains("T ["))
}

// ============================================================================
// TRANSFAC format tests
// ============================================================================

///|
test "transfac_motif_creation" {
  let pwm : Array[Array[Double]] = [
    [1.0, 2.0, 3.0, 4.0],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let motif = @src.TransfacMotif::new(
    "T00001",
    "M00001",
    "AP-1",
    "Activator protein 1",
    pwm,
    ["A", "C", "G", "T"],
    [],
  )
  assert_eq(motif.accession, "T00001")
  assert_eq(motif.id, "M00001")
  assert_eq(motif.name, "AP-1")
  assert_eq(motif.description, "Activator protein 1")
}

///|
test "transfac_parse_simple" {
  let transfac_content =
    "AC T00001\n" +
    "ID M00001\n" +
    "NA AP-1\n" +
    "DE Activator protein 1\n" +
    "P0      A      C      G      T\n" +
    "01      1      2      3      4\n" +
    "02      4      3      2      1\n" +
    "XX\n" +
    "//\n"
  let motifs = @src.parse_transfac(transfac_content)
  assert_eq(motifs.length(), 1)
  assert_eq(motifs[0].accession, "T00001")
  assert_eq(motifs[0].name, "AP-1")
  assert_eq(motifs[0].pwm.length(), 2)
}

///|
test "transfac_pwm_to_transfac_format" {
  let pwm : Array[Array[Double]] = [
    [1.0, 2.0, 3.0, 4.0],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let transfac_str = @src.pwm_to_transfac(pwm, "TestMotif", "T00001")
  assert_true(transfac_str.contains("AC T00001"))
  assert_true(transfac_str.contains("NA TestMotif"))
  assert_true(transfac_str.contains("P0"))
  assert_true(transfac_str.contains("//"))
}

// ============================================================================
// Advanced motif comparison tests
// ============================================================================

///|
test "motif_alignment_result_creation" {
  let result = @src.MotifAlignmentResult::new(0.85, 2, "+", 8)
  assert_eq(result.score, 0.85)
  assert_eq(result.offset, 2)
  assert_eq(result.strand, "+")
  assert_eq(result.aligned_length, 8)
}

///|
test "optimal_motif_alignment_identical" {
  let pwm1 : Array[Array[Double]] = [
    [0.25, 0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let pwm2 : Array[Array[Double]] = [
    [0.25, 0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let alignment = @src.optimal_motif_alignment(pwm1, pwm2)
  assert_true(alignment.score > 0.0)
  assert_eq(alignment.aligned_length, 2)
}

///|
test "optimal_motif_alignment_different_lengths" {
  let pwm1 : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.1, 0.1, 0.8],
    [0.8, 0.1, 0.1, 0.0],
  ]
  let pwm2 : Array[Array[Double]] = [
    [0.0, 0.1, 0.1, 0.8],
    [0.8, 0.1, 0.1, 0.0],
  ]
  let alignment = @src.optimal_motif_alignment(pwm1, pwm2)
  assert_true(alignment.score >= 0.0)
  assert_true(alignment.aligned_length > 0)
  assert_true(alignment.offset >= -1 && alignment.offset <= 2)
}

///|
test "motif_kl_divergence" {
  let pwm1 : Array[Array[Double]] = [
    [0.25, 0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let pwm2 : Array[Array[Double]] = [
    [0.25, 0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let kl = @src.motif_kl_divergence(pwm1, pwm2)
  assert_eq(kl, 0.0)
}

///|
test "motif_js_divergence" {
  let pwm1 : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.1, 0.1, 0.8],
  ]
  let pwm2 : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.1, 0.1, 0.8],
  ]
  let js = @src.motif_js_divergence(pwm1, pwm2)
  assert_eq(js, 0.0)
}

///|
test "motif_js_divergence_different" {
  let pwm1 : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
  ]
  let pwm2 : Array[Array[Double]] = [
    [0.0, 0.1, 0.1, 0.8],
  ]
  let js = @src.motif_js_divergence(pwm1, pwm2)
  assert_true(js >= 0.0)
}

// ============================================================================
// Motif clustering tests
// ============================================================================

///|
test "motif_cluster_creation" {
  let cluster = @src.MotifCluster::new(0, ["motif1", "motif2"], "motif1", 0.85)
  assert_eq(cluster.cluster_id, 0)
  assert_eq(cluster.motifs.length(), 2)
  assert_eq(cluster.representative, "motif1")
  assert_eq(cluster.avg_similarity, 0.85)
}

///|
test "cluster_motifs_single" {
  let names = ["motif1"]
  let pwm1 : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.1, 0.1, 0.8],
  ]
  let pwms = [pwm1]
  let clusters = @src.cluster_motifs(names, pwms, 0.5)
  assert_eq(clusters.length(), 1)
  assert_eq(clusters[0].motifs.length(), 1)
}

// ============================================================================
// Pseudocount and utility tests
// ============================================================================

///|
test "optimal_pseudocount" {
  let pc1 = @src.optimal_pseudocount(10, 8)
  let pc2 = @src.optimal_pseudocount(100, 8)
  assert_true(pc1 > 0.0)
  assert_true(pc2 > 0.0)
  assert_true(pc1 > pc2) // More sequences = lower pseudocount
}

///|
test "motif_gc_content" {
  let pwm : Array[Array[Double]] = [
    [0.1, 0.4, 0.4, 0.1],
    [0.1, 0.4, 0.4, 0.1],
  ]
  let gc = @src.motif_gc_content_adv(pwm)
  assert_true(gc > 0.7 && gc < 0.9)
}

///|
test "shuffle_motif_columns" {
  let pwm : Array[Array[Double]] = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.1, 0.1, 0.8],
    [0.8, 0.1, 0.1, 0.0],
  ]
  let shuffled = @src.shuffle_motif_columns(pwm)
  assert_eq(shuffled.length(), 3)
  assert_eq(shuffled[0].length(), 4)
}

///|
test "randomize_motif" {
  let randomized = @src.randomize_motif(5, 4)
  assert_eq(randomized.length(), 5)
  assert_eq(randomized[0].length(), 4)
  // Check that each row sums to ~1
  let mut sum = 0.0
  for i = 0; i < 4; i = i + 1 {
    sum = sum + randomized[0][i]
  }
  assert_true(sum > 0.99 && sum < 1.01)
}

// ============================================================================
// Empty input tests
// ============================================================================

///|
test "jaspar_parse_empty" {
  let motifs = @src.parse_jaspar_pfm("")
  assert_eq(motifs.length(), 0)
}

///|
test "transfac_parse_empty" {
  let motifs = @src.parse_transfac("")
  assert_eq(motifs.length(), 0)
}

///|
test "optimal_alignment_empty" {
  let empty : Array[Array[Double]] = []
  let result = @src.optimal_motif_alignment(empty, empty)
  assert_eq(result.score, 0.0)
  assert_eq(result.aligned_length, 0)
}

///|
test "kl_divergence_empty" {
  let empty : Array[Array[Double]] = []
  let kl = @src.motif_kl_divergence(empty, empty)
  assert_eq(kl, 0.0)
}

///|
test "js_divergence_empty" {
  let empty : Array[Array[Double]] = []
  let js = @src.motif_js_divergence(empty, empty)
  assert_eq(js, 0.0)
}

///|
test "pwm_to_jaspar_empty" {
  let empty : Array[Array[Double]] = []
  let result = @src.pwm_to_jaspar(empty, "test", "test")
  assert_eq(result, "")
}

///|
test "pwm_to_transfac_empty" {
  let empty : Array[Array[Double]] = []
  let result = @src.pwm_to_transfac(empty, "test", "test")
  assert_eq(result, "")
}

///|
test "cluster_motifs_empty" {
  let empty_names : Array[String] = []
  let empty_pwms : Array[Array[Array[Double]]] = []
  let clusters = @src.cluster_motifs(empty_names, empty_pwms, 0.5)
  assert_eq(clusters.length(), 0)
}

///|
test "motif_gc_content_empty" {
  let empty : Array[Array[Double]] = []
  let gc = @src.motif_gc_content_adv(empty)
  assert_eq(gc, 0.0)
}

///|
test "shuffle_motif_empty" {
  let empty : Array[Array[Double]] = []
  let result = @src.shuffle_motif_columns(empty)
  assert_eq(result.length(), 0)
}

///|
test "pfm_to_pwm_empty" {
  let empty : Array[Array[Double]] = []
  let pwm = @src.jaspar_pfm_to_pwm(empty, 0.5)
  assert_eq(pwm.length(), 0)
}
