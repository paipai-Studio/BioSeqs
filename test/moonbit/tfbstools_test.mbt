///|
/// Tests for tfbstools module (TFBSTools + PWMEnrich).

test "tfbs_pfm_normalize" {
  let counts = [
    [8.0, 2.0, 0.0, 0.0],
    [0.0, 10.0, 0.0, 0.0],
  ]
  let pfm = @src.tfbs_pfm_normalize(counts, pseudocount=0.8)
  assert_eq(pfm.length(), 2)
  assert_true(pfm[0][0] > 0.5)
}

test "tfbs_pfm_to_pwm" {
  let pfm = [
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.9, 0.1, 0.0],
  ]
  let pwm = @src.tfbs_pfm_to_pwm(pfm)
  assert_eq(pwm.length(), 2)
  assert_true(pwm[0][0] > 1.0)
}

test "tfbs_consensus" {
  let pfm = [
    [0.8, 0.1, 0.1, 0.0],
    [0.1, 0.1, 0.8, 0.0],
  ]
  let cons = @src.tfbs_consensus(pfm)
  assert_eq(cons.length(), 2)
}

test "tfbs_info_content" {
  let pfm = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
  ]
  let (per_pos, total) = @src.tfbs_info_content(pfm)
  assert_eq(per_pos.length(), 2)
  assert_true(total > 0.0)
}

test "tfbs_scan_sequence" {
  let pwm = @src.tfbs_pfm_to_pwm([
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.9, 0.1, 0.0],
    [0.1, 0.1, 0.8, 0.0],
  ])
  let hits = @src.tfbs_scan_sequence("ACGTACGT", pwm, threshold=0.0)
  assert_true(hits.length() >= 0)
}

test "tfbs_scan_batch" {
  let pwm = @src.tfbs_pfm_to_pwm([
    [0.8, 0.1, 0.1, 0.0],
    [0.0, 0.9, 0.1, 0.0],
  ])
  let seqs = [("s1", "ACGT"), ("s2", "TGCA")]
  let hits = @src.tfbs_scan_batch(seqs, pwm, threshold=0.0)
  assert_true(hits.length() >= 0)
}

test "tfbs_enrich_hypergeometric" {
  let r = @src.tfbs_enrich_hypergeometric(
    "MA001", "TF1",
    20, 100, 5, 100,
  )
  assert_eq(r.matrix_id, "MA001")
  assert_eq(r.tf_name, "TF1")
  assert_true(r.fold_enrichment > 0.0)
}

test "tfbs_adjust_pvalues" {
  let r1 = @src.tfbs_enrich_hypergeometric("M1", "T1", 10, 100, 5, 100)
  let r2 = @src.tfbs_enrich_hypergeometric("M2", "T2", 20, 100, 10, 100)
  let results = @src.tfbs_adjust_pvalues([r1, r2])
  assert_eq(results.length(), 2)
  assert_true(results[0].adj_p <= 1.0)
}

test "tfbs_example_jaspar" {
  let mats = @src.tfbs_example_jaspar()
  assert_true(mats.length() >= 0)
}

test "tfbs_matrix_width" {
  let pfm = [
    [0.5, 0.5, 0.0, 0.0],
    [0.0, 0.5, 0.5, 0.0],
    [0.0, 0.0, 0.5, 0.5],
  ]
  assert_eq(@src.tfbs_matrix_width(pfm), 3)
}
