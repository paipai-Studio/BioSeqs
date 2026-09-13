///|
/// Tests for motifmatchr — PWM motif matching.

///|
test "motifmatchr_revcomp basic" {
  assert_eq(@src.motifmatchr_revcomp("ATGC"), "GCAT")
  assert_eq(@src.motifmatchr_revcomp("ACGT"), "ACGT")
  assert_eq(@src.motifmatchr_revcomp("AAAA"), "TTTT")
  assert_eq(@src.motifmatchr_revcomp(""), "")
}

///|
test "motifmatchr_revcomp mixed case passthrough" {
  assert_eq(@src.motifmatchr_revcomp("ATN"), "NAT")
}

///|
test "motifmatchr_pwm_from_counts basic" {
  let pfm = [
    [10, 0, 0, 0],
    [0, 10, 0, 0],
    [0, 0, 10, 0],
    [0, 0, 0, 10],
  ]
  let pwm = @src.motif_pwm_from_counts(pfm)
  assert_eq(pwm.length(), 4)
  assert_eq(pwm[0].length(), 4)
  let col_sum = pwm[0][0] + pwm[1][0] + pwm[2][0] + pwm[3][0]
  assert_true((col_sum - 1.0).abs() < 1.0e-9)
}

///|
test "motifmatchr_pssm_from_pwm high info motif" {
  let pfm = [
    [10, 0, 0, 0],
    [0, 10, 0, 0],
    [0, 0, 10, 0],
    [0, 0, 0, 10],
  ]
  let pwm = @src.motif_pwm_from_counts(pfm)
  let pssm = @src.motif_pssm_from_pwm(pwm)
  assert_eq(pssm.length(), 4)
  assert_true(pssm[0][0] > 1.0)
  assert_true(pssm[1][0] < -1.0)
}

///|
test "motifmatchr_match_motifs finds ACGT" {
  let motif = @src.motifmatchr_sample_motif()
  let seqs = @src.motifmatchr_sample_sequences()
  let results = @src.match_motifs([motif], seqs, cutoff=(-5.0))
  assert_eq(results.length(), 1)
  let r = results[0]
  assert_eq(r.seq_names.length(), 4)
  assert_true(r.total_hits() > 0)
}

///|
test "motifmatchr_match_motifs hit_positions" {
  let motif = @src.motifmatchr_sample_motif()
  let seqs = ["GGGACGTGGG"]
  let results = @src.match_motifs([motif], seqs, cutoff=(-5.0))
  let positions = results[0].hit_positions()
  assert_true(positions.length() > 0)
}

///|
test "motifmatchr_match_motifs empty sequences" {
  let motif = @src.motifmatchr_sample_motif()
  let results = @src.match_motifs([motif], [], cutoff=0.0)
  assert_eq(results.length(), 1)
  assert_eq(results[0].counts.length(), 0)
}

///|
test "motifmatchr_match_motifs short sequence no hit" {
  let motif = @src.motifmatchr_sample_motif()
  let results = @src.match_motifs([motif], ["AT"], cutoff=0.0)
  assert_eq(results[0].total_hits(), 0)
}

///|
test "motifmatchr_match_motifs custom seq_names" {
  let motif = @src.motifmatchr_sample_motif()
  let results = @src.match_motifs(
    [motif],
    ["ACGTACGT"],
    seq_names=["myseq"],
    cutoff=(-5.0),
  )
  assert_eq(results[0].seq_names[0], "myseq")
}

///|
test "motifmatchr_best_score positive for matching seq" {
  let motif = @src.motifmatchr_sample_motif()
  let results = @src.match_motifs([motif], ["ACGTACGT"], cutoff=(-5.0))
  assert_true(results[0].best_score() > 0.0)
}

///|
test "motifmatchr_motif_len" {
  let motif = @src.motifmatchr_sample_motif()
  assert_eq(motif.len(), 4)
}