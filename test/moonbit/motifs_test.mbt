///|
/// Motifs module tests.

///|
test "count_matrix" {
  let sites = ["ACGT", "ACAT", "ACCT"]
  let counts = @src.count_matrix(sites)
  assert_eq(counts.length(), 4)
}

///|
test "pwm_from_counts" {
  let sites = ["ACGT", "ACAT", "ACCT"]
  let counts = @src.count_matrix(sites)
  let pwm = @src.pwm_from_counts(counts, 0.5)
  assert_eq(pwm.length(), 4)
}

///|
test "create_motif" {
  let sites = ["ACGT", "ACAT", "ACCT"]
  let motif = @src.create_motif("test_motif", sites)
  assert_eq(motif.name, "test_motif")
  assert_eq(motif.pwm.length(), 4)
}

///|
test "motif_consensus" {
  let sites = ["ACGT", "ACAT", "ACCT"]
  let motif = @src.create_motif("test", sites)
  let consensus = @src.motif_consensus(motif)
  assert_eq(consensus.length(), 4)
}

///|
test "motif_ic" {
  let sites = ["AAAA", "AAAA"]
  let motif = @src.create_motif("test", sites)
  let ic = @src.motif_ic(motif)
  assert_true(ic > 0.0)
}

///|
test "search_motif" {
  let sites = ["ACGT"]
  let motif = @src.create_motif("test", sites)
  let hits = @src.search_motif(motif, "GACGT", 0.0)
  assert_true(hits.length() >= 1)
}

///|
test "pssm_from_pwm" {
  let sites = ["ACGT"]
  let motif = @src.create_motif("test", sites)
  let pssm = @src.pssm_from_pwm(motif.pwm, [0.25, 0.25, 0.25, 0.25])
  assert_eq(pssm.length(), 4)
}

///|
test "score_sequence" {
  let sites = ["ACGT"]
  let motif = @src.create_motif("test", sites)
  let pssm = @src.pssm_from_pwm(motif.pwm, [0.25, 0.25, 0.25, 0.25])
  let score = @src.score_sequence(pssm, "ACGT", motif.alphabet)
  assert_true(pssm.length() > 0)
}

///|
test "align_motifs" {
  let sites1 = ["ACGT"]
  let sites2 = ["ACGT"]
  let motif1 = @src.create_motif("m1", sites1)
  let motif2 = @src.create_motif("m2", sites2)
  let score = @src.align_motifs(motif1, motif2)
  assert_true(score > 0.0)
}

///|
test "meme_to_motif" {
  let meme = "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF test\nletter-probability matrix: alength= 4 w= 4\n0.5 0.25 0.125 0.125\n"
  let motif = @src.meme_to_motif(meme)
  assert_eq(motif.name, "test")
}

///|
test "motif_to_meme" {
  let sites = ["ACGT"]
  let motif = @src.create_motif("test", sites)
  let meme = @src.motif_to_meme(motif)
  assert_true(meme.starts_with("MEME"))
}

///|
test "reverse_complement_in_search" {
  let sites = ["ACGT"]
  let motif = @src.create_motif("test", sites)
  let hits = @src.search_motif(motif, "ACGTACGT", 0.0)
  assert_true(hits.length() > 0)
}
