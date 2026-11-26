///|
/// Tests for biotite_application_rnaalifold.mbt.
test "rnaalifold command" {
  let cmd = @src.biotite_rnaalifold_command("aln.clustal")
  assert_eq(cmd, "RNAalifold aln.clustal")
}

///|
test "rnaalifold consensus stem loop" {
  let r = @src.biotite_rnaalifold_local(["GGGAAACCC", "GGGAAACCC"])
  assert_eq(r.status, "computed-local")
  assert_eq(r.consensus_structure.length(), 9)
  // columns 0/8, 1/7, 2/6 form G-C pairs -> outer stem closed
  assert_eq(r.consensus_structure, "(((...)))")
  assert_eq(r.conservation_scores.length(), 9)
  // fully conserved columns have conservation 1.0
  assert_true(@moonbit.assert_num_close(r.conservation_scores[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.conservation_scores[4], 1.0, 1.0e-9))
  assert_true(r.command.contains("RNAalifold"))
}

///|
test "rnaalifold AU pairs" {
  let r = @src.biotite_rnaalifold_local(["AAAATTTT", "AAAATTTT"])
  assert_eq(r.consensus_structure, "((....))")
}

///|
test "rnaalifold GU wobble pairs" {
  // the outer pairs (0,8) and (1,7) are G-U wobbles (score 1); only
  // (2,6) is a stronger G-C — all three still close the stem, and the
  // U-run interior cannot pair
  let r = @src.biotite_rnaalifold_local(["GGGUUUCUU", "GGGUUUCUU"])
  assert_eq(r.consensus_structure, "(((...)))")
}

///|
test "rnaalifold non-pairing columns" {
  // all-A columns can never pair
  let r = @src.biotite_rnaalifold_local(["AAAAAAAAA", "AAAAAAAAA"])
  assert_eq(r.consensus_structure, ".........")
  // A/C endpoints cannot pair either
  let r2 = @src.biotite_rnaalifold_local(["AAAACCCC", "AAAACCCC"])
  assert_eq(r2.consensus_structure, "........")
}

///|
test "rnaalifold gaps excluded from majority" {
  // second sequence is all gaps; consensus comes from the first
  let r = @src.biotite_rnaalifold_local(["GGGAAACCC", "---------"])
  assert_eq(r.consensus_structure.length(), 9)
  // conservation of each column is 0.5 (one of two sequences contributes)
  assert_true(@moonbit.assert_num_close(r.conservation_scores[0], 0.5, 1.0e-9))
}

///|
test "rnaalifold covariance lowers pair score" {
  // one sequence pairs G-C, the other has A-C mismatch at those columns
  let r = @src.biotite_rnaalifold_local(["GGGAAACCC", "AGGAAACCC"])
  // majority G/C but only half the sequences actually pair
  assert_false(r.consensus_structure.has_prefix("(((..)))"))
}

///|
test "rnaalifold empty alignment" {
  let r = @src.biotite_rnaalifold_local([])
  assert_eq(r.consensus_structure, "")
  assert_eq(r.conservation_scores.length(), 0)
}

///|
test "rnaalifold min loop enforced" {
  // short hairpin with 2 nt loop: min_loop=3 forbids the closing pair
  let r = @src.biotite_rnaalifold_local(["GGAAACCC", "GGAAACCC"], min_loop=3)
  // outer pair (0,7) has interior span 6; the largest possible structure
  // contains pairs with interior >= 4, which is fine, but pairs (0,3)-
  // style with tiny loops are forbidden
  assert_eq(r.consensus_structure, "((....))")
}

///|
test "rnaalifold command custom bin" {
  let cmd = @src.biotite_rnaalifold_command(
    "aln.clustal",
    bin_path="/usr/local/bin/RNAalifold",
  )
  assert_eq(cmd, "/usr/local/bin/RNAalifold aln.clustal")
}
