///|
/// Tests for biotite_application_tantan.mbt.
test "biotite_tantan_masks_low_complexity" {
  let masked = @src.biotite_tantan_mask("AAAAAAAAAAAAAAAAAAAA")
  assert_eq(masked, "NNNNNNNNNNNNNNNNNNNN")
}

///|
test "biotite_tantan_keeps_complex_sequence" {
  let sequence = "ACGTACGTACGTACGTACGT"
  assert_eq(@src.biotite_tantan_mask(sequence), sequence)
}

///|
test "biotite_tantan_entropy_profile" {
  let profile = @src.biotite_tantan_entropy_profile("AAAAAAAAAA")
  assert_eq(profile.length(), 10)
  for h in profile {
    assert_true(@moonbit.assert_num_close(h, 0.0, 1.0e-12))
  }
}

///|
test "biotite_tantan_command" {
  let command = @src.biotite_tantan_command("in.fa", "out.fa", window=30)
  assert_true(command.contains("tantan"))
  assert_true(command.contains("-in in.fa"))
  assert_true(command.contains("-out out.fa"))
  assert_true(command.contains("-window 30"))
}
