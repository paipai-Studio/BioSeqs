///|
/// Tests for biotite_infotheory3.mbt — substitution matrices and entropy.
test "biotite3_blosum_from_alignment" {
  let matrix = @src.biotite3_build_blosum(["AC", "AC"], ["A", "C"], 1.0)
  assert_true(@src.biotite3_subst_matrix_get_score(matrix, "A", "A") >= 0.0)
  assert_true(@src.biotite3_subst_matrix_get_score(matrix, "C", "C") >= 0.0)
  assert_true(@src.biotite3_subst_matrix_get_score(matrix, "A", "C") < 0.0)
  assert_true(@src.biotite3_subst_matrix_to_string(matrix).length() > 0)
}

///|
test "biotite3_pam_wrapper" {
  let matrix = @src.biotite3_build_pam(["AC", "AC"], ["A", "C"], 1.0)
  assert_eq(matrix.name, "PAM")
}

///|
test "biotite3_shannon_entropy" {
  let conserved = @src.biotite3_shannon_entropy(["AA", "AA"], 0, ["A", "C"])
  assert_true(@moonbit.assert_num_close(conserved, 0.0, 1.0e-9))
  let variable = @src.biotite3_shannon_entropy(["AC", "AG"], 1, ["C", "G"])
  assert_true(@moonbit.assert_num_close(variable, 1.0, 1.0e-9))
  let profile = @src.biotite3_shannon_entropy_profile(["AC", "AG"], [
    "A", "C", "G",
  ])
  assert_eq(profile.length(), 2)
}

///|
test "biotite3_mutual_and_relative_entropy" {
  let mi = @src.biotite3_mutual_information(["AC", "AC"], 0, 1, ["A", "C"])
  assert_true(@moonbit.assert_num_close(mi, 0.0, 1.0e-9))
  let re = @src.biotite3_relative_entropy(["AC", "AG"], 1, ["C", "G"], [
    0.5, 0.5,
  ])
  assert_true(@moonbit.assert_num_close(re, 0.0, 1.0e-9))
  let matrix = @src.biotite3_mutual_information_matrix(["AC", "AG"], [
    "A", "C", "G",
  ])
  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0][0], 0.0)
}
