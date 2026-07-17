///|
test "subsmat_blosum62_basic" {
  let matrix = @src.blosum62()
  assert_eq(matrix.name(), "BLOSUM62")
  assert_eq(matrix.get_score("A", "A"), 4)
  assert_eq(matrix.get_score("R", "R"), 5)
  assert_eq(matrix.get_score("A", "R"), -1)
  assert_eq(matrix.get_score("W", "W"), 11)
  assert_eq(matrix.get_score("C", "C"), 9)
}

///|
test "subsmat_blosum62_case_insensitive" {
  let matrix = @src.blosum62()
  assert_eq(matrix.get_score_case_insensitive("a", "A"), 4)
  assert_eq(matrix.get_score_case_insensitive("r", "R"), 5)
  assert_eq(matrix.get_score_case_insensitive("W", "w"), 11)
}

///|
test "subsmat_blosum45" {
  let matrix = @src.blosum45()
  assert_eq(matrix.name(), "BLOSUM45")
  assert_eq(matrix.get_score("A", "A"), 5)
  assert_eq(matrix.get_score("R", "R"), 7)
}

///|
test "subsmat_pam250" {
  let matrix = @src.pam250()
  assert_eq(matrix.name(), "PAM250")
  assert_eq(matrix.get_score("A", "A"), 2)
  assert_eq(matrix.get_score("W", "W"), 17)
  assert_eq(matrix.get_score("C", "C"), 12)
}

///|
test "subsmat_pam30" {
  let matrix = @src.pam30()
  assert_eq(matrix.name(), "PAM30")
  assert_eq(matrix.get_score("A", "A"), 9)
  assert_eq(matrix.get_score("R", "R"), 10)
}

///|
test "subsmat_get_matrix" {
  let m62 = @src.get_matrix(@src.matrix_blosum62())
  let m45 = @src.get_matrix(@src.matrix_blosum45())
  let p250 = @src.get_matrix(@src.matrix_pam250())
  let p30 = @src.get_matrix(@src.matrix_pam30())
  assert_eq(m62.name(), "BLOSUM62")
  assert_eq(m45.name(), "BLOSUM45")
  assert_eq(p250.name(), "PAM250")
  assert_eq(p30.name(), "PAM30")
}

///|
test "subsmat_contains" {
  let matrix = @src.blosum62()
  assert_true(matrix.contains("A"))
  assert_true(matrix.contains("Z") == false)
  assert_true(matrix.contains("a"))
}

///|
test "subsmat_standard_alphabet" {
  let alphabet = @src.standard_aa_alphabet()
  assert_eq(alphabet.length(), 20)
}

///|
test "subsmat_parse_matrix" {
  let content = "A R N\n" + "A 5 -1 0\n" + "R -1 6 1\n" + "N 0 1 5\n"
  let matrix_opt = @src.parse_matrix(content)
  assert_true(matrix_opt.is_some())
  let matrix = matrix_opt.unwrap()
  assert_eq(matrix.get_score("A", "A"), 5)
  assert_eq(matrix.get_score("R", "N"), 1)
}

///|
test "subsmat_smith_waterman_blosum62" {
  let seq1 = "ACDEFG"
  let seq2 = "ACDE"
  let result = @src.smith_waterman_blosum62(seq1, seq2, -1.0)
  assert_true(result.score() > 0.0)
}

///|
test "subsmat_needleman_wunsch_blosum62" {
  let seq1 = "ACDE"
  let seq2 = "ACDF"
  let result = @src.needleman_wunsch_blosum62(seq1, seq2, -1.0)
  assert_eq(result.aligned1().length(), result.aligned2().length())
}

///|
test "subsmat_smith_waterman_pam250" {
  let seq1 = "MVK"
  let seq2 = "MVL"
  let result = @src.smith_waterman_pam250(seq1, seq2, -2.0)
  assert_true(result.score() > 0.0)
}

///|
test "subsmat_to_string" {
  let matrix = @src.blosum62()
  let s = matrix.to_string()
  assert_true(s.length() > 0)
}
