///|
/// JASPAR PFM format parsing module tests.

///|
test "parse_jaspar_empty" {
  let result = @src.parse_jaspar("")
  assert_eq(result.length(), 0)
}

///|
test "parse_jaspar_single_matrix" {
  let content = ">MA0001.1\tFOXA1\n" +
    "    A  [    10    15    20     5    10 ]\n" +
    "    C  [    20    10     5    15    20 ]\n" +
    "    G  [    30    25    30    35    25 ]\n" +
    "    T  [    40    50    45    45    45 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  assert_eq(result[0].id, "MA0001.1")
  assert_eq(result[0].name, "FOXA1")
  assert_eq(result[0].width, 5)
  assert_eq(result[0].pfm.length(), 4) // 4 rows for A, C, G, T
}

///|
test "parse_jaspar_multiple_matrices" {
  let content = ">MA0001.1\tFOXA1\n" +
    "    A  [    10    15 ]\n" +
    "    C  [    20    10 ]\n" +
    "    G  [    30    25 ]\n" +
    "    T  [    40    50 ]\n" +
    "    >MA0002.1\tCTCF\n" +
    "    A  [    20    30 ]\n" +
    "    C  [    10    20 ]\n" +
    "    G  [    40    35 ]\n" +
    "    T  [    30    15 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 2)
  assert_eq(result[0].name, "FOXA1")
  assert_eq(result[1].name, "CTCF")
}

///|
test "jaspar_matrix_pfm_values" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10    20 ]\n" +
    "    C  [    30    40 ]\n" +
    "    G  [    50    60 ]\n" +
    "    T  [    70    80 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  assert_eq(matrix.pfm[0][0], 10) // A, pos 0
  assert_eq(matrix.pfm[3][1], 80) // T, pos 1
}

///|
test "jaspar_matrix_consensus" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10     5     1 ]\n" +
    "    C  [     5    30     2 ]\n" +
    "    G  [     1     2    40 ]\n" +
    "    T  [     2     1     3 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  let consensus = matrix.consensus()
  // Each position should have the base with highest count
  assert_true(consensus.length() > 0)
}

///|
test "jaspar_matrix_pwm" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10    10 ]\n" +
    "    C  [    10    10 ]\n" +
    "    G  [    10    10 ]\n" +
    "    T  [    10    10 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  let pwm = matrix.to_pwm()
  assert_eq(pwm.length(), 4)
  assert_eq(pwm[0].length(), 2)
}

///|
test "jaspar_matrix_information_content" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    40     0     0     0 ]\n" +
    "    C  [     0    40     0     0 ]\n" +
    "    G  [     0     0    40     0 ]\n" +
    "    T  [     0     0     0    40 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  let ic = matrix.total_information()
  assert_true(ic > 3.0) // Should be high for perfect consensus
}

///|
test "jaspar_matrix_total_counts" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10    20 ]\n" +
    "    C  [    30    40 ]\n" +
    "    G  [    50    60 ]\n" +
    "    T  [    70    80 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  assert_eq(matrix.total_counts, 160) // 10+30+50+70 = 160
}

///|
test "jaspar_alphabet_detection" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10     5 ]\n" +
    "    C  [    20    15 ]\n" +
    "    G  [    30    25 ]\n" +
    "    T  [    40    55 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  assert_eq(matrix.alphabet.length(), 4)
}

///|
test "jaspar_matrix_scan_sequence" {
  let content = ">MA0001.1\tTEST\n" +
    "    A  [    10     0     0 ]\n" +
    "    C  [     0    10     0 ]\n" +
    "    G  [     0     0    10 ]\n" +
    "    T  [     0     0     0 ]\n"
  let result = @src.parse_jaspar(content)
  assert_eq(result.length(), 1)
  let matrix = result[0]
  let hits = matrix.scan_sequence("ACGTACGT")
  assert_true(hits.length() > 0)
}
