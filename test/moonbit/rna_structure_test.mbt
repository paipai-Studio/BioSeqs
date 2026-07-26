///|
/// Tests for Bio.SeqUtils - RNA Secondary Structure Prediction

/// Helper: convert String to Array[UInt16]
fn str_to_u16_array(s : String) -> Array[UInt16] {
  let arr : Array[UInt16] = Array::new()
  for i = 0; i < s.length(); i = i + 1 {
    arr.push(s.unsafe_get(i))
  }
  arr
}

/// Helper: convert Char to UInt16
fn char_to_u16(c : Char) -> UInt16 {
  c.to_int().to_uint16()
}

///|
test "rna_predictor_new" {
  let p = @src.rna_predictor_new()
  assert_eq(p.loop_min_length, 4)
  assert_eq(p.stacking_energy, -2.0)
  assert_eq(p.loop_penalty, 0.5)
  assert_eq(p.mismatch_penalty, 1.0)
}

///|
test "rna_predictor_with_params" {
  let p = @src.rna_predictor_with_params(3, -3.0, 0.3, 0.5)
  assert_eq(p.loop_min_length, 3)
  assert_eq(p.stacking_energy, -3.0)
}

///|
test "is_complementary_au" {
  assert_true(@src.is_complementary(char_to_u16('A'), char_to_u16('U')))
  assert_true(@src.is_complementary(char_to_u16('U'), char_to_u16('A')))
}

///|
test "is_complementary_gc" {
  assert_true(@src.is_complementary(char_to_u16('G'), char_to_u16('C')))
  assert_true(@src.is_complementary(char_to_u16('C'), char_to_u16('G')))
}

///|
test "is_complementary_wobble" {
  assert_true(@src.is_complementary(char_to_u16('G'), char_to_u16('U')))
  assert_true(@src.is_complementary(char_to_u16('U'), char_to_u16('G')))
}

///|
test "is_complementary_not" {
  assert_true(!@src.is_complementary(char_to_u16('A'), char_to_u16('C')))
  assert_true(!@src.is_complementary(char_to_u16('T'), char_to_u16('G')))
}

///|
test "compute_loop_energy" {
  let e1 = @src.compute_loop_energy(4)
  let e2 = @src.compute_loop_energy(10)
  assert_true(e2 > e1)
}

///|
test "nussinov_matrix_simple" {
  let p = @src.rna_predictor_new()
  let seq = str_to_u16_array("ACGU")
  let matrix = @src.nussinov_matrix(seq, p)
  assert_eq(matrix.length(), 4)
  assert_eq(matrix[0].length(), 4)
}

///|
test "nussinov_traceback" {
  let p = @src.rna_predictor_new()
  let seq = str_to_u16_array("GCGCAAAAGCGC")
  let matrix = @src.nussinov_matrix(seq, p)
  let pairs = @src.nussinov_traceback(matrix, seq, p)
  assert_true(pairs.length() >= 0)
}

///|
test "rna_predict_simple" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_secondary_structure(p, seq)
  assert_eq(result.sequence, seq)
  assert_eq(result.structure.length(), seq.length())
}

///|
test "rna_predict_with_mfe" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_with_mfe(p, seq)
  assert_eq(result.sequence, seq)
  assert_true(result.mfe <= 0.0)
}

///|
test "rna_predict_empty" {
  let p = @src.rna_predictor_new()
  let result = @src.rna_predict_secondary_structure(p, "")
  assert_eq(result.sequence, "")
}

///|
test "rna_predict_single_base" {
  let p = @src.rna_predictor_new()
  let result = @src.rna_predict_secondary_structure(p, "A")
  assert_eq(result.structure, ".")
}

///|
test "structure_to_dot_bracket" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_secondary_structure(p, seq)
  let dot = @src.structure_to_dot_bracket(result)
  assert_eq(dot.length(), seq.length())
  for c in dot {
    let ch = c.to_int().to_uint16()
    assert_true(ch == char_to_u16('(') || ch == char_to_u16(')') || ch == char_to_u16('.'))
  }
}

///|
test "get_base_pairs" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_secondary_structure(p, seq)
  let pairs = @src.get_base_pairs(result)
  for pair in pairs {
    let (i, j) = pair
    assert_true(i >= 0 && i < seq.length())
    assert_true(j >= 0 && j < seq.length())
  }
}

///|
test "get_mfe" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_secondary_structure(p, seq)
  let mfe = @src.get_mfe(result)
  assert_true(mfe <= 0.0)
}

///|
test "identify_loop_regions" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let base_pairs = [(0, 11), (1, 10), (2, 9), (3, 8)]
  let regions = @src.identify_loop_regions(str_to_u16_array(seq), base_pairs, p)
  assert_true(regions.length() >= 0)
}

///|
test "mfe_to_dot_bracket" {
  let pairs = [(0, 5), (1, 4)]
  let dot = @src.mfe_to_dot_bracket(7, pairs)
  assert_eq(dot.length(), 7)
}

///|
test "valid_dot_bracket" {
  let p = @src.rna_predictor_new()
  let seq = "GCGCAAAAGCGC"
  let result = @src.rna_predict_secondary_structure(p, seq)
  let dot = result.structure
  let mut open_count = 0
  for c in dot {
    if c == '(' { open_count = open_count + 1 }
    if c == ')' { open_count = open_count - 1 }
    assert_true(open_count >= 0)
  }
  assert_eq(open_count, 0)
}