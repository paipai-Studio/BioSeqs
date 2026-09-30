///|
test "vienna partition function basic" {
  // Short sequence with one possible pair
  let q = @src.vienna_partition_function("ACGU", 310.15, 1)
  // Q should be 4×4 matrix
  assert_eq(q.length(), 4)
  // Q[0][3] = 1 (all unpaired) + exp(-E(A,U)/kT) (A-U pair)
  // A-U energy = -1, kT = 0.00198717 * 310.15 ≈ 0.616
  // boltzmann = exp(1/0.616) ≈ 5.07
  // Q[0][3] ≈ 1 + 5.07 ≈ 6.07
  assert_true(q[0][3] > 5.0)
  assert_true(q[0][3] < 10.0)
}

///|
test "vienna partition function no pairs" {
  // Sequence with no possible pairs
  let q = @src.vienna_partition_function("AAAA", 310.15, 1)
  // Q[0][3] = 1 (all unpaired, no pairs possible)
  assert_true((q[0][3] - 1.0).abs() < 1.0e-9)
}

///|
test "vienna bp probability basic" {
  let prob = @src.vienna_bp_probability("GCAU", 310.15, 1)
  assert_eq(prob.length(), 4)
  // G-C and A-U can pair
  // P[0][3] should be > 0 (G-C pair)
  assert_true(prob[0][3] > 0.0)
  // P[i][j] should equal P[j][i]
  assert_true((prob[0][3] - prob[3][0]).abs() < 1.0e-9)
}

///|
test "vienna bp probability no pairs" {
  let prob = @src.vienna_bp_probability("AAAA", 310.15, 1)
  // No pairs possible → all probabilities 0
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(prob[i][j] < 1.0e-15)
    }
  }
}

///|
test "vienna mea basic" {
  // Create a simple bp probability matrix
  // 4 bases, strong P(0,3) and P(1,2)
  let prob : Array[Array[Double]] = Array::makei(4, fn(_) { Array::make(4, 0.0) })
  prob[0][3] = 0.9
  prob[3][0] = 0.9
  prob[1][2] = 0.8
  prob[2][1] = 0.8
  let structure = @src.vienna_mea(prob, 2.0)
  assert_eq(structure.length(), 4)
  // Should produce "((..))" or similar with 2 pairs
  let arr = structure.to_array()
  let mut n_open = 0
  let mut n_close = 0
  for c in arr {
    if c == '(' { n_open = n_open + 1 }
    if c == ')' { n_close = n_close + 1 }
  }
  assert_eq(n_open, 2)
  assert_eq(n_close, 2)
}

///|
test "vienna mea empty" {
  let structure = @src.vienna_mea([], 2.0)
  assert_eq(structure, "")
}

///|
test "vienna mountain plot" {
  // "((..))" → heights: 1, 2, 2, 2, 2, 1
  let plot = @src.vienna_mountain_plot("((..))")
  assert_eq(plot.length(), 6)
  assert_eq(plot[0], (1, 1))
  assert_eq(plot[1], (2, 2))
  assert_eq(plot[2], (3, 2))
  assert_eq(plot[3], (4, 2))
  assert_eq(plot[4], (5, 2))
  assert_eq(plot[5], (6, 1))
}

///|
test "vienna mountain plot unpaired" {
  let plot = @src.vienna_mountain_plot("....")
  assert_eq(plot.length(), 4)
  for i = 0; i < 4; i = i + 1 {
    assert_eq(plot[i].1, 0)
  }
}

///|
test "vienna bp distance identical" {
  let d = @src.vienna_bp_distance("((..))", "((..))")
  assert_eq(d, 0)
}

///|
test "vienna bp distance different" {
  // S1 = "((..))" has pairs (0,5) and (1,4)
  // S2 = ".(..)." has pair (1,4)
  // Distance = 1 (pair (0,5) in S1 but not S2)
  let d = @src.vienna_bp_distance("((..))", ".(..).")
  assert_eq(d, 1)
}

///|
test "vienna bp distance no pairs" {
  let d = @src.vienna_bp_distance("......", "......")
  assert_eq(d, 0)
}

///|
test "vienna loop types hairpin" {
  // "(...)" → 1 hairpin, 0 others
  let (h, i, b, m, e) = @src.vienna_loop_types("(...)")
  assert_eq(h, 1)
  assert_eq(i, 0)
  assert_eq(b, 0)
  assert_eq(m, 0)
  assert_eq(e, 1)
}

///|
test "vienna loop types nested" {
  // "((..))" → 1 hairpin (inner), 1 interior (outer)
  let (h, i, b, m, e) = @src.vienna_loop_types("((..))")
  assert_eq(h, 1)
  assert_eq(i, 1)
  assert_eq(b, 0)
  assert_eq(m, 0)
  assert_eq(e, 1)
}

///|
test "vienna loop types external" {
  // "(..).." → 1 hairpin in external loop, 1 external pair
  let (h, _i, _b, _m, e) = @src.vienna_loop_types("(..)..")
  assert_eq(h, 1)
  assert_eq(e, 1)
}

///|
test "vienna structure energy" {
  // "GCAU" with "(..)" → G-C pair = -2.0, hairpin loop = 1 unpaired * 0.5 = 0.5
  // Total = -2.0 + 0.5 = -1.5
  let e = @src.vienna_structure_energy("GCAU", "(..)")
  // G-U wobble pair = -0.5, hairpin loop (2 unpaired) = 2 * 0.5 = 1.0
  // Total = -0.5 + 1.0 = 0.5
  assert_true((e - 0.5).abs() < 1.0e-9)
}

///|
test "vienna structure energy no pairs" {
  let e = @src.vienna_structure_energy("GCAU", "....")
  assert_true(e.abs() < 1.0e-9)
}

///|
test "vienna tree representation" {
  // "((..))" → "((UU))"
  let t = @src.vienna_tree_representation("((..))")
  assert_eq(t, "((UU))")
}

///|
test "vienna tree representation unpaired" {
  let t = @src.vienna_tree_representation("....")
  assert_eq(t, "UUUU")
}

///|
test "vienna tree representation mixed" {
  // "(..)." → "(UU)U"
  let t = @src.vienna_tree_representation("(..).")
  assert_eq(t, "(UU)U")
}

