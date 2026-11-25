///|
/// Tests for biotite_align_xdrop.mbt — seed extension with X-drop.
test "biotite_xdrop_exact_seed" {
  let matrix : Array[Array[Double]] = [[2.0, -1.0], [-1.0, 2.0]]
  let result = @src.biotite_align_xdrop(
    [0, 1, 0, 1],
    [0, 1, 0, 1],
    0,
    0,
    2,
    matrix,
  )
  assert_true(result.score > 0.0)
  assert_eq(result.start1, 0)
  assert_eq(result.start2, 0)
  assert_true(result.alignment.length() >= 2)
}

///|
test "biotite_xdrop_mismatch_still_aligns_seed" {
  let matrix : Array[Array[Double]] = [[2.0, -1.0], [-1.0, 2.0]]
  let result = @src.biotite_align_xdrop([0, 0, 1], [0, 0, 0], 0, 0, 2, matrix)
  assert_eq(result.alignment.length() >= 2, true)
}

///|
test "biotite_xdrop_small_threshold_truncates" {
  let matrix : Array[Array[Double]] = [[2.0, -1.0], [-1.0, 2.0]]
  let loose = @src.biotite_align_xdrop(
    [0, 1, 0, 1],
    [0, 1, 0, 1],
    0,
    0,
    1,
    matrix,
    x_drop=100.0,
  )
  let tight = @src.biotite_align_xdrop(
    [0, 1, 0, 1],
    [0, 1, 0, 1],
    0,
    0,
    1,
    matrix,
    x_drop=0.5,
  )
  assert_true(tight.alignment.length() <= loose.alignment.length())
}

///|
test "biotite_xdrop_single_element" {
  let matrix : Array[Array[Double]] = [[1.0]]
  let result = @src.biotite_align_xdrop([0], [0], 0, 0, 1, matrix)
  assert_eq(result.alignment.length(), 1)
}
