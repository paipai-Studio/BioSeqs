///|
/// Tests for biotite_permutation.mbt — ordering helpers.
test "biotite_permutation_identity" {
  let p = @src.BiotitePermutation::identity(3)
  assert_true(p.is_identity())
  let applied = p.apply([10, 20, 30])
  assert_eq(applied.length(), 3)
  assert_eq(applied[0], 10)
  assert_eq(applied[2], 30)
}

///|
test "biotite_permutation_from_order_and_inverse" {
  let p = @src.BiotitePermutation::from_order([2, 0, 1])
  let applied = p.apply(["a", "b", "c"])
  assert_eq(applied[0], "c")
  assert_eq(applied[1], "a")
  assert_eq(applied[2], "b")
  let roundtrip = p.inverse().apply(applied)
  assert_eq(roundtrip[0], "a")
  assert_eq(roundtrip[1], "b")
  assert_eq(roundtrip[2], "c")
}

///|
test "biotite_permutation_compose_inverse_is_identity" {
  let p = @src.BiotitePermutation::from_order([2, 0, 1])
  let composed = p.compose(p.inverse())
  assert_true(composed.is_identity())
}

///|
test "biotite_permutation_transpositions" {
  let p = @src.BiotitePermutation::from_order([2, 0, 1])
  assert_eq(p.to_transpositions().length(), 2)
  let swaps = @src.BiotitePermutation::identity(3)
  assert_eq(swaps.to_transpositions().length(), 0)
}

///|
test "biotite_permutation_from_distance_matrix" {
  let p = @src.BiotitePermutation::from_distance_matrix([
    [0.0, 0.1, 0.9],
    [0.1, 0.0, 0.8],
    [0.9, 0.8, 0.0],
  ])
  assert_eq(p.apply([0, 1, 2]).length(), 3)
  assert_eq(
    @src.BiotitePermutation::from_distance_matrix([]).apply([1]).length(),
    0,
  )
}
