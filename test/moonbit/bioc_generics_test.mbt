///|
/// Tests for BiocGenerics module.

test "na_int" {
  assert_eq(@src.na_int(), -2147483648)
}

test "na_double" {
  let na = @src.na_double()
  assert_true(na.is_nan())
}

test "is_na_int" {
  assert_true(@src.is_na_int(-2147483648))
  assert_false(@src.is_na_int(0))
}

test "is_na_double" {
  assert_true(@src.is_na_double(@double.not_a_number))
  assert_false(@src.is_na_double(0.0))
}

test "any_na_int" {
  assert_true(@src.any_na_int([1, 2, -2147483648, 4]))
  assert_false(@src.any_na_int([1, 2, 3, 4]))
}

test "na_omit_int" {
  let result = @src.na_omit_int([1, -2147483648, 3])
  assert_eq(result.length(), 2)
}

test "which_int" {
  let result = @src.which_int([1, 2, 1, 3], 1)
  assert_eq(result.length(), 2)
}

test "which_bool" {
  let result = @src.which_bool([true, false, true])
  assert_eq(result.length(), 2)
}

test "order_int" {
  let indices = @src.order_int([3, 1, 2])
  assert_eq(indices[0], 1)
}

test "sort_int" {
  let result = @src.sort_int([3, 1, 2])
  assert_eq(result[0], 1)
}

test "unique_int" {
  let result = @src.unique_int([1, 2, 1, 3])
  assert_eq(result.length(), 3)
}

test "unique_string" {
  let result = @src.unique_string(["a", "b", "a", "c"])
  assert_eq(result.length(), 3)
}

test "match_int" {
  let result = @src.match_int([1, 3], [1, 2, 4])
  assert_eq(result[0], 0)
}

test "match_string" {
  let result = @src.match_string(["a", "c"], ["a", "b"])
  assert_eq(result[0], 0)
}

test "intersect_int" {
  let result = @src.intersect_int([1, 2, 3], [2, 3, 4])
  assert_eq(result.length(), 2)
}

test "intersect_string" {
  let result = @src.intersect_string(["a", "b"], ["b", "c"])
  assert_eq(result.length(), 1)
}

test "union_int" {
  let result = @src.union_int([1, 2], [2, 3])
  assert_eq(result.length(), 3)
}

test "setdiff_int" {
  let result = @src.setdiff_int([1, 2, 3], [2])
  assert_eq(result.length(), 2)
}

test "table_int" {
  let result = @src.table_int([1, 2, 1, 3, 1])
  assert_eq(result.get("1").unwrap(), 3)
}

test "table_string" {
  let result = @src.table_string(["a", "b", "a"])
  assert_eq(result.get("a").unwrap(), 2)
}

test "rev_int" {
  let result = @src.rev_int([1, 2, 3])
  assert_eq(result[0], 3)
}

test "rev_string" {
  let result = @src.rev_string(["a", "b", "c"])
  assert_eq(result[0], "c")
}

test "rep_int" {
  let result = @src.rep_int(5, 3)
  assert_eq(result.length(), 3)
}

test "seq_int" {
  let result = @src.seq_int(1, 5)
  assert_eq(result.length(), 5)
}

test "seq_double" {
  let result = @src.seq_double(1.0, 3.0, by=1.0)
  assert_eq(result.length(), 3)
}

test "head_int" {
  let result = @src.head_int([1, 2, 3, 4, 5], 3)
  assert_eq(result.length(), 3)
}

test "tail_int" {
  let result = @src.tail_int([1, 2, 3, 4, 5], 2)
  assert_eq(result.length(), 2)
}

test "length_int" {
  assert_eq(@src.length_int([1, 2, 3]), 3)
}

test "cbind" {
  let result = @src.cbind([[1.0, 2.0], [3.0, 4.0]])
  assert_eq(result.length(), 2)
}

test "rbind" {
  let result = @src.rbind([[1.0, 2.0], [3.0, 4.0]])
  assert_eq(result.length(), 2)
}
