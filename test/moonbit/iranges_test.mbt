///|
test "iranges_create" {
  let ir = @src.iranges([1, 5], [10, 15])

  assert_eq(@src.ir_start(ir).length(), 2)
  assert_eq(@src.ir_end(ir).length(), 2)
  assert_eq(@src.ir_width(ir).length(), 2)
}

///|
test "iranges_shift" {
  let ir = @src.iranges([1, 5], [10, 15])
  let shifted = @src.ir_shift(ir, 3)

  assert_eq(shifted.starts[0], 4)
  assert_eq(shifted.ends[0], 13)
}

///|
test "iranges_resize_start" {
  let ir = @src.iranges([1, 5], [10, 15])
  let resized = @src.ir_resize(ir, 3, "start")

  assert_eq(resized.widths[0], 3)
}

///|
test "iranges_resize_end" {
  let ir = @src.iranges([1, 5], [10, 15])
  let resized = @src.ir_resize(ir, 3, "end")

  assert_eq(resized.widths[0], 3)
}

///|
test "iranges_reduce" {
  let ir = @src.iranges([1, 3, 8], [5, 6, 10])
  let reduced = @src.ir_reduce(ir, 0)

  assert_eq(reduced.starts.length(), 2)
}

///|
test "iranges_intersect" {
  let ir1 = @src.iranges([1, 8], [10, 15])
  let ir2 = @src.iranges([5, 12], [12, 20])
  let intersect = @src.ir_intersect(ir1, ir2)

  assert_eq(intersect.starts.length(), 3)
}

///|
test "iranges_union" {
  let ir1 = @src.iranges([1, 15], [10, 20])
  let ir2 = @src.iranges([5, 25], [12, 30])
  let union = @src.ir_union(ir1, ir2)

  assert_eq(union.starts.length(), 3)
}

///|
test "iranges_setdiff" {
  let ir1 = @src.iranges([1, 15], [20, 25])
  let ir2 = @src.iranges([5, 18], [12, 22])
  let diff = @src.ir_setdiff(ir1, ir2)

  assert_eq(diff.starts.length(), 0)
}

///|
test "iranges_overlaps" {
  let ir1 = @src.iranges([1, 20], [10, 30])
  let ir2 = @src.iranges([5, 25], [15, 35])
  let overlaps = @src.ir_overlaps(ir1, ir2)

  assert_eq(overlaps.length(), 2)
}

///|
test "iranges_count_overlaps" {
  let query = @src.iranges([1, 15], [10, 20])
  let subject = @src.iranges([5, 18], [12, 22])
  let counts = @src.ir_count_overlaps(query, subject)

  assert_eq(counts.length(), 2)
}

///|
test "iranges_single" {
  let ir = @src.iranges_single(1, 10)

  assert_eq(@src.ir_start(ir)[0], 1)
  assert_eq(@src.ir_end(ir)[0], 10)
}

///|
test "iranges_narrow" {
  let ir = @src.iranges([1, 5], [10, 15])
  let narrowed = @src.ir_narrow(ir, 2, 5)

  assert_eq(narrowed.widths[0], 4)
}

///|
test "iranges_flank" {
  let ir = @src.iranges([10, 20], [15, 25])
  let flanked = @src.ir_flank(ir, 2, true, false)

  assert_eq(flanked.starts[0], 8)
}

///|
test "iranges_disjoin" {
  let ir = @src.iranges([1, 3, 8], [5, 6, 10])
  let disjoined = @src.ir_disjoin(ir)

  assert_eq(disjoined.starts.length(), 4)
}
