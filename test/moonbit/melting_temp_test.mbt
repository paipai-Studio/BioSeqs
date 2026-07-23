///|
test "mt_wallace_basic" {
  let seq = "AGCT"
  let result = @src.mt_wallace(seq)
  assert_eq(result.method, "wallace")
  assert_eq(result.tm, 12.0)
}

///|
test "mt_wallace_gc_rich" {
  let seq = "GGCC"
  let result = @src.mt_wallace(seq)
  assert_eq(result.tm, 16.0)
}

///|
test "mt_wallace_at_rich" {
  let seq = "AATT"
  let result = @src.mt_wallace(seq)
  assert_eq(result.tm, 8.0)
}

///|
test "mt_gc_content_basic" {
  let seq = "AGCTGATCGATCGATCGATCGATCGATCGATCGATCG"
  let result = @src.mt_gc_content(seq)
  assert_eq(result.method, "gc_content")
  assert_true(result.tm > 0.0)
}

///|
test "mt_gc_content_short" {
  let seq = "AGCT"
  let result = @src.mt_gc_content(seq)
  assert_eq(result.method, "gc_content")
}

///|
test "recommend_tm_method_short" {
  let seq = "AGCT"
  let method = @src.recommend_tm_method(seq)
  assert_eq(method, "wallace")
}

///|
test "recommend_tm_method_long" {
  let seq = "AGCTGATCGATCGATCGATCGATCGATCGATCGATCG"
  let method = @src.recommend_tm_method(seq)
  assert_eq(method, "gc_content")
}

///|
test "create_example_dna" {
  let seq = @src.create_example_dna()
  assert_true(seq.length() > 0)
}