test "pseudoknot_order_nested" {
  // Simple nested structure: (((...))) — no pseudoknot
  let order = @src.pseudoknot_order("(((...)))")
  assert_eq(order, 0)
}

test "pseudoknot_order_simple_pk" {
  // Simple pseudoknot: ((...[[...))...))
  // Pairs: (0,7), (1,6) and (4,11), (5,10)
  // These cross: 0 < 4 < 7 < 11
  let order = @src.pseudoknot_order("((..[[...))..]]")
  assert_eq(order, 1)
}

test "pseudoknot_order_no_pairs" {
  let order = @src.pseudoknot_order("....")
  assert_eq(order, 0)
}

test "pseudoknot_order_levels" {
  // Nested: all order 0
  let levels = @src.pseudoknot_order_levels("(((...)))")
  for l in levels {
    assert_eq(l, 0)
  }
}
