///|
/// Tests for Suffix Array and Suffix Tree.

///|
test "suffix_array_banana" {
  let sa = @src.SuffixArray::new("banana")
  let indices = sa.sa()
  assert_eq(indices.length(), 7)
  assert_true(sa.get_sorted_suffix(0) == "$")
}

///|
test "suffix_array_sorted_suffixes" {
  let sa = @src.SuffixArray::new("banana")
  let sorted = sa.get_sorted_suffix(1)
  assert_true(sorted.length() > 0)
}

///|
test "lcp_array_banana" {
  let sa = @src.SuffixArray::new("banana")
  let lcp = @src.LCPArray::new(sa)
  assert_eq(lcp.len(), 7)
  let max_lcp = lcp.max_lcp()
  assert_true(max_lcp >= 0)
}

///|
test "suffix_array_find_lower" {
  let sa = @src.SuffixArray::new("banana")
  let lower = sa.find_lower("ana")
  assert_true(lower >= 0)
}

///|
test "suffix_array_find_upper" {
  let sa = @src.SuffixArray::new("banana")
  let upper = sa.find_upper("ana")
  assert_true(upper >= 0)
}

///|
test "suffix_array_count" {
  let sa = @src.SuffixArray::new("banana")
  assert_eq(sa.count("ana"), 2)
  assert_eq(sa.count("an"), 2)
  assert_eq(sa.count("b"), 1)
  assert_eq(sa.count("xyz"), 0)
}

///|
test "suffix_array_locate" {
  let sa = @src.SuffixArray::new("banana")
  let pos = sa.locate("ana")
  assert_eq(pos.length(), 2)
}

///|
test "suffix_array_empty" {
  let sa = @src.SuffixArray::new("")
  assert_eq(sa.len(), 1)
}

///|
test "lcp_array_basic" {
  let sa = @src.SuffixArray::new("banana")
  let lcp = @src.LCPArray::new(sa)
  assert_eq(lcp.len(), 7)
  let max_lcp = lcp.max_lcp()
  assert_true(max_lcp > 0)
}

///|
test "suffix_tree_basic" {
  let st = @src.SuffixTree::new("banana")
  assert_eq(st.text(), "banana")
  assert_true(st.size() > 0)
}

///|
test "suffix_tree_contains" {
  let st = @src.SuffixTree::new("banana")
  assert_true(st.contains("ana"))
  assert_true(st.contains("ban"))
  assert_true(st.contains("nana"))
  assert_false(st.contains("xyz"))
}

///|
test "suffix_tree_count" {
  let st = @src.SuffixTree::new("banana")
  assert_eq(st.count("ana"), 2)
  assert_eq(st.count("an"), 2)
  assert_eq(st.count("b"), 1)
  assert_eq(st.count("xyz"), 0)
}

///|
test "suffix_tree_locate" {
  let st = @src.SuffixTree::new("banana")
  let pos = st.locate("ana")
  assert_eq(pos.length(), 2)
}

///|
test "suffix_tree_longest_repeated" {
  let st = @src.SuffixTree::new("banana")
  let lrs = st.longest_repeated_substring()
  assert_true(lrs.length() > 0)
}

///|
test "suffix_tree_to_string" {
  let st = @src.SuffixTree::new("abc")
  let str = st.to_string()
  assert_true(str.length() > 0)
}

///|
test "suffix_tree_to_dot" {
  let st = @src.SuffixTree::new("abc")
  let dot = st.to_dot()
  assert_true(dot.length() > 0)
}

///|
test "suffix_tree_root" {
  let st = @src.SuffixTree::new("abc")
  let root = st.root()
  assert_eq(root, "root")
}

///|
test "suffix_array_dna" {
  let sa = @src.SuffixArray::new("ACGTACGT")
  assert_eq(sa.count("ACGT"), 2)
  let pos = sa.locate("ACGT")
  assert_eq(pos.length(), 2)
}

///|
test "suffix_tree_dna" {
  let st = @src.SuffixTree::new("ACGTACGT")
  assert_eq(st.count("ACGT"), 2)
  assert_true(st.contains("CGTA"))
  assert_false(st.contains("CGTT"))
}

///|
test "suffix_array_single_char" {
  let sa = @src.SuffixArray::new("a")
  assert_eq(sa.len(), 2)
  assert_eq(sa.count("a"), 1)
}

///|
test "suffix_tree_single_char" {
  let st = @src.SuffixTree::new("a")
  assert_true(st.contains("a"))
  assert_eq(st.count("a"), 1)
}

///|
test "suffix_tree_small_string" {
  let st = @src.SuffixTree::new("ab")
  assert_true(st.contains("a"))
  assert_true(st.contains("b"))
  assert_true(st.contains("ab"))
}

///|
test "suffix_tree_aaa" {
  let st = @src.SuffixTree::new("aaa")
  assert_true(st.contains("a"))
  assert_eq(st.count("a"), 3)
}
