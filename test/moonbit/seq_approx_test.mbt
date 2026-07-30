///|
test "count_mismatches" {
  assert_eq(@src.count_mismatches("ACGT", "ACGT"), 0)
  assert_eq(@src.count_mismatches("ACGT", "ACGA"), 1)
  assert_eq(@src.count_mismatches("ACGT", "TTTT"), 3)
  assert_eq(@src.count_mismatches("ACGT", "ACGTAC"), 0) // min length
  assert_eq(@src.count_mismatches("ACGT", "acgt"), 0) // case insensitive
}

///|
test "approx_search_exact" {
  let m = @src.approx_search("ACGT", "AACGTG", 0)
  assert_true(m.is_match)
  assert_eq(m.start, 1)
  assert_eq(m.end, 5)
  assert_eq(m.mismatches, 0)
  assert_eq(m.score, 4)
}

///|
test "approx_search_with_mismatches" {
  let m = @src.approx_search("ACGT", "AACGAG", 1)
  assert_true(m.is_match)
  assert_eq(m.mismatches, 1)
}

///|
test "approx_search_not_found" {
  let m = @src.approx_search("ACGT", "AAAAAA", 0)
  assert_true(!m.is_match)
}

///|
test "approx_search_empty" {
  let m = @src.approx_search("", "ACGT", 0)
  assert_true(!m.is_match)
  let m2 = @src.approx_search("ACGT", "", 0)
  assert_true(!m2.is_match)
}

///|
test "approx_match_accessors" {
  let m = @src.approx_match_new("pat", "query", 2, 5, 1, 3, true)
  assert_eq(m.get_pattern(), "pat")
  assert_eq(m.get_query(), "query")
  assert_eq(m.get_start(), 2)
  assert_eq(m.get_end(), 5)
  assert_eq(m.get_mismatches(), 1)
  assert_eq(m.get_score(), 3)
  assert_true(m.is_match())
}

///|
test "levenshtein_distance_identical" {
  assert_eq(@src.levenshtein_distance("ACGT", "ACGT"), 0)
}

///|
test "levenshtein_distance_substitution" {
  assert_eq(@src.levenshtein_distance("ACGT", "ACGA"), 1)
  assert_eq(@src.levenshtein_distance("ACGT", "XXGT"), 2)
}

///|
test "levenshtein_distance_insertion" {
  assert_eq(@src.levenshtein_distance("ACGT", "ACGAT"), 1)
}

///|
test "levenshtein_distance_deletion" {
  assert_eq(@src.levenshtein_distance("ACGT", "AGT"), 1)
}

///|
test "levenshtein_distance_empty" {
  assert_eq(@src.levenshtein_distance("", "ACGT"), 4)
  assert_eq(@src.levenshtein_distance("ACGT", ""), 4)
  assert_eq(@src.levenshtein_distance("", ""), 0)
}

///|
test "levenshtein_distance_case_insensitive" {
  assert_eq(@src.levenshtein_distance("ACGT", "acgt"), 0)
}

///|
test "approx_search_with_indels" {
  let m = @src.approx_search_with_indels("ACGT", "ACXGT", 1)
  assert_true(m.is_match)
  assert_eq(m.mismatches, 1)
}

///|
test "approx_search_with_indels_insertion" {
  let m = @src.approx_search_with_indels("ACGT", "ACGXT", 1)
  assert_true(m.is_match)
}

///|
test "approx_search_with_indels_deletion" {
  let m = @src.approx_search_with_indels("ACGT", "AGT", 1)
  assert_true(m.is_match)
}

///|
test "approx_search_with_indels_not_found" {
  let m = @src.approx_search_with_indels("AAAA", "TCGTCG", 1)
  assert_true(!m.is_match)
}

///|
test "approx_search_with_indels_empty" {
  let m = @src.approx_search_with_indels("", "ACGT", 0)
  assert_true(m.is_match)
  let m2 = @src.approx_search_with_indels("ACGT", "", 0)
  assert_true(!m2.is_match)
}

///|
test "approx_find_all" {
  let results = @src.approx_find_all("ACG", "AACGACGTCG", 0)
  assert_eq(results.length(), 2) // positions 1 and 5

  let results2 = @src.approx_find_all("ACG", "AACGACGTCG", 1)
  assert_true(results2.length() >= 2)
}

///|
test "approx_find_all_empty" {
  let results = @src.approx_find_all("ACGT", "A", 0)
  assert_eq(results.length(), 0)
}

///|
test "approx_word_search_found" {
  let w = @src.approx_word_search("ACGT", "AACGTG", 0)
  assert_true(w.is_found)
  assert_eq(w.edit_distance, 0)
  assert_eq(w.match_start, 1)
}

///|
test "approx_word_search_with_errors" {
  let w = @src.approx_word_search("ACGT", "AACXGTG", 1)
  assert_true(w.is_found)
}

///|
test "approx_word_search_not_found" {
  let w = @src.approx_word_search("ZZZZ", "ACGTACGT", 0)
  assert_true(!w.is_found)
}

///|
test "approx_word_search_empty_word" {
  let w = @src.approx_word_search("", "ACGT", 0)
  assert_true(!w.is_found)
}

///|
test "approx_best_match" {
  let m = @src.approx_best_match("ACGT", "ACXGTACGT", 0.25)
  assert_true(m.is_match)
}

///|
test "approx_match_to_string" {
  let m = @src.approx_match_new("pat", "query", 0, 4, 0, 4, true)
  let s = m.to_string()
  assert_true(s.length() > 0)
}
