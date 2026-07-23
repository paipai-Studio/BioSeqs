///|
/// Tests for Prosite module.

test "PrositePattern creation" {
  let pattern = @src.PrositePattern::new("PS00001", "ATP-binding", "[AG]-X(4)-G-K-[ST]")
  assert_eq(pattern.accession, "PS00001")
  assert_eq(pattern.name, "ATP-binding")
  assert_eq(pattern.pattern, "[AG]-X(4)-G-K-[ST]")
}

test "PrositeMatch creation" {
  let match_ = @src.PrositeMatch::new("PS00001", "ATP-binding", 1, 10, "AGXXXXGKST")
  assert_eq(match_.pattern_accession, "PS00001")
  assert_eq(match_.start, 1)
  assert_eq(match_.end, 10)
}

test "PrositeEntry creation" {
  let entry = @src.PrositeEntry::new("PS00001", "ATP/GTP-binding")
  assert_eq(entry.accession, "PS00001")
  assert_eq(entry.name, "ATP/GTP-binding")
}

test "prosite_pattern_to_regex" {
  let regex = @src.prosite_pattern_to_regex("[AG]-X(4)-G-K-[ST]")
  assert_true(regex.length() > 0)
}

test "prosite_search simple pattern" {
  let matches = @src.prosite_search("ST", "AASTKKST")
  assert_true(matches.length() >= 2)
}

test "prosite_search bracket pattern" {
  let matches = @src.prosite_search("[ST]", "ASTCG")
  assert_true(matches.length() >= 2)
}

test "prosite_search with wildcards" {
  let matches = @src.prosite_search("X(2)", "AAAA")
  assert_true(matches.length() >= 3)
}

test "prosite_scan" {
  let patterns = @src.prosite_create_example_patterns()
  let matches = @src.prosite_scan("AASTKKST", patterns)
  assert_true(matches.length() >= 1)
}

test "prosite_get_pattern" {
  let pattern = @src.prosite_get_pattern("PS00001")
  assert_true(pattern is Some(_))
}

test "prosite_calculate_score" {
  let match_ = @src.PrositeMatch::new("PS00001", "Test", 1, 5, "AAAAA")
  let pattern = @src.PrositePattern::new("PS00001", "Test", "AAAAA")
  let score = @src.prosite_calculate_score(match_, pattern)
  assert_true(score >= 90.0)
}

test "prosite_create_example_patterns" {
  let patterns = @src.prosite_create_example_patterns()
  assert_eq(patterns.length(), 3)
}

test "prosite_create_example_entry" {
  let entry = @src.prosite_create_example_entry()
  assert_eq(entry.accession, "PS00001")
}