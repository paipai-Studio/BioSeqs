///|
/// Tests for the trie (prefix tree) data structure module (trie.mbt).
/// Inspired by Biopython's Bio.trie and Bio.triefind modules.

// ============================================================================
// Trie Construction and Basic Operations
// ============================================================================

///|
test "trie_new_creates_empty_trie" {
  let trie = @src.Trie::new()
  assert_eq(@src.trie_size(trie), 0)
  assert_eq(@src.trie_keys(trie).length(), 0)
  assert_eq(@src.trie_items(trie).length(), 0)
}

///|
test "trie_insert_single_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  assert_eq(@src.trie_size(trie), 1)
  assert_eq(@src.trie_has_key(trie, "GAATTC"), true)
}

///|
test "trie_insert_multiple_keys" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  @src.trie_insert(trie, "AAGCTT", "HindIII")
  assert_eq(@src.trie_size(trie), 3)
}

///|
test "trie_insert_overwrite_existing" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GAATTC", "EcoRI-HF")
  // Size should not increase when overwriting
  assert_eq(@src.trie_size(trie), 1)
  // Value should be updated
  assert_eq(@src.trie_lookup(trie, "GAATTC").unwrap(), "EcoRI-HF")
}

///|
test "trie_insert_empty_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "", "empty")
  assert_eq(@src.trie_size(trie), 1)
  assert_eq(@src.trie_lookup(trie, "").unwrap(), "empty")
}

// ============================================================================
// Lookup Operations
// ============================================================================

///|
test "trie_lookup_existing_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let val = @src.trie_lookup(trie, "GAATTC")
  assert_eq(val.unwrap(), "EcoRI")
}

///|
test "trie_lookup_nonexistent_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let val = @src.trie_lookup(trie, "GGATCC")
  assert_eq(val.is_some(), false)
}

///|
test "trie_has_key_existing" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  assert_eq(@src.trie_has_key(trie, "GAATTC"), true)
}

///|
test "trie_has_key_nonexistent" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  assert_eq(@src.trie_has_key(trie, "GGATCC"), false)
}

///|
test "trie_has_key_prefix_only" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  // "GAATT" is a prefix but not a key
  assert_eq(@src.trie_has_key(trie, "GAATT"), false)
}

// ============================================================================
// Prefix Operations
// ============================================================================

///|
test "trie_has_prefix_existing" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  assert_eq(@src.trie_has_prefix(trie, "GAATT"), true)
  assert_eq(@src.trie_has_prefix(trie, "GAA"), true)
  assert_eq(@src.trie_has_prefix(trie, "G"), true)
}

///|
test "trie_has_prefix_nonexistent" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  assert_eq(@src.trie_has_prefix(trie, "GGAT"), false)
}

///|
test "trie_with_prefix_finds_all" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "ATCGATCGAT", "P1")
  @src.trie_insert(trie, "ATCGATCGAC", "P2")
  @src.trie_insert(trie, "GCTAGCTAGC", "P3")
  let results = @src.trie_with_prefix(trie, "ATCGATCGA")
  assert_eq(results.length(), 2)
}

///|
test "trie_with_prefix_no_match" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let results = @src.trie_with_prefix(trie, "GGAT")
  assert_eq(results.length(), 0)
}

///|
test "trie_with_prefix_empty_prefix_returns_all" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  let results = @src.trie_with_prefix(trie, "")
  assert_eq(results.length(), 2)
}

// ============================================================================
// Keys and Items
// ============================================================================

///|
test "trie_keys_returns_all_keys" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  @src.trie_insert(trie, "AAGCTT", "HindIII")
  let keys = @src.trie_keys(trie)
  assert_eq(keys.length(), 3)
}

///|
test "trie_items_returns_all_pairs" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  let items = @src.trie_items(trie)
  assert_eq(items.length(), 2)
  // Verify each item has correct key and value
  let mut found_ecori = false
  let mut found_bamhi = false
  for item in items {
    if item.0 == "GAATTC" {
      assert_eq(item.1, "EcoRI")
      found_ecori = true
    }
    if item.0 == "GGATCC" {
      assert_eq(item.1, "BamHI")
      found_bamhi = true
    }
  }
  assert_eq(found_ecori, true)
  assert_eq(found_bamhi, true)
}

// ============================================================================
// Remove and Clear
// ============================================================================

///|
test "trie_remove_existing_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  let removed = @src.trie_remove(trie, "GAATTC")
  assert_eq(removed, true)
  assert_eq(@src.trie_size(trie), 1)
  assert_eq(@src.trie_has_key(trie, "GAATTC"), false)
}

///|
test "trie_remove_nonexistent_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let removed = @src.trie_remove(trie, "GGATCC")
  assert_eq(removed, false)
  assert_eq(@src.trie_size(trie), 1)
}

///|
test "trie_remove_prefix_not_key" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let removed = @src.trie_remove(trie, "GAATT")
  assert_eq(removed, false)
  assert_eq(@src.trie_size(trie), 1)
}

///|
test "trie_clear_removes_all" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GGATCC", "BamHI")
  @src.trie_clear(trie)
  assert_eq(@src.trie_size(trie), 0)
  assert_eq(@src.trie_keys(trie).length(), 0)
}

// ============================================================================
// Triefind - Matching trie keys in strings
// ============================================================================

///|
test "triefind_match_longest_prefix" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GAATT", "ShortSite")
  let result = @src.triefind_match(trie, "GAATTCGGG")
  assert_eq(result.is_some(), true)
  let (key, val) = result.unwrap()
  assert_eq(key, "GAATTC")
  assert_eq(val, "EcoRI")
}

///|
test "triefind_match_no_match" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let result = @src.triefind_match(trie, "GGGGGGG")
  assert_eq(result.is_some(), false)
}

///|
test "triefind_match_empty_text" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let result = @src.triefind_match(trie, "")
  assert_eq(result.is_some(), false)
}

///|
test "triefind_match_all_multiple_matches" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  @src.trie_insert(trie, "GAATT", "ShortSite")
  @src.trie_insert(trie, "GAA", "VeryShort")
  let results = @src.triefind_match_all(trie, "GAATTCGGG")
  assert_eq(results.length(), 3)
}

///|
test "triefind_match_all_no_match" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let results = @src.triefind_match_all(trie, "TTTTTT")
  assert_eq(results.length(), 0)
}

///|
test "triefind_find_anywhere_in_string" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let results = @src.triefind_find(trie, "GGGGGAATTCGGG")
  assert_eq(results.length(), 1)
  let (key, val, start, end) = results[0]
  assert_eq(key, "GAATTC")
  assert_eq(val, "EcoRI")
  assert_eq(start, 4)
  assert_eq(end, 10)
}

///|
test "triefind_find_multiple_occurrences" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let results = @src.triefind_find(trie, "GAATTCGAATTC")
  assert_eq(results.length(), 2)
  assert_eq(results[0].2, 0)
  assert_eq(results[1].2, 6)
}

///|
test "triefind_find_no_match" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "GAATTC", "EcoRI")
  let results = @src.triefind_find(trie, "GGGGGGGGG")
  assert_eq(results.length(), 0)
}

///|
test "triefind_find_words_with_boundaries" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "EcoRI", "enzyme1")
  @src.trie_insert(trie, "BamHI", "enzyme2")
  let results = @src.triefind_find_words(trie, "Use EcoRI and BamHI for cloning")
  assert_eq(results.length(), 2)
}

///|
test "triefind_find_words_no_word_boundary" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "EcoRI", "enzyme1")
  // "EcoRI" is part of "EcoRII" without word boundary after "EcoRI"
  let results = @src.triefind_find_words(trie, "Use EcoRII for cloning")
  assert_eq(results.length(), 0)
}

// ============================================================================
// Sample Data
// ============================================================================

///|
test "trie_sample_restriction_sites" {
  let trie = @src.trie_sample_restriction_sites()
  assert_eq(@src.trie_size(trie), 7)
  assert_eq(@src.trie_has_key(trie, "GAATTC"), true)
  assert_eq(@src.trie_has_key(trie, "GGATCC"), true)
  assert_eq(@src.trie_has_key(trie, "AAGCTT"), true)
  assert_eq(@src.trie_has_key(trie, "CTGCAG"), true)
  assert_eq(@src.trie_has_key(trie, "GTCGAC"), true)
  assert_eq(@src.trie_has_key(trie, "CCCGGG"), true)
  assert_eq(@src.trie_has_key(trie, "GCGGCCGC"), true)
}

///|
test "trie_sample_primers" {
  let trie = @src.trie_sample_primers()
  assert_eq(@src.trie_size(trie), 4)
  assert_eq(@src.trie_lookup(trie, "ATCGATCGAT").unwrap(), "Primer_Fwd_1")
  assert_eq(@src.trie_lookup(trie, "ATCGATCGAC").unwrap(), "Primer_Fwd_2")
  assert_eq(@src.trie_lookup(trie, "GCTAGCTAGC").unwrap(), "Primer_Rev_1")
  assert_eq(@src.trie_lookup(trie, "TTTTGGGGCC").unwrap(), "Primer_Rev_2")
}

///|
test "trie_restriction_site_scan_sequence" {
  let trie = @src.trie_sample_restriction_sites()
  // Sequence with EcoRI site at position 4
  let results = @src.triefind_find(trie, "GGGGGAATTCGGGG")
  assert_eq(results.length(), 1)
  assert_eq(results[0].0, "GAATTC")
  assert_eq(results[0].2, 4)
}

///|
test "trie_prefix_overlapping_keys" {
  let trie = @src.Trie::new()
  @src.trie_insert(trie, "ATCG", "short")
  @src.trie_insert(trie, "ATCGATCG", "long")
  assert_eq(@src.trie_has_key(trie, "ATCG"), true)
  assert_eq(@src.trie_has_key(trie, "ATCGATCG"), true)
  assert_eq(@src.trie_has_prefix(trie, "ATC"), true)
  let prefix_results = @src.trie_with_prefix(trie, "ATCG")
  assert_eq(prefix_results.length(), 2)
}

///|
test "trie_size_reflects_operations" {
  let trie = @src.Trie::new()
  assert_eq(@src.trie_size(trie), 0)
  @src.trie_insert(trie, "AAA", "1")
  assert_eq(@src.trie_size(trie), 1)
  @src.trie_insert(trie, "BBB", "2")
  assert_eq(@src.trie_size(trie), 2)
  assert_eq(@src.trie_remove(trie, "AAA"), true)
  assert_eq(@src.trie_size(trie), 1)
  @src.trie_clear(trie)
  assert_eq(@src.trie_size(trie), 0)
}
